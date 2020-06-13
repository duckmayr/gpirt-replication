##### Setup -----
## Install version of gpirt with linear mean function
temp_lib <- tempdir()
library(devtools)
install_github("duckmayr/gpirt@e4586a1", lib = temp_lib)
## Load required packages
library(gpirt, lib.loc = temp_lib) ## For GPIRT sampling
library(ltm)                       ## For 2PL estimation
library(pROC)                      ## For area under the ROC
library(dplyr)                     ## For data summarisation
library(tidyr)                     ## For data reshaping
library(KernSmoothIRT)             ## For kernel-smoothed IRT estimation
## Ensure needed directories exist
source("code/fix_directories.R")
fix_directories()
## Source in functions for adaptive testing
Rcpp::sourceCpp("code/gpirt-cat.cpp")
source("code/cat-helpers.R")


##### Prepare the data -----
## Read in training responses and test responses for CAT
## (created in code/prepare-data.R)
responses <- readRDS(file = "data/npi-responses.rds")
cat_responses <- readRDS(file = "data/npi-cat-responses.rds")
## Put training responses in {0, 1, NA} for ltm
ltm_responses <- responses
ltm_responses[ltm_responses == -1] <- 0
## For KernSmoothIRT, we recode responses so that a "1" response is always
## a narcissistic response
key <- c(1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1,
         0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1)
ksIRTresponses <- ifelse(responses == -1, 0, responses)
ksIRTresponses <- sapply(seq_along(key), function(j) {
    if ( key[j] ) {
        return(ksIRTresponses[ , j])
    }
    return(1 - ksIRTresponses[ , j])
})


##### Get GPIRT estimates -----
## Set seed for reproducibility
set.seed(271828)
## Generate samples
chain <- gpirtMCMC(responses, sample_iterations = 2500, burn_iterations = 0,
                   store_fstar = TRUE)
## Save raw output for later inspection if desired
saveRDS(chain, file = "model-output/NPI-chain.rds")
## Get IRF estimates
theta_grid  <- seq(-5, 5, 0.01)
beta  <- chain[["beta"]]
fstar <- chain[["fstar"]]
theta <- chain[["theta"]]
f     <- chain[["f"]]
irf_estimates <- matrix(NA_real_, ncol = ncol(fstar), nrow = nrow(fstar))
iter <- 0
total <- prod(dim(irf_estimates))
for ( j in 1:ncol(irf_estimates) ) {
    for ( i in 1:nrow(irf_estimates) ) {
        iter <- iter + 1
        # cat("\r\t\t", iter, "/", total) ## Uncomment for progress updates
        ## Since the model converges quickly, we just use every sample
        ## drawn except the first variable initialization sample
        mu <- beta[1, j, -1] + theta_grid[i] * beta[2, j, -1]
        irf_estimates[i, j] <- mean(plogis(mu + fstar[i, j, -1]))
    }
}
saveRDS(irf_estimates, file = "data/npi-irf-estimates.rds")


##### Active Learning Experiment -----
## Get theta estimates for held out respondents using all items
theta_grid <- seq(-5, 5, 0.01) ## Set up the grid of theta values
prior <- dnorm(theta_grid)     ## Get the prior on those values
lp    <- log(prior)            ## as well as the log prior
theta <- sapply(1:nrow(cat_responses), function(i) {
    ## and use all items to update the prior for the test respondents
    P <- update_posterior(lp, irf_estimates, cat_responses[i, ], theta_grid)
    return(get_estimate(P, theta_grid))
})

## Compare RMSE for 16 item simulation
n <- nrow(cat_responses)
m <- ncol(cat_responses)
cat_ests  <- numeric(n)
ames_ests <- numeric(n)
rand_ests <- numeric(n)
ames <- c(30, 9, 23, 14, 12, 34, 35, 24, 7, 40, 21, 13, 32, 4, 20, 39)
ames_irfs <- irf_estimates[ , ames]
set.seed(1138)
for ( i in 1:n ) {
    data_i       <- cat_responses[i, ]
    ames_data    <- data_i[ames]
    rand         <- sample(1:ncol(irf_estimates), size = 16)
    rand_irfs    <- irf_estimates[ , rand]
    rand_data    <- data_i[rand]
    cat_post     <- gpirt_cat(16, irf_estimates, theta_grid, prior, data_i)
    ames_post    <- update_posterior(lp, ames_irfs, ames_data, theta_grid)
    rand_post    <- update_posterior(lp, rand_irfs, rand_data, theta_grid)
    cat_ests[i]  <- get_estimate(cat_post,  theta_grid)
    ames_ests[i] <- get_estimate(ames_post, theta_grid)
    rand_ests[i] <- get_estimate(rand_post, theta_grid)
}
RMSEs <- data.frame(CAT    = rmse(theta, cat_ests),
                    Fixed  = rmse(theta, ames_ests),
                    Random = rmse(theta, rand_ests),
                    check.names = FALSE)
RMSEs
(RMSEs$Random - RMSEs$Fixed) / RMSEs$Random
(RMSEs$Random - RMSEs$CAT) / RMSEs$Random


##### Held-out Experiment -----
## Rember to run code/gplvm.py *first*!
## We use python's GPy package for the GPLVM replicates,
## so those results need to exist before we can read them into R and compare
## Now let's run the replicates for the models we have in R,
## GPIRT, kernel-smoothed IRT, and traditional 2PL IRT
## Generate a seed for each replicate
set.seed(123)
seeds <- sample(x = 1:1e6, size = 20)
## Create the data frame to store results
results <- expand.grid(iter = 1:20,
                       gpirt_accuracy = 0,
                       gpirt_auc = 0,
                       ltm_accuracy = 0,
                       ltm_auc = 0,
                       ks_accuracy = 0,
                       ks_auc = 0,
                       gplvm_accuracy = 0,
                       gplvm_auc = 0,
                       gpirt_loglik = 0,
                       ltm_loglik = 0,
                       ks_loglik = 0,
                       gplvm_loglik = 0)
## This is where we'll write results to disk
outname <- "model-output/npi-heldout-results.csv"
for ( iter in 1:20 ) {
    ## Uncomment lines with cat() calls for verbose progress updates
    cat("Starting iteration", iter, ":\n")
    # cat("Setting up data...\n")
    ## Determine which responses to treat as held out,
    ## and create a new training set with those responses treated as missing
    set.seed(seeds[iter])
    altered_responses <- responses
    ks_altered_responses <- ksIRTresponses
    idx <- sample(x = seq_along(responses), size = 0.2 * length(responses))
    altered_responses[idx] <- NA_integer_
    ks_altered_responses[idx] <- NA_integer_
    actual <- ifelse(responses[idx] == 1, 1, 0)
    ## Fit the GPIRT, KS IRT, and 2PL models to the new training data
    # cat("Running gpirtMCMC()...\n")
    samples <- gpirtMCMC(altered_responses, 1000, 0)
    # outname <- paste0("model-output/npi-iter-", iter, ".rds")
    # saveRDS(samples, file = outname) ## Uncomment this & above to save output
    # cat("Running ltm()...\n")
    altered_responses2 <- as.matrix(ifelse(altered_responses == 1, 1, 0))
    ltm_fit <- ltm(altered_responses2 ~ z1, IRT.param = TRUE)
    ## cat("Running kernel-smoothed IRT...\n")
    ksIRT_fit <- ksIRT(ks_altered_responses, 1, 1)
    Ptmp <- subjOCC(ksIRT_fit, "MLTheta")
    P <- matrix(NA_real_, nrow = nrow(ksIRTresponses), ncol = ncol(ksIRTresponses))
    for ( j in 1:ncol(P) ) {
        rnumber <- "if"(ksIRTresponses[1,j] == 1, 1, 2)
        P[ , j] <- Ptmp[[j]][rnumber, ]
    }
    ## Load the GPLVM results and held out indices
    ## (the -1 and +1 things you'll see here account for python being 0-indexed
    ## while R is 1-indexed)
    f_base <- paste0("model-output/GPLVM-iter-", iter-1, "-")
    gplvm_mu <- as.matrix(read.csv(paste0(f_base, "mean.csv"), header = FALSE))
    gplvm_s2 <- as.matrix(read.csv(paste0(f_base, "var.csv"),  header = FALSE))
    ii <- as.matrix(read.csv(paste0(f_base, "idx.csv"),  header = FALSE))+1
    gplvm_actual <- ifelse(responses[ii] == 1, 1, 0)
    ## Calculate log likelihood, AUC, and classification accuracy of responses
    ## treated as missing
    # cat("Getting fit statistics...\n")
    ltm_probs   <- fitted(ltm_fit, type = "conditional-probabilities")
    ltm_class   <- ifelse(ltm_probs > 0.5, 1, 0)
    ltm_probs2  <- ifelse(responses == 1, ltm_probs, 1 - ltm_probs)
    ltm_log_lik <- sum(log(ltm_probs2[idx]))
    theta <- samples[["theta"]]
    beta  <- samples[["beta"]]
    f <- samples[["f"]]
    gpirt_log_lik <- 0
    ks_log_lik  <- 0
    gpirt_class <- numeric(length(idx))
    gpirt_probs <- numeric(length(idx))
    ks_actual   <- numeric(length(idx))
    ks_probs    <- numeric(length(idx))
    ks_class    <- numeric(length(idx))
    for ( elem in 1:length(idx) ) {
        k <- idx[elem]
        arr.ind <- arrayInd(k, dim(responses))
        i <- arr.ind[1]
        j <- arr.ind[2]
        ks_actual[elem] <- ksIRTresponses[i, j]
        ks_probs[elem]  <- P[i, j]
        ks_class[elem]  <- "if"(P[i, j] > 0.5, 1, 0)
        ks_log_lik <- ks_log_lik + log("if"(ks_actual[elem] == 1, P[i, j], 1 - P[i, j]))
        mu <- beta[1, j, -1] + theta[-1, i] * beta[2, j, -1]
        gpirt_prob <- mean(plogis(f[i, j, -1] + mu), na.rm = TRUE)
        gpirt_class[elem] <- "if"(gpirt_prob > 0.5, 1, 0)
        gpirt_probs[elem] <- gpirt_prob
        gpirt_prob2 <- "if"(responses[k] == 1, gpirt_prob, 1 - gpirt_prob)
        gpirt_log_lik <- gpirt_log_lik + log(gpirt_prob2)
    }
    idx0 <- idx
    idx <- which(results$iter == iter)
    results$gpirt_loglik[idx]   <- gpirt_log_lik
    results$gpirt_auc[idx]      <- auc(roc(actual, gpirt_probs))
    results$gpirt_accuracy[idx] <- mean(actual == gpirt_class)
    results$ks_loglik[idx]      <- ks_log_lik
    results$ks_auc[idx]         <- auc(roc(ks_actual, ks_probs))
    results$ks_accuracy[idx]    <- mean(ks_actual == ks_class)
    results$gplvm_loglik[idx]   <- sum(log(dnorm(responses[ii], gplvm_mu[ii],
                                                 sqrt(gplvm_s2[ii]))))
    results$gplvm_auc[idx]      <- auc(roc(gplvm_actual, gplvm_mu[ii]))
    results$gplvm_accuracy[idx] <- mean(ifelse(gplvm_mu[ii] > 0, 1, -1)
                                        == responses[ii])
    results$ltm_loglik[idx]     <- ltm_log_lik
    results$ltm_auc[idx]        <- auc(roc(actual, ltm_probs[idx0]))
    results$ltm_accuracy[idx]   <- mean(actual == ltm_class[idx0])
    write.csv(results, file = outname, row.names = FALSE)
    # cat("\n\n")
    rm(samples)
    Sys.sleep(1)
    gc()
    Sys.sleep(1)
    gc()
}

## Mean fit statistics across replicates
results %>%
    pivot_longer(
        cols = contains("_"),
        names_to = c("Model", "Statistic"),
        names_pattern = "(.*)_(.*)",
        values_to = "Value"
    ) %>%
    mutate(Value = ifelse(Statistic == "loglik", Value / 16000, Value)) %>%
    group_by(Model, Statistic) %>%
    summarise(Mean = mean(Value)) %>% ## Everything after this line just
    pivot_wider(names_from = Statistic, values_from = Mean) %>% ## arranges
    select(Model, loglik, auc, accuracy) %>% ## the data to look like the paper
    arrange(factor(Model, levels = c("gpirt", "ltm", "gplvm", "ks")))

## t-test
t.test(results$gpirt_loglik / 16000,
       results$ks_loglik / 16000,
       paired = TRUE)
t.test(results$gpirt_loglik / 16000,
       results$gplvm_loglik / 16000,
       paired = TRUE)
t.test(results$gpirt_loglik / 16000,
       results$ltm_loglik / 16000,
       paired = TRUE)
