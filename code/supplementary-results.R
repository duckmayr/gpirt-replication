##### Setup -----
## Install version of gpirt with quadratic mean function
temp_lib <- tempdir()
library(devtools)
install_github("duckmayr/gpirt@44a08da", lib = temp_lib, upgrade = FALSE)
## Load required packages
library(gpirt, lib.loc = temp_lib) ## For GPIRT sampling
library(dplyr)                     ## For data wrangling
library(tidyr)                     ## For data wrangling
library(pROC)                      ## For area under the ROC
library(emIRT)                     ## For Bayesian IRT EM estimates
## Source in functions for adaptive testing
Rcpp::sourceCpp("code/gpirt-cat.cpp")
source("code/cat-helpers.R")


##### Prepare the data -----
## Read in Voteview data on members (retrieved January 4, 2020)
members <- read.csv("data/H116_members.csv")
## Read in response matrix (created in code/prepare-data.R)
responses <- readRDS(file = "data/H116-responses.rds")


##### Held-out experiment -----
## Figure out which indices are already NA
NA_indices <- which(is.na(responses))
## Prepare the vector of indices that we can make NA
OK_indices <- which(!is.na(responses))
## Generate seeds for each run
set.seed(123)
seeds <- sample(x = 1:1e6, size = 20)
seeds[11] <- seeds[11] + 1 ## Had some problems (member w/ no votes)
## Create an object to hold the results
results <- expand.grid(p = c(0.1, 0.2),
                       iter = 1:20,
                       irt_ll = 0,
                       gpirt_ll = 0,
                       irt_accuracy = 0,
                       gpirt_accuracy = 0,
                       irt_auc = 0,
                       gpirt_auc = 0)
for ( iter in 1:20 ) {
    for ( p in c(0.1, 0.2) ) {
        cat("Starting iteration", iter, "for p = ", p, ":\n")
        
        cat("\tSetting up data...\n")
        set.seed(seeds[iter])
        altered_responses <- responses
        idx <- sample(x = OK_indices, size = p * length(OK_indices))
        altered_responses[idx] <- NA_integer_
        
        cat("\tRunning gpirtMCMC()...\n")
        ## Generate 1000 posterior samples for GPIRT
        samples <- gpirtMCMC(altered_responses, 1000, 0,
                             mean_function = "quadratic")
        ## Save raw samples in case we need to analyze them more later
        # outname <- paste0("model-output/house-p-", p, "-iter-", iter, ".rds")
        # saveRDS(samples, file = outname)
        
        cat("\tRunning binIRT()...\n")
        ## Set up responses, priors, and start values
        response_mat <- altered_responses
        class(response_mat) <- "matrix"
        rc <- rollcall(response_mat, nay = -1)
        rc2 <- convertRC(rc)
        prior <- makePriors(rc2$n, rc2$m, 1)
        s <- getStarts(rc2$n, rc2$m, 1)
        ## Now we can fit binIRT(); you may need to change the # of threads
        irt_fit <- binIRT(rc2, s, prior, .control = list(threads = 6))
        
        cat("\tGetting fit statistics...\n")
        actual <- ifelse(responses[idx] == -1, 0, responses[idx])
        ## Bayesian IRT
        preds  <- t(sapply(1:rc2$n, function(j) {
            sapply(1:rc2$m, function(i) {
                with(irt_fit$means, pnorm(beta[i, 1] + x[j, 1] * beta[i, 2]))
            })
        }))
        probs  <- ifelse(responses == 1, preds, 1 - preds)
        irt_ll <- sum(log(probs[idx]))
        irt_class <- ifelse(preds > 0.5, 1, 0)
        irt_accuracy <- mean(irt_class[idx] == actual, na.rm = TRUE)
        irt_auc <- auc(roc(ifelse(responses[idx] == -1, 0, responses[idx]), preds[idx]))
        ## GPIRT
        theta <- samples[["theta"]]
        beta  <- samples[["beta"]]
        f <- samples[["f"]]
        gpirt_ll <- 0
        preds <- numeric(length(idx))
        for ( K in 1:length(idx) ) {
            k <- idx[K]
            arr.ind <- arrayInd(k, dim(responses))
            i <- arr.ind[1]
            j <- arr.ind[2]
            mu <- beta[1, j, ] + theta[ , i] * beta[2, j, ] +
                theta[ , i]^2 * beta[3, j, ]
            g <- try(f[i, j, ] + mu, silent = TRUE)
            if ( inherits(g, "try-error") ) {
                print(c(i, j, k))
                print(str(f[i, j, ]))
                print(str(mu))
            }
            gpirt_prob <- mean(plogis(g), na.rm = TRUE)
            preds[K]   <- gpirt_prob
            gpirt_prob <- "if"(responses[k] == 1, gpirt_prob, 1 - gpirt_prob)
            gpirt_ll <- gpirt_ll + log(gpirt_prob)
        }
        gpirt_class <- ifelse(preds > 0.5, 1, 0)
        gpirt_accuracy <- mean(gpirt_class == actual, na.rm = TRUE)
        gpirt_auc <- auc(roc(ifelse(responses[idx] == -1, 0, responses[idx]), preds))
        index <- try(which(results$p == p & results$iter == iter))
        results[index, ] <- c(p, iter,
                              irt_ll, gpirt_ll,
                              irt_accuracy, gpirt_accuracy,
                              irt_auc, gpirt_auc)
        write.csv(results, file = "model-output/house-held-out-results.csv",
                  row.names = FALSE)
        cat("\n\n")
    }
}
## Now we can compare the results
results %>%
    filter(!is.na(gpirt_ll) & gpirt_ll != 0) %>%
    pivot_longer(cols = contains("_"),
                 names_to = c("Model", ".value"),
                 names_sep = "_") %>%
    mutate(Model = ifelse(Model == "irt", "2PL", "GPIRT")) %>%
    mutate(N = ifelse(p == 0.1, 26136, 52273)) %>%
    mutate(Mean_LL_divby_N = ll / N) %>%
    select(-iter, -N) %>%
    group_by(p, Model) %>%
    summarize_all(mean)

## t-test
t.test(results$gpirt_ll[results$p == 0.1] / 26136,
       results$irt_ll[results$p == 0.1] / 26136,
       paired = TRUE)
t.test(results$gpirt_ll[results$p == 0.2] / 52273,
       results$irt_ll[results$p == 0.2] / 52273,
       paired = TRUE)


##### Active learning experiment -----
## Select 20% to hold out
set.seed(123)
withheld_members <- sample(rownames(responses), size = 0.2 * nrow(responses))
kept_members <- setdiff(rownames(responses), withheld_members)
training_data <- responses[kept_members, ]
## Find the items we can analyze with only 80% of members
is_lopsided <- function(x, cutoff = 0.01) {
    tab  <- table(x)            ## Tabulate the votes
    prop <- min(tab) / sum(tab) ## Get proportion of minority votes
    return(prop < cutoff)       ## Return TRUE if that proporion is < cutoff
}
lopsided_votes <- which(apply(training_data, 2, is_lopsided))
training_data <- training_data[ , -lopsided_votes]
codes <- list(yea = 1, nay = -1, missing = NA)
training_data <- as.response_matrix(training_data, response_codes = codes)
analyzed_items <- colnames(training_data)
test_data <- responses[withheld_members, analyzed_items]
set.seed(456)
cat("Getting posterior samples using training data...\n")
samples <- gpirtMCMC(training_data,
                     sample_iterations = 1000,
                     burn_iterations = 0,
                     mean_function = "quadratic",
                     store_fstar = TRUE)
# saveRDS(samples, file = "model-output/house-cat-training-samples.rds")
theta_grid <- seq(-5, 5, 0.01)
n <- length(theta_grid)
m <- ncol(samples[["fstar"]])
irf_estimates <- matrix(NA_real_, ncol = m, nrow = n)
theta <- samples[["theta"]]
beta <- samples[["beta"]]
fstar <- samples[["fstar"]]
fstar_elem <- 1
mu_elem <- 1
for ( j in 1:ncol(irf_estimates) ) {
    for ( i in 1:nrow(irf_estimates) ) {
        iter <- iter + 1
        mu <- beta[1, j, ] + theta_grid[i] * beta[2, j, ]
        irf_estimates[i, j] <- mean(plogis(mu + fstar[i, j, ]), na.rm = TRUE)
    }
}
saveRDS(irf_estimates, file = "model-output/house-cat-irf-estimates.rds")
## Get theta estimates for held out respondents using all items
cat("Comparing CAT and random vs. full...\n")
update_post <- function(log_prior, IRFs, responses, theta_grid) {
    m <- length(responses)
    for ( j in 1:m ) {
        if ( is.na(responses[j]) ) {
            next()
        } else if ( responses[j] == 1 ) {
            log_prior <- log_prior + log(IRFs[ , j])
        } else { 
            log_prior <- log_prior + log(1 - IRFs[ , j])
        }
    }
    return(normalize(theta_grid, exp(log_prior)))
}
theta_grid <- seq(-5, 5, 0.01)
prior <- dnorm(theta_grid)
theta <- sapply(1:nrow(test_data), function(i) {
    P <- update_post(log(prior), irf_estimates, test_data[i, ], theta_grid)
    return(get_estimate(P, theta_grid))
})
## Compare RMSE for 20 items
n <- nrow(test_data)
m <- ncol(test_data)
cat_estimates <- numeric(n)
random_estimates <- numeric(n)
set.seed(789)
for ( i in 1:n ) {
    responses <- test_data[i, ]
    random <- sample(1:ncol(irf_estimates), size = 10)
    random_irfs <- irf_estimates[ , random]
    random_responses <- responses[random]
    cat_post <- gpirt_cat(10, irf_estimates, theta_grid, prior, responses)
    random_post <- update_post(log(prior), random_irfs, random_responses, theta_grid)
    cat_estimates[i] <- get_estimate(cat_post, theta_grid)
    random_estimates[i] <- get_estimate(random_post, theta_grid)
}
RMSEs <- data.frame(CAT = rmse(theta, cat_estimates),
                    Random = rmse(theta, random_estimates),
                    check.names = FALSE)
RMSEs
(RMSEs$Random - RMSEs$CAT) / RMSEs$Random
estimates <- data.frame(CAT = cat_estimates,
                        Random = random_estimates,
                        Full = theta)
write.csv(RMSEs, file = "data/house-cat-RMSEs.csv", row.names = FALSE)
write.csv(estimates, file = "data/house-cat-estimates.csv", row.names = FALSE)
