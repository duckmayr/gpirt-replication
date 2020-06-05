##### Setup -----
## Install version of gpirt with quadratic mean function
temp_lib <- tempdir()
library(devtools)
install_github("duckmayr/gpirt@44a08da", lib = temp_lib)
## Load required packages
library(gpirt, lib.loc = temp_lib) ## For GPIRT sampling
library(dplyr)                     ## For data wrangling
library(tidyr)                     ## For data wrangling
library(emIRT)                     ## For Bayesian IRT EM estimates
## Load custom IRF plotting function
source("code/plot_irf.R")
## Ensure needed directories exist
source("code/fix_directories.R")
fix_directories()


##### Prepare the data -----
## Read in Voteview data on members (retrieved January 4, 2020)
members <- read.csv("data/H116_members.csv")
## Read in response matrix (created in code/prepare-data.R)
responses <- readRDS(file = "data/H116-responses.rds")
## Match up party codes and NOMINATE scores for plotting
bionames      <- with(members, bioname[match(rownames(responses), icpsr)])
party_codes   <- with(members, party_code[match(rownames(responses), icpsr)])
nominate_dim1 <- with(members, nominate_dim1[match(rownames(responses), icpsr)])


##### Sample posterior -----
## Set seed for reproducibility
set.seed(42)
## Sample posterior; the model converges rather quickly, so about 1,000
## samples in practice have been sufficient. Note we use a quadratic mean
## function as we anticipate some non-monotonic items.
## Note: This can take around 3.5 hours
chain <- gpirtMCMC(responses,
                   sample_iterations = 1000,
                   burn_iterations = 0,
                   beta_proposal_sd = 1,
                   mean_function = "quadratic",
                   store_fstar = TRUE)
## Save raw results for further inspection if desired
saveRDS(chain, file = "raw-output/H116-chain.rds")


##### Estimate GPIRT quantities of interest -----
## Estimate latent traits, mean function parameters, probability of observed
## responses (i.e. {sigma(f_i(theta))}), and IRFs
## First for convenience we will create separate objects for the types of draws
theta <- chain[["theta"]]
beta  <- chain[["beta"]]
f     <- chain[["f"]]
fstar <- chain[["fstar"]]
## Now we will need to deal with reflection
theta <- -theta
beta[2, , ] <- -beta[2, , ]
fstar <- fstar[nrow(fstar):1, , ]
## We will need to construct the dense theta* grid
theta_grid  <- seq(-5, 5, 0.01)
## This is a convenience function to get the mean of the sigmoid of f draws
estimate_prob <- function(f_draws) mean(plogis(f_draws))
## These are bookkeeping variables
n <- nrow(responses)
m <- ncol(responses)
## Now we can get the latent trait and mean function parameter estimates
theta_estimates <- apply(theta,  2, mean)
beta_estimates  <- apply(beta, 1:2, mean)
## And the IRF estimates ({pi_i})
irf_estimates <- matrix(NA_real_, ncol = m, nrow = length(theta_grid))
for ( j in 1:ncol(irf_estimates) ) {
    for ( i in 1:nrow(irf_estimates) ) {
        mu <- beta[1, j, ] + theta_grid[i] * beta[2, j, ]
        + theta_grid[i]^2 * beta[3, j, ]
        irf_estimates[i, j] <- estimate_prob(fstar[i, j, ] + mu)
    }
}
## As well as the predicted probabilities
response_prob_estimates <- matrix(NA_real_, ncol = m, nrow = n)
for ( j in 1:ncol(response_prob_estimates) ) {
    for ( i in 1:nrow(response_prob_estimates) ) {
        mu <- beta[1, j, ] + theta[ , i] * beta[2, j, ]
        + theta[ , i]^2 * beta[3, j, ]
        response_prob_estimates[i, j] <- estimate_prob(mu + f[i, j, ])
    }
}
## And save it all to disk
save(theta_estimates, beta_estimates, response_prob_estimates, irf_estimates,
     file = paste0("data/H116-estimates.RData"))


##### Compare GPIRT Ideology & DW-NOMINATE Dimension 1 Ideology (Figure 2) -----
## We will plot Republican ideal points in red & Democrat ideal points in blue,
## except for members of the "Squad," who will be plotted in purple
pal    <- c(R = "#cd2626af", D = "#104e8baf", S = "#cc79a7af")
squad  <- grepl("OCASIO|TLAIB|PRESSLEY|OMAR", bionames)
ptcols <- ifelse(squad, pal["S"], ifelse(party_codes == 200, pal["R"], pal["D"]))
## We will plot Republican ideal points with circles & Democrat ideal points
## with squares, except Squad" members, who will be plotted with triangles
pttype <- ifelse(squad, 17, ifelse(party_codes == 200, 19, 15))
## Now we can make the plot
pdf("plots/gpirt-vs-nominate.pdf", width = 8, height = 5)
opar <- par(mar = c(3, 3, 1, 1) + 0.1, cex = 1.5)
plot(nominate_dim1, theta_estimates, pch = pttype, col = ptcols,
     xlab = "", ylab = "", xaxt = "n", yaxt = "n")
axis(side = 1, tick = FALSE, line = -0.75)
axis(side = 2, tick = FALSE, line = -0.75)
mtext(side = 1, line = 1.5, text = "DW-NOMINATE Dimension 1 Ideology", cex = 1.5)
mtext(side = 2, line = 1.5, text = "GPIRT Ideology", cex = 1.5)
legend("topleft", pch = c(19, 15, 17), col = pal, bty = "n",
       legend = c("Republicans", "Democrats", "\"The Squad\""))
par(opar)
dev.off()


##### Plot IRFs from Figure 3 -----
## Get 2PL estimates
## First we need to restructure the responses to how emIRT expects them
mat_responses  <- matrix(responses, dimnames = dimnames(responses),
                         nrow = nrow(responses), ncol = ncol(responses))
pscl_responses <- pscl::rollcall(mat_responses, nay = -1)
rc   <- convertRC(pscl_responses)
## Now we need to set up the priors and start values
## (and so we need to set the seed)
set.seed(63130)
p    <- makePriors(rc$n, rc$m, 1)
s    <- getStarts(rc$n, rc$m, 1)
ctrl <- list(threads = 4, verbose = TRUE)
irt_fit <- binIRT(.rc = rc, .starts = s, .priors = p, .control = ctrl)
## Subset GPIRT IRF estimates
idx        <- which(theta_grid >= -2.5 & theta_grid <= 2.5)
theta_grid <- theta_grid[idx]
irfs       <- irf_estimates[idx, ]
## Extract 2PL estimates (again, need to deal with reflection)
irt_theta  <- irt_fit$means$x
irt_beta   <- irt_fit$means$beta
irt_theta  <- -irt_theta
irt_beta[ , 2] <- -irt_beta[ , 2]
irt_grid   <- seq(floor(min(irt_theta)), ceiling(max(irt_theta)))
## Find members of "The Squad"
squad_icpsrs <- c(21955, 21975, 21950, 21949)
squad_idx <- which(rownames(responses) %in% squad_icpsrs)
## HJRES31 -- Shutdown avoidance vote from 2019/02/14 discussed in GGUM draft
## GPIRT IRF
pdf("plots/hjres31-IRF.pdf", height = 5, width = 8)
plot_irf(theta_grid, irfs[ , 70], responses[ , 70], theta_estimates,
         highlight = squad_idx)
dev.off()
## 2PL IRF
pdf("plots/hjres31-IRF-2PL.pdf", height = 5, width = 8)
plot_irf(irt_grid, irfs[ , 70], responses[ , 70], irt_theta,
         beta = irt_beta[70, ], highlight = squad_idx)
dev.off()
## SAFE ACT -- Dems' election security bill
pdf("plots/safe-act-IRF.pdf", height = 5, width = 8)
plot_irf(theta_grid, irfs[ , 384], responses[ , 384], theta_estimates,
         highlight = squad_idx)
dev.off()
pdf("plots/safe-act-IRF-2PL.pdf", height = 5, width = 8)
plot_irf(irt_grid, irfs[ , 384], responses[ , 384], irt_theta,
         beta = irt_beta[384, ], highlight = squad_idx)
dev.off()
## Cyprus amendment
pdf("plots/cyprus-IRF.pdf", height = 5, width = 8)
plot_irf(theta_grid, irfs[ , 406], responses[ , 406], theta_estimates,
         highlight = squad_idx)
dev.off()
pdf("plots/cyprus-IRF-2PL.pdf", height = 5, width = 8)
plot_irf(irt_grid, irfs[ , 406], responses[ , 406], irt_theta,
         beta = irt_beta[406, ], highlight = squad_idx)
dev.off()
