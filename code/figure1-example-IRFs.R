##### Setup -----
## Load required packages and 
library(mvtnorm) ## For drawing from a multivariate normal for the GPs
library(bggum)   ## For GGUM response probabilities
## Source in custom functions and ensure needed directories exist
source("code/covSEiso.R") ## For the squared exponential covariance function
source("code/fix_directories.R") ## To conditionally make needed directories
fix_directories()


##### Generate standard IRFs -----
theta  <- seq(-3, 3, 0.01)
twoPL  <- plogis(2.0 * theta)
fourPL <- 0.1 + 0.8 * plogis(2.0 * theta)
ggum   <- ggumProbability(rep(1, length(theta)), theta, 2, 0, c(0, -1))
lpem   <- plogis(3 * theta)^0.3
pal    <- okabe_ito()[c(7, 1, 3, 6)]
pdf("plots/traditional-models.pdf", height = 5, width = 8)
opar <- par(mar = c(3, 3, 1, 1) + 0.1, cex = 1.5)
plot(theta, twoPL, type = "l", lwd = 2, col = pal[1],
     xlab = "", ylab = "", xaxt = "n", yaxt = "n")
lines(theta, fourPL, lty = 2, lwd = 2, col = pal[2])
lines(theta, ggum,   lty = 3, lwd = 2, col = pal[3])
lines(theta, lpem,   lty = 4, lwd = 2, col = pal[4])
axis(side = 1, tick = FALSE, line = -0.75)
axis(side = 2, tick = FALSE, line = -0.75)
mtext(side = 1, line = 1.5, text = expression(theta), cex = 1.5)
mtext(side = 2, line = 1.5, text = "Pr(y = 1)", cex = 1.5)
guide <- c("2PL", "4PL", "GGUM", "LPEF")
legend("topleft", bty = "n", lty = 1:4, lwd = 2, col = pal, legend = guide)
par(opar)
dev.off()


##### Generate GPIRT IRFs -----
K <- covSEiso(theta)
set.seed(123)
GPs <- t(rmvnorm(3, sigma = K))
zero_mean      <- plogis(GPs[ , 1])
linear_mean    <- plogis(GPs[ , 2] + theta)
quadratic_mean <- plogis(GPs[ , 3] - (1/max(theta)) * theta^2)
pdf("plots/gp-irfs.pdf", height = 5, width = 8)
opar <- par(mar = c(3, 3, 1, 1) + 0.1, cex = 1.5)
plot(theta, zero_mean, type = "l", ylim = 0:1, lwd = 2, col = pal[1],
     xlab = "", ylab = "", xaxt = "n", yaxt = "n")
lines(theta, linear_mean, lty = 2, lwd = 2, col = pal[2])
lines(theta, quadratic_mean, lty = 3, lwd = 2, col = pal[3])
axis(side = 1, tick = FALSE, line = -0.75)
axis(side = 2, tick = FALSE, line = -0.75)
mtext(side = 1, line = 1.5, text = expression(theta), cex = 1.5)
mtext(side = 2, line = 1.5, text = "Pr(y = 1)", cex = 1.5)
par(opar)
dev.off()
