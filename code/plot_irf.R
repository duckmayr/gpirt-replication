#' Plots Item Response Functions
#' 
#' @param theta_grid A numeric vector giving the theta values to plot over
#' @param irf_estimate A numeric vector giving one vote probabilities at the
#'     values of theta_grid. If not provided, can be constructed if beta is
#'     not NULL (then 2PL IRFs are created)
#' @param responses An integer vector of actual responses
#' @param theta_estimate A numeric vector giving the GPIRT theta estimates
#' @param beta (Optional). A numeric vector of 2PL item parameter estimates
#' @param highlight (Optional). An integer vector of respondent indices to
#'     highlight in the rug
plot_irf <- function(theta_grid, irf_estimate, responses, theta_estimate,
                     beta = NULL, highlight = NULL) {
    bigrug <- function(...) rug(..., lwd = 3, ticksize = 0.04)
    opar <- par(mar = c(3, 3, 1, 1) + 0.1, cex = 1.5)
    on.exit(par(opar))
    if ( !is.null(beta) ) {
        irf_estimate <- plogis(beta[1] + theta_grid * beta[2])
    }
    plot(theta_grid, irf_estimate, type = "l", ylim = 0:1,
         xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    axis(side = 1, tick = FALSE, line = -0.75)
    axis(side = 2, tick = FALSE, line = -0.75)
    mtext(side = 1, line = 1.5, text = expression(theta), cex = 1.5)
    mtext(side = 2, line = 1.5, text = "Pr(y = 1)", cex = 1.5)
    yeas <- which(responses ==  1)
    nays <- which(responses == -1)
    if ( !is.null(highlight) ) {
        highlighted_yeas <- yeas[which(yeas %in% highlight)]
        highlighted_nays <- nays[which(nays %in% highlight)]
        yeas <- setdiff(yeas, highlighted_yeas)
        nays <- setdiff(nays, highlighted_nays)
        rug(theta_estimate[yeas], side = 3)
        rug(theta_estimate[nays], side = 1)
        bigrug(theta_estimate[highlighted_yeas], side = 3, col = "red")
        bigrug(theta_estimate[highlighted_nays], side = 1, col = "red")
    } else {
        rug(theta_estimate[yeas], side = 3)
        rug(theta_estimate[nays], side = 1)
    }
}
