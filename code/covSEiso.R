#' The squared expoential covariance, or radial basis, function
#' 
#' @param x1 A numeric vector, the first variable
#' @param x2 A numeric vector, the second variable
#' @param sf A numeric vector of length 1, the scale factor
#' @param ell A numeric vector of length 1, the length scale
#' 
#' @return A numeric matrix whose i, j entry gives the covariance between
#'     f(x1[i]) and f(x2[j])
covSEiso <- function(x1, x2 = x1, sf = 1.0, ell = 1.0) {
    sf  <- sf^2
    z   <- -0.5 * (1 / (ell^2))
    n   <- length(x1)
    m   <- length(x2)
    res <- matrix(0, nrow = n, ncol = m)
    for ( j in 1:m ) {
        for ( i in 1:n ) {
            res[i, j] <- sf * exp(z * (x1[i] - x2[j])^2)
        }
    }
    return(res)
}
