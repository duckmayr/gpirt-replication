## update_post() updates the posterior over theta given a (log) prior, IRFs,
## responses, and the grid over which the prior was computed. Note the columns
## of the IRFs must comport w/ response elements & the rows of the IRFs must
## comport with both the elements of log_prior and of theta_grid
update_posterior <- function(log_prior, IRFs, responses, theta_grid) {
    m <- length(responses)
    for ( j in 1:m ) {
        if ( responses[j] == 1 ) {
            log_prior <- log_prior + log(IRFs[ , j])
        } else {
            log_prior <- log_prior + log(1 - IRFs[ , j])
        }
    }
    return(normalize(theta_grid, exp(log_prior)))
}
## get_estimate() takes a value grid and posterior density and returns the mean
get_estimate <- function(post, theta_grid) {
    return(reimann_sum(theta_grid, theta_grid * post))
}
## rmse() calculates the root mean squared error
rmse <- function(actual, expected) sqrt(mean((actual - expected)^2))
