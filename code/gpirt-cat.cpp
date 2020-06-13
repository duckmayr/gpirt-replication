#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double reimann_sum(const arma::vec& x, const arma::vec& y) {
    return arma::sum(arma::diff(x) % y.head(y.n_elem - 1));
}

// [[Rcpp::export]]
arma::vec normalize(const arma::vec& x, const arma::vec& p) {
    double area = reimann_sum(x, p);
    return p / area;
}

// [[Rcpp::export]]
arma::vec binary_entropy(const arma::vec& p) {
    return -(p % arma::log2(p) ) - ((1 - p) % arma::log2(1 - p));
}

double binary_entropy(double p) {
    return -(p * std::log2(p)) - ((1 - p) * std::log2(1 - p));
}

// [[Rcpp::export]]
double mutual_information(const arma::vec& theta_grid,
                          const arma::vec& theta_prior,
                          const arma::vec& f_bar) {
    arma::vec post = theta_prior % f_bar;
    double H1 = binary_entropy(reimann_sum(theta_grid, post));
    double H2 = reimann_sum(theta_grid, binary_entropy(f_bar));
    return H1 - H2;
}

arma::vec mutual_information_mat(const arma::vec& theta_grid,
                                 const arma::vec& theta_prior,
                                 const arma::mat& f_bar) {
    arma::uword m = f_bar.n_cols;
    arma::vec res(m);
    for ( arma::uword j = 0; j < m; ++j ) {
        res[j] = mutual_information(theta_grid, theta_prior, f_bar.col(j));
    }
    return res;
}

//' GPIRT Adaptive Testing Simulation with Pre-Recorded Responses
//' 
//' This function uses pre-recorded responses to show which questions would
//' have been asked, and what the latent trait estimates would have been, after
//' presenting a certain number of items using our adaptive testing scheme.
//' 
//' @param n_questions An integer vector of length one giving the number of
//'     items to present the respondent
//' @param f_bars A numeric matrix giving the estimated IRFs (denoted pi_i
//'     in the paper)
//' @param theta_grid A numeric vector giving the grid of theta values to
//'     consider (denoted theta^* in the paper)
//' @param theta_prior A numeric vector giving the prior belief at each value
//'     in theta_grid (can be unnormalized)
//' @param responses A numeric or integer vector of responses to the items
//'     corresponding to the columns in f_bars; elements should be one of 1,
//'     -1, or NA
//' 
//' @return A numeric vector giving the posterior over theta at theta_grid,
//'     with attributes "items", giving the items our adaptive scheme
//'     presented the respondent, "choices", the respondent's responses to
//'     those items, and "estimates", the estimated latent trait after each
//'     response received.
// [[Rcpp::export]]
Rcpp::NumericVector gpirt_cat(const int n_questions,
                              const arma::mat& f_bars,
                              const arma::vec& theta_grid,
                              const arma::vec& theta_prior,
                              const arma::vec& responses) {
    // Get a normalized prior
    arma::vec theta_post = normalize(theta_grid, theta_prior);
    // We use log probabilities for increased precision
    arma::vec log_theta_post = arma::log(theta_post);
    arma::mat log_f = arma::log(f_bars);
    // Set some bookkeeping variables:
    arma::uword m = f_bars.n_cols;              // number of items
    Rcpp::IntegerVector items(n_questions);     // items we asked
    Rcpp::IntegerVector choices(n_questions);   // responses given
    Rcpp::NumericVector estimates(n_questions); // theta estimate after each Q
    // This keeps track of the set of items we haven't asked yet:
    arma::uvec item_set = arma::regspace<arma::uvec>(0, m-1);
    // We will present 'n_questions' number of items; for each one,
    for ( arma::uword t = 0; t < n_questions; ++t ) {
        // Get the mutual information for each unasked question
        arma::vec MI = mutual_information_mat(theta_grid, theta_post,
                                              f_bars.cols(item_set));
        // Find the one with the highest information
        arma::uword idx = MI.index_max();
        // Record that this is the question we'll ask
        items[t] = item_set[idx];
        // Remove this item from the unasked item set
        item_set.shed_row(idx);
        // Get the response to this item
        int choice = responses[items[t]];
        // Record the choice
        choices[t] = choice;
        // Update our posterior
        if ( choice == 1 ) {
            log_theta_post += arma::log(f_bars.col(items[t]));
        } else if ( choice != -2147483648 ) { // Deal with NAs
            log_theta_post += arma::log(1.0 - f_bars.col(items[t]));
        }
        // Normalize it
        theta_post = normalize(theta_grid, arma::exp(log_theta_post));
        // And get our updated estimate
        estimates[t] = reimann_sum(theta_grid, theta_grid % theta_post);
    }
    // Return the normalized posterior
    // (we convert to Rcpp::NumericVector to record some of the bookkeeping
    //  variables in attributes)
    Rcpp::NumericVector res = Rcpp::wrap(theta_post);
    res.attr("items") = items;
    res.attr("choices") = choices;
    res.attr("estimates") = estimates;
    return res;
}
