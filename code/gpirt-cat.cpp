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

// [[Rcpp::export]]
Rcpp::NumericVector gpirt_cat(const int n_questions,
                              const arma::mat& f_bars,
                              const arma::vec& theta_grid,
                              const arma::vec& theta_prior,
                              const arma::vec& responses) {
    arma::vec theta_post = normalize(theta_grid, theta_prior);
    arma::vec log_theta_post = arma::log(theta_post);
    arma::mat log_f = arma::log(f_bars);
    arma::uword m = f_bars.n_cols;
    Rcpp::IntegerVector items(n_questions);
    Rcpp::IntegerVector choices(n_questions);
    Rcpp::NumericVector estimates(n_questions);
    arma::uvec item_set = arma::regspace<arma::uvec>(0, m-1);
    for ( arma::uword t = 0; t < n_questions; ++t ) {
        arma::vec MI = mutual_information_mat(theta_grid, theta_post,
                                              f_bars.cols(item_set));
        arma::uword idx = MI.index_max();
        items[t] = item_set[idx];
        item_set.shed_row(idx);
        int choice = responses[items[t]];
        choices[t] = choice;
        if ( choice == 1 ) {
            log_theta_post += arma::log(f_bars.col(items[t]));
        } else if ( choice != -2147483648 ) { // Deal with NAs
            log_theta_post += arma::log(1.0 - f_bars.col(items[t]));
        }
        theta_post = normalize(theta_grid, arma::exp(log_theta_post));
        estimates[t] = reimann_sum(theta_grid, theta_grid % theta_post);
    }
    Rcpp::NumericVector res = Rcpp::wrap(theta_post);
    res.attr("items") = items;
    res.attr("choices") = choices;
    res.attr("estimates") = estimates;
    return res;
}
