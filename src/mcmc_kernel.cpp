#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List rcpp_mcmc_kernel_internal(int M, int burn_in, arma::vec time, arma::vec delta,
                                     arma::vec removal, bool is_hybrid, double t_max, double r_star) {
  int n = time.n_elem;
  arma::vec alpha_chain(M);
  arma::vec beta_chain(M);

  // Initial stable chain configuration profiles
  alpha_chain(0) = 1.5;
  beta_chain(0) = 2.0;

  double curr_alpha = alpha_chain(0);
  double curr_beta = beta_chain(0);

  for(int m = 1; m < M; ++m) {
    // MH Step for alpha (Log-Normal proposal transition matrix space)
    double prop_alpha = exp(log(curr_alpha) + R::rnorm(0, 0.05));
    double log_acc_a = 0.0;

    for(int i = 0; i < n; ++i) {
      log_acc_a += delta(i) * (log(prop_alpha) - log(curr_alpha) + (prop_alpha - curr_alpha) * log(time(i) / curr_beta)) -
        (1.0 + removal(i)) * (pow(time(i) / curr_beta, prop_alpha) - pow(time(i) / curr_beta, curr_alpha));
    }
    // Structural boundary adjustment injection: Evaluates trailing survival constraints [S(T_max)]^R_star
    if(is_hybrid) {
      log_acc_a -= r_star * (pow(t_max / curr_beta, prop_alpha) - pow(t_max / curr_beta, curr_alpha));
    }

    if(log(R::runif(0, 1)) < log_acc_a) {
      curr_alpha = prop_alpha;
    }
    alpha_chain(m) = curr_alpha;

    // MH Step for beta (Log-Normal proposal transition matrix space)
    double prop_beta = exp(log(curr_beta) + R::rnorm(0, 0.05));
    double log_acc_b = 0.0;

    for(int i = 0; i < n; ++i) {
      log_acc_b += delta(i) * (curr_alpha * log(curr_beta) - curr_alpha * log(prop_beta)) -
        (1.0 + removal(i)) * (pow(time(i) / prop_beta, curr_alpha) - pow(time(i) / curr_beta, curr_alpha));
    }
    if(is_hybrid) {
      log_acc_b -= r_star * (pow(t_max / prop_beta, curr_alpha) - pow(t_max / curr_beta, curr_alpha));
    }

    if(log(R::runif(0, 1)) < log_acc_b) {
      curr_beta = prop_beta;
    }
    beta_chain(m) = curr_beta;
  }

  // Slicing parameters post burn-in window constraints
  arma::vec post_alpha = alpha_chain.subvec(burn_in, M - 1);
  arma::vec post_beta = beta_chain.subvec(burn_in, M - 1);

  return Rcpp::List::create(Rcpp::Named("alpha") = post_alpha, Rcpp::Named("beta") = post_beta);
}
