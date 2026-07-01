.weibull_haz <- function(t, alpha, beta) { (alpha / beta) * (t / beta)^(alpha - 1) }

.weibull_cum_haz <- function(t, alpha, beta) { (t / beta)^alpha }

.compute_weibull_score <- function(time, delta, removal, theta) {
  alpha <- theta[1]; beta <- theta[2]
  d_alpha <- sum(delta / alpha + delta * log(time / beta) -
                   (1 + removal) * (time / beta)^alpha * log(time / beta))
  d_beta <- sum(-delta * alpha / beta + (1 + removal) * (alpha / beta) * (time / beta)^alpha)
  return(c(d_alpha, d_beta))
}

.compute_weibull_info <- function(time, delta, removal, theta) {
  alpha <- theta[1]; beta <- theta[2]
  d_alpha2 <- sum(-delta / alpha^2 - (1 + removal) * (time / beta)^alpha * (log(time / beta))^2)
  d_beta2  <- sum(delta * alpha / beta^2 - (1 + removal) * (alpha * (alpha + 1) / beta^2) * (time / beta)^alpha)
  d_alphabeta <- sum(-delta / beta + (1 + removal) * (1 / beta) * (time / beta)^alpha * (1 + alpha * log(time / beta)))

  mat <- matrix(c(-d_alpha2, -d_alphabeta, -d_alphabeta, -d_beta2), nrow = 2, ncol = 2)
  if (any(is.nan(mat)) || any(is.infinite(mat))) mat <- diag(2) * 1e-4
  return(mat)
}

#' Fit PCS-CLR Parameter Estimation Models
#'
#' @param time Vector of observed lifetime intervals.
#' @param delta Attrition risk indicator matrix.
#' @param removal Active dynamic progressive censorship counts.
#' @param method Selection between "RK4" paths or "Bayes_Kernel" sampler loops.
#' @param control Optional configuration list parameters.
#' @param is_hybrid Boolean indicator for hybrid censoring regimes.
#' @param t_max Truncation time milestone boundary constraint.
#' @param r_star Escaped survival components unfailed beyond cutoff point.
#' @param ... Additional structural arguments passed to internal optimization engines.
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib pcsclr, .registration = TRUE
#' @export
fit_pcsclr <- function(time, delta, removal, method = c("RK4", "Bayes_Kernel"),
                       control = list(), is_hybrid = FALSE, t_max = 3.0, r_star = 0, ...) {
  method <- match.arg(method)

  if (method == "RK4") {
    h <- if (is.null(control$h)) 0.005 else control$h
    max_iter <- if (is.null(control$max_iter)) 1000 else control$max_iter
    tol <- if (is.null(control$tol)) 1e-6 else control$tol

    theta <- c(1.5, 2.0)

    for (i in 1:max_iter) {
      score_v <- .compute_weibull_score(time, delta, removal, theta)
      inf_mat <- .compute_weibull_info(time, delta, removal, theta)
      inv_inf <- tryCatch(solve(inf_mat), error = function(e) diag(2) * 0.01)

      k1 <- h * (inv_inf %*% score_v)
      k2 <- h * (tryCatch(solve(.compute_weibull_info(time, delta, removal, theta + k1/2)), error = function(e) diag(2)*0.01) %*% .compute_weibull_score(time, delta, removal, theta + k1/2))
      k3 <- h * (tryCatch(solve(.compute_weibull_info(time, delta, removal, theta + k2/2)), error = function(e) diag(2)*0.01) %*% .compute_weibull_score(time, delta, removal, theta + k2/2))
      k4 <- h * (tryCatch(solve(.compute_weibull_info(time, delta, removal, theta + k3)), error = function(e) diag(2)*0.01) %*% .compute_weibull_score(time, delta, removal, theta + k3))

      theta_new <- theta + (k1 + 2*k2 + 2*k3 + k4) / 6
      theta_new <- pmax(theta_new, 1e-4)

      if (max(abs(theta_new - theta)) < tol) break
      theta <- theta_new
    }

    out <- list(estimates = as.vector(theta), method = "RK4", iterations = i)
  } else {
    M <- if (is.null(control$M)) 10000 else control$M
    burn_in <- if (is.null(control$burn_in)) 1000 else control$burn_in

    chains <- rcpp_mcmc_kernel_internal(M, burn_in, time, delta, removal, is_hybrid, t_max, r_star)
    out <- list(estimates = c(mean(chains$alpha), mean(chains$beta)),
                chains = chains, method = "Bayes_Kernel")
  }

  class(out) <- "pcsclr_fit"
  return(out)
}
