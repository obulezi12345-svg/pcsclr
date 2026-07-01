#' @export
print.pcsclr_fit <- function(x, ...) {
  cat("========================================================\n")
  cat(" PCS-CLR Estimation Framework Output Object             \n")
  cat("========================================================\n")
  cat("Optimization Method Engine : ", x$method, "\n")
  cat("Estimated Shape Parameter (\u03b1) : ", round(x$estimates[1], 4), "\n")
  cat("Estimated Scale Parameter (\u03b2) : ", round(x$estimates[2], 4), "\n")
  cat("--------------------------------------------------------\n")
}

#' @export
summary.pcsclr_fit <- function(object, ...) {
  print(object)
  if (object$method == "Bayes_Kernel") {
    cat("MCMC Chain Summary statistics:\n")
    cat("  \u03b1 SD: ", round(sd(object$chains$alpha), 4), " | 95% HPD: [", 
        round(quantile(object$chains$alpha, 0.025), 4), ",", round(quantile(object$chains$alpha, 0.975), 4), "]\n")
    cat("  \u03b2 SD: ", round(sd(object$chains$beta), 4), " | 95% HPD: [", 
        round(quantile(object$chains$beta, 0.025), 4), ",", round(quantile(object$chains$beta, 0.975), 4), "]\n")
  }
}

#' @importFrom graphics par plot abline
#' @importFrom stats density quantile sd rbinom runif
#' @export
plot.pcsclr_fit <- function(x, ...) {
  if (x$method != "Bayes_Kernel") {
    stop("Multi-trace graphic generation is reserved exclusively for Bayesian sampler objects.")
  }
  oldpar <- par(mfrow = c(2, 2))
  on.exit(par(oldpar))
  
  plot(x$chains$alpha, type = "l", col = "steelblue", main = "Trace plot: Shape (\u03b1)", xlab = "Iteration", ylab = "")
  plot(x$chains$beta, type = "l", col = "darkorange", main = "Trace plot: Scale (\u03b2)", xlab = "Iteration", ylab = "")
  
  plot(density(x$chains$alpha), col = "steelblue", lwd = 2, main = "Posterior Profile: \u03b1")
  abline(v = x$estimates[1], col = "red", lty = 2)
  plot(density(x$chains$beta), col = "darkorange", lwd = 2, main = "Posterior Profile: \u03b2")
  abline(v = x$estimates[2], col = "red", lty = 2)
}