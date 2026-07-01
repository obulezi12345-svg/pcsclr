#' Simulate Datasets Under PCS-CLR Configurations
#' @param n Total initial sample population count.
#' @param target Scheduled tracking truncation parameter target threshold.
#' @param alpha Shape vector coefficient.
#' @param beta Scale distribution tracker.
#' @param p Attrition or competitive latent risk probability.
#' @param scheme Setting options: "baseline" or "hybrid".
#' @param cs_layout Distribution setups: "CS_I" or "CS_II".
#' @param T_max Explicit global temporal limitation window boundary.
#' @export
sim_pcsclr <- function(n, target, alpha = 1.8, beta = 2.5, p = 0.15, 
                       scheme = c("baseline", "hybrid"), cs_layout = c("CS_I", "CS_II"), 
                       T_max = 3.0) {
  scheme <- match.arg(scheme)
  cs_layout <- match.arg(cs_layout)
  
  delta <- rbinom(n, 1, 1 - p)
  true_lifetimes <- beta * (-log(runif(n)))^(1 / alpha)
  
  order_idx <- order(true_lifetimes)
  T_sorted  <- true_lifetimes[order_idx]
  D_sorted  <- delta[order_idx]
  
  R_vector <- rep(0, n)
  T_obs    <- rep(0, target)
  final_delta <- rep(0, target)
  
  rem_pool <- n - target
  for (i in 1:target) {
    T_obs[i] <- T_sorted[i]
    final_delta[i] <- D_sorted[i]
    
    if (final_delta[i] == 0) {
      R_vector[i] <- 0
    } else {
      if (cs_layout == "CS_I" && rem_pool > 0) {
        R_vector[i] <- sample(0:rem_pool, 1)
        rem_pool <- rem_pool - R_vector[i]
      } else if (cs_layout == "CS_II" && i == target) {
        R_vector[i] <- rem_pool
        rem_pool <- 0
      }
    }
  }
  
  if (scheme == "hybrid") {
    cutoff_idx <- which(T_obs > T_max)
    if (length(cutoff_idx) > 0) {
      valid_idx <- 1:(cutoff_idx[1] - 1)
      R_star <- rem_pool + sum(R_vector[cutoff_idx[1]:length(R_vector)])
      
      T_obs <- T_obs[valid_idx]
      final_delta <- final_delta[valid_idx]
      R_vector <- R_vector[valid_idx]
      
      out <- list(T_obs = T_obs, delta = final_delta, R_vector = R_vector, 
                  scheme = scheme, is_truncated = TRUE, T_max = T_max, R_star = R_star)
      class(out) <- "pcsclr_data"
      return(out)
    }
  }
  
  out <- list(T_obs = T_obs, delta = final_delta, R_vector = R_vector, scheme = scheme, is_truncated = FALSE)
  class(out) <- "pcsclr_data"
  return(out)
}