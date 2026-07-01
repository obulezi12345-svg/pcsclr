# pcsclr: Progressive Censoring Schemes with Competitive Latent-Risk

The `pcsclr` package provides high-performance computing utilities for parameter estimation, simulation, and survival analysis under Progressive Censoring Schemes (PCS) in the presence of Competitive Latent Risks (CLR). 

The package is designed to implement the methodologies developed in accompanying research manuscripts, providing both fast numerical maximum likelihood estimation (MLE) and robust Bayesian Markov Chain Monte Carlo (MCMC) samplers.

# Features

* Numerical Optimization: Implements a 4th-order Runge-Kutta (RK4) path optimization framework to maximize log-likelihood structures and evaluate the observed Fisher Information matrix.
* High-Speed Bayesian Samplers: Features an accelerated C++ MCMC engine built via `RcppArmadillo` to execute random-walk Metropolis-Hastings update loops for Weibull baseline parameters under standard and hybrid progressive censoring setups.
* Data Simulation: Provides robust structures to simulate complex lifetime trajectories matching custom progressive removal patterns and multi-cause latent risk dropouts.



# Installation

You can install the development version of `pcsclr` directly from GitHub using the `remotes` package. Run the following commands in your R console:

```R
# Install remotes if you haven't already
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install pcsclr from the cloud repository
remotes::install_github("https://github.com/obulezi12345-svg/pcsclr")
