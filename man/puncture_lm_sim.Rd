% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/puncture_lm_sim.R
\name{puncture_lm_sim}
\alias{puncture_lm_sim}
\title{Simulation Framework for Comparing Standard and Puncture Methods}
\usage{
puncture_lm_sim(
  n = 100,
  sim_iter = 1000,
  beta_gen,
  gen_ivs,
  .formula = "Y ~ X1 + X2",
  func = stats::lm,
  combine = mean,
  alpha = 0.05,
  ...
)
}
\arguments{
\item{n}{Positive integer specifying the sample size for each simulation iteration.
Default is 100.}

\item{sim_iter}{Positive integer specifying the number of Monte Carlo iterations.
Default is 1000.}

\item{beta_gen}{Function that returns a list of true population parameters:
first element is the intercept (beta0), followed by slope coefficients
(beta1, beta2, etc.).}

\item{gen_ivs}{Function that takes sample size n as input and returns a list
containing the predictor variables and error term (epsilon).}

\item{.formula}{Character string specifying the model formula. Default is
"Y ~ X1 + X2".}

\item{func}{Function used to fit the model (e.g., stats::lm). Default is
stats::lm.}

\item{combine}{Function used to combine results from multiple bootstrap
iterations in the puncture method. Default is mean.}

\item{alpha}{Numeric value between 0 and 1 specifying the significance level
for hypothesis tests and confidence intervals. Default is 0.05.}

\item{...}{Additional arguments passed to the puncture function.}
}
\value{
A list containing two elements:
\describe{
\item{results}{A list with detailed results for both methods, including:
\itemize{
\item beta1_estimates: Point estimates
\item beta1_SE: Standard errors
\item beta1_bias: Bias from true parameter
\item beta1_CI_low: Lower confidence interval bounds
\item beta1_CI_high: Upper confidence interval bounds
\item beta1_coverage: Indicator for CI coverage of true parameter
\item beta1_rejected: Indicator for null hypothesis rejection
\item converged: Convergence indicators
}
}
\item{summary}{A data frame comparing both methods on:
\itemize{
\item Average bias
\item Average standard error
\item Coverage probability
\item Average confidence interval width
\item Rejection rate
\item Convergence rate
}
}
}
}
\description{
Conducts Monte Carlo simulations to compare the performance of standard linear
regression against the puncture method. Generates data under specified conditions,
applies both methods, and computes various performance metrics including bias,
coverage, and rejection rates.
}
\details{
The function implements the following simulation procedure:
\enumerate{
\item Generates true parameters using beta_gen()
\item Generates predictor variables and error term using gen_ivs()
\item Constructs response variable using linear model
\item Fits standard linear model
\item Applies puncture method
\item Computes and stores performance metrics
\item Repeats for specified number of iterations
\item Summarizes results across all iterations
}
}
\examples{
sim_results <- puncture_lm_sim(n = 200, sim_iter = 5,
beta_gen = function() {
  list(1,             # beta0 (intercept)
       1,             # beta1
       3)             # beta2
},
gen_ivs = function(n) {
  X1 <- stats::rnorm(n, 0, 20)
  X2 <- 0.6 * X1 + sqrt(1 - 0.6^2) * stats::rnorm(n)
  epsilon <- stats::rnorm(n, 0, abs(X1))
  list(X1 = X1, X2 = X2, epsilon = epsilon)
},
m = 3, b = 3, remove.collinear=FALSE) |> suppressWarnings()

}
