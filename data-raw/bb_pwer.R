## code to prepare `bb_pwer` dataset goes here

# Simulation framework modified from:
# https://stats.stackexchange.com/questions/1866/how-to-simulate-a-custom-power-analysis-of-an-lm-model-using-r
a = 0  #estimated intercept
b = 0.5  #desired slope
c = 0  #covariate slope
sd = 1  #residual sd
mu1 <- mu2 <- 0 # predictor means
s1 <- 1 # sd X1
s2 <- 3 # sd X2
correl <- 0.1 # correlation between X1 and X2
Sigma <- matrix(c(s1^2, s1*s2*correl, s1*s2*correl, s2^2), 2, 2) # var-cov
nsim = 1000  #1000 simulations
pval = numeric(nsim)  #placeholder for the second for loop output
Nvec = seq(25, 100, by = 1)  #vector for the range of sample sizes to be tested
power.N = numeric(length(Nvec))   #create placeholder for first for loop output
for (j in 1:length(Nvec)) {
  N = Nvec[j]
  bndat <- MASS::mvrnorm(N, mu = c(mu1, mu2), Sigma = Sigma) # generate bivariate normal data
  X1 <- bndat[,1] # extract X1
  X2 <- bndat[,2] # extract X2
  for (i in 1:nsim) {   #for this value of N, create random error nsim times
    y_det = a + b * X1 + c*X2
    y = stats::rnorm(N, mean = y_det, sd = sd)
    m = stats:::lm(y ~ X1 + X2)
    pval[i] = coef(summary(m))["X1", "Pr(>|t|)"]  #all the p values
  }  #cycle through all N values
  power.N[j] = sum(pval < 0.05)/nsim  #the proportion of correct p-values (i.e the power)
}
power.N
plot(Nvec, power.N)
bb_pwer <- data.frame(N = Nvec,
                      Pwr = power.N)
usethis::use_data(bb_pwer, overwrite = TRUE)
