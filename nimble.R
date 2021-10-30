library(nimble)  # https://r-nimble.org/ -- great package for Bayesian inference
library(dplyr)

modcode <- nimbleCode({
  # Priors on all parameters
  n ~ dunif(1, 10)
  theta_r ~ dunif(0.0003, 0.3)
  theta_s ~ dunif(0.303, 0.8)
  alpha ~ dunif(0.00011, 0.1)

  # Prior distribution on the residual standard deviation. Probably wouldn't
  # mess with this one -- this is a pretty standard uninformative prior, and
  # the Gamma distribution ensures conjugacy, which makes the sampler more
  # efficient.
  # NOTE: Tau is precision, inverse of variance
  rtau ~ dgamma(0.1, 0.1)
  rsd <- 1/sqrt(rtau)

  # Likelihood -- Gaussian.
  for (i in 1:N) {
    model[i] <- theta_r + (theta_s - theta_r) / (1 + (alpha*x[i])^n)^(1 - 1/n)
    y[i] ~ dnorm(model[i], sd = rsd)
  }
})

fit_model <- function(sdata) {
  nimbledat <- list(N = nrow(sdata), x = sdata$h_cm, y = sdata$VWC)
  samples <- nimbleMCMC(code = modcode, constants = nimbledat,
                        nchains = 3, niter = 50000,
                        samplesAsCodaMCMC = TRUE,
                        monitors = c("theta_r", "theta_s", "alpha", "n", "rsd"))
  samples
}

lu_2008_2013_data <- read.csv("data/Lu_2008_2013_Data.csv")
soil_all_param1 <- fit_model(lu_2008_2013_data %>% filter(soil == "I"))
soil_all_param2 <- fit_model(lu_2008_2013_data %>% filter(soil == "II"))
soil_all_param3 <- fit_model(lu_2008_2013_data %>% filter(soil == "III"))
soil_all_param4 <- fit_model(lu_2008_2013_data %>% filter(soil == "IV"))
soil_all_param5 <- fit_model(lu_2008_2013_data %>% filter(soil == "V"))
soil_all_param6 <- fit_model(lu_2008_2013_data %>% filter(soil == "VI"))
soil_all_param7 <- fit_model(lu_2008_2013_data %>% filter(soil == "VII"))
soil_all_param8 <- fit_model(lu_2008_2013_data %>% filter(soil == "VIII"))

all_results <- list(
  I = soil_all_param1,
  II = soil_all_param2,
  III = soil_all_param3,
  IV = soil_all_param4,
  V = soil_all_param5,
  VI = soil_all_param6,
  VII = soil_all_param7,
  VIII = soil_all_param8
)

saveRDS(all_results, "nimble_results.rds")
