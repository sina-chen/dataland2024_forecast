#-------------------------------------------------------------------------------
#
# Fit forecasting model for provoince-level elections 1984-2023
#
# Author: Sina Chen
#
#
#-------------------------------------------------------------------------------


# Libraries ---------------------------------------------------------------

{
  library(stringr)
  library(rstan)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  library(dplyr)
  library(rv)
}

gc(verbose = T)


# Data --------------------------------------------------------------------

setwd("your working directory")

# 2024 poll scenario
sc <- "E"

# load data
load(paste0("forecast_data_", sc, ".RData"))

# check model input
sapply(stan_dat, length)
sapply(stan_dat, range)


# Fit stan model ----------------------------------------------------------

resStan <- stan(file = "forecast_model_sandbox.stan", 
                data = stan_dat,
                chains = 3, iter = 500,
                seed = 8423,
                control = list(adapt_delta = 0.95, max_treedepth = 12)
) 

gc(verbose=T)


# Extract predictions -----------------------------------------------------

# convert simulations to random variable
postrv_dataland <- as.rv(resStan)
rm(resStan)
gc(verbose=T)

# extract party support for every day
theta <- postrv_dataland$theta

# filter 2024 party supporet estimates
theta2024 <- theta[(dim(theta)[1]-(length(
  which(str_detect(election_data$election, "2024")) == T)-1)):dim(theta)[1],,]

# save simulation results
saveRDS(theta2024, paste0("theta_", sc, ".rds"))





