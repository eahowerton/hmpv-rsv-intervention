library(rstan)
library(readr)
library(dplyr)
library(bayesplot)

#### PREP DATA
source("R/scotland_parameters.R")

# read in scotland data, all age groups, summarized by week
scotland_by_wk <- read_rds("data/derived_data/scotland/scotland_by_wk_overall.rds")

scotland_by_wk_pre <- scotland_by_wk %>%
  filter(wk_collected >= start_wk_pre, 
         wk_collected <= end_wk_pre) 

#### FIT SEIRS MODEL -----------------------------------------------------------
# assumptions: 
# 1. same waning rate for RSV and HMPV
# 2. same seasonal forcing (a and p)
# 3. HMPV transmission rate differs by some factor r, and is also forced by
#    observed RSV incidence by some factor c
model_SIRS <- stan_model("R/1_fit_models/SIRS_coupled.stan")

# constrain sampling of initial conditions to ensure S0 + E0 + I0 < 1
initfun = function(...){
  list(
    S0_mpv = rbeta(1, 20, 80),
    I0_mpv = rbeta(1, 20, 80),
    S0_rsv = rbeta(1, 20, 80),
    I0_rsv = rbeta(1, 20, 80)
  )
}

# run model
fit_scotland_SIRS <- sampling(
  model_SIRS,
  data = list(
    N = nrow(scotland_by_wk_pre %>% filter(pathogen == "mpv")),
    cases_mpv = scotland_by_wk_pre %>% filter(pathogen == "mpv") %>% pull(detections), 
    cases_rsv = scotland_by_wk_pre %>% filter(pathogen == "rsv") %>% pull(detections), 
    birth = scotland_by_wk_pre %>% filter(pathogen == "mpv") %>% pull(Births.registered), 
    tests_mpv = scotland_by_wk_pre %>% filter(pathogen == "mpv") %>% pull(n_tests_ma_scaled),
    tests_rsv = scotland_by_wk_pre %>% filter(pathogen == "rsv") %>% pull(n_tests_ma_scaled),
    pop = scotland_N,
    gamma = 1/(1/scotland_sigma + 1/scotland_gamma),
    mu = scotland_birthrate,
    omega = scotland_omega
  ),
  seed = 7, 
  iter = 15000,
  chain = 4, 
  init = initfun, 
  cores = 4,
  control = list(max_treedepth = 12)
)

mcmc_pairs(fit_scotland_SIRS, pars = c("b", "r", "c", "a", "p")) 

# save output
saveRDS(fit_scotland_SIRS, "data/derived_data/scotland/fit_scotland_SIRS.rds")

