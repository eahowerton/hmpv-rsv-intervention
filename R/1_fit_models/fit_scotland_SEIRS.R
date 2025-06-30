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

#### FIT COUPLED SEIRS MODEL ---------------------------------------------------
# assumptions: 
# 1. same waning rate for RSV and HMPV
# 2. same seasonal forcing (a and p)
# 3. HMPV transmission rate differs by some factor r, and is also forced by
#    observed RSV incidence by some factor c
model_SEIRS <- stan_model("R/1_fit_models/SEIRS_coupled.stan")

# constrain sampling of initial conditions to ensure S0 + E0 + I0 < 1
initfun = function(...){
  list(
    S0_mpv = rbeta(1, 20, 80),
    E0_mpv = rbeta(1, 20, 80),
    I0_mpv = rbeta(1, 20, 80),
    S0_rsv = rbeta(1, 20, 80),
    E0_rsv = rbeta(1, 20, 80),
    I0_rsv = rbeta(1, 20, 80)
  )
}

# don't save every variable
pars_to_exclude = c("Sout_rsv", "Eout_rsv", "Iout_rsv", "Rout_rsv", 
                    "Sdeath_rsv", "Edeath_rsv", "Ideath_rsv", "Rdeath_rsv", 
                    "S_rsv", "E_rsv", "I_rsv", "R_rsv", "beta_rsv", "foi_rsv", #"Ipred_rsv",
                    "Sout_mpv", "Eout_mpv", "Iout_mpv", "Rout_mpv", 
                    "Sdeath_mpv", "Edeath_mpv", "Ideath_mpv", "Rdeath_mpv", 
                    "S_mpv", "E_mpv", "I_mpv", "R_mpv", "beta_mpv", "foi_mpv", #"Ipred_mpv", 
                    "birth", "max_RSV")

# run model
fit_scotland_SEIRS <- sampling(
  model_SEIRS,
  data = list(
    N = nrow(scotland_by_wk_pre %>% filter(pathogen == "mpv")),
    cases_mpv = scotland_by_wk_pre %>% filter(pathogen == "mpv") %>% pull(detections), 
    cases_rsv = scotland_by_wk_pre %>% filter(pathogen == "rsv") %>% pull(detections), 
    birth = scotland_by_wk_pre %>% filter(pathogen == "mpv") %>% pull(Births.registered), 
    tests_mpv = scotland_by_wk_pre %>% filter(pathogen == "mpv") %>% pull(n_tests_ma_scaled),
    tests_rsv = scotland_by_wk_pre %>% filter(pathogen == "rsv") %>% pull(n_tests_ma_scaled),
    pop = scotland_N,
    sigma = scotland_sigma,
    gamma = scotland_gamma,
    mu = scotland_birthrate,
    omega = scotland_omega
  ),
  seed = 7, 
  iter = 15000,
  chain = 4, 
  init = initfun, 
  cores = 4,
  control=list(max_treedepth = 12),
  include = FALSE,        # actually, this means 'exclude' the parameters listed below
  pars = pars_to_exclude
)

mcmc_pairs(fit_scotland_SEIRS, pars = c("b", "r", "c", "a", "p")) 

# save output
saveRDS(fit_scotland_SEIRS, "data/derived_data/scotland/fit_scotland_SEIRS.rds")

#### FIT INDEPENDENT SEIRS MODEL -----------------------------------------------
model_SEIRS_ind <- stan_model("R/1_fit_models/SEIRS_independent.stan")

# constrain sampling of initial conditions to ensure S0 + E0 + I0 < 1
initfun_ind = function(...){
  list(
    S0 = rbeta(1, 20, 80),
    E0 = rbeta(1, 20, 80),
    I0 = rbeta(1, 20, 80)
  )
}

pars_to_exclude_ind = c("Sout", "Eout", "Iout", "Rout", 
  "Sdeath", "Edeath", "Ideath", "Rdeath", 
  "S", "E", "I", "R", "Ipred", "beta", "foi",
  "birth", "max")

# run model
fit_scotland_SEIRS_mpv <- sampling(
  model_SEIRS_ind,
  data = list(
    N = nrow(scotland_by_wk_pre %>% filter(pathogen == "mpv")),
    cases = scotland_by_wk_pre %>% filter(pathogen == "mpv") %>% pull(detections), 
    birth = scotland_by_wk_pre %>% filter(pathogen == "mpv") %>% pull(Births.registered), 
    tests = scotland_by_wk_pre %>% filter(pathogen == "mpv") %>% pull(n_tests_ma_scaled),
    pop = scotland_N,
    sigma = scotland_sigma,
    gamma = scotland_gamma,
    mu = scotland_birthrate,
    omega = scotland_omega
  ),
  seed = 7, 
  iter = 15000,
  chain = 4, 
  init = initfun_ind, 
  cores = 4,
  control=list(max_treedepth = 12, 
               adapt_delta = 0.7), 
  include = FALSE,        # actually, this means 'exclude' the parameters listed below
  pars = pars_to_exclude_ind
)

mcmc_pairs(fit_scotland_SEIRS_mpv, pars = c("b", "a", "p", "S0", "E0", "I0", "rho")) 

saveRDS(fit_scotland_SEIRS_mpv, "data/derived_data/scotland/fit_scotland_SEIRS_mpv.rds")

# run model
fit_scotland_SEIRS_rsv <- sampling(
  model_SEIRS_ind,
  data = list(
    N = nrow(scotland_by_wk_pre %>% filter(pathogen == "rsv")),
    cases = scotland_by_wk_pre %>% filter(pathogen == "rsv") %>% pull(detections), 
    birth = scotland_by_wk_pre %>% filter(pathogen == "rsv") %>% pull(Births.registered), 
    tests = scotland_by_wk_pre %>% filter(pathogen == "rsv") %>% pull(n_tests_ma_scaled),
    pop = scotland_N,
    sigma = scotland_sigma,
    gamma = scotland_gamma,
    mu = scotland_birthrate,
    omega = scotland_omega
  ),
  seed = 7, 
  iter = 15000,
  chain = 4, 
  init = initfun_ind, 
  cores = 4,
  control=list(max_treedepth = 12), 
  include = FALSE,        # actually, this means 'exclude' the parameters listed below
  pars = pars_to_exclude_ind
)

mcmc_pairs(fit_scotland_SEIRS_rsv, pars = c("b", "a", "p", "S0", "E0", "I0", "rho")) 

saveRDS(fit_scotland_SEIRS_rsv, "data/derived_data/scotland/fit_scotland_SEIRS_rsv.rds")

