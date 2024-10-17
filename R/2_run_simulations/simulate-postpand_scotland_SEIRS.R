library(bbmle)
library(dplyr)
library(reshape2)

#### SETUP PARAMETERS ----------------------------------------------------------
# load simulation functions
source("R/2_run_simulations/SEIRS-model.R")

# load model fit
fit_scotland_SEIRS <- readRDS("data/derived_data/scotland/fit_scotland_SEIRS.rds")

# load scotland base parameters
source("R/scotland_parameters.R")

# load google mobility data
source("R/2_run_simulations/google_mobility_functions.R")

# google mobility data
npi = get_google_mobility_data("GB", FALSE)

scotland_by_wk_full <- scotland_by_wk %>%
  filter(wk_collected >= start_wk_pre)

# create data.frame with rsv pars
# use last week as starting value for simulation
posterior_rsv = as.data.frame(fit_scotland_SEIRS)[,c("S0_rsv", "E0_rsv", "I0_rsv", 
                                                     "rho_rsv", "b", "a", "p")]
colnames(posterior_rsv) = c("S0", "E0", "I0", "rho", "b", "a", "p")
posterior_rsv = cbind(posterior_rsv, r = 1, sigma = scotland_sigma, omega = scotland_omega,
                      gamma = scotland_gamma, mu = scotland_birthrate)

# mpv pars
posterior_mpv = as.data.frame(fit_scotland_SEIRS)[,c("S0_mpv", "E0_mpv","I0_mpv",
                                                     "rho_mpv", "b", "a", "p", 
                                                     "r", "c")]
colnames(posterior_mpv) = c("S0", "E0", "I0", "rho", "b", "a", "p", "r", "c")
posterior_mpv = cbind(posterior_mpv, sigma = scotland_sigma,  omega = scotland_omega,
                      gamma = scotland_gamma, mu = scotland_birthrate)

# randomly select a subset of simulations
n_samp <- 5000
samp <- sample(1:nrow(posterior_mpv), n_samp)
posterior_rsv <- posterior_rsv[samp, ]
posterior_mpv <- posterior_mpv[samp, ]


#### FIND BEST NPI SCALAR ------------------------------------------------------
# first find the npi scalar (in increment of 0.05) that gives best fit to
# out of sample data
fit_npi_scalar <- function(npi_scale, obs_rsv, obs_mpv, scotland_by_wk_full, posterior_rsv, posterior_mpv, 
                           scotland_N, end_wk_pre, start_wk_pre){
  # print(npi_scale)
  npi_scaled <- scotland_by_wk_full %>%
    filter(pathogen == "mpv") %>%
    left_join(npi, by = join_by(wk_collected)) %>%
    mutate(mean_mob = ifelse(is.na(mean_mob), 1, 1+(mean_mob*npi_scale/100))) %>%
    pull(mean_mob)
  pps <- SEIRS_force_discrete_both(
    N = length(npi_scaled), parms_rsv = posterior_rsv, parms_mpv = posterior_mpv,
    pop = scotland_N, npi = npi_scaled,
    testing_scalar_mpv = scotland_by_wk_full %>% filter(pathogen == "mpv") %>% pull(n_tests_ma_scaled),
    testing_scalar_rsv = scotland_by_wk_full %>% filter(pathogen == "rsv") %>% pull(n_tests_ma_scaled),
    rsvscaling_cuttoff = (end_wk_pre - start_wk_pre)/7 + 1
  )
  pps <- pps %>% 
    filter(variable == "Ipred") %>%
    left_join(t_to_wk_scotland_full, by = join_by(t))
  # get estimated detections from model
  pred_rsv <- pps %>% 
    filter(wk_collected > end_wk_pre, pathogen == "rsv") %>%
    left_join(obs_rsv, by = join_by(wk_collected)) %>%
    mutate(ll = dpois(x = detections, lambda = value, log = TRUE)) # calculate ll for each
  pred_mpv <- pps %>% 
    filter(wk_collected > end_wk_pre, pathogen == "mpv") %>%
    left_join(obs_mpv, by = join_by(wk_collected)) %>%
    mutate(ll = dpois(x = detections, lambda = value, log = TRUE)) # calculate ll for each
  # sum log-likelihood
  ll_rsv <- sum(pred_rsv$ll)
  ll_mpv <- sum(pred_mpv$ll)
  return(-(ll_rsv + ll_mpv))
}

npi_scalar2 <- unique(round(sort(c(seq(0.6, 1.3, 0.1), seq(0.9, 1.2, 0.05))),3))
npi_nll <- rep(NA, length.out = length(npi_scalar2))
for(i in 1:length(npi_scalar2)){
  print(npi_scalar2[i])
  npi_nll[i] = fit_npi_scalar(npi_scale = npi_scalar2[i], 
                              obs_rsv = scotland_by_wk %>% 
                                filter(wk_collected >= start_wk_post, pathogen == "rsv") %>% 
                                select(wk_collected, detections), 
                              obs_mpv = scotland_by_wk %>% 
                                filter(wk_collected >= start_wk_post, pathogen == "mpv") %>% 
                                select(wk_collected, detections), 
                              scotland_by_wk_full = scotland_by_wk_full, 
                              posterior_rsv = posterior_rsv, 
                              posterior_mpv = posterior_mpv, 
                              scotland_N = scotland_N, 
                              end_wk_pre = end_wk_pre, 
                              start_wk_pre = start_wk_pre)
}

# plot(npi_scalar2, npi_nll)
# lines(npi_scalar2, npi_nll)

npi_scalar <- npi_scalar2[which(npi_nll == min(npi_nll))]

#### SIMULATE POST-PANDEMIC REBOUND --------------------------------------------

# simulate forward with npi_scalar
npi_scaled <- scotland_by_wk_full %>%
  filter(pathogen == "mpv") %>%
  left_join(npi) %>%
  mutate(mean_mob = ifelse(is.na(mean_mob), 1, 1+(mean_mob*npi_scalar/100))) %>%
  pull(mean_mob)
postpand_sims <- SEIRS_force_discrete_both(
  N = length(npi_scaled), parms_rsv = posterior_rsv, parms_mpv = posterior_mpv,
  pop = scotland_N, npi = npi_scaled,
  testing_scalar_mpv = scotland_by_wk_full %>% filter(pathogen == "mpv") %>% pull(n_tests_ma_scaled),
  testing_scalar_rsv = scotland_by_wk_full %>% filter(pathogen == "rsv") %>% pull(n_tests_ma_scaled),
  rsvscaling_cuttoff = (end_wk_pre - start_wk_pre)/7 + 1
)

postpand_sims <- postpand_sims %>%
  mutate(draw_id = samp[draw_id])

# saveRDS(postpand_sims, "data/derived_data/scotland/postpand-sims_scotland_SEIRS.rds")