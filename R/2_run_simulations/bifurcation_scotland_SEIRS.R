library(dplyr)
library(reshape2)

#### SETUP PARAMETERS ----------------------------------------------------------
# load simulation functions
source("R/2_run_simulations/SEIRS-model.R")

# load model fit
fit_scotland_SEIRS <- readRDS("data/derived_data/scotland/fit_scotland_SEIRS.rds")

# load scotland base parameters
source("R/scotland_parameters.R")

# extract median parameters
fitted_pars = apply(as.array(fit_scotland_SEIRS)[,,1:14], 3, median)

# create vector with rsv pars
fitted_pars_rsv = fitted_pars[c("S0_rsv", "E0_rsv", "I0_rsv", "rho_rsv", "b", "a", "p")]
names(fitted_pars_rsv) = c("S0", "E0", "I0", "rho", "b", "a", "p")
fitted_pars_rsv = c(fitted_pars_rsv, r = 1, sigma = scotland_sigma, 
                    gamma = scotland_gamma, mu = scotland_birthrate, 
                    omega = scotland_omega)

# mpv pars
fitted_pars_mpv = fitted_pars[c("S0_mpv", "E0_mpv", "I0_mpv", "rho_mpv", "b", "a", "p", "r", "c")]
names(fitted_pars_mpv) = c("S0", "E0", "I0", "rho", "b", "a", "p", "r", "c")
fitted_pars_mpv = c(fitted_pars_mpv, sigma = scotland_sigma, 
                    gamma = scotland_gamma, mu = scotland_birthrate, 
                    omega = scotland_omega)

n_wks = 52*40

c_slices = sort(c(0, fitted_pars["c"]*0.9, fitted_pars["c"], fitted_pars["c"]*1.1))
a_range = sort(c(seq(0.01, 0.99, 0.001), fitted_pars["a"]))

bifur_sims = vector("list", length(c_slices))
bifur_points = vector("list", length(c_slices))
for(j in 1:length(c_slices)){
  print(c_slices[j])
  bifur_sims_ctmp = vector("list", length(a_range))
  bifur_points_ctmp = vector("list", length(a_range))
  for(i in 1:length(a_range)){
    tst_pars_rsv = fitted_pars_rsv
    tst_pars_rsv["a"] = a_range[i]
    tst_pars_mpv = fitted_pars_mpv
    tst_pars_mpv["a"] = a_range[i]
    tst_pars_mpv["c"] = c_slices[j]
    tst_pars_rsv <- matrix(tst_pars_rsv, nrow = 1)
    colnames(tst_pars_rsv) <- names(fitted_pars_rsv)
    tst_pars_mpv <- matrix(tst_pars_mpv, nrow = 1)
    colnames(tst_pars_mpv) <- names(fitted_pars_mpv)
    bifur_sims_ctmp[[i]] = SEIRS_force_discrete_both(
      N = n_wks, parms_rsv = tst_pars_rsv, parms_mpv = tst_pars_mpv, pop = scotland_N, 
      npi = 1, testing_scalar_mpv = rep(1, n_wks),
      testing_scalar_rsv = rep(1, n_wks),
      rsvscaling_cuttoff = (end_wk_pre - start_wk_pre)/7 + 1)
    bifur_points_ctmp[[i]] = bifur_sims_ctmp[[i]] %>%
      filter(t > 30*52, t%%52 == 12) %>%
      mutate(a = a_range[i], 
             c = c_slices[j])
  }
  bifur_sims[[j]] = bifur_sims_ctmp
  bifur_points[[j]] = bifur_points_ctmp
}

write_rds(bifur_sims, "data/derived_data/scotland/bifur-sims_scotland_SEIRS.rds")
write_rds(bifur_points, "data/derived_data/scotland/bifur-points_scotland_SEIRS.rds")

