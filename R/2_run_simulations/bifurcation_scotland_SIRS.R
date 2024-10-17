#### SETUP PARAMETERS ----------------------------------------------------------
# load simulation functions
source("R/2_run_simulations/SIRS-model.R")

# load model fit
fit_scotland_SIRS <- read_rds("data/derived_data/scotland/fit_scotland_SIRS.rds")

# load scotland base parameters
source("R/scotland_parameters.R")

# extract median parameters
fitted_pars = summary(fit_scotland_SIRS)$summary[1:14,6]

# create vector with rsv pars
fitted_pars_rsv = fitted_pars[c("S0_rsv", "I0_rsv", "rho_rsv", "b", "a", "p", "omega")]
names(fitted_pars_rsv) = c("S0", "I0", "rho", "b", "a", "p", "omega")
fitted_pars_rsv = c(fitted_pars_rsv, r = 1, gamma = 1/(1/scotland_sigma + 1/scotland_gamma), mu = scotland_birthrate)

# mpv pars
fitted_pars_mpv = fitted_pars[c("S0_mpv", "I0_mpv", "rho_mpv", "b", "a", "p", "omega", "r", "c")]
names(fitted_pars_mpv) = c("S0", "I0", "rho", "b", "a", "p", "omega", "r", "c")
fitted_pars_mpv = c(fitted_pars_mpv, gamma = 1/(1/scotland_sigma + 1/scotland_gamma), mu = scotland_birthrate)

n_wks = 52*40

c_slices = sort(c(c_slices, fitted_pars["c"]))
a_range = sort(c(a_range, fitted_pars["a"]))

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
    bifur_sims_ctmp[[i]] = SIRS_force_discrete_both(
      N = n_wks, parms_rsv = tst_pars_rsv, parms_mpv = tst_pars_mpv, pop = scotland_N, 
      npi = 1, testing_scalar_mpv = rep(1, n_wks),
      testing_scalar_rsv = rep(1, n_wks))
    bifur_points_ctmp[[i]] = bifur_sims_ctmp[[i]] %>%
      filter(t > 30*52, t%%52 == 12) %>%
      mutate(a = a_range[i], 
             c = c_slices[j])
  }
  bifur_sims[[j]] = bifur_sims_ctmp
  bifur_points[[j]] = bifur_points_ctmp
}

write_rds(bifur_sims, "data/derived_data/scotland/bifur-sims_scotland_SIRS.rds")
write_rds(bifur_points, "data/derived_data/scotland/bifur-points_scotland_SIRS.rds")

