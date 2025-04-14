library(dplyr)
library(reshape2)
library(ggplot2)
library(cowplot)

#### SETUP PARAMETERS ----------------------------------------------------------
# load simulation functions
source("R/2_run_simulations/SEIRS-model.R")

# load model fit
fit_scotland_SEIRS <- readRDS("data/derived_data/scotland/fit_scotland_SEIRS.rds")

# load scotland base parameters
source("R/scotland_parameters.R")

# create data.frame with rsv pars
# use last week as starting value for simulation
posterior_rsv = as.data.frame(fit_scotland_SEIRS)[,c("S0_rsv", "E0_rsv", "I0_rsv", 
                                                     "rho_rsv", "b", "a", "p")]
colnames(posterior_rsv) = c("S0", "E0", "I0", "rho", "b", "a", "p")
posterior_rsv = cbind(posterior_rsv, r = 1, c = 0, sigma = scotland_sigma, omega = scotland_omega,
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


n_wks = 52*40

bifur_sims = SEIRS_force_discrete_both(
  N = n_wks, parms_rsv = posterior_rsv, parms_mpv = posterior_mpv, pop = scotland_N, 
  npi = 1, testing_scalar_mpv = rep(1, n_wks),
  testing_scalar_rsv = rep(1, n_wks),
  rsvscaling_cuttoff = (end_wk_pre - start_wk_pre)/7 + 1)

posterior_mpv_noc = posterior_mpv
posterior_mpv_noc["c"] = 0
bifur_sims_noc = SEIRS_force_discrete_both(
  N = n_wks, parms_rsv = posterior_rsv, parms_mpv = posterior_mpv_noc, pop = scotland_N, 
  npi = 1, testing_scalar_mpv = rep(1, n_wks),
  testing_scalar_rsv = rep(1, n_wks),
  rsvscaling_cuttoff = (end_wk_pre - start_wk_pre)/7 + 1)

t_to_wk_scotland_bifur <- data.frame(wk_collected = as.Date(start_wk_pre + (0:n_wks)*7),
                                     t = 1:(n_wks+1))

diff_sims <- left_join(bifur_sims %>% filter(pathogen == "mpv", variable == "Ipred") %>% rename(fit = value), 
                       bifur_sims_noc %>% filter(pathogen == "mpv", variable == "Ipred") %>% 
                         rename(vacc = value)) %>%
  left_join(t_to_wk_scotland_bifur) %>%
  # filter(year(wk_collected) >= 2036 & year(wk_collected) < 2046) %>%
  mutate(year = lubridate::epiyear(wk_collected),
         epiwk_std = lubridate::epiweek(wk_collected),
         seas = ifelse(epiwk_std <= 40, paste(year-1, year, sep = "-"), paste(year, year + 1, sep = "-"))) %>%
  mutate(epiwk = ifelse(epiwk_std <= 40, epiwk_std + (52-40), epiwk_std-40)) %>%
  filter(seas %in% paste(2035:2044, 2036:2045, sep = "-"))

diff_sims <- diff_sims  %>%
  summarize(cum_fit = sum(fit), 
            cum_vacc = sum(vacc), 
            max_fit = max(fit), 
            max_vacc = max(vacc), 
            max_fit_wk = epiwk[fit == max(fit)], 
            max_vacc_wk = epiwk[vacc == max(vacc)], .by = c( "seas", "draw_id")) %>% #
  mutate(cum_diff = cum_vacc - cum_fit, 
         cum_diff_pct = (cum_vacc - cum_fit)/cum_fit, 
         max_diff = max_vacc - max_fit, 
         max_diff_pct = (max_vacc - max_fit)/max_fit, 
         max_wk_diff = max_vacc_wk - max_fit_wk) %>%
  reshape2::melt(c("seas", "draw_id")) %>%
  summarize(mn = median(value), .by = c("variable", "draw_id"))

diff_sims %>% 
  summarize(
    lwr = quantile(mn, 0.05), 
    med = median(mn), 
    upr = quantile(mn, 0.95), 
    .by = c("variable")
  )

#### PLOT PROBABILITY OF DOUBLING ----------------------------------------------
sims_summary <- bifur_sims %>% 
  filter(variable == "Ipred", t > 52*11, pathogen == "mpv") %>%
  summarize(cum_inc = sum(value), 
            peak_inc = max(value), .by = c("draw_id", "pathogen"))

sims_noint_summary <- bifur_sims_noc %>% 
  filter(variable == "Ipred", t > 52*11, pathogen == "mpv") %>%
  summarize(cum_inc_noint = sum(value), 
            peak_inc_noint = max(value), .by = c("draw_id", "pathogen"))


summ <- left_join(sims_summary, sims_noint_summary) %>% 
  mutate(cum_abs_change = (cum_inc_noint - cum_inc), 
         cum_rel_change = (cum_inc_noint - cum_inc)/cum_inc, 
         peak_abs_change = (peak_inc_noint - peak_inc), 
         peak_rel_change = (peak_inc_noint - peak_inc)/peak_inc, 
  )
  
reference = bifur_sims %>%
  filter(variable == "Ipred", t > 52*11) %>% 
  summarize(cum_inc = sum(value), 
            peak_inc = max(value), .by = c("draw_id", "pathogen")) %>%
  reshape2::melt(c("draw_id", "pathogen")) %>%
  reshape2::dcast(draw_id + variable ~ pathogen) %>% 
  mutate(pct_change = (rsv-mpv)/mpv) %>% 
  summarize(mean = mean(pct_change), .by = c("variable"))



p1 = data.frame(pct_increase = seq(0,2, length.out = 50), 
                probability = sapply(seq(0,2, length.out = 50), function(i){nrow(summ %>% filter(peak_rel_change > i))/nrow(summ)})) %>% 
  ggplot(aes(x = pct_increase, y = probability)) + 
  geom_line() +
  scale_x_continuous(labels = scales::percent, name = "% increase in annual hMPV peak magnitude") + 
  scale_y_continuous(name = "poseterior probability") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank())

p2 = data.frame(pct_increase = seq(0,0.03,length.out = 50), 
                probability = sapply(seq(0,0.03,length.out = 50), function(i){nrow(summ %>% filter(cum_rel_change > i))/nrow(summ)})) %>% 
  ggplot(aes(x = pct_increase, y = probability)) + 
  geom_line() +
  scale_x_continuous(labels = scales::percent, name = "% increase in annual hMPV burden") + 
  scale_y_continuous(name = "poseterior probability") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank())

plot_grid(p2, p1, nrow = 1, labels = LETTERS[1:2])
ggsave("figures/probability_of_hmpv_change.pdf", width = 8, height = 3)  
