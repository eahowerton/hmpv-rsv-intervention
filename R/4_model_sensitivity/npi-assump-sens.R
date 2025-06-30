library(bbmle)
library(dplyr)
library(reshape2)
library(ggplot2)
library(cowplot)

#### FUNCTIONS -----------------------------------------------------------------
# first find the npi scalar (in increment of 0.05) that gives best fit to
# out of sample data
fit_npi_scalar <- function(npi_scale_rsv, npi_scale_mpv, obs_rsv, obs_mpv, 
                           scotland_by_wk_full, posterior_rsv, posterior_mpv, 
                           scotland_N, end_wk_pre, start_wk_pre, summarize_across_drawids = TRUE){
  # print(npi_scale)
  npi_scaled <- scotland_by_wk_full %>%
    filter(pathogen == "mpv") %>%
    left_join(npi, by = join_by(wk_collected)) %>%
    mutate(mean_mob_rsv = ifelse(is.na(mean_mob), 1, 1+(mean_mob*npi_scale_rsv/100)), 
           mean_mob_mpv = ifelse(is.na(mean_mob), 1, 1+(mean_mob*npi_scale_mpv/100))) 
  pps <- SEIRS_force_discrete_both_sepnpi(
    N = nrow(npi_scaled), parms_rsv = posterior_rsv, parms_mpv = posterior_mpv,
    pop = scotland_N, npi_rsv = npi_scaled %>% pull(mean_mob_rsv), npi_mpv = npi_scaled %>% pull(mean_mob_mpv),
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
    mutate(ll = dpois(x = detections, lambda = value, log = TRUE), # calculate ll for each
           period = ifelse(wk_collected < as.Date("2020-03-01"), "pre", ifelse(wk_collected >= as.Date("2021-01-01"), "post", NA))
    ) 
  pred_mpv <- pps %>% 
    filter(wk_collected > end_wk_pre, pathogen == "mpv") %>%
    left_join(obs_mpv, by = join_by(wk_collected)) %>%
    mutate(ll = dpois(x = detections, lambda = value, log = TRUE), # calculate ll for each
           period = ifelse(wk_collected < as.Date("2020-03-01"), "pre", ifelse(wk_collected >= as.Date("2021-01-01"), "post", NA))
    ) 
  # sum log-likelihood
  summ = bind_rows(pred_rsv, pred_mpv) %>%
    summarize(nll = -sum(ll), 
              c = cor(log(value+1), log(detections+1)),
              .by = c("draw_id", "period", "pathogen")) %>%
    bind_rows(bind_rows(pred_rsv, pred_mpv) %>%
                summarize(nll = -sum(ll), 
                          c = cor(log(value+1), log(detections+1)),
                          .by = c("draw_id", "pathogen")) %>%
                mutate(period = "all")) %>%
    bind_rows(bind_rows(pred_rsv, pred_mpv) %>%
                summarize(nll = -sum(ll), 
                          c = cor(log(value+1), log(detections+1)),
                          .by = c("draw_id", "period")) %>%
                mutate(pathogen = "both")) %>%
    bind_rows(bind_rows(pred_rsv, pred_mpv) %>%
                summarize(nll = -sum(ll), 
                          c = cor(log(value+1), log(detections+1)),
                          .by = c("draw_id")) %>%
                mutate(period = "all", pathogen = "both")) %>%
    filter(!is.na(period))
  # browser()
  if(any(is.na(summ$c) | is.na(summ$nll))){browser()}
    if(summarize_across_drawids){
      summ = summ %>%
        reshape2::melt(c("period", "pathogen", "draw_id"), variable.name = "metric") %>%
        summarize(Q5 = quantile(value, 0.05),
                  Q25 = quantile(value, 0.25),
                  Q50 = median(value),
                  Q75 = quantile(value, 0.75),
                  Q95 = quantile(value, 0.95),
                  tot = sum(value),
                  .by = c("period", "pathogen", "metric"))
# 
#         summarize(nll = mean(nll),
#                   c = mean(c), .by = c("period", "pathogen"))
    }
  return(summ)  
} 

#### SETUP PARAMETERS ----------------------------------------------------------
source("R/2_run_simulations/SEIRS-model.R")

# load scotland base parameters
source("R/scotland_parameters.R")

# load google mobility data
source("R/2_run_simulations/google_mobility_functions.R")

# google mobility data
npi = get_google_mobility_data("GB", FALSE)

scotland_by_wk_full <- scotland_by_wk %>%
  filter(wk_collected >= start_wk_pre)

# define range of npi scalars to test
npi_scalar <- unique(round(sort(c(seq(0.4, 1.6, 0.1))),3))
npi_scalar_both <- expand.grid(rsv = npi_scalar, mpv = npi_scalar)

set.seed(100)
n_samp <- 5000
samp <- sample(1:(7500*4), n_samp)

# plot of dynamics under different scalars
expand.grid(wk_collected = scotland_by_wk_full %>% filter(wk_collected >= as.Date("2020-01-01"), wk_collected <= as.Date("2023-03-01")) %>% pull(wk_collected), 
            npi_scalar = npi_scalar[which(round(npi_scalar %% 0.2, 4) != 0.1)]) %>%
  left_join(npi) %>%
  mutate(mean_mob = ifelse(is.na(mean_mob), 1, 1+(mean_mob*npi_scalar/100))) %>%
  ggplot(aes(x = wk_collected, y = mean_mob, color = as.factor(npi_scalar))) + 
  geom_line() +
  coord_cartesian(ylim = c(0, 1.05)) + 
  theme_bw() +
  theme(legend.position = "bottom")

#### FULL NPI ANALYSIS FOR COUPLED MODEL ---------------------------------------
# load model fit
fit_scotland_SEIRS <- readRDS("data/derived_data/scotland/fit_scotland_SEIRS.rds")

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
posterior_rsv <- posterior_rsv[samp, ]
posterior_mpv <- posterior_mpv[samp, ]


npi_nll <- list()
for(i in 1:nrow(npi_scalar_both)){
  print(npi_scalar_both[i, ])
  npi_nll[[i]] = fit_npi_scalar(npi_scale_rsv = npi_scalar_both[i, "rsv"], 
                                npi_scale_mpv = npi_scalar_both[i, "mpv"],
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

npi_nll_coupled_winteraction <- bind_rows(npi_nll, .id = "id") %>%
  mutate(id = as.integer(id)) %>%
  left_join(npi_scalar_both %>% mutate(id = seq_len(n())))
saveRDS(npi_nll_coupled_winteraction, "data/derived_data/scotland/npi_sens_coupled_winteraction.rds")

#### FULL NPI ANALYSIS FOR COUPLED MODEL WITHOUT INTERACTION -------------------
posterior_mpv_noc = posterior_mpv
posterior_mpv_noc$c = 0 # remove interaction term from coupled model

npi_nll <- list()
for(i in 1:nrow(npi_scalar_both)){
  print(npi_scalar_both[i, ])
  npi_nll[[i]] = fit_npi_scalar(npi_scale_rsv = npi_scalar_both[i, "rsv"], 
                                npi_scale_mpv = npi_scalar_both[i, "mpv"],
                                obs_rsv = scotland_by_wk %>% 
                                  filter(wk_collected >= start_wk_post, pathogen == "rsv") %>% 
                                  select(wk_collected, detections), 
                                obs_mpv = scotland_by_wk %>% 
                                  filter(wk_collected >= start_wk_post, pathogen == "mpv") %>% 
                                  select(wk_collected, detections), 
                                scotland_by_wk_full = scotland_by_wk_full, 
                                posterior_rsv = posterior_rsv, 
                                posterior_mpv = posterior_mpv_noc, 
                                scotland_N = scotland_N, 
                                end_wk_pre = end_wk_pre, 
                                start_wk_pre = start_wk_pre)
}

npi_nll_coupled_nointeraction <- bind_rows(npi_nll, .id = "id") %>%
  mutate(id = as.integer(id)) %>%
  left_join(npi_scalar_both %>% mutate(id = seq_len(n())))
saveRDS(npi_nll_coupled_nointeraction, "data/derived_data/scotland/npi_sens_coupled_nointeraction.rds")

#### FULL NPI ANALYSIS FOR INDEPENDENT MODEL  ----------------------------------
# load model fit
fit_scotland_SEIRS_rsv <- readRDS("data/derived_data/scotland/fit_scotland_SEIRS_rsv.rds")
fit_scotland_SEIRS_mpv <- readRDS("data/derived_data/scotland/fit_scotland_SEIRS_mpv.rds")

posterior_rsv_ind = as.data.frame(fit_scotland_SEIRS_rsv)[,c("S0", "E0", "I0", 
                                                         "rho", "b", "a", "p")]
posterior_rsv_ind = cbind(posterior_rsv_ind, r = 1, sigma = scotland_sigma, omega = scotland_omega,
                      gamma = scotland_gamma, mu = scotland_birthrate)

# mpv pars
posterior_mpv_ind = as.data.frame(fit_scotland_SEIRS_mpv)[,c("S0", "E0","I0",
                                                         "rho", "b", "a", "p")]
posterior_mpv_ind = cbind(posterior_mpv_ind, sigma = scotland_sigma,  omega = scotland_omega,
                      gamma = scotland_gamma, mu = scotland_birthrate, r = 1, c = 0) # r = 1, c = 0 because independent

# randomly select a subset of simulations
posterior_rsv_ind <- posterior_rsv_ind[samp, ]
posterior_mpv_ind <- posterior_mpv_ind[samp, ]

npi_nll <- list()
for(i in 1:nrow(npi_scalar_both)){
  print(npi_scalar_both[i, ])
  npi_nll[[i]] = fit_npi_scalar(npi_scale_rsv = npi_scalar_both[i, "rsv"], 
                                npi_scale_mpv = npi_scalar_both[i, "mpv"],
                                obs_rsv = scotland_by_wk %>% 
                                  filter(wk_collected >= start_wk_post, pathogen == "rsv") %>% 
                                  select(wk_collected, detections), 
                                obs_mpv = scotland_by_wk %>% 
                                  filter(wk_collected >= start_wk_post, pathogen == "mpv") %>% 
                                  select(wk_collected, detections), 
                                scotland_by_wk_full = scotland_by_wk_full, 
                                posterior_rsv = posterior_rsv_ind, 
                                posterior_mpv = posterior_mpv_ind, 
                                scotland_N = scotland_N, 
                                end_wk_pre = end_wk_pre, 
                                start_wk_pre = start_wk_pre)
}

npi_nll_independent <- bind_rows(npi_nll, .id = "id") %>%
  mutate(id = as.integer(id)) %>%
  left_join(npi_scalar_both %>% mutate(id = seq_len(n())))
saveRDS(npi_nll_independent, "data/derived_data/scotland/npi_sens_coupled_independent.rds")

#### PLOT NLL RESULTS ----------------------------------------------------------

df = bind_rows(npi_nll_coupled_winteraction %>% mutate(model = "coupled, with interaction"), 
               npi_nll_coupled_nointeraction %>% mutate(model = "coupled, without interaction")) %>%
  bind_rows(npi_nll_independent %>% mutate(model = "independent")) %>%
  filter(metric == "nll", period == "all") %>%
  mutate(value = tot/n_samp) # determine which value to plot
df_opt_per_path = df %>% 
  mutate(min_value = min(value), .by = c("model", "pathogen", "metric")) %>% 
  filter(value == min_value, period == "all")
df_opt = df %>% 
  mutate(min_value = min(value), .by = c("pathogen", "metric")) %>% 
  filter(value == min_value, period == "all")
df_opt_per_path_sameNPI = df %>% 
  filter(rsv == mpv) %>%
  mutate(min_value = min(value), .by = c("model", "pathogen", "metric")) %>% 
  filter(value == min_value, period == "all")
df_opt_sameNPI = df %>% 
  filter(rsv == mpv) %>%
  mutate(min_value = min(value), .by = c("pathogen", "metric")) %>% 
  filter(value == min_value, period == "all")

## likelihood ratio test
# df_lrt = df %>%
#   filter(metric == "nll", period == "all") %>%
#   # select(i, pathogen, rsv, mpv, model) %>%
#   left_join(df_opt %>% select(pathogen, min_value) %>% unique()) %>%
#   mutate(lr_stat = pchisq(-2*(-value+min_value), 2))

p1 = ggplot(data = df %>% filter(pathogen == "rsv"), aes(x = rsv, y = mpv, fill = value)) + 
  geom_tile() + 
  geom_point(data = df_opt_per_path %>% filter(pathogen == "rsv"), shape = 20, color = "black", size = 1) + 
  geom_point(data = df_opt_sameNPI %>% filter(pathogen == "rsv"), shape = 1, color = "black", size = 3) + 
  geom_point(data = df_opt %>% filter(pathogen == "rsv"), shape = 0, color = "black", size = 5) + 
  ggtitle("negative log-likelihood") + 
  facet_grid(cols = vars(pathogen), rows = vars(model), switch = "y") +
  scale_fill_distiller(palette = "Reds", direction = 0, name = "RSV") + 
  scale_x_continuous(expand = c(0,0), name = " ") + 
  scale_y_continuous(expand = c(0,0), name = "hMPV NPI scalar") + 
  theme(legend.key.width = unit(0.75, "cm"), 
        legend.position = "bottom", 
        strip.background = element_blank(), 
        strip.placement = "outside", 
        strip.text.x = element_blank())
p2 = ggplot(data = df %>% filter(pathogen == "mpv"), aes(x = rsv, y = mpv, fill = value)) + 
  geom_tile() + 
  geom_point(data = df_opt_per_path %>% filter(pathogen == "mpv"), shape = 20, color = "black", size = 1) + 
  geom_point(data = df_opt_sameNPI %>% filter(pathogen == "mpv"), shape = 1, color = "black", size = 3) + 
  geom_point(data = df_opt %>% filter(pathogen == "mpv"), shape = 0, color = "black", size = 5) + 
  facet_grid(cols = vars(pathogen), rows = vars(model)) +
  scale_fill_distiller(palette = "Blues", direction = 0, name = "hMPV") + 
  scale_x_continuous(expand = c(0,0), name = "RSV NPI scalar") + 
  scale_y_continuous(expand = c(0,0), name = "hMPV NPI scalar") + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), 
        legend.key.width = unit(0.75, "cm"), 
        legend.position = "bottom", 
        strip.background = element_blank(), 
        strip.text = element_blank())
p3 = ggplot(data = df %>% filter(pathogen == "both"),
            aes(x = rsv, y = mpv, fill = value)) + 
  geom_tile() + 
  geom_point(data = df_opt_per_path %>% filter(pathogen == "both"), shape = 20, color = "black", size = 1) + 
  geom_point(data = df_opt_sameNPI %>% filter(pathogen == "both"), shape = 1, color = "black", size = 3) + 
  geom_point(data = df_opt %>% filter(pathogen == "both"), shape = 0, color = "black", size = 5) + 
  facet_grid(cols = vars(pathogen), rows = vars(model)) +
  scale_fill_distiller(palette = "Purples", direction = 0, name = "both") + 
  scale_x_continuous(expand = c(0,0), name = "") + 
  scale_y_continuous(expand = c(0,0), name = "hMPV NPI scalar") + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), 
        legend.key.width = unit(0.75, "cm"), 
        legend.position = "bottom", 
        strip.background = element_blank(), 
        strip.text = element_blank())
p = cowplot::plot_grid(p1, p2, p3, nrow = 1, rel_widths = c(0.4, 0.3, 0.3), align = "h", axis = "tb")
ggsave("figures/model_comparison_allNPI_nll.pdf", p, width = 8, height = 8)

#### PLOT CORRELATION RESULTS --------------------------------------------------
df = bind_rows(npi_nll_coupled_winteraction %>% mutate(model = "coupled, with interaction"), 
               npi_nll_coupled_nointeraction %>% mutate(model = "coupled, without interaction")) %>%
  bind_rows(npi_nll_independent %>% mutate(model = "independent")) %>%
  filter(metric == "c", period == "all") %>%
  mutate(value = tot/n_samp, 
         pathogen = factor(pathogen, levels = c("rsv", "mpv", "both"))) # determine which value to plot
df_opt_per_path = df %>% 
  mutate(max_value = max(value), .by = c("model", "pathogen", "metric")) %>% 
  filter(value == max_value, period == "all")
df_opt = df %>% 
  mutate(max_value = max(value), .by = c("pathogen", "metric")) %>% 
  filter(value == max_value, period == "all")
df_opt_per_path_sameNPI = df %>% 
  filter(rsv == mpv) %>%
  mutate(max_value = max(value), .by = c("model", "pathogen", "metric")) %>% 
  filter(value == max_value, period == "all")
df_opt_sameNPI = df %>% 
  filter(rsv == mpv) %>%
  mutate(max_value = max(value), .by = c("pathogen", "metric")) %>% 
  filter(value == max_value, period == "all")

pathogen_labs = c("RSV", "hMPV", "both")
names(pathogen_labs) = c("rsv", "mpv", "both")

ggplot(data = df, aes(x = rsv, y = mpv, fill = value)) + 
  geom_tile() + 
  geom_point(data = df_opt_per_path, shape = 20, color = "black", size = 1) + 
  # geom_point(data = df_opt_per_path_sameNPI %>% filter(pathogen == "rsv"), shape = 8, color = "red", size = 2) + 
  geom_point(data = df_opt_sameNPI, shape = 1, color = "black", size = 3) + 
  geom_point(data = df_opt, shape = 0, color = "black", size = 5) + 
  facet_grid(cols = vars(pathogen), rows = vars(model), switch = "y", labeller = labeller(pathogen = pathogen_labs)) +
  scale_fill_viridis_c(name = "correlation") +
  # scale_fill_distiller(palette = "Reds", name = "RSV") + 
  scale_x_continuous(expand = c(0,0), name = "RSV NPI scalar") + 
  scale_y_continuous(expand = c(0,0), name = "hMPV NPI scalar") + 
  theme(legend.key.width = unit(1, "cm"), 
        legend.position = "bottom",
        strip.background = element_blank(), 
        strip.placement = "outside")
ggsave("figures/model_comparison_allNPI_corr.pdf", width = 8, height = 8)


### DEMONSTRATE DIFFERENCES BETWEEN MODEL WITH AND WITHOUT INTERACTION ---------
npi_scalar_rsv = 1.1
npi_scalar_mpv_high = 1.1
npi_scalar_mpv_low = 0.6

npi_scaled_rsv <- scotland_by_wk_full %>%
  filter(pathogen == "mpv") %>%
  left_join(npi) %>%
  mutate(mean_mob = ifelse(is.na(mean_mob), 1, 1+(mean_mob*npi_scalar_rsv/100))) %>%
  pull(mean_mob)
npi_scaled_mpv_high <- scotland_by_wk_full %>%
  filter(pathogen == "mpv") %>%
  left_join(npi) %>%
  mutate(mean_mob = ifelse(is.na(mean_mob), 1, 1+(mean_mob*npi_scalar_mpv_high/100))) %>%
  pull(mean_mob)
npi_scaled_mpv_low <- scotland_by_wk_full %>%
  filter(pathogen == "mpv") %>%
  left_join(npi) %>%
  mutate(mean_mob = ifelse(is.na(mean_mob), 1, 1+(mean_mob*npi_scalar_mpv_low/100))) %>%
  pull(mean_mob)

postpand_sims <- SEIRS_force_discrete_both_sepnpi(
  N = length(npi_scaled_rsv), parms_rsv = posterior_rsv, parms_mpv = posterior_mpv,
  pop = scotland_N, npi_rsv = npi_scaled_rsv, npi_mpv = npi_scaled_mpv_high, 
  testing_scalar_mpv = scotland_by_wk_full %>% filter(pathogen == "mpv") %>% pull(n_tests_ma_scaled),
  testing_scalar_rsv = scotland_by_wk_full %>% filter(pathogen == "rsv") %>% pull(n_tests_ma_scaled),
  rsvscaling_cuttoff = (end_wk_pre - start_wk_pre)/7 + 1
)

postpand_sims_noc_high <- SEIRS_force_discrete_both_sepnpi(
  N = length(npi_scaled_rsv), parms_rsv = posterior_rsv, parms_mpv = posterior_mpv_noc,
  pop = scotland_N, npi_rsv = npi_scaled_rsv, npi_mpv = npi_scaled_mpv_high, 
  testing_scalar_mpv = scotland_by_wk_full %>% filter(pathogen == "mpv") %>% pull(n_tests_ma_scaled),
  testing_scalar_rsv = scotland_by_wk_full %>% filter(pathogen == "rsv") %>% pull(n_tests_ma_scaled),
  rsvscaling_cuttoff = (end_wk_pre - start_wk_pre)/7 + 1
)

postpand_sims_noc_low <- SEIRS_force_discrete_both_sepnpi(
  N = length(npi_scaled_rsv), parms_rsv = posterior_rsv, parms_mpv = posterior_mpv_noc,
  pop = scotland_N, npi_rsv = npi_scaled_rsv, npi_mpv = npi_scaled_mpv_low, 
  testing_scalar_mpv = scotland_by_wk_full %>% filter(pathogen == "mpv") %>% pull(n_tests_ma_scaled),
  testing_scalar_rsv = scotland_by_wk_full %>% filter(pathogen == "rsv") %>% pull(n_tests_ma_scaled),
  rsvscaling_cuttoff = (end_wk_pre - start_wk_pre)/7 + 1
)

plot_df = bind_rows(
  postpand_sims %>% 
    filter(variable == "Ipred") %>% 
    mutate(model = "with interaction\nhMPV NPI = 1.1"), 
  postpand_sims_noc_high %>% 
    filter(variable == "Ipred") %>% 
    mutate(model = "without interaction\nhMPV NPI = 1.1"), 
  postpand_sims_noc_low %>% 
    filter(variable == "Ipred") %>% 
    mutate(model = "without interaction\nhMPV NPI = 0.6")
) %>%
  left_join(t_to_wk_scotland_full) %>% 
  filter(wk_collected >= start_wk_post, draw_id %in% sample(1:n_samp, 100))

p2 = ggplot(data = plot_df %>%  filter(wk_collected >= start_wk_post, pathogen == "mpv"), 
             aes(x = wk_collected, y = value)) +
  geom_line(aes(group = draw_id), color = "black", alpha = 0.02) + 
  geom_point(data = scotland_by_wk_full %>% filter(wk_collected >= start_wk_post, pathogen == "mpv"), 
             aes(y = detections, color = pathogen), shape = 21, color = "#D95F02") + 
  facet_grid(rows = vars(model), cols = vars(pathogen), switch = "y", labeller = labeller(pathogen = pathogen_labs)) + 
  labs(y = "hMPV detections") +
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        legend.position = "none", 
        panel.grid.minor = element_blank(), 
        strip.background = element_blank(), 
        strip.placement = "outside",
        strip.text.x = element_blank()
        )

ll_check =  bind_rows(
  postpand_sims %>% 
    filter(variable == "Ipred") %>% 
    mutate(model = "with interaction\nhMPV NPI = 1.1"), 
  postpand_sims_noc_high %>% 
    filter(variable == "Ipred") %>% 
    mutate(model = "without interaction\nhMPV NPI = 1.1"), 
  postpand_sims_noc_low %>% 
    filter(variable == "Ipred") %>% 
    mutate(model = "without interaction\nhMPV NPI = 0.6")
) %>%
  left_join(t_to_wk_scotland_full) %>% 
  filter(wk_collected >= start_wk_post) %>%
  left_join(scotland_by_wk_full) %>% 
  mutate(ll = dpois(detections, value, log = TRUE))

p3 = ll_check %>%
  mutate(year = lubridate::year(wk_collected)) %>%
  filter(pathogen == "mpv", wk_collected >= start_wk_post) %>%
  summarize(ll = -sum(ll)/n_samp, .by = c("model", "wk_collected")) %>%
  mutate(cum_ll = cumsum(ll), .by = c("model")) %>%
  ggplot(aes(x = wk_collected, y = cum_ll, color = model)) + 
  geom_line(size = 0.8) + 
  labs(y = "cumulative negative log likelihood") + 
  scale_color_manual(values = c("#56B4E9","#009E73", "#F5C710")) +
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        panel.grid.minor = element_blank())

plot_grid(p, 
          plot_grid(p2, p3, ncol = 1, labels = c("B", "C"), 
                    align = "v", axis = "lr"), rel_widths = c(0.65, 0.35), 
          nrow = 1, labels = c("A", NA))  
ggsave("figures/model_comparison_allNPI_nll_full.pdf", width = 14, height = 8)


