library(dplyr)
library(reshape2)
library(ggplot2)

#### POST-PANDEMIC SIMULATIONS -------------------------------------------------
# load simulation functions
source("R/4_model-sensitivity/lags-SEIRS-model.R")

# load scotland base parameters
source("R/scotland_parameters.R")

# load google mobility data
source("R/2_run_simulations/google_mobility_functions.R")

# google mobility data
npi = get_google_mobility_data("GB", FALSE)

# load fits and simulate
fits_to_include = 1:8
postpand_sim = vector("list", length(fits_to_include))
for(z in 1:length(fits_to_include)){
  print(z)
  browser()
  f <- paste0("data/derived_data/scotland/fit_scotland_SEIRS_NEW_lag",fits_to_include[z],".rds")
  fit <- read_rds(f)
  postpand_sim[[z]] = postpandemic_analysis(
    fit = fit, samp = 1:100, observations = scotland_by_wk_full,
    mobility_df = npi, k = k_tst[z], pop = scotland_N, start_wk_pre = start_wk_pre, 
    end_wk_pre = end_wk_pre, t_to_wk_full = t_to_wk_scotland_full)
}

#### PLOT RESULTS --------------------------------------------------------------
# plot fits across lags
p = grid.arrange(grobs = lapply(1:length(postpand_sim), function(i){postpand_sim[[i]]$plot_sims + labs(subtitle = paste0(LETTERS[i], ") ", k_tst[i]-1, " week lag"))}), ncol = 1)
ggsave("figures/lagged_fits.pdf", p, width = 8, height = 12)

# plot out of sample performance
oos_labs = c("full period (2018-2024)", "pre-lockdown (2018-2020)", "post-lockdown (2021-2024)")
names(oos_labs) = c("full", "pre-pand", "post-pand")

performance_labs = c("out of sample log likelihood", "out of sample correlation")
names(performance_labs) = c("ll", "cor")

bind_rows(lapply(postpand_sim, function(i){i$performance})) %>% 
  mutate(pathogen = factor(pathogen, levels = c("rsv", "mpv")), 
         oos_period = factor(oos_period, levels = c("full", "pre-pand", "post-pand"))) %>%
  # add log likelihood (because it should be maximized like correlation)
  melt(c("pathogen", "oos_period", "variable", "k"), variable.name = "quantile") %>%
  mutate(value = ifelse(variable == "nll", -value, value), 
         variable = ifelse(variable == "nll", "ll", "cor")) %>%
  dcast(pathogen + oos_period + k + variable ~ quantile) %>%
  ggplot(aes(x = k-1, color = pathogen)) + 
  geom_segment(aes(xend = k-1, y = Q5, yend = Q95), linewidth = 0.5) + 
  geom_segment(aes(xend = k-1, y = Q25, yend = Q75), linewidth = 0.9) + 
  geom_point(aes(y = Q50)) + 
  facet_grid(cols = vars(oos_period), rows = vars(variable), 
             labeller = labeller(oos_period = oos_labs, 
                                 variable = performance_labs),
             scales = "free", switch = "y") + 
  labs(x = "lag in RSV forcing") +
  scale_color_brewer(palette = "Dark2") + 
  theme_bw() + 
  theme(axis.title.y = element_blank(), 
        legend.position = "bottom", 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(), 
        strip.placement = "outside")
ggsave("figures/lagged_postpand_performance.pdf", width = 8, height = 6)
  
  