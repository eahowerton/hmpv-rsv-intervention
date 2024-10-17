library(ggplot2)
library(cowplot)
library(gridExtra)
library(readr)
library(tidyr)
library(scales)
library(dplyr)
library(reshape2)

# load simulations
# postpand_sims <- readRDS("data/derived_data/scotland/postpand-sims_scotland_SEIRS.rds")
source("R/2_run_simulations/simulate-postpand_scotland_SEIRS.R")
plot_samps <- sample(samp,100)
source("R/2_run_simulations/simulate-postpand_scotland_SEIRS-ind.R")
plot_samps_ind <- sample(samp_ind,100)

# load scotland base parameters and data
source("R/scotland_parameters.R")

# load google mobility data
source("R/2_run_simulations/google_mobility_functions.R")
npi = get_google_mobility_data("GB", FALSE)

# create data.frame with npi scalars
df_plot_npis <- expand.grid(npi_scalar = npi_scalar, # npi_scalar optimized in simulate-postpand_scotland_SEIRS.R
                            wk_collected = seq.Date(start_wk_post, end_wk_post, by = 7)) %>%
  left_join(npi) %>%
  mutate(mean_mob = ifelse(is.na(mean_mob), 1, 1+(mean_mob*npi_scalar/100))) %>%
  mutate(plot_npi = 50 + mean_mob*10, 
         pathogen = "mpv")

facet_var_lvls = c("S-rsv", "S-mpv", "I-rsv", "I-mpv", "Ipred-rsv", "Ipred-mpv")

facet_var_labs = c("RSV susceptible", "HMPV susceptible", "RSV infected", "HMPV infected", "RSV detections", "HMPV detections")
names(facet_var_labs) = facet_var_lvls

postpand_sims %>% 
  # filter(npi_scalar == npi_scalar) %>%
  left_join(t_to_wk_scotland_full) %>%
  filter(draw_id %in% plot_samps, 
         variable %in% c("S", "I", "Ipred")) %>%
  mutate(facet_var = factor(paste(variable, pathogen, sep = "-"), levels = facet_var_lvls)) %>%
  # mutate(variable = factor(variable, levels = c("S", "I", "Ipred"))) %>%
  ggplot(aes(x = wk_collected, color = pathogen)) + 
  geom_vline(xintercept = start_wk_post) +
  geom_line(data = postpand_sims_ind %>% 
              left_join(t_to_wk_scotland_full) %>%
              filter(draw_id %in% plot_samps_ind, 
                     variable %in% c("Ipred"), 
                     wk_collected >= start_wk_post) %>%
              mutate(facet_var = factor(paste(variable, pathogen, sep = "-"), levels = facet_var_lvls)), 
            aes(y = value, group = draw_id), alpha = 0.2, color = "gray30") +
  geom_line(aes(y = value, group = draw_id), alpha = 0.2) +
  geom_point(data = scotland_by_wk %>%
               filter(wk_collected > start_wk_pre) %>%
               mutate(facet_var = paste("Ipred", pathogen, sep = "-")) %>%
               mutate(facet_var = factor(facet_var, levels = facet_var_lvls)),
             aes(y = detections), shape = 21, fill = "white", alpha = 0.6, size = 0.9) +
  geom_line(data = df_plot_npis %>% 
              filter(npi_scalar == npi_scalar) %>%
              mutate(facet_var = factor("S-rsv", levels = facet_var_lvls), 
                     npi_plot = 3e6-4e5 + mean_mob*4e5),
            aes(y = npi_plot), color = "darkgray") + 
  geom_line(data = df_plot_npis %>% 
              filter(npi_scalar == npi_scalar) %>%
              mutate(facet_var = factor("S-mpv", levels = facet_var_lvls), 
                     npi_plot = 3e6-4e5 + mean_mob*4e5),
            aes(y = npi_plot), color = "darkgray") + 
  geom_text(data = data.frame(y = 3e6 + 2e5, 
                              x = as.Date(rep(c(mean(c(start_wk_pre, end_wk_pre)), 
                                    mean(c(start_wk_post, end_wk_post))),2), "%Y-%m-%d"), 
                              label = rep(c("fit pre-pandemic period", "simulate post-pandemic period"), 2), 
                              facet_var = factor(c("S-mpv", "S-mpv","S-rsv","S-rsv"), levels = facet_var_lvls)), 
            aes(x = x, y = y, label = label), vjust = 1, color = "black", size = 2) + 
  facet_wrap(vars(facet_var), labeller = labeller(facet_var = facet_var_labs), scales = "free", ncol = 2) + 
  scale_color_brewer(palette = "Dark2", direction = -1) + 
  scale_x_date(date_breaks = "2 years", 
               date_labels = "%Y", 
               expand = c(0,0)) +
  scale_y_continuous(labels = comma) + 
  theme_bw() + 
  theme(axis.title = element_blank(), 
        legend.position = "none", 
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank())

ggsave("figures/full-timeseries_SEIRS_scotland.pdf", width = 10, height = 6)

