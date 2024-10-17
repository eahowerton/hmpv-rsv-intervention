library(ggplot2)
library(cowplot)
library(gridExtra)
library(readr)
library(tidyr)
library(scales)

source("R/scotland_parameters.R")

bifur_points <- read_rds("data/derived_data/scotland/bifur-points_scotland_SIRS.rds")
bifur_sims <- read_rds("data/derived_data/scotland/bifur-sims_scotland_SIRS.rds")

fit_scotland_SIRS <- read_rds("data/derived_data/scotland/fit_scotland_SIRS.rds")
fitted_pars = summary(fit_scotland_SIRS)$summary[1:14,6]

c_slices = sort(c(c_slices, fitted_pars["c"]))
a_range = sort(c(a_range, fitted_pars["a"]))

a_ex = sort(c(fitted_pars["a"], 0.3,0.45, 0.8))

c_names = c("c = 0 (vaccination)", paste("c =", round(fitted_pars["c"],2) , "(fitted)"))

p1 = bind_rows(bifur_points) %>%
  select(a, c, t, I_rsv) %>%
  melt(c("a", "c", "t")) %>%
  filter(c == 0) %>%
  ggplot(aes(x = a, y = value/scotland_N)) + 
  geom_point() + 
  geom_vline(xintercept = a_ex, linetype = "dotted", size = 0.8) + 
  labs(x = "amplitude of seasonal forcing", 
       y = "RSV prevalence") + 
  theme_bw()

p2 = bind_rows(bifur_points) %>%
  select(a, c, t, I_mpv) %>%
  melt(c("a", "c", "t")) %>%
  filter(c %in% c(0, fitted_pars["c"])) %>%
  mutate(c_name = factor(ifelse(c == 0, c_names[1], c_names[2]), levels = c_names)) %>%
  ggplot(aes(x = a, y = value/scotland_N, color = as.factor(c_name))) + 
  geom_point() +
  geom_vline(xintercept = a_ex, linetype = "dotted", size = 0.8) + 
  labs(x = "amplitude of seasonal forcing", 
       y = "HMPV prevalence") + 
  scale_color_manual(values = c("gray", "black"))+
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank())

t_to_wk_scotland_bifur <- data.frame(wk_collected = as.Date(start_wk_pre + (0:nrow(bifur_sims[[1]][[1]]))*7),
                                     t = 1:(nrow(bifur_sims[[1]][[1]])+1))

aexplts_I <- bind_rows(
  bifur_sims[[which(c_slices == fitted_pars["c"])]][sapply(a_ex, function(i){which(a_range == i)})], 
  .id = "draw_id") %>%
  mutate(c =fitted_pars["c"] ) %>%
  bind_rows(
    bind_rows(
      bifur_sims[[which(c_slices == 0)]][sapply(a_ex, function(i){which(a_range == i)})], 
      .id = "draw_id") %>%
      mutate(c = 0)) %>%
  mutate(a = a_ex[as.integer(draw_id)]) %>%
  select(t, a, c, I_rsv, I_mpv) %>% #, Ipred_rsv, Ipred_mpv
  melt(c("t", "a", "c")) %>%
  tidyr::separate(variable, into = c("variable", "pathogen")) %>%
  left_join(t_to_wk_scotland_bifur) %>%
  filter(year(wk_collected) > 2036) %>%
  mutate(c_name = ifelse(c == 0, "c = 0 (vaccination)", paste("c =", round(c,2), "(fitted)"))) %>%
  ggplot(aes(x = wk_collected, y = value, color = pathogen)) + 
  geom_line(aes(linetype = c_name)) + 
  facet_wrap(vars(paste("a = ", round(a,2))), ncol = 1, scales = "free") + 
  scale_color_brewer(palette = "Dark2", 
                     labels = c("HMPV", "RSV"),
                     name = "pathogen") + 
  scale_linetype_manual(values = c("solid", "dashed"), 
                        name = "cross protection scenario") + 
  scale_x_date(date_breaks = "year") +
  scale_y_continuous(labels = comma, name = "prevalence") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.minor = element_blank(), 
        strip.background = element_blank())

plot_grid(plot_grid(p1, p2, ncol = 1), aexplts_I)
ggsave("figures/bifurcation_SIRS_scotland_I.pdf", width = 14, height = 8)


# repeat for Ipred

aexplts_Ipred <- bind_rows(
  bifur_sims[[which(c_slices == fitted_pars["c"])]][sapply(a_ex, function(i){which(a_range == i)})], 
  .id = "draw_id") %>%
  mutate(c =fitted_pars["c"] ) %>%
  bind_rows(
    bind_rows(
      bifur_sims[[which(c_slices == 0)]][sapply(a_ex, function(i){which(a_range == i)})], 
      .id = "draw_id") %>%
      mutate(c = 0)) %>%
  mutate(a = a_ex[as.integer(draw_id)]) %>%
  select(t, a, c, Ipred_rsv, Ipred_mpv) %>% #
  melt(c("t", "a", "c")) %>%
  tidyr::separate(variable, into = c("variable", "pathogen")) %>%
  left_join(t_to_wk_scotland_bifur) %>%
  filter(year(wk_collected) > 2036) %>%
  mutate(c_name = ifelse(c == 0, "c = 0 (vaccination)", paste("c =", round(c,2), "(fitted)"))) %>%
  ggplot(aes(x = wk_collected, y = value, color = pathogen)) + 
  geom_line(aes(linetype = c_name)) + 
  facet_wrap(vars(paste("a = ", round(a,2))), ncol = 1, scales = "free") + 
  scale_color_brewer(palette = "Dark2", 
                     labels = c("HMPV", "RSV"),
                     name = "pathogen") + 
  scale_linetype_manual(values = c("solid", "dashed"), 
                        name = "cross protection scenario") + 
  scale_x_date(date_breaks = "year") +
  scale_y_continuous(labels = comma, name = "detections") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.minor = element_blank(), 
        strip.background = element_blank())

plot_grid(plot_grid(p1, p2, ncol = 1), aexplts_Ipred)
ggsave("figures/bifurcation_SIRS_scotland_Ipred.pdf", width = 14, height = 8)

