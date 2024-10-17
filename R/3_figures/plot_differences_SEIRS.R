library(ggplot2)
library(cowplot)
library(gridExtra)
library(readr)
library(tidyr)
library(scales)

source("R/scotland_parameters.R")

bifur_points <- readRDS("data/derived_data/scotland/bifur-points_scotland_SEIRS.rds")
bifur_sims <- readRDS("data/derived_data/scotland/bifur-sims_scotland_SEIRS.rds")

fit_scotland_SEIRS <- readRDS("data/derived_data/scotland/fit_scotland_SEIRS.rds")
fitted_pars = apply(as.array(fit_scotland_SEIRS)[,,1:14], 3, median)

c_slices = sort(c(0, fitted_pars["c"]*0.9, fitted_pars["c"], fitted_pars["c"]*1.1))
a_range = sort(c(seq(0.01, 0.99, 0.001), fitted_pars["a"]))

t_to_wk_scotland_bifur <- data.frame(wk_collected = as.Date(start_wk_pre + (0:nrow(bifur_sims[[1]][[1]]))*7),
                                     t = 1:(nrow(bifur_sims[[1]][[1]])+1))

# calculate the change in vaccination for different levels of a
diff_df <- bind_rows(
  bifur_sims[[which(c_slices == fitted_pars["c"])]], 
  .id = "draw_id") %>%
  mutate(c = fitted_pars["c"]) %>%
  bind_rows(
    bind_rows(
      bifur_sims[[which(c_slices == 0)]], 
      .id = "draw_id") %>%
      mutate(c = 0)) %>%
  bind_rows(
    bind_rows(
      bifur_sims[[which(c_slices == fitted_pars["c"]*0.9)]], 
      .id = "draw_id") %>%
      mutate(c = fitted_pars["c"]*0.9)) %>%
  bind_rows(
    bind_rows(
      bifur_sims[[which(c_slices == fitted_pars["c"]*1.1)]], 
      .id = "draw_id") %>%
      mutate(c = fitted_pars["c"]*1.1)) %>%
  mutate(a = a_range[as.integer(draw_id)]) %>%
  filter(pathogen == "mpv", variable == "Ipred") %>%
  select(t, a, c, value) %>% 
  mutate(c = ifelse(c == 0, "vacc", paste0("fit", round(c,2)))) %>%
  dcast(t + a ~ c, value.var = "value") %>%
  left_join(t_to_wk_scotland_bifur) %>%
  # filter(year(wk_collected) >= 2036 & year(wk_collected) < 2046) %>%
  mutate(year = lubridate::epiyear(wk_collected),
         epiwk_std = lubridate::epiweek(wk_collected),
         seas = ifelse(epiwk_std <= 40, paste(year-1, year, sep = "-"), paste(year, year + 1, sep = "-"))) %>%
  mutate(epiwk = ifelse(epiwk_std <= 40, epiwk_std + (52-40), epiwk_std-40)) %>%
  filter(seas %in% paste(2035:2044, 2036:2045, sep = "-"))

diff_df <- diff_df  %>%
  melt(c("t", "a", "vacc", "wk_collected", "year", "epiwk_std", "seas", "epiwk")) %>%
  rename(fit = value, 
         cfit = variable) %>%
  summarize(cum_fit = sum(fit), 
            cum_vacc = sum(vacc), 
            max_fit = max(fit), 
            max_vacc = max(vacc), 
            max_fit_wk = epiwk[fit == max(fit)], 
            max_vacc_wk = epiwk[vacc == max(vacc)], .by = c("a", "seas", "cfit")) %>% #
  mutate(cum_diff = cum_vacc - cum_fit, 
         cum_diff_pct = (cum_vacc - cum_fit)/cum_fit, 
         max_diff = max_vacc - max_fit, 
         max_diff_pct = (max_vacc - max_fit)/max_fit, 
         max_wk_diff = max_vacc_wk - max_fit_wk) %>%
  reshape2::melt(c("a", "seas", "cfit")) %>%
  summarize(mn = median(value), .by = c("a", "cfit", "variable"))
  

labs_diffs <- c("difference in annual HMPV burden (reported cases)", "relative change in annual burden", "difference in annual HMPV peak magnitude (reported cases)", "relative change in peak magnitude (reported cases)", "difference in annual HMPV peak timing (weeks)")
names(labs_diffs) <- c("cum_diff", "cum_diff_pct", "max_diff", "max_diff_pct", "max_wk_diff")

labs_abs <- c("annual HMPV burden (reported cases)", "annual HMPV peak magnitude (reported cases)", "annual HMPV peak timing (weeks)")
names(labs_abs) <- c("cum", "max", "max_wk")

p1 = diff_df %>%
  filter(a < 0.5) %>%
  mutate(c = as.double(substr(cfit, 4, nchar(as.character(cfit))))) %>%
  mutate(cfit = paste("c =", round(c, 2))) %>%
  mutate(cfit = ifelse(c == round(fitted_pars["c"],2), paste(cfit, "(fitted)"), cfit)) %>%
  filter(variable %in% c("cum_fit", "cum_vacc", "max_fit", "max_vacc", "max_fit_wk", "max_vacc_wk")) %>%
  mutate(cfit = ifelse(variable %in% c("cum_vacc", "max_vacc", "max_vacc_wk"), "c=0 (vaccination)", cfit),
         variable = ifelse(variable %in% c("max_fit_wk", "max_vacc_wk"), "max_wk", substr(variable,1,3))) %>%
  ggplot(aes(x = a)) + 
  geom_vline(xintercept = fitted_pars_mpv["a"], linetype = "dashed") + 
  geom_line(aes(y = mn, color = cfit), linewidth = 0.8) + 
  facet_wrap(vars(variable), labeller = labeller(variable = labs_abs), scales = "free", ncol = 1) + 
  labs(x = "amplitude of seasonal forcing") + 
  scale_color_manual(values = c(RColorBrewer::brewer.pal(3, "Set1"),"black")) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw() + 
  theme(axis.title.y = element_blank(), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

p2 = diff_df %>%
  filter(a < 0.5) %>%
  mutate(c = as.double(substr(cfit, 4, nchar(as.character(cfit))))) %>%
  mutate(cfit = paste("c =", round(c, 2))) %>%
  mutate(cfit = ifelse(c == round(fitted_pars["c"],2), paste(cfit, "(fitted)"), cfit)) %>%
  filter(variable %in% c("cum_diff", "max_diff", "max_wk_diff")) %>%
  ggplot(aes(x = a)) + 
  geom_vline(xintercept = fitted_pars_mpv["a"], linetype = "dashed") + 
  geom_line(aes(y = mn, color = cfit), linewidth = 0.8) + 
  facet_wrap(vars(variable), labeller = labeller(variable = labs_diffs), scales = "free", ncol = 1) + 
  labs(x = "amplitude of seasonal forcing") + 
  scale_color_manual(values = c(RColorBrewer::brewer.pal(3, "Set1"),"black")) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw() + 
  theme(axis.title.y = element_blank(), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        strip.background = element_blank())

l <- cowplot::get_plot_component(p2, 'guide-box-bottom', return_all = TRUE)

plot_grid(
  plot_grid(p1 + theme(legend.position = "none"), 
            p2 + theme(legend.position = "none"), 
            rel_widths = c(0.45, 0.45), nrow = 1, labels = c("A", "B")), 
  l, rel_heights = c(0.9, 0.1), ncol = 1)



ggsave("figures/vaccddifferences_SEIRS_scotland_Ipred.pdf", width = 8, height = 10)
  
