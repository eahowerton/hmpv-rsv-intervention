library(ggplot2)
library(cowplot)
library(gridExtra)
library(readr)
library(tidyr)
library(scales)

source("R/2_run_simulations/simulate-postpand_scotland_SEIRS.R")
source("R/2_run_simulations/simulate-postpand_scotland_SEIRS-ind.R")


c = postpand_sims_ind %>%
  filter(variable == "Ipred") %>%
  left_join(t_to_wk_scotland_full) %>%
  left_join(scotland_by_wk %>%
              filter(wk_collected > start_wk_pre) %>%
              select(wk_collected, pathogen, detections)) %>%
  mutate(model = "ind") %>%
  bind_rows(postpand_sims %>%
              filter(variable == "Ipred") %>%
              left_join(t_to_wk_scotland_full) %>%
              left_join(scotland_by_wk %>%
                          filter(wk_collected > start_wk_pre) %>%
                          select(wk_collected, pathogen, detections)) %>%
              mutate(model = "coupled")) 

c_all = c %>%
  filter(!is.na(detections), wk_collected >= start_wk_post) %>% #wk_collected >= as.Date("2020-03-01")
  summarize(c = cor(log(value+1), log(detections+1)), .by = c("draw_id", "model", "pathogen"))

c_pre = c %>%
  filter(!is.na(detections), wk_collected >= start_wk_post,wk_collected < as.Date("2020-03-01")) %>% #wk_collected >= as.Date("2020-03-01")
  summarize(c = cor(log(value+1), log(detections+1)), .by = c("draw_id", "model", "pathogen"))

c_post = c %>%
  filter(!is.na(detections), wk_collected >= as.Date("2021-01-01")) %>% 
  summarize(c = cor(log(value+1), log(detections+1)), .by = c("draw_id", "model", "pathogen"))


## nll 
nll_all <- c %>% 
  filter(wk_collected >= start_wk_post, !is.na(detections)) %>%
  mutate(ll = dpois(x = detections, lambda = value, log = TRUE)) %>% # calculate ll for each
  summarize(nll = -sum(ll), .by = c("draw_id", "model", "pathogen"))

nll_pre <- c %>% 
  filter(wk_collected >= start_wk_post,wk_collected < as.Date("2020-03-01"), !is.na(detections)) %>%
  mutate(ll = dpois(x = detections, lambda = value, log = TRUE)) %>% # calculate ll for each
  summarize(nll = -sum(ll), .by = c("draw_id", "model", "pathogen"))

nll_post <- c %>% 
  filter(wk_collected >= as.Date("2021-01-01"), !is.na(detections)) %>%
  mutate(ll = dpois(x = detections, lambda = value, log = TRUE)) %>% # calculate ll for each
  summarize(nll = -sum(ll), .by = c("draw_id", "model", "pathogen"))


lbs = c("correlation between posterior prediction and observation", 
        "log likelihood")
names(lbs) = c("c", "ll")


# try adding pathogen
p1 = bind_rows(nll_all %>% mutate(period = "all"), 
               nll_pre %>% mutate(period = "pre"),
               nll_post %>% mutate(period = "post")) %>%
  mutate(period = factor(period, levels = c("all", "pre", "post")), 
         model = factor(model, levels = rev(c("coupled", "ind"))),
         pathogen = factor(pathogen, levels = c("rsv", "mpv")),
         ll = -nll) %>% 
  select(-nll) %>%
  summarize(Q5 = quantile(ll, 0.05), 
            Q25 = quantile(ll, 0.25), 
            Q50 = median(ll), 
            Q75 = quantile(ll, 0.75), 
            Q95 = quantile(ll, 0.95), 
            .by = c("model", "period", "pathogen")) %>%
  mutate(y = -as.numeric(period) + ifelse(pathogen == "rsv", 0.1,-0.1) + ifelse(model == "coupled", 0.05, -0.05)) %>%
  ggplot(aes(color = pathogen)) + 
  geom_segment(aes(x = Q5, xend = Q95, y = y, yend = y, group = model)) +
  geom_segment(aes(x = Q25, xend = Q75, y = y, yend = y, group = model), linewidth = 1.2) + 
  geom_point(aes(x = Q50, y = y, shape = model), size = 1.8) + 
  scale_color_brewer(palette = "Dark2", labels = c("RSV", "hMPV")) + 
  scale_shape_manual(labels = c("independent model", "coupled model"), 
                     values = c(3,16)) +
  scale_x_continuous(name = "out of sample log likelihood") + 
  scale_y_continuous(breaks = -c(1,2,3), 
                     labels = c("full period\n(2018-2024)", 
                                "pre-pandemic\n(2018-2020)", 
                                "post-pandemic\n(2021-2024)")) +
  theme_bw() + 
  theme(axis.title.y = element_blank(), 
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        strip.background = element_blank(),
        strip.placement = "outside")


p2 = bind_rows(c_all %>% mutate(period = "all"), 
               c_pre %>% mutate(period = "pre"),
               c_post %>% mutate(period = "post")) %>%
  mutate(period = factor(period, levels = c("all", "pre", "post")), 
         model = factor(model, levels = rev(c("coupled", "ind"))),
         pathogen = factor(pathogen, levels = c("rsv", "mpv"))) %>% 
  summarize(Q5 = quantile(c, 0.05), 
            Q25 = quantile(c, 0.25), 
            Q50 = median(c), 
            Q75 = quantile(c, 0.75), 
            Q95 = quantile(c, 0.95), 
            .by = c("model", "period", "pathogen")) %>%
  mutate(y = -as.numeric(period) + ifelse(pathogen == "rsv", 0.1,-0.1) + ifelse(model == "coupled", 0.05, -0.05)) %>%
  ggplot(aes(color = pathogen)) + 
  geom_segment(aes(x = Q5, xend = Q95, y = y, yend = y)) +
  geom_segment(aes(x = Q25, xend = Q75, y = y, yend = y), linewidth = 1.2) + 
  geom_point(aes(x = Q50, y = y, shape = model), size = 1.8) + 
  scale_color_brewer(palette = "Dark2", labels = c("RSV", "hMPV")) + 
  scale_shape_manual(labels = c("independent model", "coupled model"), 
                       values = c(3,16)) +
  scale_x_continuous(name = "out of sample correlation") + 
  scale_y_continuous(breaks = -c(1,2,3), 
                     labels = c("full period\n(2018-2024)", 
                                "pre-pandemic\n(2018-2020)", 
                                "post-pandemic\n(2021-2024)")) +
  theme_bw() + 
  theme(axis.title.y = element_blank(), 
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        strip.background = element_blank(),
        strip.placement = "outside")

l <- cowplot::get_plot_component(p1, 'guide-box-bottom', return_all = TRUE)

plot_grid(
  plot_grid(p1 + theme(legend.position = "none"), 
            p2 + theme(legend.position = "none"), 
            nrow = 1, labels = LETTERS[1:2]), 
  l, ncol = 1, rel_heights = c(0.9, 0.1)
)

ggsave("figures/model_comparison.pdf", width = 8, height = 3)




