library(ggplot2)
library(cowplot)
library(gridExtra)
library(readr)
library(tidyr)
library(scales)
library(lubridate)

source("R/scotland_parameters.R")

bifur_points <- read_rds("data/derived_data/scotland/bifur-points_scotland_SEIRS.rds")
bifur_sims <- read_rds("data/derived_data/scotland/bifur-sims_scotland_SEIRS.rds")

fit_scotland_SEIRS <- read_rds("data/derived_data/scotland/fit_scotland_SEIRS.rds")
fitted_pars = apply(as.array(fit_scotland_SEIRS)[,,1:14], 3, median)

c_slices = sort(c(c_slices, fitted_pars["c"]))
a_range = sort(c(a_range, fitted_pars["a"]))

a_ex = as.vector(sort(c(0.05, fitted_pars["a"], 0.8)))

c_names = c("c = 0 (vaccination)", paste("c =", round(fitted_pars["c"],2) , "(fitted)"))


find_seasonality <- function(pts, thrsh = 0.0001){
    pts_sort = sort(pts)
    d = diff(pts_sort)
    l = length(d[d>thrsh])
    return(ifelse(l == 0, "annual", ifelse(l == 1, "biennial", "3-year +")))
}

 ### PLOT BIFURCATION FOR I ----------------------------------------------------

seas = bind_rows(bifur_points) %>%
  filter(variable == "I") %>%
  filter(c == fitted_pars["c"]) %>%
  summarize(seas = find_seasonality(value/scotland_N, thrsh = 0.0005), .by = c("a", "pathogen", "draw_id")) %>%
  mutate(seas = factor(seas, levels = c("annual", "biennial", "3-year +")))

rsv_region_yval = max(bind_rows(bifur_points) %>% filter(c == 0) %>% 
                        filter(variable == "I", pathogen == "rsv") %>% 
                        pull(value))/scotland_N+0.01
mpv_region_yval = max(bind_rows(bifur_points) %>% filter(c == fitted_pars["c"]) %>% 
                        filter(variable == "I", pathogen == "mpv") %>% 
                        pull(value))/scotland_N+0.01

# plot bifurcation for RSV
p1 = bind_rows(bifur_points) %>%
  filter(variable == "I", pathogen == "rsv") %>%
  filter(c == 0) %>%
  ggplot(aes(x = a, y = value/scotland_N)) + 
  geom_point() + 
  geom_vline(xintercept = a_ex, linetype = "dotted", size = 0.8) + 
  geom_point(data = seas %>% filter(pathogen == "rsv") %>% mutate(y = rsv_region_yval),
             aes(x = a, y = y, color = as.factor(seas)), size = 8, shape = "|") +
  labs(x = "amplitude of seasonal forcing", y = "RSV prevalence") + 
  scale_color_brewer(palette = "Purples", direction = -1) + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_bw() + 
  theme(legend.position = "none", 
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

# plot bifurcation for HMPV
p2 = bind_rows(bifur_points) %>%
  filter(variable == "I", pathogen == "mpv") %>%
  filter(c == 0 ) %>%
  mutate(c_name = factor(ifelse(c == 0, c_names[1], c_names[2]), levels = c_names)) %>%
  arrange(c_name) %>%
  ggplot(aes(x = a, y = value/scotland_N)) + 
  geom_point(aes(color = c_names[1])) +
  geom_point(data = bind_rows(bifur_points) %>%
               filter(variable == "I", pathogen == "mpv") %>%
               select(a, c, t, value) %>%
               melt(c("a", "c", "t")) %>%
               filter(c == fitted_pars["c"]) %>%
               mutate(c_name = factor(ifelse(c == 0, c_names[1], c_names[2]), levels = c_names)) %>%
               arrange(c_name), aes(color = c_names[2])) +
  geom_vline(xintercept = a_ex, linetype = "dotted", size = 0.8) + 
  geom_point(data = seas %>% filter(pathogen == "mpv") %>% mutate(y = mpv_region_yval),
             aes(x = a, y = y, color = as.factor(seas)), size = 8, shape = "|") +
  guides(color=guide_legend(nrow = 1)) + 
  labs(x = "amplitude of seasonal forcing", 
       y = "HMPV prevalence") + 
  scale_color_manual(values = c(rev(RColorBrewer::brewer.pal(3, "Purples")),"black", "gray"))+
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

# plot example timeseries
t_to_wk_scotland_bifur <- data.frame(wk_collected = as.Date(start_wk_pre + (0:nrow(bifur_sims[[1]][[1]]))*7),
                                     t = 1:(nrow(bifur_sims[[1]][[1]])+1))

a_labs = paste("a =", round(a_ex,2))
a_labs[which(a_ex == fitted_pars["a"])] = paste(a_labs[which(a_ex == fitted_pars["a"])], "(fitted)")
names(a_labs) = a_ex

aexplts_I <- bind_rows(
  bifur_sims[[which(c_slices == fitted_pars["c"])]][sapply(a_ex, function(i){which(round(a_range,5) == round(i,5))})], 
  .id = "draw_id") %>%
    mutate(c =fitted_pars["c"] ) %>%
    bind_rows(
      bind_rows(
        bifur_sims[[which(c_slices == 0)]][sapply(a_ex, function(i){which(round(a_range,5) == round(i,5))})], 
        .id = "draw_id") %>%
        mutate(c = 0)) %>%
  mutate(a = a_ex[as.integer(draw_id)]) %>%
  filter(variable == "Ipred") %>%
  left_join(t_to_wk_scotland_bifur) %>%
  filter(year(wk_collected) >= 2036, year(wk_collected) < 2046) %>%
  mutate(c_name = ifelse(c == 0, "c = 0 (vaccination)", paste("c =", round(c,2), "(fitted)"))) %>%
  ggplot(aes(x = wk_collected, y = value, color = pathogen)) + 
  geom_line(aes(linetype = c_name)) + 
  facet_wrap(vars(a), labeller = labeller(a = a_labs), ncol = 1, scales = "free") + 
  scale_color_brewer(palette = "Dark2", 
                     direction = -1,
                     labels = c("HMPV", "RSV"),
                     name = "pathogen") + 
  scale_linetype_manual(values = c("solid", "dashed"), 
                        name = "cross protection scenario") + 
  scale_x_date(date_breaks = "year") +
  scale_y_continuous(labels = comma, name = "detections") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.minor = element_blank(), 
        strip.background = element_blank())

l <- cowplot::get_plot_component(p2, 'guide-box-bottom', return_all = TRUE)

plot_grid(
  plot_grid(p1 + theme(legend.position = "none"), 
            p2 + theme(legend.position = "none"), 
            l, rel_heights = c(0.45, 0.45, 0.1), ncol = 1, labels = c("A", "B")), 
  aexplts_I, labels = c(NA, "C"))

ggsave("figures/bifurcation_SEIRS_scotland.pdf", width = 14, height = 8)
