library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

theme_set(theme_bw(base_size = 11))

#### LOAD DATA AND SETUP -------------------------------------------------------
scotland_by_wk <- read_rds("data/derived_data/scotland/scotland_by_wk_overall.rds")
canada_by_wk <- read_rds("data/derived_data/canada/canada_by_wk_overall.rds")
load("/Users/emilyhowerton/Documents/GitHub/seasonal-resp-viruses/timeseries-data/data/data_korea_ari.rda")

data_all_locs <- bind_rows(
  scotland_by_wk %>%
    filter(wk_collected > as.Date("2006-01-01")) %>%
    select(wk_collected, pathogen, detections) %>%
    mutate(region = "Scotland", 
           wk_collected = ceiling_date(wk_collected, "week")), 
  canada_by_wk %>%
    filter(!(region %in% c("CANADA", "Territories"))) %>%
    select(region, wk_collected, pos_rsv, pos_mpv) %>%
    melt(c('region', "wk_collected"), value.name = "detections") %>%
    separate(variable, into = c("variable", "pathogen")) %>%
    select(-variable) %>%
    mutate(region = paste0(ifelse(substr(region,1,3) == "Pro", region, paste(region, "Region")), ", Canada"), 
           wk_collected = ceiling_date(wk_collected, "week")), 
  data_korea_ari %>%
    filter(key %in% c("Human metapneumovirus", "RSV")) %>%
    mutate(region = "Korea", 
           pathogen = ifelse(key == "RSV", "rsv", "mpv"),
           wk_collected = ceiling_date(as.Date(paste(year, week, 1, sep="-"), "%Y-%U-%u"), "week")) %>%
    rename(detections = cases) %>%
    select(region, wk_collected, pathogen, detections)
)

biennial_regions <- c("Province of British Columbia, Canada", "Prairies Region, Canada")

region_lvls <- c(
  "Scotland",  "Atlantic Region, Canada","Province of British Columbia, Canada",
  "Prairies Region, Canada", "Province of Ontario, Canada", "Province of Québec, Canada", "Korea"
)


#### CALCULATE CROSS CORRELATIONS ----------------------------------------------
rgs = unique(data_all_locs$region)
lag_corr_df <- list(length(rgs))
for(i in 1:length(rgs)){
  tmp_mpv = data_all_locs %>% 
    filter(region == rgs[i], pathogen == "mpv", wk_collected <= "2020-03-01") %>% 
    pull(detections)
  tmp_rsv = data_all_locs %>% 
    filter(region == rgs[i], pathogen == "rsv", wk_collected <= "2020-03-01") %>% 
    pull(detections)
  lag_corr <- ccf(x = tmp_rsv, y = tmp_mpv, na.action = na.pass)
  lag_corr_df[[i]] <- data.frame(lag = lag_corr$lag[,,1], 
                                 n.used = lag_corr$n.used, 
                                 corr = unlist(lag_corr$acf[,,1]), 
                                 region = rgs[i])
}

max_corr <- bind_rows(lag_corr_df) %>%
  mutate(max_corr = max(corr), .by = c("region")) %>%
  filter(corr == max_corr) %>%
  select(-max_corr) %>%
  mutate(lwr = tanh(atanh(corr) + qnorm(p=0.025)*(1/sqrt(n.used-3))),
         upr = tanh(atanh(corr) + qnorm(p=0.975)*(1/sqrt(n.used-3))), 
         biennial_flag = ifelse(region %in% biennial_regions, TRUE, FALSE))

# for supplement
bind_rows(lag_corr_df) %>%
  filter(!(region %in% biennial_regions)) %>%
  mutate(upper = qnorm((1 + 0.95)/2)/sqrt(n.used), 
         lower = -qnorm((1 + 0.95)/2)/sqrt(n.used),
         sig = ifelse(corr > upper | corr < lower, TRUE, FALSE), 
         region = factor(region, levels = c(region_lvls[2:6], region_lvls[c(1,7)]))) %>%
  ggplot(aes(x = lag, y = corr)) +
  geom_hline(aes(yintercept = 0), linewidth = 0.3) +
  geom_hline(aes(yintercept = lower), linetype = "dashed", linewidth = 0.3) +
  geom_hline(aes(yintercept = upper), linetype = "dashed", linewidth = 0.3) +
  geom_segment(aes(x = lag, xend = lag, y = 0, yend = corr), linewidth = 0.5) + 
  facet_wrap(vars(region)) + 
  labs(x = "RSV-HPMV lag (weeks)", y = "correlation") + 
  theme_bw() + 
  theme(legend.position = "none", 
        panel.grid = element_blank(), 
        strip.background = element_blank())
ggsave("figures/cross-correlation_allregions.pdf", width = 7, height = 4)


####  PLOT ALL TIMESERIES ------------------------------------------------------

lag_text <- max_corr %>%
  mutate(x = as.Date(ifelse(region == "Scotland", "2007-01-01", 
                            ifelse(region == "Korea", "2023-07-01", "2014-01-01"))), 
         y = Inf, 
         txt = ifelse(region %in% biennial_regions, "biennial", 
                      paste0("annual\n", -lag, " week lag")), 
         hjust = ifelse(region == "Korea", 1, 0)) %>%
  mutate(region = factor(region, levels = region_lvls))

txt_size = 3

p1_scotland = data_all_locs %>%
  filter(detections > 0, 
         region %in% c("Scotland")) %>%
  mutate(region = factor(region, levels = region_lvls), 
         pathogen = factor(pathogen, levels = c("rsv", "mpv"))) %>%
  ggplot(aes(x = wk_collected, y = detections, color = pathogen)) + 
  geom_line() + 
  geom_text(data = lag_text %>% filter(region %in% c("Scotland")),
            aes(x = x, y = y, label = txt, hjust = hjust), color = "black",  size = txt_size, vjust = 1) + 
  facet_wrap(vars(region), nrow = 1, scales = "free") +
  scale_color_brewer(palette = "Dark2", labels = c("RSV", "hMPV")) + 
  scale_x_date(expand = c(0,0)) +
  scale_y_sqrt(expand = c(0,0)) + 
  theme(axis.title.x = element_blank(), 
        legend.margin = margin(rep(0.1,4)),
        legend.position = c(0.91, 0.79),
        legend.title = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank())
p1_scotland

p1_biennial = data_all_locs %>%
  filter(detections > 0, 
         region %in% c("Atlantic Region, Canada", "Province of Ontario, Canada")) %>%
  mutate(region = factor(region, levels = region_lvls), 
         pathogen = factor(pathogen, levels = c("rsv", "mpv"))) %>%
  ggplot(aes(x = wk_collected, y = detections, color = pathogen)) + 
  geom_line() + 
  geom_text(data = lag_text %>% filter(region %in% c("Atlantic Region, Canada", "Province of Ontario, Canada")),
            aes(x = x, y = y, label = txt, hjust = hjust), color = "black",  size = txt_size, vjust = 1) + 
  # geom_point(size = 0.6) + 
  facet_wrap(vars(region), nrow = 1, scales = "free") + #, labeller = labeller(region = region_labs)
  scale_color_brewer(palette = "Dark2", labels = c("RSV", "hMPV")) + 
  scale_x_date(expand = c(0,0)) +
  scale_y_sqrt(expand = c(0,0)) + 
  theme(axis.title.x = element_blank(), 
        legend.position = "none", 
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

p1_rest <- data_all_locs %>%
  filter(detections > 0, 
         region %in% c("Province of British Columbia, Canada", "Prairies Region, Canada",
                        "Province of Québec, Canada", "Korea")) %>%
  mutate(region = factor(region, levels = region_lvls), 
         pathogen = factor(pathogen, levels = c("rsv", "mpv"))) %>%
  ggplot(aes(x = wk_collected, y = detections, color = pathogen)) + 
  geom_line() + 
  geom_text(data = lag_text %>% filter(!(region %in% c("Scotland", "Atlantic Region, Canada", "Province of Ontario, Canada"))),
            aes(x = x, y = y, label = txt, hjust = hjust), color = "black", size = txt_size, vjust = 1) + 
  facet_wrap(vars(region), nrow = 2, scales = "free") +
  scale_color_brewer(palette = "Dark2", labels = c("RSV", "hMPV")) + 
  scale_x_date(expand = c(0,0)) +
  scale_y_sqrt(expand = c(0,0)) + 
  theme(axis.title.x = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

p1 = plot_grid(
  plot_grid(p1_scotland, p1_biennial, ncol = 1), 
  p1_rest + theme(axis.title.y = element_blank(), legend.position = "none"), 
  nrow = 1
)

ggsave("figures/all-timeseries.pdf", width = 10, height = 4)


#### PLOT TESTING TRENDS -------------------------------------------------------
path_labs = c("RSV", "hMPV")
names(path_labs) = c("rsv", "mpv")

scotland_by_wk %>% 
  filter(wk_collected > start_wk_pre) %>%
  mutate(inc_flag = ifelse(n_tests == 0, "F", "T"), 
         pathogen = factor(pathogen, levels = c("rsv", "mpv"))) %>%
  ggplot(aes(x = wk_collected)) + 
  geom_point(aes(y = n_tests), shape = 21) + 
  geom_line(aes(y = n_tests_ma)) + 
  facet_wrap(vars(pathogen), scales = "free", ncol = 1, labeller = labeller(pathogen = path_labs)) + 
  scale_x_date(expand = c(0,0)) +
  scale_y_continuous(name = "number of tests") + 
  theme(axis.title.x = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave("figures/testing-timeseries.pdf", width = 10, height = 6)
