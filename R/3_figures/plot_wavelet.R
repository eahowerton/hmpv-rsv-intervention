library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(biwavelet)
library(lubridate)
library(reshape2)
library(scales)


#### FUNCTIONS -----------------------------------------------------------------
implement_wavelet <- function(detections_rsv, detections_mpv, 
                              sqrt_norm_transform = TRUE){
  if(sqrt_norm_transform){
    # square root transform and normalize
    detections_rsv_sqrt <- sqrt(detections_rsv)
    detections_rsv_trans <- (detections_rsv_sqrt - mean(detections_rsv_sqrt))/
      sd(detections_rsv_sqrt)
    detections_mpv_sqrt <-sqrt(detections_mpv)
    detections_mpv_trans <- (detections_mpv_sqrt - mean(detections_mpv_sqrt))/
      sd(detections_mpv_sqrt)
  }
  else{
    detections_rsv_trans <- detections_rsv
    detections_mpv_trans <- detections_mpv
  }
  # prepare data for wavelet transform
  cases_rsv <- cbind(1:length(detections_rsv), detections_rsv_trans)
  cases_mpv <- cbind(1:length(detections_mpv), detections_mpv_trans)
  # perform cross wavelet transform (XWT)
  xwt_result <- xwt(cases_rsv, cases_mpv)
  # perform wavelet coherence analysis (WTC)
  wtc_result <- wtc(cases_rsv, cases_mpv)
  return(list(xwt = xwt_result, wtc = wtc_result))
}

plot_wavelet <- function(wavelet_result, title, dates, ...){
  # customize x-axis ticks
  start_date = min(dates)
  end_date = max(dates)
  start_year <- as.numeric(format(start_date, "%Y")) + 1
  end_year <- as.numeric(format(end_date, "%Y"))
  custom_ticks <- seq(start_year, end_year, 1)
  # place tick on first week of each year
  tick_weeks <- sapply(custom_ticks, function(i){min(which(lubridate::year(dates) == i))}) 
  plot(wavelet_result, plot.coi = TRUE, plot.phase = TRUE, # plot.cb = TRUE,
       xaxt = "n", xlab = "", ylab = "period (weeks)", ...)
  axis(1, at = tick_weeks, labels = custom_ticks)
  mtext(title, side = 3, line = 0, at = 0, adj = 0)
}

#### LOAD DATA -----------------------------------------------------------------

scotland_by_wk <- readRDS("data/derived_data/scotland/scotland_by_wk_overall.rds")
scotland_by_wk <- read_rds("data/derived_data/scotland/scotland_by_wk_overall.rds")
canada_by_wk <- read_rds("data/derived_data/canada/canada_by_wk_overall.rds")
load("/Users/emilyhowerton/Documents/GitHub/seasonal-resp-viruses/timeseries-data/data/data_korea_ari.rda")

## prep scotland data
source("R/scotland_parameters.R")
scotland_by_wk_pre <- scotland_by_wk %>%
  filter(wk_collected >= start_wk_pre, 
         wk_collected <= as.Date("2020-03-01")) 

## prep canada data 
canada_by_wk_pre <- canada_by_wk %>%
  filter(!(region %in% c("CANADA", "Territories")), 
         wk_collected <= as.Date("2020-03-01")) %>%
  select(region, wk_collected, pos_rsv, pos_mpv) %>%
  melt(c('region', "wk_collected"), value.name = "detections") %>%
  separate(variable, into = c("variable", "pathogen")) %>%
  select(-variable) %>%
  mutate(region = paste0(ifelse(substr(region,1,3) == "Pro", region, paste(region, "Region")), ", Canada"), 
         wk_collected = ceiling_date(wk_collected, "week"), 
         detections = ifelse(is.na(detections), 0, detections))

## prep korea data
korea_by_wk_pre <- data_korea_ari %>%
  filter(key %in% c("Human metapneumovirus", "RSV")) %>%
  mutate(region = "Korea", 
         pathogen = ifelse(key == "RSV", "rsv", "mpv"),
         wk_collected = ceiling_date(as.Date(paste(year, week, 1, sep="-"), "%Y-%U-%u"), "week")) %>%
  rename(detections = cases) %>%
  select(region, wk_collected, pathogen, detections) %>%
  filter(!is.na(wk_collected), !is.na(detections), wk_collected <= as.Date("2020-03-01"))

#### IMPLEMENT WAVELET ANALYSIS ------------------------------------------------
wavelet_scotland <- implement_wavelet(
  detections_rsv = scotland_by_wk_pre %>% filter(pathogen == "rsv") %>% pull(detections), 
  detections_mpv = scotland_by_wk_pre %>% filter(pathogen == "mpv") %>% pull(detections), 
  sqrt_norm_transform = TRUE)

wavelet_canada_bc <- implement_wavelet(
  detections_rsv = canada_by_wk_pre %>% 
    filter(pathogen == "rsv", region == "Province of British Columbia, Canada") %>% 
    pull(detections), 
  detections_mpv = canada_by_wk_pre %>% 
    filter(pathogen == "mpv", region == "Province of British Columbia, Canada") %>% 
    pull(detections), 
  sqrt_norm_transform = TRUE)

wavelet_canada_ontario <- implement_wavelet(
  detections_rsv = canada_by_wk_pre %>% 
    filter(pathogen == "rsv", region == "Province of Ontario, Canada") %>% 
    pull(detections), 
  detections_mpv = canada_by_wk_pre %>% 
    filter(pathogen == "mpv", region == "Province of Ontario, Canada") %>% 
    pull(detections), 
  sqrt_norm_transform = TRUE)

wavelet_canada_atlantic <- implement_wavelet(
  detections_rsv = canada_by_wk_pre %>% 
    filter(pathogen == "rsv", region == "Atlantic Region, Canada") %>% 
    pull(detections), 
  detections_mpv = canada_by_wk_pre %>% 
    filter(pathogen == "mpv", region == "Atlantic Region, Canada") %>% 
    pull(detections), 
  sqrt_norm_transform = TRUE)

wavelet_canada_quebec <- implement_wavelet(
  detections_rsv = canada_by_wk_pre %>% 
    filter(pathogen == "rsv", region == "Province of Québec, Canada") %>% 
    pull(detections), 
  detections_mpv = canada_by_wk_pre %>% 
    filter(pathogen == "mpv", region == "Province of Québec, Canada") %>% 
    pull(detections), 
  sqrt_norm_transform = TRUE)

wavelet_canada_praries <- implement_wavelet(
  detections_rsv = canada_by_wk_pre %>% 
    filter(pathogen == "rsv", region == "Prairies Region, Canada") %>% 
    pull(detections), 
  detections_mpv = canada_by_wk_pre %>% 
    filter(pathogen == "mpv", region == "Prairies Region, Canada") %>% 
    pull(detections), 
  sqrt_norm_transform = TRUE)

wavelet_korea <- implement_wavelet(
  detections_rsv = korea_by_wk_pre %>% 
    filter(pathogen == "rsv") %>% 
    pull(detections), 
  detections_mpv = korea_by_wk_pre %>% 
    filter(pathogen == "mpv") %>% 
    pull(detections), 
  sqrt_norm_transform = TRUE)

#### PLOT RESULTS --------------------------------------------------------------

## standardize colors and periods across all plots
wavelet_list <- list(wavelet_canada_bc$xwt, wavelet_canada_praries$xwt, 
                     wavelet_canada_atlantic$xwt, wavelet_canada_ontario$xwt,
                     wavelet_canada_quebec$xwt, wavelet_scotland$xwt, wavelet_korea$xwt)

# zlim chosen as the largest range of values across all wavelet analyses
# where zlim = range(c(-1, 1) * max(log2(x$power.corr)/(x$d1.sigma * x$d2.sigma)))
maxzval_all <- lapply(wavelet_list, function(i){max(log2(i$power.corr)/(i$d1.sigma * i$d2.sigma))})
zlim =  range(c(-1, 1) * max(unlist(maxzval_all)))

# find ylim across all wavelets (exclude Korea, because time series is shorter)
# where ylim = range(log2(x$period)) in plot.biwavelet()
ylim_all <- lapply(wavelet_list[-length(wavelet_list)], 
                   function(i){range(i$period)})
# choose smallest possible bounds
ylim = rev(c(max(unlist(lapply(ylim_all, function(i){i[1]}))), 
         min(unlist(lapply(ylim_all, function(i){i[2]})))))

# pdf("figures/cross-wavelet.pdf", width = 15, height = 6)
pdf("figures/cross-wavelet.pdf", width = 8, height = 10)
par(mfrow = c(4,2), mar = c(3, 4, 2, 1), oma = c(0, 0, 0, 0))
plot_wavelet(wavelet_scotland$xwt, title = "Scotland", 
             dates = sort(unique(scotland_by_wk_pre$wk_collected)), 
             ylim = rev(ylim), zlim = zlim, lwd.sig = 2)
plot_wavelet(wavelet_korea$xwt, title = "Korea", 
             dates = unique(korea_by_wk_pre$wk_collected), 
             ylim = rev(ylim), zlim = zlim, lwd.sig = 2)
plot_wavelet(wavelet_canada_bc$xwt, title = "British Columbia, Canada", 
             dates = unique(canada_by_wk_pre$wk_collected), 
             ylim = rev(ylim), zlim = zlim, lwd.sig = 2)
plot_wavelet(wavelet_canada_praries$xwt, title = "Prairies Region, Canada", 
             dates = unique(canada_by_wk_pre$wk_collected), 
             ylim = rev(ylim), zlim = zlim, lwd.sig = 2)
plot_wavelet(wavelet_canada_atlantic$xwt, title = "Atlantic Region, Canada", 
             dates = unique(canada_by_wk_pre$wk_collected), 
             ylim = rev(ylim), zlim = zlim, lwd.sig = 2)
plot_wavelet(wavelet_canada_ontario$xwt, title = "Ontario, Canada", 
             dates = unique(canada_by_wk_pre$wk_collected), 
             ylim = rev(ylim), zlim = zlim, lwd.sig = 2)
plot_wavelet(wavelet_canada_quebec$xwt, title = "Quebec, Canada", 
             dates = unique(canada_by_wk_pre$wk_collected), 
             ylim = rev(ylim), zlim = zlim, lwd.sig = 2)
dev.off()

