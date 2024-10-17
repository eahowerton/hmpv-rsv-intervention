library(dplyr)
library(readr)

# read in scotland data, all age groups, summarized by week
scotland_by_wk <- read_rds("data/derived_data/scotland/scotland_by_wk_overall.rds")

# parameters for model fit
start_wk_pre <- as.Date("2006-10-09") #min(scotland_by_wk$wk_collected)
end_wk_pre <- as.Date("2018-06-18") #2018-06-18

# parameters for post-pandemic simulations
start_wk_post <- as.Date("2018-06-25")
end_wk_post <- as.Date("2024-03-04")

# demographic parameters
scotland_N <- 5.436e6 # 2022 population 
scotland_sigma <-  1/(1/7)
scotland_gamma <- 1/(6/7)
scotland_omega <- 0.51/52 # (from White et al.)
scotland_birthrate <- mean(scotland_by_wk %>% 
                             filter(wk_collected >= start_wk_pre, 
                                    wk_collected <= end_wk_pre) %>%
                             pull(Births.registered))/scotland_N

# data frames to convert simulation weeks to calendar weeks
t_to_wk_scotland <- data.frame(wk_collected = seq.Date(start_wk_pre, end_wk_pre, by = 7),
                               t = 1:length(seq.Date(start_wk_pre, end_wk_pre, by = 7)))

t_to_wk_scotland_post <- data.frame(wk_collected = seq.Date(start_wk_post, end_wk_post, by = 7),
                                    t = 1:length(seq.Date(start_wk_post, end_wk_post, by = 7)))

t_to_wk_scotland_full <- data.frame(wk_collected = seq.Date(start_wk_pre, end_wk_post, by = 7),
                                   t = 1:length(seq.Date(start_wk_pre, end_wk_post, by = 7)))

