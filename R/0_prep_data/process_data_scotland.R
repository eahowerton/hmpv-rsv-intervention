# process scotland data
library(readxl)
library(dplyr)
library(reshape2)
library(purrr)
library(tseries)
library(zoo)
library(readr)

#### UTILS ------------------------------------------------------------

# labels
labs_health_board <- c(
  "Ayrshire and\nArran", "Argyll and Bute", "Borders", "Dumfries and\nGalloway", 
  "Eire", "England", "Fife", "Forensic Medicine Glasgow", "Forth valley", 
  "Greater Glasgow\nand Clyde ", "Golden Jubilee", "Grampian", "Highland", 
  "Lanarkshire", "Lothian", "Tayside", "Western Isles", "Armed Services"
)
names(labs_health_board) <- c(
  "AA", "ABC", "BOR", "DG", "EIRE", "ENG", "FIFE", "FM", "FV", "GG", "GOLDJ",
  "GR", "HIG", "LAN", "LOT", "TAY", "WI", "SERVIC")

labs_pathogen <- c("RSV", "HMPV")
names(labs_pathogen) <- c("rsv", "mpv")

# health boards with meaningful numbers of cases
major_health_boards <- c("AA", "DG", "FV", "GG", "LAN")

#### COLLECT AND TIDY DATA -----------------------------------------------------

scotland_data1 <- read_excel("data/raw_data/scotland/rsvhump_011006-040814_RG.xlsx", 
                             sheet = "rsvhump_011023-04082014")
scotland_data2 <- read_excel("data/raw_data/scotland/Amended_rsvhump_040814-280224_RGv2.xlsx", 
                             sheet = "rsvhump_040814-280224")

# combine files
scotland_data <- bind_rows(
  scotland_data1 %>%
    select(-"...13") %>%
    rename(all_of(c(lab_spec_id = "Lab No/Spec No", 
                    age = "Age", 
                    sex = "Sex", 
                    date_collected = "Date collectedC", 
                    date_received = "Date Received", 
                    loc_code = "Location code", 
                    loc_expand = "Location expansion", 
                    health_board = "Health board", 
                    result_rsv = "RSV PCR", 
                    ct_rsv = "ct rsv", 
                    result_mpv = "HUMP", 
                    ct_mpv = "Ct"))), 
  scotland_data2 %>%
    rename(all_of(c(lab_spec_id = "Lab No/Spec No", 
                    age = "Age", 
                    sex = "Sex", 
                    date_collected = "DC", 
                    date_received = "DR", 
                    loc_code = "LOC", 
                    loc_expand = "Location", 
                    health_board = "DRLG", 
                    result_rsv = "RSV PCR", 
                    ct_rsv = "Ct RSV", 
                    result_mpv = "HUMPCR", 
                    ct_mpv = "Ct."))) 
)

# standardize result codes
# ****NOTE: decide what to do with "DL" class (currently counting this
#           as a positive, but might deserve more nuance)
rsv_result_codes <- data.frame(
  original_code_rsv = c("A", "B", "E", "N", "P", 
                        "D", "DL", "I", "IN", "N", "NA", "RD"), 
  original_description_rsv = c(
    "RSV subtype A", "RSV subtype B", "Equivocal", 
    "Not detected by pcr","Detected by pcr", 
    "Detected by pcr", "Detected at low levels by pcr", "Insufficient",
    "Inhibited", "Not detected by pcr", "Not available", "Request declined"), 
  file = c(rep(2, 5), rep(1,7)), 
  new_code_rsv = c(
    "D", "D", "I", "ND", "D", 
    "D", "D", "I", "I", "ND", "NA", "NA"), 
  new_description_rsv = c(
    "detected", "detected", "insufficient", "not detected", "detected", 
    "detected", "detected", "insufficient", "insufficient", "not detected",
    "not tested", "not tested")
)

mpv_result_codes <- data.frame(
  original_code_mpv = c(
    "E", "Insuff", "LP", "N", "P", 
    "D", "DL", "I","IN", "N", "RD"), 
  original_description_mpv = c(
    "Equivocal", "Insufficient", "Detected at low levels by pcr",
    "Not detected by pcr", "Detected by pcr", 
    "Detected by pcr", "Detected at low levels by pcr", 
    "Insufficient", "Inhibited", "Not detected by pcr", "Request declined"), 
  file = c(rep(2,5), rep(1, 6)), 
  new_code_mpv = c(
    "I", "I", "D", "ND", "D", "D", "D", "I", "I", "ND", "NA"), 
  new_description_mpv = c(
    "insufficient", "insufficient", "detected", "not detected", 
    "detected", "detected", "detected", "insufficient", "insufficient", 
    "not detected", "not tested")
)

scotland_data <- scotland_data %>%
  left_join(rsv_result_codes %>% select(-original_description_rsv, -file) %>% unique(), 
            by = c("result_rsv" = "original_code_rsv")) %>% 
  left_join(mpv_result_codes %>% select(-original_description_mpv, -file) %>% unique(), 
            by = c("result_mpv" = "original_code_mpv"))

# standardize health board codes
health_board_codes <- data.frame(
  original_code = c(
    "AA", "AC", "BOR", "DG", "EIRE", "ENG", "FIFE", "FM", "FV", "GGHB",
    "GOLDJ", "GRA", "HIG", "LAN", "LOT", "TAY", "WI", 
    "AAHB", "ABHGHB", "DGHB", "FFHB", "FVHB", "GGCHB", "GOLJUB", "GRHB",
    "HGHB", "LNHB", "LOHB", "SERVIC", "TYHB", "WIHB"), 
  file = c(rep(1,17), rep(2, 14)), 
  new_code = c(
    "AA", "AB", "BOR", "DG", "EIRE", "ENG", "FIFE", "FM", "FV", "GG",
    "GOLDJ", "GR", "HIG", "LAN", "LOT", "TAY", "WI", 
    "AA", "AB", "DG", "FIFE", "FV", "GG", "GOLDJ", "GR", "HIG", "LAN",
    "LOT", "TAY", "WI", "SERVIC")
)

scotland_data <- scotland_data %>% 
  left_join(health_board_codes %>% select(-file),
            by = c("health_board" = "original_code")) %>%
  rename(health_board_new = new_code)

# convert age to numeric, reformat dates
scotland_data <- scotland_data %>%
  mutate(age_cont = ifelse(substr(age, nchar(age), nchar(age)) == "d", 
                           as.integer(substr(age,1,nchar(age)-1))/365, 
                           ifelse(substr(age,nchar(age), nchar(age)) == "w", 
                                  as.integer(substr(age, 1, nchar(age)-1))*7/365, 
                                  as.numeric(age))), 
         date_collected = as.Date(date_collected), 
         date_received = as.Date(date_received), 
         wk_collected = lubridate::ceiling_date(date_collected, unit = "weeks", week_start = 1))

# add some age bins
# note: if length(bin_cutoffs) = n, length(bin_names) = n + 1
convert_age_to_bin <- function(x, bin_cutoffs, bin_names){
  bins_fit <- which(bin_cutoffs > x)
  if(is_empty(bins_fit)){return(bin_names[length(bin_names)])}
  else{return(bin_names[min(bins_fit)])}
}

age_bin_cutoffs <- c(1, 5, 65)
age_bins <- c("<1 yr", "1-5 yr", "5-65 yr", ">65 yr")

age_bin_cutoffs_expanded <- c(seq(9,41, 8)*7/365,1:10, seq(20,100,10))
age_bins_expanded <- c("<2 mo",
                       paste0(seq(2,10,2),"-", seq(4,12,2), " mo"),
                       paste0(c(1:9,seq(10,90,10)), "-", c(2:10,seq(20,100,10)), "yr"), ">100 yr")


scotland_data <- scotland_data %>%
  mutate(age_bin_expanded = unlist(purrr::map(age_cont, convert_age_to_bin, 
                                              bin_cutoffs = age_bin_cutoffs_expanded, 
                                              bin_names = age_bins_expanded)), 
         age_bin = unlist(purrr::map(age_cont, convert_age_to_bin, 
                                     bin_cutoffs = age_bin_cutoffs, 
                                     bin_names = age_bins)))

# create long version
scotland_data_long <- scotland_data %>%
  melt(c("lab_spec_id", "age", "age_cont", "age_bin", "age_bin_expanded", "sex", 
         "date_collected", "wk_collected", "date_received", "loc_code", 
         "loc_expand", "health_board", "health_board_new"
  )) %>%
  mutate(pathogen = substr(variable, nchar(as.character(variable))-2, nchar(as.character(variable))), 
         variable = substr(variable, 1, nchar(as.character(variable))-4))

#### SUMMARIZE BY WEEK ---------------------------------------------------------

scotland_births <- read.csv("data/raw_data/scotland/scotland-weekly-births.csv")
scotland_births$Week.beginning <- as.Date(scotland_births$Week.beginning, 
                                          format = "%m/%d/%y")

scotland_by_wk <- scotland_data_long %>%
  filter(variable == "new_code", 
         value == "D") %>%
  summarize(detections = n(), .by = c("wk_collected", "pathogen"))

all_wks <- expand.grid(wk_collected = seq(min(scotland_by_wk$wk_collected), 
                                          max(scotland_by_wk$wk_collected), 
                                          by = 7), 
                       pathogen = c("rsv", "mpv"))

scotland_by_wk <- left_join(all_wks, scotland_by_wk) %>% 
  mutate(detections = ifelse(is.na(detections), 0, detections))

# merge births with incidence data
scotland_by_wk <- scotland_by_wk %>%  
  left_join(scotland_births %>% select(Week.beginning, Births.registered), 
            by = c("wk_collected" = "Week.beginning")) %>%
  # some weekly births are missing, so set to 0 if NA
  mutate(Births.registered = ifelse(is.na(Births.registered), 10, Births.registered))

# add week
scotland_by_wk <- scotland_by_wk %>%
  mutate(wk = lubridate::week(wk_collected)) %>%
  mutate(wk = ifelse(wk == 53, 52, wk)) %>% # remove wk 53 
  mutate(bi_wk = ceiling(wk/2), 
         year = as.numeric(as.factor(lubridate::year(wk_collected)))) %>%
  mutate(seas = ifelse(wk <= 26, year-1, year))

# % positive
all_wks_ages <- expand.grid(wk_collected = seq(min(scotland_by_wk$wk_collected), 
                                          max(scotland_by_wk$wk_collected), 
                                          by = 7), 
                       pathogen = c("rsv", "mpv"), 
                       age_bin = c(age_bins, "overall"))

pct_pos <- bind_rows(scotland_data_long %>%
                       filter(variable == "new_code", 
                              value == "D", 
                              date_collected > "2006-07-01") %>%
                       summarize(detections = n(), .by = c("wk_collected", "pathogen")) %>%
                       mutate(age_bin = "overall"),
                     scotland_data_long %>%
                       filter(variable == "new_code", 
                              value == "D", 
                              date_collected > "2006-07-01") %>%
                       summarize(detections = n(), .by = c("wk_collected", "pathogen", "age_bin"))) %>% 
  left_join(bind_rows(scotland_data_long %>% 
                        filter(variable == "new_code", value %in% c("ND", "D")) %>% # note: assuming insufficient and rows with no values (NA) were not tested
                        summarize(n_tests = n(), .by = c("wk_collected", "age_bin", "pathogen")), 
                      scotland_data_long %>% 
                        filter(variable == "new_code", value %in% c("ND", "D")) %>%
                        summarize(n_tests = n(), .by = c("wk_collected", "pathogen")) %>%
                        mutate(age_bin = "overall"))) %>%
  right_join(all_wks_ages) %>% 
  # filter(wk_collected >= as.Date("2006-10-01")) %>%
  arrange(wk_collected) %>%
  mutate(detections = ifelse(is.na(detections),0, detections), 
         n_tests = ifelse(is.na(n_tests), 0, n_tests)) %>%
  mutate(pct_pos = ifelse(n_tests == 0, 0, detections/n_tests), 
         age_bin = factor(age_bin, levels = c("overall", age_bins))) %>%
  # smooth number of tests across 52 weeks, center the smoothed time series
  mutate(n_tests_ma = rollmean(n_tests, 52, na.pad = TRUE), .by = c("pathogen", "age_bin")) %>%
  # correct for lower/upper bounds that are not included in rollmean()
  mutate(n_tests_ma = ifelse(is.na(n_tests_ma) & wk_collected < as.Date("2009-01-01"), 
                             min(n_tests_ma[wk_collected < as.Date("2009-01-01")], na.rm = TRUE), 
                             ifelse(is.na(n_tests_ma) & wk_collected > as.Date("2023-01-01"), 
                                    min(n_tests_ma[wk_collected > as.Date("2023-01-01")], na.rm = TRUE), n_tests_ma)),
         .by = c("pathogen", "age_bin")) %>% # NOTE: min might not be the first value but seems to work
  mutate(n_tests_ma_scaled = n_tests_ma/max(n_tests_ma), .by = c("pathogen", "age_bin")) %>%
  # smooth positivity % across 9 weeks, center the smoothed time series
  mutate(pct_pos_smooth = c(pct_pos[1:4], rollmean(pct_pos, 9), pct_pos[(n()-3):n()]), .by = c("pathogen", "age_bin")) %>%
  mutate(pct_pos_rel = (pct_pos+0.00001)/max(pct_pos+0.00001), 
         pct_pos_rel_smooth = (pct_pos_smooth+0.00001)/max(pct_pos_smooth+0.00001), .by = c("pathogen", "age_bin"))

scotland_by_wk_overall = left_join(scotland_by_wk, pct_pos %>% filter(age_bin == "overall"), 
                                   by = c("wk_collected", "pathogen", "detections")) %>%
  select(-age_bin)

# save output
write_rds(scotland_by_wk_overall, "data/derived_data/scotland/scotland_by_wk_overall.rds")

