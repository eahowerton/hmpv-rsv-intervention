#### REGION MAPPING ------------------------------------------------------------
region_map <- c(
  "CANADA", rep("Atlantic",8), "Province of Québec", "Province of Ontario", 
  "Prairies", "Prairies", "Province of Alberta", 
  rep("Province of British Columbia",2), rep("Territories", 7)
)
names(region_map) = c(
  "Canada", 
  # atlantic
  "Newfoundland and Labrador", "Newfoundland.and.Labrador.3", 
  "Prince Edward Island","Prince.Edward.Island", "Nova Scotia", "Nova.Scotia", 
  "New Brunswick", "New.Brunswick", 
  # quebec
  "Quebec", 
  # ontario
  "Ontario", 
  # prairies
  "Saskatchewan", "Manitoba",  # note: alberta is in prairies per CA DPH, but excluding here because of data issue
  # alberta
  "Alberta", 
  # british columbia
  "British Columbia", "British.Columbia", 
  # territories
  "Yukon", "Northwest Territories", "Northwest Territories 5","Northwest.Territories.3", 
  "Nunavut", "Nunavut 5", "Nunavut.3"
)


#### DETECTIONS ----------------------------------------------------------------
canada_by_wk <- readr::read_rds("data/derived_data/canada/data_canada_scrape.rds")

# exclue alberta from prairies region
canada_by_wk_prairies <- canada_by_wk %>%
  filter(region == "Prairies" & !(reporting_lab %in% c("Province of Alberta", "Prairies"))) %>%
  summarize(rsv_test = sum(rsv_test, na.rm = TRUE), 
            rsv_pos = sum(rsv_pos, na.rm = TRUE), 
            mpv_test = sum(mpv_test, na.rm = TRUE), 
            mpv_pos = sum(mpv_pos, na.rm = TRUE), 
            .by = c("wk_collected", "region")) %>%
  mutate(region_flag = 1, reporting_lab = "Prairies")

canada_by_wk <- canada_by_wk %>%
  filter(region_flag == 1, region != "Prairies") %>%
  bind_rows(canada_by_wk_prairies)

# fill in missing weeks
canada_by_wk <- canada_by_wk %>%
  right_join(expand.grid(wk_collected = seq(min(canada_by_wk$wk_collected), 
                                           max(canada_by_wk$wk_collected), by = 7), 
                         region = unique(canada_by_wk$region)))

# reshape
canada_by_wk <- canada_by_wk %>% 
  select(-reporting_lab, -region_flag) %>%
  melt(c("region", "wk_collected", "rsv_pos", "mpv_pos"), value.name = "n_tests") %>%
  arrange(wk_collected) %>%
  mutate(n_tests_ma = zoo::rollapply(n_tests, 52, fill = NA, FUN=function(x){ifelse(all(is.na(x)), NA,  mean(x, na.rm=TRUE))}),
         .by = c("variable", "region")) %>%
  mutate(n_tests_ma = ifelse(is.na(n_tests_ma) & wk_collected < as.Date("2015-06-01"), 
                             min(n_tests_ma[wk_collected < as.Date("2015-06-01")], na.rm = TRUE), 
                             ifelse(is.na(n_tests_ma) & wk_collected > as.Date("2023-01-01"), 
                                    min(n_tests_ma[wk_collected > as.Date("2023-01-01")], na.rm = TRUE), n_tests_ma)),
         .by = c("variable", "region")) %>% # NOTE: min might not be the first value (but hack for now), seems to work
  mutate(n_tests_ma_scaled = n_tests_ma/max(n_tests_ma), .by = c("variable", "region")) %>%
  mutate(variable = substr(variable, 1,3)) %>%
  pivot_wider(id_cols = c("region", "wk_collected", "rsv_pos", "mpv_pos"), 
              names_from = variable, 
              values_from = c(n_tests, n_tests_ma, n_tests_ma_scaled)) %>%
  rename(pos_rsv = rsv_pos, pos_mpv = mpv_pos) %>%
  mutate(year = lubridate::year(wk_collected))

# save
write_rds(canada_by_wk, "data/derived_data/canada/canada_by_wk_overall.rds")

