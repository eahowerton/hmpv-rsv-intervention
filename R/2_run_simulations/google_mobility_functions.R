# raw files from https://www.google.com/covid19/mobility/

#' function to extract google mobility data and manipulate into single multiplier
#' 
#' @param country_name character indicating country abbreviation for desired 
#' mobility data (`"CA"` = Canada, `"GB"` = United Kingdom, `"KR"` = Korea)
#' @param regional logical indicating whether regional mobility data should be 
#' returned, defaults to FALSE
get_google_mobility_data <- function(country_abbreviation, regional = FALSE){
  mobility_2020 <- read.csv(paste0("data/raw_data/google-mobility/2020_",
                                   country_abbreviation,"_Region_Mobility_Report.csv"))
  mobility_2021 <- read.csv(paste0("data/raw_data/google-mobility/2021_",
                                   country_abbreviation,"_Region_Mobility_Report.csv"))
  mobility_2022 <- read.csv(paste0("data/raw_data/google-mobility/2022_",
                                   country_abbreviation,"_Region_Mobility_Report.csv"))
  all_mobility <- bind_rows(mobility_2020, mobility_2021, mobility_2022)
  all_mobility_transformed <- transform_google_mobility_data(all_mobility)
  return(all_mobility_transformed)   
}

#' transform multiple mobility metrics into a single multiplier
#' by taking mean of retail and recreation, grocery and pharmacy, transit, and 
#' workplaces % change from baseline
transform_google_mobility_data <- function(raw_mobility, regional = FALSE){
  if(regional){
    raw_mobility_transform <- raw_mobility %>%
      filter(sub_region_1 != "")
  }
  else{
    raw_mobility_transform <- raw_mobility %>%
      filter(sub_region_1 == "") 
  }
  raw_mobility_transform <- raw_mobility_transform %>%
    mutate(region = sub_region_1) %>%
    mutate(mean_mob = (retail_and_recreation_percent_change_from_baseline +
                         grocery_and_pharmacy_percent_change_from_baseline + 
                         transit_stations_percent_change_from_baseline + 
                         workplaces_percent_change_from_baseline)/4, 
           .by = c(region, date)) %>%
    mutate(wk_collected = lubridate::ceiling_date(as.Date(date), unit = "weeks", week_start = 1)) %>%
    summarize(mean_mob = mean(mean_mob, na.rm = TRUE), 
              .by = c(region, wk_collected))
  return(raw_mobility_transform)
}
