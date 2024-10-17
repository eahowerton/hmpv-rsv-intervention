library(xml2)
library(rvest)
library(purrr)
library(dplyr)

clean_up_table_vals <- function(table_df, wk_collected){
  if(any(grepl("Reporting Laboratory",unlist(table_df[,1])))){
    k = which(grepl("Reporting Laboratory",unlist(table_df[,1])))
    colnames(table_df) = paste(colnames(table_df),table_df[k,])
    table_df = table_df[-k,]}
  if(any(grepl("Delays",unlist(table_df[,1])))){
    k = which(grepl("Delays",unlist(table_df[,1])))
    table_df = table_df[-k,]}
  if(any(grepl("*Territories include",unlist(table_df[,1])))){
    k = which(grepl("*Territories include",unlist(table_df[,1])))
    table_df = table_df[-k,]}
  if(any(grepl("Specimens from YT",unlist(table_df[,1])))){
    k = which(grepl("Specimens from YT",unlist(table_df[,1])))
    table_df = table_df[-k,]}
  na_strings <- c("N/A", "1N.A.","Not Available", "not available", "Not available", "N.A.",
                  "N.R.","N.R" , "N.C.", "N.D.", "Not Tested", "not tested", "Not tested", "non testé")
  table_df <- table_df %>% 
    mutate(across(where(is.character), 
                  ~replace(.x, .x %in% na_strings, NA)))
  table_df[,2:ncol(table_df)] <- sapply(table_df[,2:ncol(table_df)], as.integer)
  table_df <- table_df %>%
    dplyr::mutate(wk_collected = wk_collected)
  return(table_df)
}


# function to pull table names
get_names <- function(n){
  as.character(strsplit(strsplit(n, "<")[[1]][2],">")[[1]][2])
}

seasons <- c(paste0("201",3:8, "-201",4:9), "2019-2020", paste0("202",0:3, "-202",1:4))

df_seasons <- cum_df_seasons <- list()
start.time = Sys.time()
for(j in 1:length(seasons)){
  if(seasons[j] == "2023-2024"){
    raw_html_level1 <- xml2::read_html("https://www.canada.ca/en/public-health/services/surveillance/respiratory-virus-detections-canada.html")}
  else{
    raw_html_level1 <- xml2::read_html(paste0("https://www.canada.ca/en/public-health/services/surveillance/respiratory-virus-detections-canada/",
                                              seasons[j],".html"))
  }
  all_links <- rvest::html_elements(raw_html_level1, "a") %>%
    rvest::html_attr('href')
  keep_link_indicator <- paste0("/en/public-health/services/surveillance/respiratory-virus-detections-canada/",seasons[j])
  surveillance_report_links <- all_links[unlist(purrr::map(all_links, function(x){grepl(keep_link_indicator,x)}))]
  if(seasons[j] == "2016-2017"){ # manually fix issue with some links in 2016-2017
    other_links <- paste0("/en/public-health/services/surveillance/respiratory-virus-detections-canada/2016-2017/week-", 19:29, 
                          "-ending-", c("may-13", "may-20", "may-27", "june-3", "june-10", "june-17", "june-24", "july-1", "july-8", "july-15", "july-22"),"-2017.html") 
    surveillance_report_links <- c(surveillance_report_links, other_links)
  }
  if(seasons[j] == "2020-2021"){ # anomaly in 2020-2021 links: after Jan. 1, 2021, links use 2021-2022
    keep_link_indicator2 <- "/en/public-health/services/surveillance/respiratory-virus-detections-canada/2021-2022"
    surveillance_report_links <- c(surveillance_report_links, 
                                   all_links[unlist(purrr::map(all_links, function(x){grepl(keep_link_indicator2,x)}))])
  }
  df_list <- cum_df_list <- list()
  for(i in 1:length(surveillance_report_links)){
    print(paste("j:",j, "i:",i))
    # 1. read HTML from website
    raw_html_level2 <- xml2::read_html(paste0("https://www.canada.ca", surveillance_report_links[i]))
    # 2. extract all tables
    tables_html <- rvest::html_elements(raw_html_level2, "table")
    # 3. convert Table 1 to data.frame
    # (I believe everything else is a summary of these data)
    tab_num = 1
    if(j == 5 & i == 49){browser();tab_num = 2} # j == 7
    if(j == 7 & i == 30 | j == 7 & i == 36){browser();print("skipping - no Table 1");next}
    if(j < 7 & !grepl("Table 1", as.character(tables_html[[tab_num]]))){if(!(j == 4 & i == 52)){browser();print("skipping - no Table 1");next}}
    table_df <- html_table(tables_html[[tab_num]])
    table_df_cum <- html_table(tables_html[[tab_num+1]])
    date_char <- strsplit(strsplit(surveillance_report_links[i], "ending-")[[1]][2], ".html")[[1]][1]
    if(is.na(date_char)){browser()}
    if(grepl("february-29", date_char)){date_char = paste0("february-28-", substr(date_char, 13,16))}
    wk_collected <- lubridate::mdy(date_char)
    # clean up table a bit
    df_list[[i]] <- clean_up_table_vals(table_df, wk_collected)
    cum_df_list[[i]] <- clean_up_table_vals(table_df_cum, wk_collected)
  }
  browser()
  df_seasons[[j]] <- df_list
  cum_df_seasons[[j]] <- cum_df_list
}

readr::write_rds(df_seasons, "timeseries-data/data/canada/data_canada_scrape_raw.rds")


colnames_lab <- c("Reporting Laboratoryr", "Reporting Laboratory", "Reporting Laboratory Reporting Laboratory")
colnames_rsvtest <- c("R.S.V.\nTest", "R.S.V. Test","RSV Tested", "R.S.V.\r\n          Test", "RSV Test")
colnames_rsvpos <- c("R.S.V.\nPos.", "R.S.V. Pos.", "RSV Positive", "R.S.V.\r\n          Pos.", "RSV Pos.")
colnames_mpvtest <- c("hMPV\nTest", "hMPV Test", "hMPV Tested","HMPV Tested", "hMPV\r\n          Test")
colnames_mpvpos <- c("hMPV\nPos.", "hMPV Pos.", "hMPV Positive", "HMPV Positive", "hMPV\r\n          Pos.")


filt_to_rsv_mpv <- function(d){
  if(any(is.null(d))){return(NA)}
  d_filt <- d[,which(colnames(d) %in% c(colnames_rsvtest, colnames_rsvpos, 
                                        colnames_mpvtest, colnames_mpvpos, 
                                        colnames_lab, "wk_collected"))]
  colnames(d_filt) <- ifelse(colnames(d_filt) %in% colnames_rsvtest, "rsv_test", 
                            ifelse(colnames(d_filt) %in% colnames_rsvpos, "rsv_pos", 
                                   ifelse(colnames(d_filt) %in% colnames_mpvtest, "mpv_test", 
                                          ifelse(colnames(d_filt) %in% colnames_mpvpos, "mpv_pos",
                                                 ifelse(colnames(d_filt) %in% colnames_lab, "reporting_lab",
                                                 colnames(d_filt))))))
  return(d_filt)
}

data_rsv_mpv <- do.call(c, df_seasons)
data_rsv_mpv <- lapply(data_rsv_mpv, filt_to_rsv_mpv)
data_rsv_mpv <- data_rsv_mpv[!sapply(data_rsv_mpv, function(i){all(is.na(i))})]
data_rsv_mpv <- bind_rows(data_rsv_mpv)

cum_data_rsv_mpv <- do.call(c, cum_df_seasons)
cum_data_rsv_mpv <- lapply(cum_data_rsv_mpv, filt_to_rsv_mpv)
cum_data_rsv_mpv <- cum_data_rsv_mpv[!sapply(cum_data_rsv_mpv, function(i){all(is.na(i))})]
cum_data_rsv_mpv <- bind_rows(cum_data_rsv_mpv)

# fix a few things
data_rsv_mpv <- data_rsv_mpv %>%
  mutate(
    reporting_lab = ifelse(reporting_lab %in% c("Canada", "CAN.A."), "CANADA", reporting_lab), 
    wk_collected = ifelse(is.na(wk_collected), as.Date("2021-08-07"), wk_collected)) %>%
  mutate(
    wk_collected = case_when( # fix errors in the dates embedded in URL text
      wk_collected == as.Date("2013-09-29") ~ as.Date("2013-09-28"),
      wk_collected == as.Date("2013-11-24") ~ as.Date("2013-11-30"),
      wk_collected == as.Date("2013-12-29") ~ as.Date("2013-12-28"),
      wk_collected == as.Date("2014-06-29") ~ as.Date("2014-06-28"),
      wk_collected == as.Date("2015-03-29") ~ as.Date("2015-03-28"),
      wk_collected == as.Date("2015-11-29") ~ as.Date("2015-11-28"), 
      wk_collected == as.Date("2016-05-29") ~ as.Date("2016-05-28"),
      wk_collected == as.Date("2017-01-29") ~ as.Date("2017-01-28"), 
      wk_collected == as.Date("2017-07-30") ~ as.Date("2017-07-29"),
      wk_collected == as.Date("2020-02-28") ~ as.Date("2020-02-29"),
      wk_collected == as.Date("2023-09-03") ~ as.Date("2023-09-02"),
      .default = as.Date(wk_collected)
    ))

cum_data_rsv_mpv <- cum_data_rsv_mpv %>%
  mutate(
    reporting_lab = ifelse(reporting_lab %in% c("Canada", "CAN.A."), "CANADA", reporting_lab), 
    wk_collected = ifelse(is.na(wk_collected), as.Date("2021-08-07"), wk_collected)) %>%
  mutate(
    wk_collected = case_when( # fix errors in the dates embedded in URL text
      wk_collected == as.Date("2013-09-29") ~ as.Date("2013-09-28"),
      wk_collected == as.Date("2013-11-24") ~ as.Date("2013-11-30"),
      wk_collected == as.Date("2013-12-29") ~ as.Date("2013-12-28"),
      wk_collected == as.Date("2014-06-29") ~ as.Date("2014-06-28"),
      wk_collected == as.Date("2015-03-29") ~ as.Date("2015-03-28"),
      wk_collected == as.Date("2015-11-29") ~ as.Date("2015-11-28"), 
      wk_collected == as.Date("2016-05-29") ~ as.Date("2016-05-28"),
      wk_collected == as.Date("2017-01-29") ~ as.Date("2017-01-28"), 
      wk_collected == as.Date("2017-07-30") ~ as.Date("2017-07-29"),
      wk_collected == as.Date("2020-02-28") ~ as.Date("2020-02-29"),
      wk_collected == as.Date("2023-09-03") ~ as.Date("2023-09-02"),
      .default = as.Date(wk_collected)
    ))


# update reporting_lab names
region_key = c(
  rep("Atlantic", 5), 
  rep("Province of Québec", 7), 
  rep("Province of Ontario", 23), 
  rep("Territories", 4), 
  rep(NA, 2), # for Saskatoon and Regina, which are subsets of Manitoba
  rep("Prairies", 4), 
  "Province of British Columbia",
  "CANADA"
)

names(region_key) = c(
  # atlantic
  "Newfoundland","Prince Edward Island","Nova Scotia","New Brunswick",
  "Atlantic",
  # quebec
  "Région Nord-Est","Québec-Chaudière-Appalaches","Centre-du-Québec",
  "Montréal-Laval","Ouest du Québec","Montérégie",
  "Province of Québec",
  # ontario
  "CHEO - Ottawa","Ottawa","EORLA","Kingston","Toronto Medical Laboratory", "Toronto",
  "C. - Toronto","Sick Kids' Hospital - Toronto","Sault Ste. Marie", "Shared Hospital Laboratory",
  "Timmins","St. Joseph's - London","London","Orillia","Thunder Bay", "Sault Area Hospital",
  "Sudbury","Hamilton","Peterborough","St. Joseph's - Hamilton", "UHN/Mount Sinai Hospital",
  "Sunnybrook & Women's College HSC", 
  "Province of Ontario",
  # territories
  "Yukon", "Northwest Territories","Nunavut",
  "Territories",
  # prairies
  "Saskatoon", "Regina", # subsets of Manitoba, so do not assign to Prairies region
  "Province of Manitoba", "Province of Saskatchewan","Province of Alberta",
  "Prairies", 
  # British Columbia
  "Province of British Columbia",
  # national
  "CANADA"
)

data_rsv_mpv <- data_rsv_mpv %>%
  mutate(reporting_lab = gsub("\\P.H.L.", "", reporting_lab)) %>%
  mutate(reporting_lab = gsub("\\P.H.O.L. - ", "", reporting_lab)) %>%
  mutate(reporting_lab = gsub("\\P.H.O.L - ", "", reporting_lab)) %>%
  mutate(reporting_lab = gsub("\\s+"," ", reporting_lab)) %>% # remove extra spaces
  mutate(reporting_lab = trimws(reporting_lab)) %>%
  mutate(reporting_lab = case_when(
    reporting_lab == "CHEO/HEE0 - Ottawa" ~ "CHEO - Ottawa", 
    reporting_lab == "UHN / Mount Sinai Hospital" ~ "UHN/Mount Sinai Hospital", 
    reporting_lab == "St. Joseph's Healthcare - Hamilton" ~ "St. Joseph's - Hamilton", 
    reporting_lab == "Province of" ~ "Province of Alberta", # error in table data, manual override
    reporting_lab == "Alberta" ~ "Province of Alberta",
    reporting_lab == "Ontario" ~ "Province of Ontario",
    reporting_lab == "Saskatchewan" ~ "Province of Saskatchewan",
    reporting_lab == "Manitoba" ~ "Province of Manitoba",
    reporting_lab == "British Columbia" ~ "Province of British Columbia",
    reporting_lab == "Québec" ~ "Province of Québec",
    reporting_lab %in% c("Sick Kids Hospital - Toronto", "Sick Kids'Hospital - Toronto") ~ "Sick Kids' Hospital - Toronto", 
    reporting_lab %in% c("Territories*","Territories/territoires") ~ "Territories", 
    .default = reporting_lab
  )) %>% 
  mutate(region = region_key[reporting_lab]) %>%
  mutate(region_flag = ifelse(reporting_lab == region, 1, 0))

readr::write_rds(data_rsv_mpv, "timeseries-data/data/canada/data_canada_scrape.rds")


cum_data_rsv_mpv <- cum_data_rsv_mpv %>%
  mutate(reporting_lab = gsub("\\P.H.L.", "", reporting_lab)) %>%
  mutate(reporting_lab = gsub("\\P.H.O.L. - ", "", reporting_lab)) %>%
  mutate(reporting_lab = gsub("\\P.H.O.L - ", "", reporting_lab)) %>%
  mutate(reporting_lab = gsub("\\s+"," ", reporting_lab)) %>% # remove extra spaces
  mutate(reporting_lab = trimws(reporting_lab)) %>%
  mutate(reporting_lab = case_when(
    reporting_lab == "CHEO/HEE0 - Ottawa" ~ "CHEO - Ottawa", 
    reporting_lab == "UHN / Mount Sinai Hospital" ~ "UHN/Mount Sinai Hospital", 
    reporting_lab == "St. Joseph's Healthcare - Hamilton" ~ "St. Joseph's - Hamilton", 
    reporting_lab == "Province of" ~ "Province of Alberta", # error in table data, manual override
    reporting_lab == "Alberta" ~ "Province of Alberta",
    reporting_lab == "Ontario" ~ "Province of Ontario",
    reporting_lab == "Saskatchewan" ~ "Province of Saskatchewan",
    reporting_lab == "Manitoba" ~ "Province of Manitoba",
    reporting_lab == "British Columbia" ~ "Province of British Columbia",
    reporting_lab == "Québec" ~ "Province of Québec",
    reporting_lab %in% c("Sick Kids Hospital - Toronto", "Sick Kids'Hospital - Toronto") ~ "Sick Kids' Hospital - Toronto", 
    reporting_lab %in% c("Territories*","Territories/territoires") ~ "Territories", 
    .default = reporting_lab
  )) %>% 
  mutate(region = region_key[reporting_lab]) %>%
  mutate(region_flag = ifelse(reporting_lab == region, 1, 0))

readr::write_rds(cum_data_rsv_mpv, "timeseries-data/data/canada/cum_data_canada_scrape.rds")


ggplot(data = data_rsv_mpv %>% 
         filter(region_flag == 1) %>%
         melt(c("region","reporting_lab", "region_flag","wk_collected")) %>%
         separate(variable, into = c("pathogen", "variable")) %>% 
         filter(variable == "pos", region != "Territories"), 
       aes(x = wk_collected, y = value, color = pathogen)) + 
  geom_line(size = 0.8) + 
  facet_wrap(vars(region, pathogen), ncol = 2, scales = "free") + 
  theme_bw()

ggplot(data = data_rsv_mpv %>% 
         filter(region_flag == 1) %>%
         melt(c("region","reporting_lab", "region_flag","wk_collected")) %>%
         separate(variable, into = c("pathogen", "variable")) %>% 
         filter(variable == "pos", region != "Territories"), 
       aes(x = wk_collected, y = value, color = pathogen)) + 
  geom_line(size = 0.8) + 
  facet_wrap(vars(region), scales = "free") + 
  theme_bw()




tmp_inc <- df_seasons[[2]]
tmp_inc <- lapply(tmp_inc, filt_to_rsv_mpv)
tmp_inc <- tmp_inc[!sapply(tmp_inc, function(i){all(is.na(i))})]
tmp_inc <- bind_rows(tmp_inc)

tmp_cum <- cum_df_seasons[[2]]
tmp_cum <- lapply(tmp_cum, filt_to_rsv_mpv)
tmp_cum <- tmp_cum[!sapply(tmp_cum, function(i){all(is.na(i))})]
tmp_cum <- bind_rows(tmp_cum)

# convert cum to inc
tmp_cum <- tmp_cum %>%
  melt(c("reporting_lab", "wk_collected")) %>%
  arrange(wk_collected) %>%
  mutate(inc_from_cum = c(0,diff(value)), .by = c("reporting_lab", "variable")) #%>%
# dcast(reporting_lab + wk_collected ~ variable, value.var = "inc") %>%
# rename(rsv_test_fromcum = rsv_test, 
#        rsv_pos_fromcum = rsv_pos, 
#        mpv_test_fromcum = mpv_test, 
#        mpv_pos_fromcum = mpv_pos)

# merge and compare
View(tmp_inc %>%
       melt(c("reporting_lab", "wk_collected")) %>%
       arrange(wk_collected) %>%
       left_join(tmp_cum %>% rename(cum_value = value)) %>%
       mutate(diff = abs(value - inc_from_cum)) %>%
       filter(diff > 0 & wk_collected != "2013-08-31"))


data_rsv_mpv %>%
  melt(c("reporting_lab", "wk_collected", "region", "region_flag")) %>%
  arrange(wk_collected) %>%
  left_join(cum_data_rsv_mpv %>% 
              melt(c("reporting_lab", "wk_collected", "region", "region_flag")) %>%
              arrange(wk_collected) %>%
              mutate(inc_from_cum = c(0,diff(value)), .by = c("reporting_lab", "variable", "region", "region_flag")) %>%
              rename(cum_value = value)) %>%
  filter(variable == "rsv_pos" & inc_from_cum > 0 & reporting_lab == "Province of Alberta") %>%
  ggplot(aes(x = wk_collected)) +
  geom_line(aes(y = value, color = "inc")) + 
  geom_line(aes(y = inc_from_cum, color = "cum")) + 
  facet_wrap(vars(reporting_lab), scales = 'free')

