library(dplyr)
library(rstan)
library(bayesplot)
library(readr)
library(ggplot2)
library(reshape2)
library(beepr)
library(gridExtra)
library(bbmle)
library(RcppRoll)


#### PREP DATA -----------------------------------------------------------------
source("R/scotland_parameters.R")
source("R/utils.R")

# read in scotland data, all age groups, summarized by week
scotland_by_wk <- read_rds("data/derived_data/scotland/scotland_by_wk_overall.rds")

scotland_by_wk_pre <- scotland_by_wk %>%
  filter(wk_collected >= start_wk_pre, 
         wk_collected <= end_wk_pre) 

#### FIT WITH LAGS -------------------------------------------------------------
# constrain sampling of initial conditions to ensure S0 + E0 + I0 < 1
initfun = function(...){
  list(
    S0_mpv = rbeta(1, 20, 80),
    E0_mpv = rbeta(1, 20, 80),
    I0_mpv = rbeta(1, 20, 80),
    S0_rsv = rbeta(1, 20, 80),
    E0_rsv = rbeta(1, 20, 80),
    I0_rsv = rbeta(1, 20, 80)
  )
}

k_tst = c(1:8)
n_samp <- 5000
samp <- sample(1:7500, n_samp)

#### FIT MODEL -----------------------------------------------------------------
model_SEIRS <- stan_model("R/1_fit_models/lags-SEIRS_RSVforce_RSVfit.stan")

# fit_scotland_SEIRS_all <- vector("list", length(k_tst))
for(i in 1:length(k_tst)){
  print(k_tst[i])
  warning(paste0("iteration: ", i," (lag of ", k_tst[i], " weeks)"))
  fit_scotland_SEIRS_all[[i]] <- sampling(
    model_SEIRS,
    data = list(
      N = nrow(scotland_by_wk_pre %>% filter(pathogen == "mpv")),
      cases_mpv = scotland_by_wk_pre %>% filter(pathogen == "mpv") %>% pull(detections),
      cases_rsv = scotland_by_wk_pre %>% filter(pathogen == "rsv") %>% pull(detections), 
      birth = scotland_by_wk_pre %>% filter(pathogen == "mpv") %>% pull(Births.registered), 
      tests_mpv = scotland_by_wk_pre %>% filter(pathogen == "mpv") %>% pull(n_tests_ma_scaled),
      tests_rsv = scotland_by_wk_pre %>% filter(pathogen == "rsv") %>% pull(n_tests_ma_scaled),
      pop = scotland_N,
      sigma = scotland_sigma,
      gamma = scotland_gamma,
      mu = scotland_birthrate,
      omega = scotland_omega,
      k = k_tst[i]
    ),
    seed = 7, 
    iter = 10000,
    chain = 4, 
    init = initfun, 
    cores = 4,
    control=list(max_treedepth = 12)
  )
  saveRDS(fit_scotland_SEIRS_all[[i]], paste0("data/derived_data/scotland/fit_scotland_SEIRS_NEW_lag",k_tst[i],".rds"))
}  
error <- names(warnings())
out <- file(paste0("data/derived_data/scotland/", "warnings_NEW.txt"))
writeLines(error, out)
close(out)

#### LOAD POSTERIORS (IF ALREADY SAVED FITS PREVIOUSLY) ------------------------
fit_scotland_SEIRS_all = vector("list", length(k_tst))
for(z in 1:length(k_tst)){
  print(z)
  f <- paste0("data/derived_data/scotland/fit_scotland_SEIRS_NEW_lag",
              k_tst[z],".rds")
  fit_scotland_SEIRS_all[[z]] <- read_rds(f)
}


#### PLOT FITS -----------------------------------------------------------------
for(i in 1:length(k_tst)){
  print(
    bind_rows(
      filter_posterior(as.data.frame(fit_scotland_SEIRS_all[[i]] ), "Ipred_mpv\\[") %>%
        mutate(pathogen = "mpv"), 
      filter_posterior(as.array(fit_scotland_SEIRS_all[[i]] ), "Ipred_rsv\\[") %>%
        mutate(pathogen = "rsv")
    ) %>%
      left_join(t_to_wk_scotland) %>%
      filter(draw_id %in% sample(1:n_samp,100)) %>%
      mutate(pathogen = factor(pathogen, levels = c("rsv", "mpv"))) %>%
      ggplot(aes(x = wk_collected, color = pathogen)) + 
      geom_point(data = scotland_by_wk_pre, 
                 aes(y = detections), shape = 21) + 
      geom_line(data = scotland_by_wk_pre, 
                aes(y = detections), alpha = 0.2) + 
      geom_line(aes(y = value, group = draw_id), alpha = 0.2) +
      ggtitle(paste0("lag = ", k_tst[i])) +
      facet_grid(rows = vars(pathogen), scales = "free", switch = "y") + 
      scale_color_brewer(palette = "Dark2") + 
      theme_bw() + 
      theme(axis.title.x = element_blank(), 
            legend.position = "none", 
            panel.grid.minor.y = element_blank(),
            strip.background = element_blank(), 
            strip.placement = "outside")
  )
}

#### PLOT PARAMETER ESTIMATES --------------------------------------------------
# parameter estimates across all lags
summarize_param_est = function(lag, fit){
  pars = c("b", "r", "a", "p", "c", "S0_rsv", "S0_mpv", "E0_rsv", "E0_mpv", 
           "I0_rsv", "I0_mpv", "max_RSV_inc", "rho_mpv", "rho_rsv")
  ext <- extract(fit, pars)
  ext[["c_prime"]] = ext$c/ext$max_RSV_inc*scotland_N/lag
  ext[["beta_mpv"]] = ext$b*ext$r
  sum <- sapply(ext, function(i){return(c(
    quantile(i, 0.05),
    quantile(i, 0.25), 
    quantile(i, 0.5), 
    quantile(i, 0.75),
    quantile(i, 0.95)))}) %>%
    as.data.frame() %>% 
    mutate(quantile = paste0("Q", c(5, 25, 50, 75, 95))) %>%
    melt(c("quantile"), variable.name = "parameter") %>%
    dcast(parameter ~ quantile, value.var = "value") %>% 
    mutate(lag = lag)
  return(sum)
}

par_summary <- lapply(1:length(fit_scotland_SEIRS_all), function(i){summarize_param_est(i, fit_scotland_SEIRS_all[[i]])}) %>%
  bind_rows()

par_labs = c("RSV transmission rate", "HMPV transmission rate", "seasonal phase", "seasonal amplitude", "HMPV relative transmissibility", 
             "effect of RSV on HMPV transmission", "RSV reporting rate", "HMPV reporting rate", "RSV initial susceptible", "RSV initial infected", "HMPV initial susceptible", "HMPV initial infected", "RSV initial exposed", "HMPV initial exposed", "per capita effect of RSV on HMPV transmission")
names(par_labs) = c("b", "beta_mpv", "p", "a", "r", "c", "rho_rsv", "rho_mpv", "S0_rsv", "I0_rsv", "S0_mpv", "I0_mpv", "E0_rsv", "E0_mpv", "c_prime")

par_summary %>% filter(!(parameter %in% c("c", "max_RSV_inc"))) %>%
  mutate(parameter = factor(parameter, levels = c("b", "r", "beta_mpv", "c_prime", "p", "a","rho_rsv", "rho_mpv", "S0_rsv", "I0_rsv", "S0_mpv", "I0_mpv", "E0_rsv", "E0_mpv"))) %>%
  ggplot(aes(x = lag-1)) + 
  geom_segment(aes(xend = lag-1, y = Q5, yend = Q95)) +
  geom_segment(aes(xend = lag-1, y = Q25, yend = Q75), size = 1.5) +
  geom_point(aes(y = Q50), size = 2.25) + 
  facet_wrap(vars(parameter), scales = "free", labeller = labeller(parameter = par_labs)) + 
  labs(x = "lag in RSV incidence (weeks)", y = "posterior distribution") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank())
ggsave("figures/lagged_param_estimates.pdf", width = 13, height = 9)

# plot transmission rate
p_beta = vector("list", length(fit_scotland_SEIRS_all))
for(i in k_tst){#1:length(fit_scotland_SEIRS_all)){
  ext <- extract(fit_scotland_SEIRS_all[[i]], paste0("beta_mpv[", 1:611, "]"))
  p_beta[[i]] <- sapply(ext, function(i){return(c(
    quantile(i, 0.05),
    quantile(i, 0.25), 
    quantile(i, 0.5), 
    quantile(i, 0.75),
    quantile(i, 0.95)))}) %>%
    as.data.frame() %>%
    mutate(quantile = paste0("Q", c(5, 25, 50, 75, 95))) %>%
    melt(c("quantile"), variable.name = "variable") %>%
    mutate(t = as.integer(gsub("\\]", "", gsub("beta_mpv\\[", "", variable)))) %>% 
    left_join(t_to_wk_scotland) %>%
    dcast(wk_collected ~ quantile, value.var = "value") %>%
    filter(wk_collected >= end_wk_pre-3*365) %>%
    ggplot(aes(x = wk_collected)) + 
    geom_ribbon(aes(ymin = Q5, ymax = Q95), alpha = 0.3) + 
    geom_ribbon(aes(ymin = Q25, ymax = Q75), alpha = 0.5) +
  geom_line(aes(y = Q50)) + 
    ggtitle(paste("lag =", k_tst[i]-1)) + 
    labs(y = "transmission rate") + 
    scale_x_date(date_labels = "%b %Y") + 
    theme_bw() + 
    theme(axis.title.x = element_blank(), 
          legend.position = "none", 
          panel.grid.minor.y = element_blank(),
          strip.background = element_blank(), 
          strip.placement = "outside")
}

p = grid.arrange(grobs = p_beta, ncol = 3)
ggsave("figures/lagged_beta.pdf", p, width = 8, height = 6)

