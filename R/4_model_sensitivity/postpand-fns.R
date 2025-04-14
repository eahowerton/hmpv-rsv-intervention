#' Function to perform entire analysis of out-of-sample postpandemic fits
#' This includes the following steps: 
#'   1. get posterior distribution for RSV and hMPV 
#'   2. simulate post-pandemic period for different npi_scalar values
#'   3. calculate negative log likelihood and choose the npi_scalar value that
#'   gives the best post-pandemic fit
#'   4. re-simulate entire time series using chosen npi_scalar
#'   5. save results
#'   6. plot outcomes
#'   7. calculate performance for each out-of-sample period and overall (nll and correlation)
#' @param fit STAN object with fitted model
#' @param samp subset of simulations to simulate forward (if NA use all simulations)
#' @param observations data.frame with observed detection data for RSV and HMPV
#' @param k number of weeks to include in the lagged sum for hMPV forcing
postpandemic_analysis = function(fit, samp, observations, mobility_df, k, 
                                 pop, start_wk_pre, end_wk_pre, t_to_wk_full){
  # 1. get posterior distribution for RSV and hMPV
  posterior = get_posterior(fit, samp)
  # 2/3. simulate post-pandemic period for different npi_scalar values and 
  # find best npi_scalar
  npi_scalar_nll = find_npi_scalar(
    posterior_rsv = posterior$rsv, posterior_mpv = posterior$mpv, 
    observations = observations, mobility_df = mobility_df, pop = pop, 
    start_wk_pre = start_wk_pre, end_wk_pre = end_wk_pre, k = k, 
    t_to_wk_full = t_to_wk_full)
  best_npi_scalar = npi_scalar_nll %>% filter(nll == min(nll)) %>% pull(npi_scale)
  # 4. re-simulate whole time series with chosen npi scalar
  npi_scaled <- get_npi_values(observations, mobility_df, best_npi_scalar)
  full_sim <- SEIRS_force_discrete_both(
    N = length(npi_scaled), parms_rsv = posterior$rsv, parms_mpv = posterior$mpv,
    testing_scalar_mpv = observations %>% filter(pathogen == "mpv") %>% pull(n_tests_ma_scaled),
    testing_scalar_rsv = observations %>% filter(pathogen == "rsv") %>% pull(n_tests_ma_scaled),
    pop = pop, npi = npi_scaled, rsvscaling_cuttoff = (end_wk_pre - start_wk_pre)/7 + 1, k = k
  )
  full_sim = full_sim %>% left_join(t_to_wk_full, by = join_by(t))
  # 5. save results
  saveRDS(full_sim, paste0("data/derived_data/scotland/postpand-sims_scotland_SEIRS_lag",k,".rds"))
  # 6. plot outcomes
  plot_samp = full_sim$draw_id[sample(1:nrow(full_sim), 100)]
  p = plot_full_sim(simulations = full_sim, observations = observations, 
                    plot_samp = plot_samp, npi_df = npi_scaled, 
                    start_wk_pre = start_wk_pre, start_wk_post = start_wk_post)
  # 7. calculate out-of-sample performance
  oos_performance = bind_rows(
    # overall performance 
    calculate_nll(full_sim %>% filter(variable == "Ipred", wk_collected > end_wk_pre), observations) %>% 
      summarize(nll = sum(nll), cor = cor(detections, value), .by = c("pathogen", "draw_id")) %>% 
      mutate(oos_period = "full"),
    # period-specific performance
    calculate_nll(full_sim %>% filter(variable == "Ipred", wk_collected > end_wk_pre), observations) %>% 
      mutate(oos_period = ifelse(wk_collected <= "2020-03-01", "pre-pand", ifelse(wk_collected >= "2021-01-01", "post-pand", "pand"))) %>%
      filter(oos_period != "pand") %>%
      summarize(nll = sum(nll), cor = cor(log(detections+1), log(value+1)), .by = c("pathogen", "oos_period", "draw_id"))
  ) %>% 
    melt(c("pathogen", "oos_period", "draw_id")) %>%
    summarize(Q5 = quantile(value, 0.05), 
              Q25 = quantile(value, 0.25), 
              Q50 = quantile(value, 0.5), 
              Q75 = quantile(value, 0.75), 
              Q95 = quantile(value, 0.95), .by = c("pathogen", "oos_period", "variable")) %>%
    mutate(k = k)
  return(list(npi_scalar_nll = npi_scalar_nll, plot_sims = p, performance = oos_performance,
              beta = full_sim %>% filter(variable == "beta")))
}


#### HELPERS -------------------------------------------------------------------
#' function to create data.frame with posterior distribution for each pathogen
#' 
#' @param fit stan object with posterior distribution
#' @param samp subset of simulations to simulate forward
get_posterior = function(fit, samp = NA){
  fit2 <- extract(fit, pars = c("S0_rsv", "E0_rsv", "I0_rsv", "rho_rsv", "b", "a", "p", 
                                "S0_mpv", "E0_mpv", "I0_mpv", "rho_mpv", "r", "c"))
  # get posterior distribution
  posterior_rsv = data.frame(
    S0 = fit2$S0_rsv, E0 = fit2$E0_rsv, I0 = fit2$I0_rsv, 
    rho = fit2$rho_rsv, b = fit2$b, a = fit2$a, p = fit2$p,
    # fixed parameters
    r = 1, sigma = scotland_sigma, omega = scotland_omega,
    gamma = scotland_gamma, mu = scotland_birthrate
  )
  posterior_mpv  = data.frame(
    S0 = fit2$S0_mpv, E0 = fit2$E0_mpv, I0 = fit2$I0_mpv,
    rho = fit2$rho_mpv, b = fit2$b, a = fit2$a, p = fit2$p,
    r = fit2$r, c = fit2$c, 
    # fixed parameters
    sigma = scotland_sigma, omega = scotland_omega,
    gamma = scotland_gamma, mu = scotland_birthrate
  )
  # randomly select a subset of simulations
  if(!any(is.na(samp))){
    posterior_rsv <- posterior_rsv[samp, ]
    posterior_mpv <- posterior_mpv[samp, ]
  }
  return(list(rsv = posterior_rsv, mpv = posterior_mpv))
}

#' @param simulations data.frame of simulated values
#' @param observations data.frame of observed detection data for RSV and HMPV
calculate_nll = function(simulations, observations){
  r = simulations %>%
    left_join(observations, by = join_by(wk_collected, pathogen)) %>%
    mutate(nll = -dpois(x = detections, lambda = value, log = TRUE)) %>%
  return(r)
}

get_npi_values = function(observations, mobility_df, npi_scalar){
  npi_scaled_base <- observations %>% 
    filter(pathogen == "mpv") %>% # can choose either pathogen, just need one
    left_join(mobility_df, by = join_by(wk_collected)) %>%
    mutate(mean_mob = ifelse(is.na(mean_mob), 1, 1+(mean_mob*npi_scalar/100))) %>%
    pull(mean_mob)
}

find_npi_scalar <- function(
    posterior_rsv, posterior_mpv, observations, mobility_df, pop, start_wk_pre, end_wk_pre, 
    k, t_to_wk_full, npi_scalar_to_test = unique(sort(c(seq(0.6, 1.3, 0.1), seq(0.9, 1.2, 0.05)))) 
                            ){
  npi_nll <- rep(NA, length.out = length(npi_scalar_to_test))
  for(i in 1:length(npi_scalar_to_test)){
    npi_scaled <- get_npi_values(observations, mobility_df, npi_scalar_to_test[i])
    pps <- SEIRS_force_discrete_both(
      N = length(npi_scaled), parms_rsv = posterior_rsv, parms_mpv = posterior_mpv,
      pop = pop, npi = npi_scaled,
      testing_scalar_mpv = observations %>% filter(pathogen == "mpv") %>% pull(n_tests_ma_scaled),
      testing_scalar_rsv = observations %>% filter(pathogen == "rsv") %>% pull(n_tests_ma_scaled),
      rsvscaling_cuttoff = (end_wk_pre - start_wk_pre)/7 + 1, k = k
    )
    pps <- pps %>% 
      filter(variable == "Ipred") %>%
      left_join(t_to_wk_full, by = join_by(t)) %>% 
      filter(wk_collected > end_wk_pre)
    npi_nll[i] = calculate_nll(pps, observations) %>% pull(nll) %>% sum()
  }
  nll_df = data.frame(npi_scale = npi_scalar_to_test, k = k,
                      nll = npi_nll)
  return(nll_df)
}

plot_full_sim <- function(simulations, observations, plot_samp, npi_df, start_wk_pre, start_wk_post){
  facet_var_lvls = c("S-rsv", "S-mpv", "I-rsv", "I-mpv", "Ipred-rsv", "Ipred-mpv")
  facet_var_labs = c("RSV susceptible", "HMPV susceptible", "RSV infected", "HMPV infected", "RSV detections", "HMPV detections")
  names(facet_var_labs) = facet_var_lvls
  p = simulations %>% 
    filter(draw_id %in% plot_samp, variable == "Ipred") %>%
    mutate(facet_var = factor(paste(variable, pathogen, sep = "-"), levels = facet_var_lvls)) %>%
    ggplot(aes(x = wk_collected, color = pathogen)) + 
    geom_vline(xintercept = start_wk_post) +
    geom_line(aes(y = value, group = draw_id), alpha = 0.2) +
    geom_point(data = observations %>%
                 filter(wk_collected > start_wk_pre) %>%
                 mutate(facet_var = paste("Ipred", pathogen, sep = "-")) %>%
                 mutate(facet_var = factor(facet_var, levels = facet_var_lvls)),
               aes(y = detections), shape = 21, fill = "white", alpha = 0.6, size = 0.9) +
    facet_wrap(vars(facet_var), labeller = labeller(facet_var = facet_var_labs), scales = "free", ncol = 2) + 
    scale_color_brewer(palette = "Dark2", direction = -1) + 
    scale_x_date(date_breaks = "2 years", date_labels = "%Y", expand = c(0,0)) +
    scale_y_continuous(labels = comma) + 
    theme_bw() + 
    theme(axis.title = element_blank(), 
          legend.position = "none", 
          panel.grid.minor.y = element_blank(),
          strip.background = element_blank())
  return(p)
}
