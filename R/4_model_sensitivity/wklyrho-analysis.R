library(rstan)
library(readr)
library(dplyr)
library(bayesplot)
library(cowplot)
library(ggplot2)

#### PREP DATA
source("R/scotland_parameters.R")

# read in scotland data, all age groups, summarized by week
scotland_by_wk <- read_rds("data/derived_data/scotland/scotland_by_wk_overall.rds")

scotland_by_wk_pre <- scotland_by_wk %>%
  filter(wk_collected >= start_wk_pre, 
         wk_collected <= end_wk_pre) 

#### FIT COUPLED SEIRS MODEL ---------------------------------------------------
# assumptions: 
# 1. same waning rate for RSV and HMPV
# 2. same seasonal forcing (a and p)
# 3. HMPV transmission rate differs by some factor r, and is also forced by
#    observed RSV incidence by some factor c
# 4. reporting rate accounts for weekly average and annual trends in testnig
model_SEIRS_wklyrho <- stan_model("R/1_fit_models/wklyrho-SEIRS_coupled.stan")

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

pathogen_labs = c("RSV", "hMPV")
names(pathogen_labs) = c("rsv", "mpv")

#### CALCULATE WEEKLY TESTING --------------------------------------------------
rolling_avg = function(wk, value, lag){
  n_wks = length(value)
  val_vec = rep(value,2)
  wk_vec = rep(wk, 2)
  mean_vec = roll_mean(val_vec, lag, align = "center", fill = NA)
  min_indx = min(which(!is.na(mean_vec)))
  return(data.frame(wk = wk_vec[min_indx:(min_indx+n_wks-1)], 
                    mean = mean_vec[min_indx:(min_indx+n_wks-1)]))
}

wkly_testing_force = scotland_by_wk_pre %>% 
  mutate(n_tests_scaled = n_tests/sum(n_tests), .by = c("year", "pathogen")) %>%
  filter(n_tests_scaled > 0) %>% 
  summarize(mean = mean(n_tests_scaled, na.rm = TRUE), .by = c("wk","pathogen")) %>%
  arrange(wk+10) %>%
  reframe(rolling_avg(wk = wk, value = mean, lag = 8), .by = c("pathogen")) %>%
  arrange(wk)
wkly_testing_force = scotland_by_wk_pre %>% 
  select(wk_collected, wk, pathogen, n_tests_ma_scaled) %>%
  left_join(wkly_testing_force %>% rename(wkly_tests = mean)) %>%
  mutate(rr = n_tests_ma_scaled * wkly_tests) %>%
  # scale for max at 1
  mutate(rr = rr/max(rr), .by = c("pathogen"))

p1 = scotland_by_wk_pre %>% 
  mutate(year = year(wk_collected)) %>%
  filter(year < 2018, year > 2006) %>%
  mutate(n_tests_scaled = n_tests/sum(n_tests), .by = c("year", "pathogen")) %>%
  filter(n_tests_scaled > 0) %>%
  ggplot(aes(x = wk, y = n_tests_scaled, color = as.factor(year))) + 
  geom_line() + 
  geom_line(data = scotland_by_wk_pre %>% 
              mutate(n_tests_scaled = n_tests/sum(n_tests), .by = c("year", "pathogen")) %>%
              filter(n_tests_scaled > 0) %>% 
              summarize(mean = mean(n_tests_scaled, na.rm = TRUE), .by = c("wk","pathogen")) %>%
              arrange(wk+10) %>%
              reframe(rolling_avg(wk = wk, value = mean, lag = 8), .by = c("pathogen")), 
            aes(y = mean), color = "black", linewidth = 1) +
  facet_wrap(vars(pathogen), labeller = labeller(pathogen = pathogen_labs), ncol = 1) + 
  labs(x = "week of year", y = "% of tests per year") + 
  scale_y_continuous(label = percent) + 
  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.position = "none", 
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
        )

p2 = wkly_testing_force %>%
  ggplot(aes(x = wk_collected, y = rr)) + 
  geom_line() + 
  labs(y = "testing scalar") + 
  facet_wrap(vars(pathogen), ncol = 1, scales = "free", labeller = labeller(pathogen = pathogen_labs)) + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank())


# run model
fit_scotland_SEIRS <- sampling(
  model_SEIRS_wklyrho,
  data = list(
    N = nrow(scotland_by_wk_pre %>% filter(pathogen == "mpv")),
    cases_mpv = scotland_by_wk_pre %>% filter(pathogen == "mpv") %>% pull(detections), 
    cases_rsv = scotland_by_wk_pre %>% filter(pathogen == "rsv") %>% pull(detections), 
    birth = scotland_by_wk_pre %>% filter(pathogen == "mpv") %>% pull(Births.registered), 
    tests_mpv = wkly_testing_force %>% filter(pathogen == "mpv") %>% pull(rr),
    tests_rsv = wkly_testing_force %>% filter(pathogen == "rsv") %>% pull(rr),
    pop = scotland_N,
    sigma = scotland_sigma,
    gamma = scotland_gamma,
    mu = scotland_birthrate,
    omega = scotland_omega
  ),
  seed = 7, 
  iter = 15000,
  chain = 4, 
  init = initfun, 
  cores = 4,
  control=list(max_treedepth = 12)
)


p3 = bind_rows(
  filter_posterior(as.data.frame(fit_scotland_SEIRS), "Ipred_mpv\\[") %>%
    mutate(pathogen = "mpv"), 
  filter_posterior(as.data.frame(fit_scotland_SEIRS), "Ipred_rsv\\[") %>%
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
  facet_wrap(vars(pathogen), scales = "free",
             labeller = labeller(pathogen = pathogen_labs), ncol = 1) + 
  scale_color_brewer(palette = "Dark2") + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        legend.position = "none", 
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank())


pars = c("b", "r", "c", "a", "p", "rho_rsv", "rho_mpv")

par_labs = c("RSV transmission\nrate", "hMPV relative\ntransmissibility", 
             "effect of RSV\non HMPV transmission","seasonal\namplitude",
             "seasonal phase",  
              "RSV\nreporting rate", "hMPV\nreporting rate")
names(par_labs) = pars


p4 = as.data.frame(apply(as.data.frame(fit_scotland_SEIRS)[,pars], 2, function(i){quantile(i, c(0.05, 0.25, 0.5, 0.75, 0.95))})) %>%
  mutate(quantile = paste0("Q", c(5, 25, 50, 75, 95))) %>%
  melt(c("quantile"), variable.name = "parameter") %>%
  dcast(parameter ~ quantile) %>%
  ggplot() + 
  geom_segment(aes(x = 1, xend = 1, y = Q5, yend = Q95)) +
  geom_segment(aes(x = 1, xend = 1, y = Q25, yend = Q75), size = 1.5) +
  geom_point(aes(x = 1, y = Q50), size = 2.25) + 
  facet_wrap(vars(parameter), scales = "free", labeller = labeller(parameter = par_labs), nrow = 2) + 
  labs(y = "posterior distribution") +
  # scale_x_discrete(expand = c(0.1, 0.1)) + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.position = "bottom", 
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

plot_grid(
  plot_grid(
    plot_grid(p1, p2, nrow = 1, labels = LETTERS[1:2]), 
    p4, ncol = 1, labels = c(NA, "C"), rel_heights = c(0.6,0.4)), 
  p3, nrow = 1, labels = c(NA, "D"), rel_widths = c(0.6, 0.4)
)
ggsave("figures/wkly-reporting-summary.pdf", width = 12, height = 6)

