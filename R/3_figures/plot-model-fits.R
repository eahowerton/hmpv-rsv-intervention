library(ggplot2)
library(RColorBrewer)

# load helper functions
source("R/utils.R")

fit_scotland_SIRS <- readRDS("data/derived_data/scotland/fit_scotland_SIRS.rds")
fit_scotland_SEIRS <- readRDS("data/derived_data/scotland/fit_scotland_SEIRS.rds")
fit_scotland_SEIRS_rsv <- readRDS("data/derived_data/scotland/fit_scotland_SEIRS_rsv.rds")
fit_scotland_SEIRS_mpv <- readRDS("data/derived_data/scotland/fit_scotland_SEIRS_mpv.rds")

# read in scotland data, all age groups, summarized by week
scotland_by_wk <- readRDS("data/derived_data/scotland/scotland_by_wk_overall.rds")

# load scotland base parameters
source("R/scotland_parameters.R")

scotland_by_wk_pre <- scotland_by_wk %>%
  filter(wk_collected >= start_wk_pre, 
         wk_collected <= end_wk_pre) 

n_samp <- dim(as.array(fit_scotland_SIRS))[1]

# PREDICTIONS VS. DETECTIONS
bind_rows(
  filter_posterior(as.data.frame(fit_scotland_SIRS), "Ipred_mpv\\[") %>%
    mutate(pathogen = "mpv"), 
  filter_posterior(as.array(fit_scotland_SIRS), "Ipred_rsv\\[") %>%
    mutate(pathogen = "rsv")
) %>%
  left_join(t_to_wk_scotland) %>%
  filter(draw_id %in% sample(1:n_samp,100)) %>%
  ggplot(aes(x = wk_collected, color = pathogen)) + 
  geom_point(data = scotland_by_wk_pre, 
             aes(y = detections), shape = 21) + 
  geom_line(data = scotland_by_wk_pre, 
            aes(y = detections), alpha = 0.2) + 
  geom_line(aes(y = value, group = draw_id), alpha = 0.2) +
  facet_grid(rows = vars(pathogen), scales = "free", switch = "y") + 
  scale_color_brewer(palette = "Dark2") + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        legend.position = "none", 
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(), 
        strip.placement = "outside")
ggsave("figures/fits_SIRS_scotland.pdf", width = 12, height = 5)

bind_rows(
  filter_posterior(as.data.frame(fit_scotland_SEIRS), "Ipred_mpv\\[") %>%
    mutate(pathogen = "mpv"), 
  filter_posterior(as.data.frame(fit_scotland_SEIRS), "Ipred_rsv\\[") %>%
    mutate(pathogen = "rsv")
) %>%
  left_join(t_to_wk_scotland) %>%
  filter(draw_id %in% sample(1:n_samp,100)) %>%
  ggplot(aes(x = wk_collected, color = pathogen)) + 
  geom_point(data = scotland_by_wk_pre, 
             aes(y = detections), shape = 21) + 
  geom_line(data = scotland_by_wk_pre, 
            aes(y = detections), alpha = 0.2) + 
  geom_line(aes(y = value, group = draw_id), alpha = 0.2) +
  facet_grid(rows = vars(pathogen), scales = "free", switch = "y") + 
  scale_color_brewer(palette = "Dark2") + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        legend.position = "none", 
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(), 
        strip.placement = "outside")
ggsave("figures/fits_SEIRS_scotland.pdf", width = 12, height = 5)


#### PLOT POSTERIOR ESTIMATES --------------------------------------------------
param_groups = data.frame(parameter = colnames(as.data.frame(fit_scotland_SEIRS)[,1:14]), 
                         group = c(1:3, 1:3, 4, 4, 5, 6, 7, 8, 9, 10))

ind_rsv_pars <- as.data.frame(apply(as.data.frame(fit_scotland_SEIRS_rsv)[,1:7], 2, function(i){quantile(i, c(0.05, 0.25, 0.5, 0.75, 0.95))}))
colnames(ind_rsv_pars) <- c("S0_rsv", "E0_rsv", "I0_rsv", "rho_rsv", "beta_rsv", "a", "p")
ind_mpv_pars <- as.data.frame(apply(as.data.frame(fit_scotland_SEIRS_mpv)[,1:7], 2, function(i){quantile(i, c(0.05, 0.25, 0.5, 0.75, 0.95))}))
colnames(ind_mpv_pars) <- c("S0_mpv", "E0_mpv", "I0_mpv", "rho_mpv", "beta_mpv", "a", "p")


fit_scotland_SIRS_long = as.data.frame(fit_scotland_SIRS) %>% 
  mutate(beta_rsv = b, beta_mpv = b*r)
fit_scotland_SEIRS_long = as.data.frame(fit_scotland_SEIRS) %>% 
  mutate(beta_rsv = b, beta_mpv = b*r)

pars = c("beta_rsv", "beta_mpv", "p", "a", "r", "c", "rho_rsv", "rho_mpv", "S0_rsv", "I0_rsv", "S0_mpv", "I0_mpv", "E0_rsv", "E0_mpv")

par_labs = c("RSV transmission rate", "HMPV transmission rate", "seasonal phase", "seasonal amplitude", "HMPV relative transmissibility", 
             "effect of RSV on HMPV transmission", "RSV reporting rate", "HMPV reporting rate", "RSV initial susceptible", "RSV initial infected", "HMPV initial susceptible", "HMPV initial infected", "RSV initial exposed", "HMPV initial exposed")
names(par_labs) = pars

bind_rows(
  as.data.frame(apply(fit_scotland_SIRS_long[,pars[1:12]], 2, function(i){quantile(i, c(0.05, 0.25, 0.5, 0.75, 0.95))})) %>%
    mutate(quantile = paste0("Q", c(5, 25, 50, 75, 95)), model = "SIRS-both"),
  as.data.frame(apply(fit_scotland_SEIRS_long[,pars], 2, function(i){quantile(i, c(0.05, 0.25, 0.5, 0.75, 0.95))})) %>%
    mutate(quantile = paste0("Q", c(5, 25, 50, 75, 95)), model = "SEIRS-both"), 
  ind_rsv_pars %>%
    mutate(quantile = paste0("Q", c(5, 25, 50, 75, 95)), model = "SEIRS-RSV"),
  ind_mpv_pars %>%
    mutate(quantile = paste0("Q", c(5, 25, 50, 75, 95)), model = "SEIRS-HMPV")
) %>%
  melt(c("quantile", "model"), variable.name = "parameter") %>%
  dcast(parameter + model ~ quantile) %>%
  mutate(parameter = factor(parameter, levels = c("beta_rsv", "rho_rsv", "S0_rsv", "E0_rsv", "I0_rsv",
                                                  "beta_mpv", "rho_mpv", "S0_mpv", "E0_mpv", "I0_mpv",
                                                  "a", "p",  "r", "c")), 
         model = factor(model, levels = c("SEIRS-both", "SIRS-both", "SEIRS-RSV", "SEIRS-HMPV"))) %>%
  filter(!is.na(Q50)) %>%
  ggplot() + 
  geom_segment(aes(x = model, xend = model, y = Q5, yend = Q95)) +
  geom_segment(aes(x = model, xend = model, y = Q25, yend = Q75), size = 1.5) +
  geom_point(aes(x = model, y = Q50), size = 2.25) + 
  facet_wrap(vars(parameter), scales = "free", ncol = 5, labeller = labeller(parameter = par_labs)) + 
  labs(y = "posterior distribution") +
  # scale_x_discrete(expand = c(0.1, 0.1)) + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        legend.position = "bottom", 
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
  
ggsave("figures/params_SIRS-SEIRS_scotland.pdf", width = 15, height = 7)

