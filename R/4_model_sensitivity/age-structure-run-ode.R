library(deSolve)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidytable)
library(epimdr2)
library(fields)
library(cowplot)
library(rstan)
library(readr)
library(dplyr)
library(bayesplot)

source("R/4_model-sensitivity/age-structure-utils.R")
source("R/scotland_parameters.R")


#### MODEL SETUP ---------------------------------------------------------------
compartments = c("S", "E", "I", "R")
start_pop = scotland_N
fert = (1-exp(-scotland_birthrate))
mort = (1-exp(-scotland_birthrate))

# load model fit
fit_scotland_SEIRS <- readRDS("data/derived_data/scotland/fit_scotland_SEIRS.rds")

# extract median parameters
fitted_pars = c(unlist(lapply(extract(fit_scotland_SEIRS, c("b", "r", "a", "p")), median)),
                c_prime = median(unlist(extract(fit_scotland_SEIRS, "c"))/unlist(extract(fit_scotland_SEIRS, "max_RSV_inc"))*scotland_N))
  
params_rsv = c(
  b = as.double(fitted_pars["b"]), r = 1, p = as.double(fitted_pars["p"]), a = as.double(fitted_pars["a"]), c = 0, 
  sigma = (1-exp(-scotland_sigma)), gamma = (1-exp(-scotland_gamma)), omega = (1-exp(-scotland_omega)))

params_mpv = c(
  b = as.double(fitted_pars["b"]), r = as.double(fitted_pars["r"]), p = as.double(fitted_pars["p"]), a = as.double(fitted_pars["a"]), 
  c = as.double(fitted_pars["c_prime"]), sigma = (1-exp(-scotland_sigma)), 
  gamma = (1-exp(-scotland_gamma)), omega = (1-exp(-scotland_omega)))

params_mpv_noint = c(
  b = as.double(fitted_pars["b"]), r = as.double(fitted_pars["r"]), 
  p = as.double(fitted_pars["p"]), a = as.double(fitted_pars["a"]), 
  c = 0, sigma = (1-exp(-scotland_sigma)), 
  gamma = (1-exp(-scotland_gamma)), omega = (1-exp(-scotland_omega)))

#### TEST WITH ONE COMPARTMENT -------------------------------------------------
age_classes1 =  c(100)
rslts1 <- run_ode(
  age_classes = age_classes1, mort = mort, fert = fert, 
  start_pop = start_pop, IC_type = "std", max_t = 20*52, 
  params_rsv = params_rsv, params_mpv = params_mpv,
  plot_flag = TRUE, plot_title = "1 age class")

#### TEST WITH TWO COMPARTMENTS ------------------------------------------------
age_classes2 = c(20*52, 40*52)
rslts2 <- run_ode(
  age_classes = age_classes2, mort = rep(mort, length(age_classes2)),
  fert = rep(fert, length(age_classes2)), start_pop = start_pop,
  IC_type = "std", max_t = 20*52, params_rsv = params_rsv, params_mpv = params_mpv,
  plot_flag = TRUE, plot_title = "2 age classes")

#### TEST WITH REALISTIC STRUCTURE ---------------------------------------------
# up to 5 years in months, up to 10 in years, and up to 70 in 5 years 
# (multiply by 4 to go from months to weeks)
age_classes = c(1:60, seq(72,120,by=12), seq(180,840,by=60))*4 
bin_width = diff(c(0, age_classes))
names(bin_width) = age_classes

start.time <- Sys.time()
rslts3 <- run_ode(
  age_classes = age_classes, mort = rep(mort, length(age_classes)),
  fert = rep(fert, length(age_classes)), start_pop = start_pop,
  IC_type = "std", max_t = 20*52, params_rsv = params_rsv, params_mpv = params_mpv,
  plot_flag = TRUE, plot_title = "realistic structure, constant WAIFW")
Sys.time() - start.time

rslts3_nointeraction <- run_ode(
  age_classes = age_classes, mort = rep(mort, length(age_classes)),
  fert = rep(fert, length(age_classes)), start_pop = start_pop,
  IC_type = "std", max_t = 20*52, params_rsv = params_rsv, params_mpv = params_mpv_noint,
  plot_flag = TRUE, plot_title = "realistic structure, constant WAIFW, no interaction")

#### TEST WITH POLYMOD CONTACTS ------------------------------------------------
# now add more realistic mixing with POLYMOD data
W = create_polymod_matrix(age_classes)

start.time <- Sys.time()
rslts4 <- run_ode(
  age_classes = age_classes, mort = rep(mort, length(age_classes)),
  fert = rep(fert, length(age_classes)), start_pop = start_pop, waifw = W,
  IC_type = "std", max_t = 20*52, params_rsv = params_rsv, params_mpv = params_mpv,
  adjust_beta_flag = TRUE, plot_flag = TRUE, plot_title = "realistic structure, POLYMOD")
Sys.time() - start.time


rslts4_nointeraction <- run_ode(
  age_classes = age_classes, mort = rep(mort, length(age_classes)),
  fert = rep(fert, length(age_classes)), start_pop = start_pop, waifw = W,
  IC_type = "std", max_t = 20*52, params_rsv = params_rsv, params_mpv = params_mpv_noint, 
  adjust_beta_flag = TRUE, plot_flag = TRUE, plot_title = "realistic structure, POLYMOD, no interaction")

# verify that beta adjustment gives the same results with and without interaction
all_results <- bind_rows(
  rslts3 %>% mutate(mixing = "constant", interaction = TRUE),
  rslts3_nointeraction  %>% mutate(mixing = "constant", interaction = FALSE),
  rslts4 %>% mutate(mixing = "polymod", interaction = TRUE),
  rslts4_nointeraction %>% mutate(mixing = "polymod", interaction = FALSE)
)
  
all_results %>%
  filter(variable == "I", time > 15*52) %>%
  summarize(value = sum(value), .by = c("time", "pathogen", "mixing", "interaction")) %>% 
  ggplot(aes(x = time, y = value, color = pathogen, linetype = mixing)) + 
  geom_line() + 
  theme_bw() + 
  facet_wrap(vars(interaction))

#### COMPARE MEAN AGES OF INFECTION --------------------------------------------
mean_age = all_results %>% 
  filter(time > 15*52) %>%
  summarize(value = sum(value), .by = c("variable", "age", "pathogen", "mixing", "interaction")) %>%
  mutate(age = as.double(age)) %>%
  left_join(data.frame(age = age_classes, bin_width = diff(c(0, age_classes)))) %>% 
  filter( variable == "I", age <= (3*52)) %>% 
  summarize(mean_age = sum(value*age)/ sum(value), .by = c("variable", "pathogen", "mixing", "interaction"))

txt = mean_age %>% 
  filter(paste(pathogen, interaction) != "rsv FALSE") %>% 
  dcast(variable + mixing ~ interaction + pathogen, value.var = "mean_age") %>%
  mutate(lab = paste0("\nRSV: ", round(TRUE_rsv/52, 2), #mean age of infection\n
                         " years \nhMPV without interaction: ", round(FALSE_mpv/52, 2), 
                      " years \nhMPV with interaction: ", round(TRUE_mpv/52, 2), " years "))


p1 = all_results %>%
  filter(variable == "I", time > 15*52, interaction) %>%
  summarize(value = sum(value), .by = c("time", "pathogen", "mixing")) %>% 
  mutate(pathogen = factor(pathogen, levels = c("rsv", "mpv")), 
         r_mixing = paste(mixing, "mixing")) %>%
  ggplot(aes(x = time/52, y = value/scotland_N, color = pathogen, linetype = r_mixing)) + 
  geom_line(size = 0.8, alpha = 0.7) + 
  labs(x = "years since start of simulation", y = "all infected") +
  scale_color_brewer(palette = "Dark2", labels = c("RSV", "hMPV")) + 
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        panel.grid.minor = element_blank())

p2 = all_results %>% 
  mutate(age = as.integer(age), 
         pathogen = factor(pathogen, levels = c("rsv", "mpv"))) %>%
  filter(time > 15*52, age <= (3*52), variable == "I", paste(pathogen, interaction) != "rsv FALSE") %>%
  summarize(value = sum(value), .by = c("variable", "age", "pathogen", "mixing", "interaction")) %>%
  mutate(value_rel = value/sum(value), .by = c("variable", "pathogen", "mixing", "interaction")) %>%
  mutate(interaction = ifelse(interaction, "with interaction", "without interaction")) %>%
  ggplot(aes(x = age/52, y = value_rel, color = pathogen)) + 
  geom_line(aes(linetype = interaction), size = 0.8) + 
  geom_text(data = txt, aes(x = Inf, y = Inf, label = lab), hjust = 1, vjust = 1, color = "black", size = 2.5) + 
  facet_grid(cols = vars(paste(mixing, "mixing"))) +
  scale_color_brewer(palette = "Dark2", labels = c("RSV", "hMPV")) + 
  scale_linetype_manual(values = c("solid", "dashed")) + 
  scale_x_continuous(name = c("age (years)")) +
  scale_y_continuous(name = c("% of infections under 3 yrs old"), 
                     labels = scales::percent) +
  theme_bw() + 
  theme(legend.key.width = unit(1, "cm"),
    legend.position = "bottom", 
        legend.title = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.background = element_blank())

plot_grid(p1, p2,  ncol = 1, labels = LETTERS[1:2], rel_heights = c(0.43, 0.57))
ggsave("figures/mean_age_of_infection.pdf", width = 6, height = 6)

