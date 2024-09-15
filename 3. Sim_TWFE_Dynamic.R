# This is Section 7.1 Figure 5 in the manuscript #

# Note 1: ggplot code and some simulation code is adapted from:
#         1) Callaway's DID quick guide: https://bcallaway11.github.io/did/articles/TWFE.html
#         2) code from Baker et al. 2022. 
#         We thank Callaway and Baker et al. for making the plotting process easier. 
#         We thank Baker et al,. for sharing the code for knowledge dissemination. 

# load packages
pkgs = c(
  'tidyverse','patchwork','fastDummies','ggthemes','did','bacondecomp',
  'kableExtra','fixest','ggplot2','readxl','readr','tidyr',
  'dplyr','stringr','lme4','RColorBrewer','broom.mixed', 'TwoWayFEWeights', 
  'DIDmultiplegt', 'here')

kwy = lapply(pkgs, library, character.only=TRUE)

# load in compustat data
comp <- read_rds(here::here("Data", "sim_tpt_data.rds"))

# set plot theme
theme_set(theme_clean() + theme(plot.background = element_blank(),
                                legend.background = element_blank()))

# estimate the fixed effects with firm and year fe
# estimate the model without exogenous variables, use "1" as a placeholder
mod <- feols(roa ~ 1 | gvkey + fyear, cluster = "incorp", data = comp)

# get sd skewness and kurtosis from mod
resid_sd <- sd(mod$residuals)
resid_skew <- skewness(mod$residuals)
resid_kurtosis <- kurtosis(mod$residuals)

# create firm year long format data
shell <- comp %>% select(gvkey, fyear)

# get the firm fe, year fe, and sd of ROA
firm_fes <- fixef(mod)$gvkey
n_firm_fes <- length(fixef(mod)$gvkey)
year_fes <- fixef(mod)$fyear
n_year_fes <- length(fixef(mod)$fyear)
sd_roa <- sd(comp$roa)

##### 1.1 No Pre-trends #####

# function to run simulation (firm fe, year fe, residuals)
run_sim <- function(i, p) {
  p()
  # get firm FE from empirical distribution (replace = TRUE: with replacement)
  sim_firm_fe <- tibble(
    gvkey = unique(shell$gvkey),
    firm_fe = sample(firm_fes, n_firm_fes, replace = TRUE),
    incorp = sample(state.abb, n_firm_fes, replace = TRUE)# state.abb: two letter 50 states
  )
  
  # get year FE
  sim_year_fe <- tibble(
    fyear = unique(shell$fyear),
    year_fe = sample(year_fes, n_year_fes, replace = TRUE)
  )
  
  # merge data, create roa
  data <- shell %>% 
    left_join(sim_firm_fe, by = "gvkey") %>% 
    left_join(sim_year_fe, by = "fyear") %>% 
    mutate(resid = sample(mod$residuals, length(mod$residuals), replace = TRUE),
           roa = firm_fe + year_fe + resid)
  
  # randomly assign 50 states (state.abb) into treatment groups
  random_states <- sample(state.abb, length(state.abb), replace = FALSE) 
  
  # add treatment effect -  staggered and Dynamic Treatment Effects
  dt <- data %>%
    mutate( 
      # three groups with exogenous shocks 1997, 2008, 2012
      group = case_when(
        incorp %in% random_states[1:17] ~ 1997,
        incorp %in% random_states[18:35] ~ 2009,
        incorp %in% random_states[36:50] ~ 2015
      ),
      # add in treatment effects - number of standard deviation of ROA added per year
      delta_base = case_when(
        fyear >= group & group == 1997 ~ 0.08*sd_roa,
        fyear >= group & group == 2009 ~ 0.05*sd_roa,
        fyear >= group & group == 2015 ~ 0.03*sd_roa, 
        TRUE ~ 0
      ),
      delta = delta_base * (fyear - group + 1), # cumulative sum for treatment effect
      treat_roa = roa + delta,# treatment effect with delta added
      treat = ifelse(fyear >= group, 1, 0), # treatment group dummy
      rel_year = fyear - group) %>% # relative year for event study
    # create year dummy and bin <-5 and >5 years
    mutate(Pre = ifelse(rel_year < -5, 1, 0),
           `rel_year_-5` = if_else(rel_year == -5, 1, 0),
           `rel_year_-4` = if_else(rel_year == -4, 1, 0),
           `rel_year_-3` = if_else(rel_year == -3, 1, 0),
           `rel_year_-2` = if_else(rel_year == -2, 1, 0),
           rel_year_0 = if_else(rel_year == 0, 1, 0),
           rel_year_1 = if_else(rel_year == 1, 1, 0),
           rel_year_2 = if_else(rel_year == 2, 1, 0),
           rel_year_3 = if_else(rel_year == 3, 1, 0),
           rel_year_4 = if_else(rel_year == 4, 1, 0),
           rel_year_5 = if_else(rel_year == 5, 1, 0),
           Post = ifelse(rel_year > 5, 1, 0))
  
  # put year dummy into vector
  indicators <- c("Pre", paste0("`", "rel_year_", c(-5:-2, 0:5), "`"), "Post")
  
  # estimate model
  mod <- feols(treat_roa ~ .[indicators] | gvkey + fyear, data = dt,
               cluster = "incorp")
  
  # export result
  broom::tidy(mod) %>% 
    # drop the binned year dummies (<-5 and >5)
    filter(!(term %in% c("Pre", "Post"))) %>% 
    mutate(t = c(-5:-2, 0:5)) %>% 
    select(t, estimate) %>% 
    bind_rows(tibble(t = -1, estimate = 0)) %>% # t = -1 omitted. Added back
    arrange(t) %>% 
    mutate(sim = i) %>% 
    mutate(true_te = map_dbl(c(-5:5), 
                             function(x) {dt %>% filter(rel_year == x) %>% 
                                 pull(delta) %>% mean}))
}

# 1000 Simulation
x <- 1:1000
options(future.globals.maxSize= 999999999)
set.seed(20230707)
plan(multisession, workers = 6)
with_progress({
  p <- progressor(steps = length(x))
  out <- future_map_dfr(
    .x = x, 
    .f = run_sim,
    p = p,
    .options = furrr_options(globals = c("mod", "shell", "firm_fes", "n_firm_fes",
                                         "year_fes", "n_year_fes", "sd_roa"),
                             packages = c("tidyverse", "fixest", "fastDummies", "broom"),
                             seed = TRUE)
  )})

# plot
p1
p1 <- out %>% 
  group_by(t) %>% 
  summarize(est = mean(estimate),
            true_effect = mean(true_te),
            lower_ci = quantile(estimate, probs = 0.025),
            upper_ci = quantile(estimate, probs = 0.975)) %>% 
  # split the errors by pre-post
  mutate(band_groups = case_when(
    t < -1 ~ "Pre",
    t >= 0 ~ "Post",
    t == -1 ~ ""
  )) %>%
  # plot
  ggplot(aes(x = t, y = est)) + 
  geom_line(aes(x = t, y = true_effect, color = "True Effect"), 
            linetype = "dashed", linewidth = 2) + 
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
              color = "lightgrey", alpha = 1/4) + 
  geom_pointrange(aes(ymin = lower_ci, ymax = upper_ci, 
                      color = "TWFE Estimated Effect"), show.legend = FALSE) + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = -0.5, linetype = "dashed") + 
  scale_x_continuous(breaks = -5:5) + 
  labs(x = "Year Relative to Exogenous Shock", y = "Estimate") +
  scale_color_manual(values = c("#0029a5","#993441")) + 
  ggtitle("TWFE Estimates") + 
  ylim(c(-0.02, 0.08)) + 
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(angle = 360, hjust = 0.5, vjust = 0.5))

##### 1.2 Actual Pre-trends #####
# Opposite scenario - real pre-trends

run_sim_2 <- function(i, p) {
  
  p()
  sim_firm_fe <- tibble(
    gvkey = unique(shell$gvkey),
    firm_fe = sample(firm_fes, n_firm_fes, replace = TRUE),
    incorp = sample(state.abb, n_firm_fes, replace = TRUE)
  )
  
  sim_year_fe <- tibble(
    fyear = unique(shell$fyear),
    year_fe = sample(year_fes, n_year_fes, replace = TRUE)
  )
  data <- shell %>% 
    left_join(sim_firm_fe, by = "gvkey") %>% 
    left_join(sim_year_fe, by = "fyear") %>% 
    mutate(resid = sample(mod$residuals, length(mod$residuals), replace = TRUE),
           roa = firm_fe + year_fe + resid)
  random_states <- sample(state.abb, length(state.abb), replace = FALSE) 
  
  dt <- data %>%
    mutate(
      group = case_when(
        incorp %in% random_states[1:17] ~ 1997,
        incorp %in% random_states[18:35] ~ 2009,
        incorp %in% random_states[36:50] ~ 2015
      ),
      delta_base = case_when(
        fyear < group & group == 1997 ~ -0.01*sd_roa,
        fyear < group & group == 2009 ~ -0.02*sd_roa,
        fyear < group & group == 2015 ~ -0.03*sd_roa, 
        TRUE ~ 0
      ),
      delta = delta_base * (group - fyear),
      treat_roa = roa + delta,
      treat = ifelse(fyear >= group, 1, 0),
      rel_year = fyear - group) %>% 
    mutate(Pre = ifelse(rel_year < -5, 1, 0),
           `rel_year_-5` = if_else(rel_year == -5, 1, 0),
           `rel_year_-4` = if_else(rel_year == -4, 1, 0),
           `rel_year_-3` = if_else(rel_year == -3, 1, 0),
           `rel_year_-2` = if_else(rel_year == -2, 1, 0),
           rel_year_0 = if_else(rel_year == 0, 1, 0),
           rel_year_1 = if_else(rel_year == 1, 1, 0),
           rel_year_2 = if_else(rel_year == 2, 1, 0),
           rel_year_3 = if_else(rel_year == 3, 1, 0),
           rel_year_4 = if_else(rel_year == 4, 1, 0),
           rel_year_5 = if_else(rel_year == 5, 1, 0),
           Post = ifelse(rel_year > 5, 1, 0))
  
  indicators <- c("Pre", paste0("`", "rel_year_", c(-5:-2, 0:5), "`"), "Post")
  
  mod <- feols(treat_roa ~ .[indicators] | gvkey + fyear, data = dt, cluster = "incorp")
  
  broom::tidy(mod) %>% 
    filter(!(term %in% c("Pre", "Post"))) %>% 
    mutate(t = c(-5:-2, 0:5)) %>% 
    select(t, estimate) %>% 
    bind_rows(tibble(t = -1, estimate = 0)) %>% 
    arrange(t) %>% 
    mutate(sim = i) %>% 
    mutate(true_te = map_dbl(c(-5:5), 
                             function(x) {dt %>% filter(rel_year == x) %>% pull(delta) %>% mean}))
}

# do 1000 simulations
x <- 1:1000
options(future.globals.maxSize= 999999999)
set.seed(20230707)
plan(multisession, workers = 6)
with_progress({
  p <- progressor(steps = length(x))
  out <- future_map_dfr(
    .x = x, 
    .f = run_sim_2,
    p = p,
    .options = furrr_options(globals = c("mod", "shell", "firm_fes", "n_firm_fes",
                                         "year_fes", "n_year_fes", "sd_roa"),
                             packages = c("tidyverse", "fixest", "fastDummies", "broom"), 
                             seed = TRUE)
  )})

# plot
p2 <- out %>% 
  group_by(t) %>% 
  summarize(est = mean(estimate),
            true_effect = mean(true_te),
            lower_ci = quantile(estimate, probs = 0.025),
            upper_ci = quantile(estimate, probs = 0.975)) %>% 
  mutate(band_groups = case_when(
    t < -1 ~ "Pre",
    t >= 0 ~ "Post",
    t == -1 ~ ""
  )) %>%
  
  ggplot(aes(x = t, y = est)) + 
  geom_line(aes(x = t, y = true_effect, color = "True Effect"), 
            linetype = "dashed", size = 2) + 
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
              color = "lightgrey", alpha = 1/4) + 
  geom_pointrange(aes(ymin = lower_ci, ymax = upper_ci, 
                      color = "TWFE Estimated Effect"), show.legend = FALSE) + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = -0.5, linetype = "dashed") + 
  scale_x_continuous(breaks = -5:5) + 
  labs(x = "Year Relative to Exogenous Shock", y = "Estimate") +
  scale_color_manual(values = c("#0029a5","#993441")) + 
  ggtitle("With Pre-Trends") + 
  ylim(c(-0.04, 0.02)) + 
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(angle = 360, hjust = 0.5, vjust = 0.5))

##### 1.3 Combine plots #####
plot <- p1 + p2
plot

ggsave(plot, filename = here::here("Figs_Tables", "bad_es_tpt_4849.png"), 
       dpi = 500,width = 8, height = 4)


ggsave(p1, filename = here::here("Figs_Tables", "TWFE_NoPretrend.png"), 
       dpi = 500,width = 6, height = 4.5)


