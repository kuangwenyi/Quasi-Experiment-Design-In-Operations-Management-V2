# This replicates Online Supplment 2 in the manuscript #
# Code partially adopted from Callaway and Sant'Anna 2021 and Baker et al. 2022

# load library
pkg = c('tidyverse','RPostgres','fixest','e1071','kableExtra', 
        'ggthemes', 'patchwork','did','furrr', 'latex2exp',
        'bacondecomp','ggforce','fastDummies','progressr')
newpkg <- pkg[!(pkg %in% installed.packages()[,1])]
install.packages(newpkg)

load_pkg = lapply(pkg, library, character.only=TRUE)

# set plot theme
theme_set(theme_clean() + theme(plot.background = element_blank(),
                                legend.background = element_blank()))

# load in compustat data
comp <- read_rds(here::here("Data", "sim_tpt_data.rds"))

# estimate the fixed effects regression of ROA on firm and year fixed effects
mod <- feols(roa ~ 1 | gvkey + fyear, cluster = "incorp", data = comp)
summary(mod)

# get the moments for the residuals from the baseline model
resid_sd <- sd(mod$residuals)
resid_skew <- skewness(mod$residuals)
resid_kurtosis <- kurtosis(mod$residuals)

# get firm and years and state of incorporation
shell <- comp %>% select(gvkey, fyear)

# get the firm and year fes, as well as the standard deviation of ROA
firm_fes <- fixef(mod)$gvkey
n_firm_fes <- length(fixef(mod)$gvkey)
year_fes <- fixef(mod)$fyear
n_year_fes <- length(fixef(mod)$fyear)
sd_roa <- sd(comp$roa)

# function to run simulation, pull firm and year FE, as well as the residuals 
# from their empirical distributions then add in treatment effects

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
                        incorp %in% random_states[1:25] ~ "T",
                        incorp %in% random_states[26:50] ~ "C"),
      delta_base = case_when(fyear >= 1998 & group == "T" ~ 0.05*sd_roa,
                             TRUE ~ 0),
      delta = delta_base * (fyear - 1998 + 1), # cumulative sum for treatment effect
      treat_roa = roa + delta,# treatment effect with delta added
      treat = ifelse(group == "T" & fyear >= 1998, 1, 0),
      rel_year = fyear - 1998) %>%
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
  mod <- feols(treat_roa ~ .[indicators], data = dt,
               cluster = "incorp")
  
  # export result
  broom::tidy(mod) %>% 
    # drop the binned year dummies (<-5 and >5)
    filter(!(term %in% c("Pre", "Post","(Intercept)"))) %>% 
    mutate(t = c(-5:-2, 0:5)) %>% 
    select(t, estimate) %>% 
    bind_rows(tibble(t = -1, estimate = 0)) %>% # t = -1 omitted. Added back
    arrange(t) %>% 
    mutate(sim = i) %>% 
    mutate(true_te = map_dbl(c(-5:5), 
                             function(x) {dt %>% filter(rel_year == x) %>% 
                                 pull(delta) %>% mean}))
}

# parallelize and do 500 simulations
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

class(out)

out %>% 
  group_by(t)# %>% 
  summarize(est = mean(estimate),
            true_effect = mean(true_te),
            lower_ci = quantile(estimate, probs = 0.025),
            upper_ci = quantile(estimate, probs = 0.975))

# plot

p <- out %>% 
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
  ggtitle("Multiple Group Multiple Time Period - Single Treatment:\nDynamic Effect") + 
  ylim(c(-0.02, 0.09)) + 
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(angle = 360, hjust = 0.5, vjust = 0.5))

ggsave(p, filename = here::here("Figs_Tables", "SingleT_Dynamic.png"), 
       dpi = 500,
       width = 6, height = 4)
