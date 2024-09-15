# This replicates Online Supplement 1 in the Manuscript #
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
  
  sim_mod <- feols(roa ~ 1 | gvkey + fyear, cluster = "incorp", data = data)
  mom <- c(sd(sim_mod$residuals), skewness(sim_mod$residuals), kurtosis(sim_mod$residuals))
  
  random_states <- sample(state.abb, length(state.abb), replace = FALSE) 
  
  # Scenario 1: One Treatment Period LEVEL Shift of effect
  data1 <- data %>% 
    mutate(
      group = case_when(
        incorp %in% random_states[1:25] ~ "T",
        incorp %in% random_states[26:50] ~ "C"),
      # treatment effects - constant half of a standard deviation of ROA
      delta = case_when(fyear >= 1998 & group == "T" ~ 0.5*sd_roa,
                        TRUE ~ 0),
      treat_roa = roa + delta,
      treat = ifelse(group == "T" & fyear >= 1998, 1, 0),
      rel_year = fyear - 1998)
  
  # Scenario 2: One treatment but treatment effect increasing over time
  data2 <- data %>% 
    mutate(
      group = case_when(
        incorp %in% random_states[1:25] ~ "T",
        incorp %in% random_states[26:50] ~ "C"),
      delta_base = case_when(fyear >= 2002 & group == "T" ~ 0.05*sd_roa,
                             TRUE ~ 0),
      delta = delta_base * (fyear - 2002 + 1),
      treat_roa = roa + delta,
      treat = ifelse(group == "T" & fyear >= 2002, 1, 0),
      rel_year = fyear - 2002)
  
  
  # make function to get estimates and treatment effects from data for k in different scenarios
  get_est <- function(k) {
        dt <- get(paste0("data", k))

    # full treatment effect
    # observation level
    full_te_1 <- dt %>% filter(treat == 1) %>% pull(delta) %>% mean()
    # firm level average
    full_te_2 <- dt %>% filter(treat == 1) %>% group_by(gvkey) %>% 
      summarize(m = mean(delta)) %>% pull(m) %>% mean()

    
    # get twfe estimates on full data
    twfe <- feols(treat_roa ~ treat | gvkey + fyear, cluster = "incorp", data = dt)$coefficients[1]
    
    tibble(sim = i, dt = k, 
           full_te_1 = full_te_1, 
           full_te_2 = full_te_2, 
           twfe = twfe)
  }
  
  # run it for our six simulations and store results in a dataset
  estimates <- map_dfr(1:2, get_est)
  
  # get moments into a tibble as well
  moments <- tibble(
    moment = 1:3,
    value = mom
  ) %>% 
    mutate(sim = i)
  
  # output a list of both objects that we want
  list(moments = moments,
       estimates = estimates)
  
}

# run 1000 simulations
x <- 1:1000
options(future.globals.maxSize= 999999999)
set.seed(20230707)
plan(multisession, workers = 6)
with_progress({
  p <- progressor(steps = length(x))
  out <- future_map(
    .x = x, 
    .f = run_sim,
    p = p,
    .options = furrr_options(globals = c("mod", "shell", "firm_fes", "n_firm_fes",
                                         "year_fes", "n_year_fes", "sd_roa"),
                             packages = c("tidyverse", "fixest", "e1071", "did", "fastDummies"),
                             seed = TRUE)
  )})

# unpack the two datasets
moments <- do.call(function(...) mapply(bind_rows, ..., SIMPLIFY=F), args = out)$moments
estimates <- do.call(function(...) mapply(bind_rows, ..., SIMPLIFY=F), args = out)$estimates


# Trend of the mean for control and treatment group
# pull firm FE from empirical distribution with replacement
set.seed(20230707)
sim_firm_fe <- tibble(
  gvkey = unique(shell$gvkey),
  firm_fe = sample(firm_fes, n_firm_fes, replace = TRUE),
  incorp = sample(state.abb, n_firm_fes, replace = TRUE)
)

# pull year FE from the empirical distribution with replacement
sim_year_fe <- tibble(
  fyear = unique(comp$fyear),
  year_fe = sample(year_fes, n_year_fes, replace = TRUE)
)

# merge in the FE to the firm/year/state observations and add in residuals from the 
# empirical distribution. ROA is the linear combination of the FEs and the residual
data <- shell %>% 
  left_join(sim_firm_fe, by = "gvkey") %>% 
  left_join(sim_year_fe, by = "fyear") %>% 
  mutate(resid = sample(mod$residuals, length(mod$residuals), replace = TRUE),
         roa = firm_fe + year_fe + resid)

# randomly assign the states into treatment groups
random_states <- sample(state.abb, length(state.abb), replace = FALSE)

#  One Treatment Period, Constant Treatment Effects
newdt <- data %>% 
  mutate(
    group = case_when(
      incorp %in% random_states[1:25] ~ "Treatment",
      incorp %in% random_states[26:50] ~ "Control"),
    delta_base = case_when(fyear >= 2001 & group == "Treatment" ~ 0.05*sd_roa,
                           TRUE ~ 0),
    delta = delta_base * (fyear - 2001 + 1),
    treat_roa = roa + delta,
    treat = ifelse(group == "Treatment" & fyear >= 2001, 1, 0))


# plot for trend
roa_means <- newdt %>% 
  ggplot(aes(x = fyear, y = treat_roa, group = gvkey)) +
  # unit specific lines
  geom_line(alpha = 1/30, color = "black") + 
  # group specific averages
  geom_line(
    data = . %>% 
      group_by(group, fyear) %>% 
      summarize(treat_roa = mean(treat_roa)),
    aes(x = fyear, y = treat_roa, group = factor(group),
        color = factor(group)), linewidth = 1) + 
  labs(x = "", y = "ROA", color = "Group") + 
  geom_vline(xintercept = 2001.5, color = '#0029a5',
             linetype = "dashed", linewidth = 1) + 
  scale_y_continuous(limits = c(-1*sd_roa, 3*sd_roa)) + 
  scale_color_manual(values = c("#0029a5","#993441")) + 
  ggtitle("Multiple Group Multiple Time Period - Single Treatment\nStatic Effect") +  
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.y = element_text(hjust = 0.5, vjust = 0.5, angle = 360))

# check plot
roa_means

#save
ggsave(roa_means, filename = here::here("Figs_Tables", "SingleT_Trend.png"), 
       dpi = 500,
       width = 6, height = 4)


# function to make distribution plot
make_dist_plot <- function(i, name) {
  estimates %>% 
    filter(dt == i) %>% 
    ggplot(aes(x = twfe)) + 
    geom_density(aes(fill = "TWFE Estimates"), alpha = 3/5) + 
    geom_vline(aes(xintercept = mean(estimates %>% filter(dt == i) %>% pull(full_te_1)),
                   color = "Observation Average"),
               linetype = "dashed", size = 1, alpha = 3/5,) +
    geom_vline(aes(xintercept = mean(estimates %>% filter(dt == i) %>% pull(full_te_2)),
                   color = "Firm Average"),
               linetype = "dashed", size = 1, alpha = 3/5) + 
    ggtitle(paste0("Multiple Group Multiple Time Period - Single Treatment\nStatic Effect")) + 
    labs(y = "", x = "") + 
    scale_color_manual(name = "", values = c("Observation Average" = "#0029a5",
                                             "Firm Average" = "#993441")) + 
    scale_fill_manual(name = "", values = c("TWFE Estimates" = "#999999")) + 
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          axis.title.y = element_text(hjust = 0.5, vjust = 0.5, angle = 360),
          legend.position =  "bottom")
}

# run function 
sim1_estimates <- make_dist_plot(2) # dt==2 from estimates 
sim1_estimates

# save
ggsave(sim1_estimates, filename = here::here("Figs_Tables", "SingleT_Static.png"), 
       dpi = 500,
       width = 6, height = 4)











