# This is Section 7.1 Figure 6 in the manuscript #

# Note 1: ggplot code and some simulation code is adapted from:
#         1) Callaway's DID quick guide: https://bcallaway11.github.io/did/articles/TWFE.html
#         2) code from Baker et al. 2022. 
#         We thank Callaway and Baker et al. for making the plotting process easier. 
#         We thank Baker et al,. for sharing the code for knowledge dissemination. 
#

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

## Last Plot - show the alternative methods and how they work over time
run_sim_es <- function(i, p) {
  
  p()
  
  # pull firm FE from empirical distribution with replacement and 
  # randomly assign state of incorporation
  sim_firm_fe <- tibble(
    gvkey = unique(shell$gvkey),
    firm_fe = sample(firm_fes, n_firm_fes, replace = TRUE),
    incorp = sample(state.abb, n_firm_fes, replace = TRUE)
  )
  
  # pull year FE from the empirical distribution with replacement
  sim_year_fe <- tibble(
    fyear = unique(shell$fyear),
    year_fe = sample(year_fes, n_year_fes, replace = TRUE)
  )
  
  # merge in the FE to the firm/year/state observations and add in residuals from the 
  # empirical distribution. ROA is the linear combination of the FEs and the residual
  data <- shell %>% 
    left_join(sim_firm_fe, by = "gvkey") %>% 
    left_join(sim_year_fe, by = "fyear") %>% 
    mutate(resid = sample(mod$residuals, length(mod$residuals), replace = TRUE),
           roa = firm_fe + year_fe + resid)
  
  # randomly assign the state of incorporation into treatment groups
  # put random states into a vector
  random_states <- sample(state.abb, length(state.abb), replace = FALSE) 
  
  # now add in the treatment effect -  Multiple Treatment Periods and Dynamic Treatment Effects
  data6 <- data %>%
    mutate(
      # figure out treatment group based on random ordering of states of incorporation
      group = case_when(
        incorp %in% random_states[1:17] ~ 1997,
        incorp %in% random_states[18:35] ~ 2009,
        incorp %in% random_states[35:50] ~ 2015
      ),
      # add in treatment effects - varying percent of standard deviation of ROA added per year
      delta_base = case_when(
        fyear >= group & group == 1997 ~ .08*sd_roa,
        fyear >= group & group == 2009 ~ .05*sd_roa,
        fyear >= group & group == 2015 ~ .03*sd_roa, 
        TRUE ~ 0
      ),
      # true treatment effect is the cumulative sum of this - dynamic trend break treatment effect
      delta = delta_base * (fyear - group + 1),
      # new ROA is the sum of the old ROA and the treatment effect
      treat_roa = roa + delta,
      # make indicator variable for obs when treatment is turned on for the TWFE regs
      treat = ifelse(fyear >= group, 1, 0),
      # make a relative-to-treatment year variable
      rel_year = fyear - group,
      # get a first treat variable for CS
      first_treat = group)  
  
  # true treatment effects
  dt <- tibble(
    sim = rep(i, 11),
    t = -5:5,
    # get the average of the imputed treatment effects for each relative time period. 
    # Drop 2015 cohort bc can't estimate.
    true_te = map_dbl(-5:5, function(x) {data6 %>% filter(rel_year == x & first_treat < 2015) %>% pull(delta) %>% mean})
  )
  
  # get CS estimates
  # first full the full set of attgts
  CS_out <- att_gt(yname = "treat_roa", 
                   data = data6,
                   gname = "first_treat",
                   idname = "gvkey", 
                   tname = "fyear", 
                   bstrap = F, 
                   cband = F,
                   est_method = "reg",
                   control_group = "notyettreated",
                   print_details = F,
                   panel = TRUE,
                   allow_unbalanced_panel = TRUE)
  
  # get cs estimates aggregated to event time
  cs <- aggte(CS_out, type = "dynamic", min_e = -5, max_e = 5, bstrap = FALSE, cband = FALSE)
  
  # add into the data
  dt$cs <- cs$att.egt
  
  ##### Stacked regressions
  # first make the stacked data
  # get the treatment cohorts
  cohorts <- data6 %>% 
    # drop never treated units, and also year 2015 when everyone is treated
    filter(!(first_treat %in% c(0, 2015))) %>% 
    pull(first_treat) %>% 
    unique()
  
  # make formula to create the sub-dataset
  getdata <- function(j) {
    
    #keep what we need
    data6 %>% 
      # keep treated units and all units not treated within -5 to 5
      filter(first_treat == j | first_treat == 0 | first_treat > j + 5) %>% 
      # keep just year -5 to 5
      filter(fyear >= j - 5 & fyear <= j + 5) %>%
      # create an indicator for the dataset
      mutate(df = j)
  }
  
  # get data stacked
  stacked_data <- map_df(cohorts, getdata) %>% 
    mutate(rel_year = if_else(df == group, rel_year, NA_real_)) %>% 
    fastDummies::dummy_cols("rel_year", ignore_na = TRUE) %>% 
    mutate(across(starts_with("rel_year_"), ~replace_na(., 0)))
  
  # get stacked value. Note that -1 was omitted. 
  stacked <- feols(treat_roa ~ `rel_year_-5` + `rel_year_-4` + `rel_year_-3` + 
                     `rel_year_-2` + rel_year_0 + rel_year_1 + rel_year_2 + rel_year_3 + 
                     rel_year_4 + rel_year_5 | gvkey^df + fyear^df, data = stacked_data)$coefficients
  
  # add in 0 for omitted -1
  stacked <- c(stacked[1:4], 0, stacked[5:10])
  
  # add in
  dt$stacked <- stacked
  
  #### Sun and Abraham
  # Drop 2015 and beyond when everyone is treated. 
  # Use less than 5 years of treatment
  sa_data <- data6 %>% 
    filter(rel_year <= 5) %>%
    filter(!first_treat == 0) %>%
    filter(fyear < 2015) 
  
  # tidy up sun abraham estimates
  sun_ab <- feols(treat_roa ~ 1 + sunab(first_treat, fyear) | gvkey + fyear, sa_data)
  sa <- tidy(sun_ab)[16:25, ] %>% pull(estimate) # make sure to extract the correct coefficients
  sa <- c(sa[1:4], 0, sa[5:10])
  dt$sa <- sa
  
  # export results
  dt
  
}

# parallelize and do 1000 simulations
x <- 1:1000
options(future.globals.maxSize=999999999)
set.seed(20040707)
plan(multisession, workers = 6)
with_progress({
  p <- progressor(steps = length(x))
  out <- future_map_dfr(
    .x = x, 
    .f = run_sim_es,
    p = p,
    .options = furrr_options(globals = c("mod", "shell", "firm_fes", "n_firm_fes",
                                         "year_fes", "n_year_fes", "sd_roa"),
                             packages = c("tidyverse", "fixest", "e1071", "did", "fastDummies"),
                             seed = TRUE)
  )})

## make plots
# function to make plot
make_es_plot <- function(name, title) {
  out %>% 
    group_by(t) %>% 
    summarize(true_effect = mean(true_te),
              est = mean({{name}}),
              lower_ci = quantile({{name}}, probs = 0.025),
              upper_ci = quantile({{name}}, 0.975)) %>% 
    # split the error bands by pre-post
    mutate(band_groups = case_when(
      t < -1 ~ "Pre",
      t >= 0 ~ "Post",
      t == -1 ~ ""
    )) %>%
    # plot
    ggplot(aes(x = t, y = est)) + 
    geom_line(aes(x = t, y = true_effect, color = "True Effect"), linetype = "dashed",linewidth = 2) + 
   # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),color = "lightgrey", alpha = 1/4) + 
    geom_point(aes(#ymin = lower_ci, ymax = upper_ci, 
                        color = "Estimated Effect"), show.legend = FALSE, 
               size = 3) + 
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = -0.5, linetype = "dashed") + 
    scale_x_continuous(breaks = -5:5) + 
    labs(x = "Year Relative to Exogenous Shock", y = "Estimate") +
    scale_color_manual(values = c("#993441","#0029a5")) + 
    ylim(c(-0.02, 0.08)) + 
    ggtitle(title) + 
    theme(legend.position = if_else(title == "Callaway and Sant'Anna Estimator", "bottom", "none"),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_text(angle = 360, hjust = 0.5, vjust = 0.5))
}

# run plots
cs_es <- make_es_plot(cs, "Callaway and Sant'Anna Estimator")
sa_es <- make_es_plot(sa, "Sun and Abraham Estimator")
stacked_es <- make_es_plot(stacked, "Stacked Regression TWFE")

# check plots
cs_es
sa_es
stacked_es

# combine all plots
total = cs_es + sa_es + stacked_es

ggsave(total, filename = here::here("Figs_Tables", "Sim_Alternatives.png"), 
       dpi = 500,
       width = 12, height = 3.5)

