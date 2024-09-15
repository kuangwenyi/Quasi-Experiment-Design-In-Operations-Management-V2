############## Replication of Wang et al. 2022 #######################
# This is Section 7.2 in the manuscript 
# Figure 7, 8, 9, 10
# Table 4, 5

# load packages
pkgs = c(
  'tidyverse','patchwork','fastDummies','ggthemes','did','bacondecomp',
  'kableExtra','fixest','ggplot2','readxl','readr','tidyr',
  'dplyr','stringr','lme4','RColorBrewer','broom.mixed', 'TwoWayFEWeights', 
  'DIDmultiplegt', 'here')

kwy = lapply(pkgs, library, character.only=TRUE)

# set theme
theme_set(theme_clean() + theme(plot.background = element_blank(),
                                legend.background = element_blank()))

#### 1. Data Prep #### 
# read in file. file can be retrieved from Wang et al. 2022.
library(haven)
wang = read_dta('OriginDestKiVolAnalysis.dta') 

# check year and event year
sort(unique(wang$year))
sort(unique(wang$EventYear))

# create OD pair ID
OD_list = data.frame(unique(wang$OrgDestAirport))
OD_list$OD_Dummy = seq(1:1958)
colnames(OD_list) = c("OrgDestAirport", "OD_Dummy")

# join back to data
wang1 = wang %>% left_join(OD_list, by= "OrgDestAirport")

# create occasion ID
occasion_list = data.frame(sort(unique(wang$year)))
occasion_list$Occasion_Dummy = seq(1:16)
colnames(occasion_list) = c("year", "Occasion_Dummy")

wang2 = wang1 %>% left_join(occasion_list, by= "year")

# create FirstTreat Occasion (EventYear variable in Wang et al. 2022)
FirstTreat = occasion_list
colnames(FirstTreat) = c("EventYear", "FirstTreat")

wang3 = wang2 %>% left_join(FirstTreat, by = "EventYear")

# create Treated variable (TreatmentAfter variable in Wang et al. 2022)
wang3$Treated = ifelse(wang3$treatment_pair == 1 & wang3$EventYear <= wang3$year, 
                       1, 
                       0)

# only use >800 miles OD pair based on Wang et al.2022
wang_final = wang3 %>% filter(airport_dist >= 800)
write.csv(wang_final, file = "wang_final.csv")


#### 2. Estimate using Callaway and Sant'Anna 2021 Estimator -----------
##### 2.1 Model 1 of Wang et al. Table 3 ---------

# replace NA with 0 following Callaway and Sant'Anna 2021 for FristTreat
unique(wang_final$FirstTreat) # check NAs
wang_final$FirstTreat[is.na(wang_final$FirstTreat)] = 0
wang_final$FirstTreat = as.numeric(wang_final$FirstTreat)

# cs estimator
cs_wang_m1 <- att_gt(yname = "LnOrganVolume",
                     data = wang_final,
                     tname = "Occasion_Dummy",
                     idname = "OD_Dummy",
                     gname = "FirstTreat",
                     xformla = ~ OD_Dummy + Occasion_Dummy, # controls from Model 1 Table 3
                     control_group = c("notyettreated","nevertreated"),
                     est_method = "ipw",
                     print_details = FALSE,
                     bstrap = TRUE,
                     cband = TRUE,
                     base_period = "varyiing",
                     clustervars = "OD_Dummy",
                     panel = TRUE,
                     anticipation = 0,
                     allow_unbalanced_panel = TRUE)

# event study type aggregation pre 3 post 3 based on Wang et al. 2022
ES_wang_m1 <- aggte(cs_wang_m1, type="dynamic", min_e = -3,  
                    max_e = 3, na.rm = TRUE)

# get the coefficients for the time indicators
summary(ES_wang_m1)


# Aggregate by groups to delineate the effects
cs_wang_group <- aggte(cs_wang_m1, type="group", min_e = -3,  
                       max_e = 3, na.rm = TRUE)
summary(cs_wang_group)

groupplot = ggdid(cs_wang_group) +
  ggtitle("Average Effect by Group of New Routes")

ggsave(groupplot, 
       filename = here::here("Figs_Tables", "Group_plot_Wang.png"),
       dpi = 500, width = 5, height = 4)


#####  2.2 Static Effect as in Table 2 of Wang et al. 2022 ---------
# post 3 years
cs_wang_static <- aggte(cs_wang_m1, type="simple", min_e = -3,  
                        max_e = 3, na.rm = TRUE)
summary(cs_wang_static)


##### 2.3 Event Study Plot ----------

cs_plot_wang <- tidy(ES_wang_m1) %>% 
  as_tibble() %>% 
  mutate(group = as.factor(case_when(
    event.time < 0 ~ 1,
    event.time >= 0 ~ 2
  ))) %>% 
  # plot
  ggplot(aes(x = event.time, y = estimate)) + 
  geom_point(aes(fill= factor(group)), shape = 21) + geom_line() + 
  scale_fill_manual(values = c("#993441", "#0029a5")) + 
  ggtitle("Impact of New Routes on Kidney Tranplantation") + 
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, 
                    color = factor(group)), 
                linetype = "longdash", show.legend = FALSE) + 
  scale_color_manual(values = c("#993441", "#0029a5"))+
  geom_hline(yintercept = 0,  linetype = "longdash", color = "gray") + 
  geom_vline(xintercept = 0,  linetype = "longdash", color = "gray") + 
  labs(y = "Treatment Effect", x = "Years Relative to New Routes Introduction",
       subtitle = "Callaway and Sant’Anna (2021) Estimator") + 
  scale_x_continuous(breaks = seq(-8, 8, by = 1)) + 
  scale_y_continuous(breaks = seq(-4, 4, by = 0.1)) + 
  theme(axis.title.y = element_text(hjust = 0.5, vjust = 0.5, angle = 90),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")

# save 
ggsave(cs_plot_wang, 
       filename = here::here("Figs_Tables", "cs_plot_Wang.png"),
       dpi = 500, width = 6, height = 4)

##### 2.4 Model 4 in Table 3 of Wang et al. 2022 ------------
cs_wang_m4 <- att_gt(yname = "LnOrganVolume",
                     data = wang_final,
                     tname = "Occasion_Dummy",
                     idname = "OD_Dummy",
                     gname = "FirstTreat",
                     xformla = ~ OD_Dummy + Occasion_Dummy + 
                       CanAge + CanHeight + CanWeight + CanGenderFemale + CanRaceWhite + 
                       CanRaceBlack + CanRaceHispanic + CanBloodA + CanBloodB +
                       CanBloodO + CanDiabetesYes + CanDiabetesUnknown + CanDialysisYes + 
                       CanDialysisYear + DonAge + DonHeigth + DonWeight + DonCreat + 
                       DonGenderFemale + DonRaceWhite + DonRaceBlack + DonRaceHispanic + 
                       DonBloodA + DonBloodB + DonBloodO + 
                       DonDiabetesYes + DonDiabetesUnknown + DonHypertensionYes + 
                       DonHypertensionUnknown + DonHepCYes + DonHepCUnknown +
                       DonCodAnoxia + DonCodCVS + DonCodTrauma + DonECDYes, # controls
                     control_group = c("notyettreated","nevertreated"),
                     est_method = "ipw",
                     print_details = FALSE,
                     bstrap = TRUE,
                     cband = TRUE,
                     base_period = "varyiing",
                     clustervars = "OD_Dummy",
                     panel = TRUE,
                     anticipation = 0,
                     allow_unbalanced_panel = TRUE)


# event study type aggregation
ES_wang_m4 <- aggte(cs_wang_m4, type="dynamic", min_e = -3,  
                    max_e = 3, na.rm = TRUE)

# get the coefficients
summary(ES_wang_m4)

#### 3. Using Stacked Regression #### 
##### 3.1 Create Clean Controls Data --------
# function to get treat-year specific cohorts for the event window
make_dt <- function(tyr) {
  wang_final %>% 
    # filter(Occasion_Dummy <= 16) %>% # drop observation after everyone is treated
    filter(FirstTreat == tyr | FirstTreat > tyr + 3) %>% 
    filter(Occasion_Dummy %>% between(tyr - 4, tyr + 3)) %>% # pre 4 post 3 (-4 will be omitted in regression)
    mutate(FirstTreat = if_else(FirstTreat == tyr, FirstTreat, NA_real_),
           rel_occasion = Occasion_Dummy - FirstTreat) %>% 
    select(LnOrganVolume, 
           rel_occasion, 
           FirstTreat,
           Treated, 
           year, 
           Occasion_Dummy, 
           OD_Dummy
           ) %>% 
    mutate(dt = as.character(tyr))
}


wang_final$EventYear
# treats
treats <- wang_final %>% 
  filter(FirstTreat != 0) %>% # get rid of 0
  filter(EventYear < 2017) %>% # less than 2017 in which year every route was treated
  pull(FirstTreat) %>% 
  unique() %>% 
  sort()


# stack the data 
stacked_data <- map_dfr(treats, make_dt) %>% 
  dummy_cols(select_columns = "rel_occasion", 
             remove_selected_columns = FALSE,
             ignore_na = TRUE) %>% 
  mutate(across(starts_with("rel_occasion_"), ~replace_na(., 0))) %>% 
  mutate(cluster = paste0(OD_Dummy, "_", dt))

# make time indicators: pre 3 post 3 
indicatorStacked <- c(paste0("`", "rel_occasion_", c(-3:-1, 0:3), "`"))


##### 3.2 Model Estimation Using Stacked Regression #####
# Model 1 in Table 3 of Wang et al. 
stack_Reg <- feols( LnOrganVolume ~ .[indicatorStacked]
                       |OD_Dummy^dt + Occasion_Dummy^dt, # include OD pair and Time Fixed Effects
                       cluster = "OD_Dummy",
                       data = stacked_data)

summary(stack_Reg)

# check the coefficients
test = broom::tidy(stack_Reg, conf.int = TRUE)

ES_Wang_Stacked

ES_Wang_Stacked <- broom::tidy(stack_Reg, conf.int = TRUE)[1:7,] %>%
  mutate(t = c(-3:-1, 0:3)) %>% 
  select(t, estimate, conf.low, conf.high) %>% 
  mutate(group = as.factor(case_when(
    t < 0 ~ 1,
    t >= 1 ~ 2
  ))) %>% 
  ggplot(aes(x = t, y = estimate)) + 
  geom_point(aes(fill= factor(group)),shape = 21) +
  scale_fill_manual(values = c("#993441", "#0029a5")) + geom_line() +
  ggtitle("Impact of New Routes on Kidney Transplantation") + 
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, color=factor(group)), 
                linetype = "longdash", show.legend = FALSE) + 
  scale_color_manual("group", breaks=c(1,2),values=c("#993441", "#0029a5"))+
  geom_hline(yintercept = 0,  linetype = "longdash", color = "gray") + 
  geom_vline(xintercept = 0,  linetype = "longdash", color = "gray") + 
  labs(y = "Point Change", x = "Years Relative to New Route Introduction",
       subtitle = "Stacked Regression") + 
  scale_x_continuous(breaks = seq(-6, 6, by = 1)) + 
  scale_y_continuous(breaks = seq(-1, 1, by = 0.1)) + 
  theme(axis.title.y = element_text(hjust = 0.5, vjust = 0.5, angle = 90),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")


ggsave(ES_Wang_Stacked, 
       filename = here::here("Figs_Tables", "ES_Wang_Stacked.png"), 
       dpi = 500, width = 6, height = 4)


# static effect estimated using Stacked Data
stack_static <- feols(LnOrganVolume ~ Treated
                    |OD_Dummy + Occasion_Dummy,
                    cluster = "OD_Dummy",
                    data = stacked_data)

summary(stack_static)


#### 4. Bacon Decomposition -------------
##### 4.1 Calculate Decomposition -----------------------

bacon_out <- bacon(LnOrganVolume ~ Treated,
                   data = wang_final,
                   id_var = "OD_Dummy",
                   time_var = "Occasion_Dummy") 

save(bacon_out, file = "bacon_out_Wang2022.rda")

##### 4.2 Plot of Bacon Decomposition ------------------

# get the total weight for each group
total_weights <- bacon_out %>% 
  group_by(type) %>% 
  summarize(weight = sum(weight))
# get the weighted average within group
group_avg <- bacon_out %>% 
  group_by(type) %>% 
  summarize(avg = weighted.mean(estimate, weight),
            weights = sum(weight))


#### Make Bacon Decomposition Plots
# Group 1: earlier v later 
EvL <- bacon_out %>% 
  filter(type == "Earlier vs Later Treated") %>% 
  ggplot(aes(x = weight, y = estimate)) + 
  geom_point(size = 3, alpha = 1/2) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = group_avg$avg[1], color = "#014182", size = 1.5) + 
  labs(x = "Weight", y = expression(widehat(delta^'DD'))) + 
  ggtitle(paste0("Early vs Later Treated \n Total Weight = ", scales::percent(total_weights$weight[1]))) + 
  scale_y_continuous(limits = c(-.12, 0.12)) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_text(angle = 360, hjust = 0.5, vjust = 0.5))

# Group 2: later VS Always
LvA <- bacon_out %>% 
  filter(type == "Later vs Always Treated") %>% 
  ggplot(aes(x = weight, y = estimate)) + 
  geom_point(size = 3, alpha = 1/2) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = group_avg$avg[2], color = "#014182", size = 1.5) + 
  labs(x = "Weight", y = expression(widehat(delta^'DD'))) + 
  scale_y_continuous(limits = c(-.12, 0.12)) + 
  ggtitle(paste0("Later vs Always Treated \n Total Weight = ", scales::percent(total_weights$weight[2]))) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_text(angle = 360, hjust = 0.5, vjust = 0.5))


# Group 3: Later v Earlier
LvE <- bacon_out %>% 
  filter(type == "Later vs Earlier Treated") %>% 
  ggplot(aes(x = weight, y = estimate)) + 
  geom_point(size = 3, alpha = 1/2) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = group_avg$avg[3], color = "#014182", size = 1.5) + 
  labs(x = "Weight", y = expression(widehat(delta^'DD'))) + 
  scale_y_continuous(limits = c(-.12, 0.12)) + 
  ggtitle(paste0("Later vs Early Treated \n Total Weight = ", scales::percent(total_weights$weight[3]))) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_text(angle = 360, hjust = 0.5, vjust = 0.5))

# Group 4: Treated VS Untreated
TvU <- bacon_out %>% 
  filter(type == "Treated vs Untreated") %>% 
  ggplot(aes(x = weight, y = estimate)) + 
  geom_point(size = 3, alpha = 1/2) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = group_avg$avg[[4]], color = "#014182", size = 1.5) + 
  labs(x = "Weight", y = expression(widehat(delta^'DD'))) + 
  scale_y_continuous(limits = c(-.12, 0.12)) + 
  ggtitle(paste0("Treated vs Untreated \n Total Weight = ", scales::percent(total_weights$weight[4]))) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_text(angle = 360, hjust = 0.5, vjust = 0.5))

# combine the figures and save
BLL_decomp_plot <- EvL + LvA + LvE + TvU
plot(BLL_decomp_plot) # check plot before saving

# save (This is Figure 11)
ggsave(BLL_decomp_plot, 
       filename = here::here("Figs_Tables", "Wang_decomp_plot.png"),
       dpi = 500, width = 8, height = 6)


#### 5. Treatment Timing Plot --------------- 

wang_timing <- wang_final %>% 
  filter(treatment_pair == 1) %>%
  select(OD_Dummy, Occasion_Dummy, FirstTreat) %>% 
  mutate(OD_n = as.factor(OD_Dummy)) %>%
  mutate(OD = fct_reorder(OD_n, rank(desc(OD_n)))) %>% 
  mutate(post = if_else(Occasion_Dummy < FirstTreat, "Pre", "Post")) %>% 
  mutate(post = factor(post, levels = c("Pre", "Post"))) %>% 
  ggplot(aes(x = Occasion_Dummy, y = OD)) + 
  geom_tile(aes(fill = as.factor(post)), alpha = 3/4) + 
  scale_fill_manual(values = c("#5e001f", "#030E4F")) + 
  ggtitle("New Routes Adoption")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom',
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.background = element_rect(color = "white")) 

plot(wang_timing) # check plot

ggsave(wang_timing, 
       filename = here::here("Figs_Tables", "Wang_Timing_plot.png"),
       dpi = 500, width = 5, height = 5)





