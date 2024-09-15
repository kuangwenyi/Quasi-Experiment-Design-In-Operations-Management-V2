# This is the code to download Compustat data for simulation 1 and simulation 2#
# Code partially adapted from Baker et al. 2022

#### SECTION 1 DOWNLOAD COMPUSTAT DATA AND SAVE ####

# load library
pkg1 = c('tidyverse','RPostgres','fixest','e1071','kableExtra', 
         'ggthemes', 'patchwork','did','furrr', 'latex2exp',
         'bacondecomp','ggforce','fastDummies','progressr')
newpkg1 <- pkg1[!(pkg1 %in% installed.packages()[,1])]
install.packages(newpkg1)

load_pkg1 = lapply(pkg1, library, character.only=TRUE)


# Connect to WRDS Server
# replace username and password using your own
#### NEED TO LOG IN TO WRDS AND VERIFY THROUGH DUO BEFORE running below
wrds <- dbConnect(Postgres(),
                  host = 'wrds-pgdata.wharton.upenn.edu',
                  port = 9737,
                  user = 'jianliang_hao',
                  password = 'JIAYOU1ke1de!',
                  dbname = 'wrds',
                  sslmode = 'require')


# download compustat data. only dowload necessary variables
# gvkey:Global Company Key
# fyear: fiscal year
# opiti: Operating Income - Total
# oibdp: Operating Income Before Depreciation
# at:Assets - Total
# indfmt:Industry Format
# datafmt:Data Format
# popsrc:Population Source
# consol:Level of Consolidation - Company Annual Descriptor
# fic:Current ISO Country Code - Incorporation
# sich:Standard Industrial Classification - Historical
# naicsh: North American Industrial Classification - Historical

comp <- tbl(wrds, sql("SELECT gvkey, fyear, oibdp, at, indfmt, datafmt, popsrc, consol, fic, sich,naicsh FROM comp.funda")) %>%
  filter(indfmt == 'INDL' & 
           datafmt == 'STD' & 
           popsrc == 'D' & 
           consol == 'C' & 
           !is.na(fyear) & 
           fyear %>% between(1988, 2019) & 
           !is.na(at) & 
           at > 0 & 
           fic == "USA") %>% 
  group_by(gvkey) %>% 
  mutate(roa = oibdp / lag(at),
         gvkey = as.numeric(gvkey)) %>% 
  filter(!is.na(roa)) %>% 
  filter(fyear >= 1989) %>% 
  ungroup() %>% 
  collect()

comp$naicsh_n = substr(comp$naicsh,start=1,stop=2)

# download comp header to bind


comp_header <- tbl(wrds, sql("SELECT * FROM crsp.comphead")) %>% 
  mutate(gvkey = as.numeric(gvkey)) %>% 
  collect()
# industry code: Bali, Engle, and Murray (2016)
comp_header$naics_n = substr(comp_header$naics,start=1,stop=2)

# merge in state of incorporation and industry information
comp <- comp %>% 
  left_join(comp_header %>% select(gvkey, incorp, naics_n)) %>% 
  # drop if state of incorporation is missing or not 50 states
  filter(!is.na(incorp) & incorp %in% state.abb) %>% 
  # clean up SIC code - use historical sic code if available, if not use header sic
  mutate(sich = coalesce(naicsh_n, naics_n)) %>% 
  # Keep transportation industry only 4000 - 4899 
  filter(( naicsh_n%in% c(48:49))) # manufacturing 2000 - 3999

# make sure that each firm has at least 20 observations
comp <- comp %>% 
  group_by(gvkey) %>% 
  add_tally() %>% 
  filter(n >= 20) %>% 
  ungroup()

# winsorize ROA at 0.99, and censor at -1
wins <- function(x) {
  case_when(
    is.na(x) ~ NA_real_,
    x < -1 ~ -1,
    x > quantile(x, 0.99, na.rm = TRUE) ~ quantile(x, 0.99, na.rm = TRUE),
    TRUE ~ x
  )
}

# winsorize ROA by years
comp <- comp %>% 
  group_by(fyear) %>% 
  mutate(roa = wins(roa)) %>% 
  arrange(gvkey, fyear) %>% 
  ungroup()

# doublecheck ROA >0.90 and -1
comp = comp[comp$roa <=0.99, ]
comp = comp[comp$roa !=-1, ]

# save

# Create a folder named Data in the file path displayed by here:here()
saveRDS(comp, here::here("Data", "sim_tpt_data.rds"))


#### SECTION 2 PLOT Average ROA ####
avgroa =comp %>% group_by(fyear) %>% summarise(avg_roa = mean(roa))

# summarize
total = comp %>% group_by (gvkey) %>% summarise(count = n_distinct(fyear))
# write_csv(test, file = "avgroa.csv")

p1= ggplot(avgroa, aes(fyear, avg_roa)) +  
  geom_line() +
  geom_vline(xintercept = 2002,
             col = "#0029a5", lwd = 0.5) +
  ggtitle("ROA of US Tranportation and Warehousing Sector:\n1989 - 2019")+
  labs(x = "Year", y = "ROA")+ 
  scale_x_continuous(breaks=seq(1989, 2019, 3))

p3= ggplot(avgroa, aes(fyear, avg_roa)) +  
  geom_line() +
  geom_vline(xintercept = c(1997,2009,2015),
             col = "#0029a5", lwd = 0.5) +
  ggtitle("ROA of US Tranportation and Warehousing Sector:\n1989 - 2019")+
  labs(x = "Year", y = "ROA")+ 
  scale_x_continuous(breaks=seq(1989, 2019, 3))

# save image
ggsave(p1, filename = here::here("Figs_Tables", "ROA_1cut.png"), 
       dpi = 500,
       width = 6, height = 4)
ggsave(p3, filename = here::here("Figs_Tables", "ROA_4849.png"), 
       dpi = 500,
       width = 6, height = 4)



  
  
  