### Purpose: Prepare HRS data for estimation of age-specific dual function rates
### Author:  S Bauldry
### Date:    Dec 15, 2022

### Load packages and set working directory
setwd("~/desktop")
library(haven)
library(tidyverse)
library(panelr)


### Load and prepare Langa-Weir cognition data (available from ICPSR)
lw_vars <- c("hhid", "pn", "cogfunction2000", "cogfunction2002",
             "cogfunction2004", "cogfunction2006", "cogfunction2008",
             "cogfunction2010", "cogfunction2012", "cogfunction2014",
             "cogfunction2016", "proxy2000", "proxy2002", "proxy2004",
             "proxy2006", "proxy2008", "proxy2010", "proxy2012", 
             "proxy2014", "proxy2016")

lw <- read_dta("langa-weir-cogfinalimp_9516wide.dta") %>%
  select(all_of( lw_vars )) %>%
  mutate(hhidpn = as.numeric( paste(hhid, pn, sep = "") ),
         r5cpx  = ifelse(proxy2000 == 5, 0, 1),
         r6cpx  = ifelse(proxy2002 == 5, 0, 1),
         r7cpx  = ifelse(proxy2004 == 5, 0, 1),
         r8cpx  = ifelse(proxy2006 == 5, 0, 1),
         r9cpx  = ifelse(proxy2008 == 5, 0, 1),
         r10cpx = ifelse(proxy2010 == 5, 0, 1),
         r11cpx = ifelse(proxy2012 == 5, 0, 1),
         r12cpx = ifelse(proxy2014 == 5, 0, 1),
         r13cpx = ifelse(proxy2016 == 5, 0, 1)) %>%
  rename(r5cgf  = cogfunction2000, r6cgf  = cogfunction2002, 
         r7cgf  = cogfunction2004, r8cgf  = cogfunction2006, 
         r9cgf  = cogfunction2008, r10cgf = cogfunction2010, 
         r11cgf = cogfunction2012, r12cgf = cogfunction2014, 
         r13cgf = cogfunction2016) %>%
  select(-c(hhid, pn, proxy2000, proxy2002, proxy2004, proxy2006,
            proxy2008, proxy2010, proxy2012, proxy2014, proxy2016))


### Load RAND HRS longitudinal file (available from RAND)
d1 <- read_dta("randhrs1992_2018v1.dta")


### Identify analysis variables
v_ti <- c("hhidpn", "ragender", "raracem", "rahispan", "rabplace", "raehsamp", "raestrat")
v_tv <- c("wtcrnh", "agey_b", "adla", "proxy")
v_00 <- paste0("r5", v_tv)
v_02 <- paste0("r6", v_tv)
v_04 <- paste0("r7", v_tv)
v_06 <- paste0("r8", v_tv)
v_08 <- paste0("r9", v_tv)
v_10 <- paste0("r10", v_tv)
v_12 <- paste0("r11", v_tv)
v_14 <- paste0("r12", v_tv)
v_16 <- paste0("r13", v_tv)
vars <- c(v_ti, v_00, v_02, v_04, v_06, v_08, v_10, v_12, v_14, v_16)


### Extract analysis variables and prepare time-invariant variables
d2 <- d1 %>%
  select(all_of( vars )) %>%
  mutate(fem = ifelse(ragender == 2, 1, 0),  # women
         fbr = ifelse(rabplace == 11, 1, 0), # foreign-born
         rce = case_when(
           rahispan == 1 ~ 1,  # Hispanic
           raracem  == 2 ~ 2,  # Black
           raracem  == 1 ~ 3,  # White
           raracem  == 3 ~ 4)) %>%
  select(-c(ragender, rahispan, rabplace, raracem))

### Merge with Langa-Weir data
d3 <- left_join(d2, lw, by = "hhidpn")


### Reshape to long form
d4 <- long_panel(d3, prefix = "r", begin = 5, end = 13, label_location = "beginning")


### Analysis sample selection
# establishing baseline sample
d5 <- d4 %>% drop_na(wtcrnh) %>% filter(agey_b >= 50, wtcrnh != 0, rce != 4)
length(d5$hhidpn)
n_distinct(d5$hhidpn)

# drop missing dual function variables
d5 <- d5 %>% drop_na(c(adla, cgf))
length(d5$hhidpn)

# drop missing demographic variables
d5 <- d5 %>% drop_na(c(fem, rce, fbr))
length(d5$hhidpn)
n_distinct(d5$hhidpn)


### Prepare age intervals and dual-function indicators
d6 <- d5 %>% 
  rename(age = agey_b) %>%
  mutate(
    # dual-function indicators (no ADLs & no dementia)
    df = ifelse(adla == 0 & cgf <= 2, 1, 0),
    
    # indicator for proxy respondent
    pxy = ifelse(proxy == 1 | cpx == 1, 1, 0),
    
    # age intervals
    agei = cut(age, breaks = c(seq(50, 85, 5), 110), right = F, labels = F)) %>%
  select(-c(proxy, cpx, adla, cgf))
  

### Save data for bootstrapping and analysis
write_csv(d6, "redfl-asdf-data.csv")
