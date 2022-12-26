### Purpose: Prepare NHIS-LMF for estimates of age-specific mortality rates
### Author:  S Bauldry
### Date:    Dec 15, 2022

### Loading packages
setwd("~/desktop")
library(ipumsr)
library(tidyverse)
library(lubridate)
library(clock)

### Read data extract
ddi <- read_ipums_ddi("nhis_00003.xml")
d1 <- read_ipums_micro(ddi)

### Prepare variables for analysis
d1 <- d1 %>%
  rename_with(tolower) %>%
  filter(age >= 50, mortelig == 1) %>%
  mutate(
    # prepare weights and covariates
    fem = ifelse(sex == 2, 1, 0), 
    rce = case_when(
      hispeth != 10 ~ 1, # Hispanic
      racea == 200 ~ 2,  # non-Hispanic Black
      racea == 100 ~ 3,  # non-Hispanic White
      racea >= 310 ~ 4), # other races
    fbr = case_when(
      usborn == 10 | usborn == 11 | usborn == 12 ~ 1,
      usborn == 20 ~ 0),
    
    # set top-coded age at 85 to missing
    age = ifelse(age < 85, age, NA),
    
    # set missing birth month to mid-year
    bmo = ifelse(birthmo < 97, birthmo, 7),
    
    # set missing values for birth year
    byr = ifelse(birthyr < 9997, birthyr, NA)
  )

### Set bottom-coded birth year to missing
for(i in 1:16) {
  bry = 1914 + i
  sry = 1998 + i
  d1$byr[d1$year == sry & d1$byr == bry] = NA
}

### Calculate dates and age on Jan 1
d1 <- d1 %>%
  mutate(
    # interview date
    intyr = ifelse(intervwyr == 9998, year, intervwyr),
    doi = date_build(intyr, intervwmo, 15),
    
    # birth date
    byr = ifelse(is.na(byr), intyr - age, byr),
    dob = date_build(byr, bmo, 15),
    
    # date of death or censoring
    dth = ifelse(mortstat == 1, 1, 0),
    dod = case_when(
      dth == 0 ~ as_date( date_build(2015, 12, 31) ),
      mortdodq == 1 ~ as_date( date_build(mortdody, 2, 15) ),
      mortdodq == 2 ~ as_date( date_build(mortdody, 5, 15) ),
      mortdodq == 3 ~ as_date( date_build(mortdody, 8, 15) ),
      mortdodq == 4 ~ as_date( date_build(mortdody, 11, 15) )),
    
    # fix 121 cases with interview date after death date
    dod = case_when(
      doi <= dod ~ dod,
      doi > dod ~ doi + 1),
    
    # age on Jan 1
    jage = as.numeric( round( (date_build(intyr, 1, 1) - dob)/365.25 ) )
  )

### Select analysis sample and variables before expanding
# baseline sample size
d2 <- d1 %>% filter(jage >= 50, mortwt != 0)
length(d1$nhispid)

# drop other races
d2 <- d2 %>% filter(rce != 4)
length(d1$nhispid)

# drop missing date of death
d2 <- d2 %>% drop_na(dod)
length(d1$nhispid)

# drop missing demographic variables
d2 <- d2 %>% drop_na(c(fem, rce, fbr))
length(d1$nhispid)

# select analysis variables
d2 <- d2 %>% select(nhispid, intyr, fem, fbr, rce, jage, doi, dod, dth, mortwt, 
                    strata, psu)

### Expand to person-year format
d3 <- d2 %>%
  mutate(ny = ifelse(intyr < year(dod), year(dod) + 1 - intyr, 1)) %>%
  uncount(ny)
length(d3$nhispid)

### Prepare person-year data for survival models
d4 <- d3 %>%
  group_by(nhispid) %>%
  mutate(
    
    # calculate current year
    cyear = row_number() + intyr - 1,
    
    # calculate age in current year
    cage = jage + row_number() - 1,
    cage = ifelse(cage > 85, 85, cage),
    
    # create indicator for died in current year
    died = ifelse(dth == 1 & cyear == max(cyear), 1, 0),
    
    # create exposure variable in current year
    cyrb = date_build(cyear, 1, 1),
    cyre = date_build(cyear, 12, 31),
    cerd = ifelse(cyear == year(dod) & died == 1, (dod - cyrb + 1)/365.25, 1),
    cerd = ifelse(cyear == intyr & died == 0, (cyre - doi + 1)/365.25, cerd)
  ) %>%
  select(nhispid, fem, fbr, rce, cage, cerd, died, mortwt, strata, psu)

### Save data for bootstrapping
write_csv(d4, "redfl-asmr-data.csv")