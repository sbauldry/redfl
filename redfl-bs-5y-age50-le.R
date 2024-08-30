### Purpose: Prepare bootstrap age-50 life expectancy estimates
### Author:  S Bauldry
### Date:    Dec 15, 2022

### load packages and set working directory
setwd("~/desktop")
library(tidyverse)
library(survival)
library(survey)


### Function for Sullivan life tables
sullivan_life_table <- function(m, d, age) {
  n            = c(diff(age), Inf)
  a            = n*0.5
  a[length(m)] = 1/m[length(m)]
  q            = (n*m)/(1 + (n - a)*m)
  q[length(m)] = 1
  p            = 1 - q
  l            = head(cumprod( c(1, p) ), -1)
  L            = (n*l*p) + (a*l*q)
  L[length(m)] = l[length(m)]/m[length(m)] 
  e            = rev(cumsum( rev(L) ))/l
  Ld           = (1 - d)*L
  ed           = rev(cumsum( rev(Ld) ))/l
  le           = cbind(e, ed)
  return(le)
}


### Function to calculate age-50 life expectancy for bootstrapped sample
le_bs <- function(i, d1, d2) {
  
  # obtain estimates of age-specific mortality rates
  if(i == 1) {
    nbs <- d1
  }
  if(i > 1) {
    nuid <- unique(d1$nhispid)
    nsam <- tibble( sample(nuid, size = length(nuid), replace = T) )
    colnames(nsam) <- "nhispid"
    nbs <- nsam %>% left_join(d1, by = "nhispid")
  }
  nhis_des <- svydesign(id = ~1, weights = ~mortwt, data = nbs)
  nm  <- svysurvreg(Surv(cerd, died) ~ cage, dist = "exponential", 
                    data = nbs, design = nhis_des)
  aslp <- predict(nm, data.frame(cage = 50:94), type = "lp")
  ashr <- exp(-aslp)
  ai <- c( sort( rep(1:7, 5) ), rep(8, 10) )
  mx <- tapply(ashr, ai, mean)
  
  # obtain estimates of age-specific dual function rates
  if(i == 1) {
    hbs <- d2
  }
  if(i > 1) {
    huid <- unique(d2$hhidpn)
    hsam <- tibble( sample(huid, size = length(huid), replace = T) )
    colnames(hsam) <- "hhidpn"
    hbs <- hsam %>% left_join(d2, by = "hhidpn")
  }
  hrs_des <- svydesign(id = ~1, weights = ~wtcrnh, data = hbs)
  m1 <- svyglm(df ~ as.factor(agei), data = hbs, family = quasibinomial, 
               design = hrs_des)
  aged <- data.frame( agei = 1:8 )
  asdf  <- cbind( aged, predict(m1, newdata = aged, type = "response"))

  # calculate age-50 total and dual-function life expectancy
  agei <- seq(50, 85, 5)
  dfp  <- asdf$response 
  e1   <- sullivan_life_table(m = mx, d = dfp, a = agei)
  
  # combine estimates
  ec <- c( e1[1,1], e1[1,2] )
    
  # return estimates
  return(ec)
}


### Read prepared NHIS-LMF and HRS data
nhis <- read_csv("redfl-asmr-data.csv")
hrs  <- read_csv("redfl-asdf-data.csv")

### Set seed for reproducing bootstrap results and number of bootstrap samples
setwd("~/dropbox/research/hlthineq/redfl/redfl-wrk")
set.seed(23730026)
nb <- 501

### Bootstrap estimates by nativity and gender
for(s in 0:1) {
  for(g in 0:1) {
    nhis_sub <- nhis %>% filter(rce == 1 & fbr == s & fem == g)
    hrs_sub  <- hrs %>% filter(rce == 1 & fbr == s & fem == g)
    est <- matrix(NA, nb, 3)
    for(i in 1:nb) {
      print( c(s, g, i) )
      be <- le_bs(i,nhis_sub, hrs_sub)
      est[i,] <- c(i, be[1], be[2])
    }
    colnames(est) <- c("bsam", "e50", "dfe50")
    est <- data.frame(est)
    fn  <- paste("redfl-5y-age50-le-ng-", paste0(s,g) , ".csv", sep = "") 
    write_csv(est, fn)
  }
}

### Bootstrap estimates by race/ethnicity and gender
for(s in 2:3) {
  for(g in 0:1) {
    nhis_sub <- nhis %>% filter(rce == s & fem == g)
    hrs_sub  <- hrs %>% filter(rce == s & fem == g)
    est <- matrix(NA, nb, 3)
    for(i in 1:nb) {
      print( c(s, g, i) )
      be <- le_bs(i,nhis_sub, hrs_sub)
      est[i,] <- c(i, be[1], be[2])
    }
    colnames(est) <- c("bsam", "e50", "dfe50")
    est <- data.frame(est)
    fn  <- paste("redfl-5y-age50-le-rg-", paste0(s,g) , ".csv", sep = "") 
    write_csv(est, fn)
  }
}
