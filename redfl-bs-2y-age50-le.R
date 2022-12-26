### Purpose: Prepare bootstrap age-50 life expectancy estimates
### Author:  S Bauldry
### Date:    Dec 15, 2022

### load packages and set working directory
setwd("~/desktop")
library(tidyverse)
library(survival)
library(survey)


### Function for life table calculations
lifetable <- function(x) {
  
  # calculate qx and lx
  for(i in 1:nrow(x)) {
    if(i < nrow(x)) {
      x[i,"qx"] = (x[i,"nx"]*x[i,"mx"])/(1 + (x[i,"nx"] - x[i,"ax"])*x[i,"mx"])
    }
    if(i == nrow(x)) {
      x[i,"qx"] = 1
    }
    
    if(i == 1) {
      x[i,"lx"] = 100000
    }
    if(i > 1) {
      x[i,"lx"] = x[i-1,"lx"]*(1 - x[i-1,"qx"])
    }
  }
  
  # calculate dx, Lx, and dLx
  for(i in 1:nrow(x)) {
    if(i < nrow(x)) {
      x[i,"dx"] = x[i,"lx"] - x[i+1,"lx"]
      x[i,"Lx"] = x[i,"nx"]*x[i+1,"lx"] + x[i,"dx"]*x[i,"ax"]
    }
    if(i == nrow(x)) {
      x[i,"dx"] = x[i,"lx"]
      x[i,"Lx"] = x[i,"lx"]/x[i,"mx"]
    }
    x[i,"dLx"] = x[i,"Lx"]*(x[i,"dfp"])
  }
  
  # calculate Tx and dTx
  for(i in nrow(x):1) {
    if(i < nrow(x)) {
      x[i,"Tx"]  = x[i+1,"Tx"] + x[i,"Lx"]
      x[i,"dTx"] = x[i+1,"dTx"] + x[i,"dLx"]
    } 
    if(i == nrow(x)) {
      x[i,"Tx"]  = x[i,"Lx"]
      x[i,"dTx"] = x[i,"dLx"]
    }
  }
  
  # calculate ex and dex
  for(i in 1:nrow(x)) {
    x[i,"ex"]  = x[i,"Tx"]/x[i,"lx"]
    x[i,"dex"] = x[i,"dTx"]/x[i,"lx"]
  }
  
  # return age-50 ex and dex
  e <- x[1, c("agei", "ex", "dex")]
  return(e)
}


### Function to calculate age-50 life expectancy for bootstrapped sample
le_est <- function(d1, d2) {
  
  # obtain estimates of age-specific mortality rates
  nuid <- unique(d1$nhispid)
  nsam <- tibble( sample(nuid, size = length(nuid), replace = T) )
  colnames(nsam) <- "nhispid"
  nbs <- nsam %>% left_join(d1, by = "nhispid")
  nhis_des <- svydesign(id = ~1, weights = ~mortwt, data = nbs)
  nm  <- svysurvreg(Surv(cerd, died) ~ cage, dist = "exponential", 
                    data = nbs, design = nhis_des)
  aslp <- predict(nm, data.frame(cage = 50:95), type = "lp")
  ashr <- exp(-aslp)
  ai <- c( sort( rep(1:18, 2) ), rep(19, 10) )
  mx <- tapply(ashr, ai, mean)
  
  # obtain estimates of age-specific dual function rates
  huid <- unique(d2$hhidpn)
  hsam <- tibble( sample(huid, size = length(huid), replace = T) )
  colnames(hsam) <- "hhidpn"
  hbs <- hsam %>% left_join(d2, by = "hhidpn")
  hrs_des <- svydesign(id = ~1, weights = ~wtcrnh, data = hbs)
  m1 <- svyglm(df ~ as.factor(agei), data = hbs, family = quasibinomial, 
               design = hrs_des)
  aged <- data.frame( agei = 1:19 )
  asdf  <- cbind( aged, predict(m1, newdata = aged, type = "response"))

  # calculate age-50 total and dual-function life expectancy
  agei <- seq(50, 86, 2)
  nx   <- c( rep(2, 18), 10 )
  ax   <- c( rep(1, 18), 5 )
  dfp  <- asdf$response  
  lt   <- tibble(agei, nx, ax, mx, dfp)
  e1   <- lifetable(lt)
  
  # combine estimates
  ec <- c( e1[[1,1]], e1[[1,2]], e1[[1,3]] )
    
  # return estimates
  return(ec)
}


### Read prepared NHIS-LMF and HRS data
nhis <- read_csv("redfl-asmr-data.csv")
hrs  <- read_csv("redfl-asdf-data.csv")

### Set seed for reproducing bootstrap results and number of bootstrap samples
setwd("~/dropbox/research/hlthineq/redfl/redfl-wrk")
set.seed(23730026)
nb <- 500

### Bootstrap estimates by nativity and gender
for(s in 0:1) {
  for(g in 0:1) {
    nhis_sub <- nhis %>% filter(rce == 1 & fbr == s & fem == g)
    hrs_sub  <- hrs %>% filter(rce == 1 & fbr == s & fem == g)
    est <- matrix(NA, nb, 3)
    for(i in 1:nb) {
      print( c(s, g, i) )
      be <- le_est(nhis_sub, hrs_sub)
      est[i,] <- c(i, be[2], be[3])
    }
    colnames(est) <- c("bsam", "e50", "dfe50")
    est <- data.frame(est)
    fn  <- paste("redfl-2y-age50-le-ng-", paste0(s,g) , ".csv", sep = "") 
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
      be <- le_est(nhis_sub, hrs_sub)
      est[i,] <- c(i, be[2], be[3])
    }
    colnames(est) <- c("bsam", "e50", "dfe50")
    est <- data.frame(est)
    fn  <- paste("redfl-2y-age50-le-rg-", paste0(s,g) , ".csv", sep = "") 
    write_csv(est, fn)
  }
}
