### Purpose: Conduct analysis of race-ethnicity and 2FLEs
### Author:  S Bauldry
### Date:    Dec 15, 2022

### Load packages and set working directory
setwd("~/desktop")
library(tidyverse)
library(ggpubr)
library(survey)
library(weights)
library(ragg)


### Read prepared age-specific dual-function and mortality data
hrs  <- read_csv("redfl-asdf-data.csv", col_types = list(agei = "f", fem = "f")) %>%
  mutate(
    ref = as.factor(case_when (
      fbr == 1 & rce == 1 ~ 1,   #FB Hispanic
      fbr == 0 & rce == 1 ~ 2,   #US Hispanic
      rce == 2            ~ 3,   #Black
      rce == 3            ~ 4))) #White

### Proportion of observations based on proxy reports in HRS
table(hrs$pxy)


### Plot age-group specific dual-function rates (Figure 1)
# function to calculate age-group specific weighted percentages and CIs
dfwp <- function(d1) {
  des <- svydesign(id = ~1, weights = ~wtcrnh, data = d1)
  df <- data.frame( expand.grid( agei = 1:8, ref = c(1, 2, 3, 4) ) )
  df$agei <- factor(df$agei)
  df$ref  <- factor(df$ref) 
  m1  <- svyglm(df ~ agei + ref + agei*ref, data = d1, family = quasibinomial, 
                design = des)
  pdf1 <- predict(m1, newdata = df, type = "response")
  pdf <- cbind(df, pdf1)
}

# dual-function rates for women and men
dfp_fem <- hrs %>%
  filter(fem == 1) %>%
  dfwp() %>%
  rename(fdfp = response, fdfse = SE) %>%
  mutate(pdf = fdfp*100,
         plb = (fdfp - 1.96*fdfse)*100,
         pub = (fdfp + 1.96*fdfse)*100,
         agei = as.numeric(as.character(agei)))

dfp_mal <- hrs %>%
  filter(fem == 0) %>%
  dfwp() %>%
  rename(mdfp = response, mdfse = SE) %>%
  mutate(pdf = mdfp*100,
         plb = (mdfp - 1.96*mdfse)*100,
         pub = (mdfp + 1.96*mdfse)*100,
         agei = as.numeric(as.character(agei)))

# racial-ethnic gaps for women and men
dfp_fem_gap <- dfp_fem %>%
  select(agei, ref, pdf) %>%
  pivot_wider(names_from = ref, values_from = pdf) %>%
  rename(dfp_fbh = `1`, dfp_ubh = `2`, dfp_b = `3`, dfp_w = `4`) %>%
  mutate(blk_w = dfp_b - dfp_w,
         ubh_w = dfp_ubh - dfp_w,
         fbh_w = dfp_fbh - dfp_w)

dfp_mal_gap <- dfp_mal %>%
  select(agei, ref, pdf) %>%
  pivot_wider(names_from = ref, values_from = pdf) %>%
  rename(dfp_fbh = `1`, dfp_ubh = `2`, dfp_b = `3`, dfp_w = `4`) %>%
  mutate(blk_w = dfp_b - dfp_w,
         ubh_w = dfp_ubh - dfp_w,
         fbh_w = dfp_fbh - dfp_w)

# function to produce figures of age-group specific dual function rates
fga <- function(dfn, tit) {
  fig <- ggplot(data = dfn, mapping = aes(x = agei, y = pdf, color = ref)) +
    geom_line() +
    geom_point() +
    geom_ribbon(mapping = aes(ymin = plb, ymax = pub), alpha = 0.2, linetype = 0) +
    scale_color_manual(values = c("1" = "blue", "2" = "dark green", "3" = "dark red", "4" = "purple")) +
    scale_x_continuous(breaks = 1:8, labels = c("50", "55", "60", "65", "70", "75", "80", "85"), name = "age interval") +
    scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100), name = "percentage dual functional") +
    labs(title = tit) +
    guides(fill = "none", color = "none") +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          text = element_text(family = "Times New Roman"))
  return(fig)
}

fgb <- function(dfn, tit) {
  fig <- ggplot(data = dfn, mapping = aes(x = agei)) +
    geom_segment(mapping = aes(x = agei, xend = agei, y = 0, yend = blk_w), color = "dark red", size = 2) +
    geom_segment(mapping = aes(x = agei + 0.2, xend = agei + 0.2, y = 0, yend = ubh_w), color = "dark green", size = 2) +
    geom_segment(mapping = aes(x = agei + 0.4, xend = agei + 0.4, y = 0, yend = fbh_w), color = "blue", size = 2) +
    scale_x_continuous(breaks = 1:8, labels = c("50", "55", "60", "65", "70", "75", "80", "85"), name = "age interval") +
    scale_y_continuous(limits = c(-30, 0), name = "difference in percentage") +
    labs(title = tit) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          text = element_text(family = "Times New Roman"))
  return(fig)
}
  
fig1a <- fga(dfp_fem, "Panel A: Women")
fig1b <- fga(dfp_mal, "Panel B: Men")
fig1c <- fgb(dfp_fem_gap, "Panel C: Women - Gaps")
fig1d <- fgb(dfp_mal_gap, "Panel D: Men - Gaps")
fig1 <- ggarrange(fig1a, fig1b, fig1c, fig1d) 
ragg::agg_tiff("redfl-fig1.tiff", units = "in", width = 6.5, height = 6.5, res = 300)
fig1
dev.off()


### Age-50 total and dual-function life expectancy (Table 1)
# function to read bootstrap data and extract results
bsres <- function(fn) {
  bsd <- read_csv(fn)
  est <- c( mean(bsd$e50), sd(bsd$e50), quantile(bsd$e50, 0.025), quantile(bsd$e50, 0.975), 
            mean(bsd$dfe50), sd(bsd$dfe50), quantile(bsd$dfe50, 0.025), quantile(bsd$dfe50, 0.975), 
            mean(bsd$dfe50)/mean(bsd$e50) )
  return(est)
}

# function to calculate gaps and confidence intervals
bsdif <- function(mt1, mt2) {
  tle_dif    <- mt1[1] - mt2[1]
  tle_dif_se <- sqrt( mt1[2]^2 + mt2[2]^2 )
  tle_dif_lb <- tle_dif - 1.96*tle_dif_se
  tle_dif_ub <- tle_dif + 1.96*tle_dif_se
  
  dfle_dif    <- mt1[5] - mt2[5]
  dfle_dif_se <- sqrt( mt1[6]^2 + mt2[6]^2 )
  dfle_dif_lb <- dfle_dif - 1.96*dfle_dif_se
  dfle_dif_ub <- dfle_dif + 1.96*dfle_dif_se
  
  dif_est <- c( tle_dif, tle_dif_lb, tle_dif_ub, dfle_dif, dfle_dif_lb, dfle_dif_ub )
  return(dif_est)
}

# estimates by gender and race-ethnicity
le_fem_fbh <- bsres("redfl-5y-age50-le-ng-11.csv")
le_fem_ubh <- bsres("redfl-5y-age50-le-ng-01.csv")
le_fem_blk <- bsres("redfl-5y-age50-le-rg-21.csv")
le_fem_wht <- bsres("redfl-5y-age50-le-rg-31.csv")

le_mal_fbh <- bsres("redfl-5y-age50-le-ng-10.csv")
le_mal_ubh <- bsres("redfl-5y-age50-le-ng-00.csv")
le_mal_blk <- bsres("redfl-5y-age50-le-rg-20.csv")
le_mal_wht <- bsres("redfl-5y-age50-le-rg-30.csv")

# Age-50 total and dual-function life expectancy
rbind(le_fem_fbh, le_fem_ubh, le_fem_blk, le_fem_wht)
rbind(le_mal_fbh, le_mal_ubh, le_mal_blk, le_mal_wht)

# Gaps in age-50 total and dual-function life expectancy
bsdif(le_fem_fbh, le_fem_wht)
bsdif(le_fem_ubh, le_fem_wht)
bsdif(le_fem_blk, le_fem_wht)

bsdif(le_mal_fbh, le_mal_wht)
bsdif(le_mal_ubh, le_mal_wht)
bsdif(le_mal_blk, le_mal_wht)



### Age-50 total and dual-function life expectancy -- 2-year intervals
le_fem_fbh_2 <- bsres("redfl-2y-age50-le-ng-11.csv")
le_fem_ubh_2 <- bsres("redfl-2y-age50-le-ng-01.csv")
le_fem_blk_2 <- bsres("redfl-2y-age50-le-rg-21.csv")
le_fem_wht_2 <- bsres("redfl-2y-age50-le-rg-31.csv")

le_mal_fbh_2 <- bsres("redfl-2y-age50-le-ng-10.csv")
le_mal_ubh_2 <- bsres("redfl-2y-age50-le-ng-00.csv")
le_mal_blk_2 <- bsres("redfl-2y-age50-le-rg-20.csv")
le_mal_wht_2 <- bsres("redfl-2y-age50-le-rg-30.csv")

# Age-50 total and dual-function life expectancy
rbind(le_fem_fbh_2, le_fem_ubh_2, le_fem_blk_2, le_fem_wht_2)
rbind(le_mal_fbh_2, le_mal_ubh_2, le_mal_blk_2, le_mal_wht_2)

diff_2y_5y_fem_fbh <- c(le_fem_fbh[1] - le_fem_fbh_2[1], le_fem_fbh[5] - le_fem_fbh_2[5])
diff_2y_5y_fem_ubh <- c(le_fem_ubh[1] - le_fem_ubh_2[1], le_fem_ubh[5] - le_fem_ubh_2[5])
diff_2y_5y_fem_blk <- c(le_fem_blk[1] - le_fem_blk_2[1], le_fem_blk[5] - le_fem_blk_2[5])
diff_2y_5y_fem_wht <- c(le_fem_wht[1] - le_fem_wht_2[1], le_fem_wht[5] - le_fem_wht_2[5])
rbind(diff_2y_5y_fem_fbh, diff_2y_5y_fem_ubh, diff_2y_5y_fem_blk, diff_2y_5y_fem_wht)

diff_2y_5y_mal_fbh <- c(le_mal_fbh[1] - le_mal_fbh_2[1], le_mal_fbh[5] - le_mal_fbh_2[5])
diff_2y_5y_mal_ubh <- c(le_mal_ubh[1] - le_mal_ubh_2[1], le_mal_ubh[5] - le_mal_ubh_2[5])
diff_2y_5y_mal_blk <- c(le_mal_blk[1] - le_mal_blk_2[1], le_mal_blk[5] - le_mal_blk_2[5])
diff_2y_5y_mal_wht <- c(le_mal_wht[1] - le_mal_wht_2[1], le_mal_wht[5] - le_mal_wht_2[5])
rbind(diff_2y_5y_mal_fbh, diff_2y_5y_mal_ubh, diff_2y_5y_mal_blk, diff_2y_5y_mal_wht)

