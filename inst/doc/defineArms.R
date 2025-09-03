## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  cache = TRUE,
  cache.path = 'cache/defineArms/',
  comment = '#>',
  dpi = 300,
  out.width = '100%'
)

## ----setup, message = FALSE, echo = FALSE-------------------------------------
library(dplyr)
library(simdata)
library(TrialSimulator)
set.seed(12345)

## ----garalbaiola--------------------------------------------------------------
rng <- function(n, pfs_rate, os_rate, psa_mean, psa_sd, corr_matrix){
  
  dist <- list()
  dist[['PFS']] <- function(x) qexp(x, rate = pfs_rate)
  dist[['OS']] <- function(x) qexp(x, rate = os_rate)
  dist[['PSA_baseline']] <- function(x) qnorm(x, mean = psa_mean, sd = psa_sd)
  dist[['PSA_year1']] <- function(x) qnorm(x, mean = psa_mean - 12, sd = psa_sd)
  dsgn = simdata::simdesign_norta(cor_target_final = corr_matrix, 
                                dist = dist, 
                                transform_initial = data.frame,
                                names_final = names(dist), 
                                seed_initial = 1)
  
  simdata::simulate_data(dsgn, n_obs = n) %>% 
    mutate(PFS = pmin(PFS, OS)) %>% 
    mutate(PFS_event = 1, OS_event = 1)
  
}

## ----adldaieafg, results='asis'-----------------------------------------------
ep1 <- endpoint(name = c('PSA_baseline', 'PSA_year1', 'OS', 'PFS'), 
                type = c('non-tte', 'non-tte', 'tte', 'tte'), 
                readout = c(PSA_baseline = 0, PSA_year1 = 1), 
                generator = rng, 
                pfs_rate = log(2)/2.5, os_rate = log(2)/4.5, 
                psa_mean = 20, psa_sd = 4, 
                corr_matrix = matrix(c(1, .6, -.5, -.4, 
                                       .6, 1, -.4, -.3, 
                                       -.5, -.4, 1, .7, 
                                       -.4, -.3, .7, 1), nrow = 4))

ep1

## ----ioagjie, results='asis'--------------------------------------------------
ep2 <- endpoint(name = 'biomarker', 
                type = 'non-tte', 
                readout = c(biomarker = 0), 
                generator = rbinom, 
                size = 1, prob = .3)
ep2

## ----dlaieafj, results='asis'-------------------------------------------------
trt <- arm(name = 'treated')
trt$add_endpoints(ep1, ep2)
trt

## ----eiaojda, results='asis'--------------------------------------------------
trt <- arm(name = 'treated', PSA_baseline > 10 & PSA_year1 > 0)
trt$add_endpoints(ep1, ep2)
trt

## ----eiaojfa------------------------------------------------------------------
## not recommended
tmp <- trt$generate_data(100)
head(tmp)

