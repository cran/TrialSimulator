## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  cache = TRUE,
  cache.path = 'cache/defineLongitudinalEndpoints/',
  comment = '#>',
  dpi = 300,
  out.width = '100%'
)

## ----setup, message = FALSE, echo=FALSE---------------------------------------
library(dplyr)
library(TrialSimulator)
set.seed(12345)

## ----adfa, echo = FALSE-------------------------------------------------------
locked_data_ <- data.frame(
  baseline = runif(6, 130, 150), 
  bp2 <- runif(6, 110, 140), 
  bp4 <- runif(6, 105, 130)
) %>% 
  mutate(
    bp_cfb2 = bp2 - baseline, 
    bp_cfb4 = bp4 - baseline
  ) %>% 
  select(baseline, bp_cfb2, bp_cfb4) %>% 
  tibble()

## ----dalfjal, echo = FALSE----------------------------------------------------
locked_data_

## ----eioajf-------------------------------------------------------------------
library(mvtnorm)
bp_generator <- function(n, bp_means, bp_vcov){
  dat <- rmvnorm(n, mean = bp_means, sigma = bp_vcov) %>% 
    as.data.frame()
  names(dat) <- c('baseline', 'bp2', 'bp4')
  dat %>% 
    mutate(
      bp_cfb2 = bp2 - baseline, 
      bp_cfb4 = bp4 - baseline
    ) %>% 
    select(baseline, bp_cfb2, bp_cfb4)
}

## ----adlfjadl-----------------------------------------------------------------
vcov1 <- matrix(
  c(2, 1.5, 1, 
    1.5, 3, 1.5, 
    1, 1.5, 4),
  nrow = 3
)

ep_in_trt1 <- endpoint(
  name = c('bp_cfb2', 'baseline', 'bp_cfb4'), 
  type = rep('non-tte', 3), 
  readout = c(baseline = 0, bp_cfb4 = 4, bp_cfb2 = 2), 
  generator = bp_generator, 
  bp_means = c(140, 125, 120), 
  bp_vcov = vcov1
)

## ----eiieala, results='asis'--------------------------------------------------
ep_in_trt1

## ----adlfeqadl, results='asis'------------------------------------------------
vcov2 <- matrix(
  c(2, 1.5, 1, 
    1.5, 3, 1, 
    1, 1, 4),
  nrow = 3
)

ep_in_trt2 <- endpoint(
  name = c('bp_cfb2', 'baseline', 'bp_cfb4'), 
  type = rep('non-tte', 3), 
  readout = c(baseline = 0, bp_cfb4 = 4, bp_cfb2 = 2), 
  generator = bp_generator, 
  bp_means = c(140, 127, 122), 
  bp_vcov = vcov2
)

ep_in_trt2

## ----daliea, results='asis'---------------------------------------------------
trt1 <- arm(name = 'treatment 1')
trt2 <- arm(name = 'treatment 2')

trt1$add_endpoints(ep_in_trt1)
trt2$add_endpoints(ep_in_trt2)

trt1

trt2

