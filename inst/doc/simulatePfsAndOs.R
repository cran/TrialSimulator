## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  cache.path = 'cache/simulatePfsAndOs/',
  comment = '#>',
  dpi = 300,
  out.width = '100%'
)

## ----setup, echo = FALSE, message = FALSE-------------------------------------
library(dplyr)
library(knitr)
library(ggplot2)
library(TrialSimulator)

## ----pplloott-----------------------------------------------------------------
knitr::include_graphics('three_state_ill_death_model.png')

## ----laiela-------------------------------------------------------------------
pars <- solveThreeStateModel(median_pfs = 5, median_os = 12, 
                             corr = seq(.55, .65, by = .05), 
                             h12 = seq(.05, .15, length.out = 50))
plot(pars)

## ----eiaoljf, results='asis'--------------------------------------------------
pfs_and_os <- endpoint(name = c('pfs', 'os'), 
                       type = c('tte', 'tte'), 
                       generator = CorrelatedPfsAndOs3, 
                       h01 = .11, h02 = .03, h12 = .10, 
                       pfs_name = 'pfs', os_name = 'os')
pfs_and_os

## ----llea---------------------------------------------------------------------
dat <- CorrelatedPfsAndOs3(n = 1e6, h01 = .11, h02 = .03, h12 = .10)
head(dat, 2)

## should be close to 0.6
with(dat, cor(pfs, os))

## should be close to 5.0
with(dat, median(pfs))

## should be close to 12.0
with(dat, median(os))
with(dat, all(pfs <= os))

