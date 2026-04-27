## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  cache.path = 'cache/simulatePfsAndOsGumbel/',
  comment = '#>',
  dpi = 300,
  out.width = '100%'
)

## ----setup, echo = FALSE, message = FALSE-------------------------------------
library(dplyr)
library(knitr)
library(ggplot2)
library(TrialSimulator)

if (!requireNamespace('copula', quietly = TRUE)) {
  knitr::opts_chunk$set(eval = FALSE)
}

## ----example-endpoint, results = 'asis'---------------------------------------
pfs_and_os <- endpoint(name = c('pfs', 'os'),
                       type = c('tte', 'tte'),
                       generator = CorrelatedPfsAndOs2,
                       median_pfs = 5,
                       median_os = 11,
                       kendall = 0.6,
                       pfs_name = 'pfs',
                       os_name = 'os')

pfs_and_os

## ----direct-generator---------------------------------------------------------
set.seed(123)
dat <- CorrelatedPfsAndOs2(n = 10000,
                           median_pfs = 5,
                           median_os = 11,
                           kendall = 0.6)

head(dat, 2)

## ----validation---------------------------------------------------------------
with(dat, median(pfs))
with(dat, median(os))
with(dat, cor(pfs, os, method = 'kendall'))
with(dat, all(pfs <= os))

