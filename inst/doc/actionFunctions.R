## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  cache = TRUE,
  cache.path = 'cache/defineActionFunctions/',
  comment = '#>',
  dpi = 300,
  out.width = '100%'
)

## ----setup, message = FALSE, echo=FALSE---------------------------------------
library(dplyr)
library(TrialSimulator)
set.seed(12345)

## ----eval=FALSE, include=TRUE-------------------------------------------------
# listener <- listener()
# #' register milestones with listener
# listener$add_milestones(interim, final)
# controller <- controller(trial, listener)
# controller$run()

