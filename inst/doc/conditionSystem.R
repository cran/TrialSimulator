## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  cache = TRUE,
  cache.path = 'cache/conditionSystem/',
  comment = '#>',
  dpi = 300,
  out.width = '100%'
)

## ----setup, message = FALSE, echo=FALSE---------------------------------------
library(TrialSimulator)

## ----oiea---------------------------------------------------------------------
## Note: the unit of time depends on the context of a trial
##       It is users responsibility to align it with trial's parameters
calendarTime(time = 6)

## ----ioea---------------------------------------------------------------------
enrollment(n = 520)

## ----uuea---------------------------------------------------------------------
## condition is based on number of event in the placebo arm
## Note: 'pfs' is the name of the endpoint, i.e., endpoints(name = 'pfs', ...)
## Note: 'pbo' is the name of the placebo arm, i.e., arm(name = 'pbo', ...)
eventNumber(endpoint = 'pfs', n = 340, arms = 'pbo')

## condition can be based on number of event in specific arms
eventNumber(endpoint = 'pfs', n = 340, arms = c('pbo', 'trt'))

## condition can be based on number of event in the trial
## Note: for example, a trial of more than two arms. 
##       Here the number of PFS events are counted on all arms in trial
eventNumber(endpoint = 'pfs', n = 340)

## ----eaiggo-------------------------------------------------------------------
## observe at least 200 OS events on at least 500 enrolled patients, or
## the trial has been running for at least 12 months
(enrollment(n = 500) & eventNumber(endpoint = 'os', n = 200)) | 
  calendarTime(time = 12)

## observe 320 OS or 480 PFS event when the trial has been running for 6 months
calendarTime(time = 6) & 
  (eventNumber('os', n = 320) | eventNumber('pfs', n = 480))

