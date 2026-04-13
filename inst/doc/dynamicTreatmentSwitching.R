## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  cache.path = 'cache/dynamicTreatmentSwitching/',
  comment = '#>',
  dpi = 300,
  out.width = '100%'
)

## ----setup, echo = FALSE, message = FALSE-------------------------------------
library(TrialSimulator)
library(mvtnorm)
library(dplyr)
library(kableExtra)

## ----eval=FALSE---------------------------------------------------------------
# treatment_allocator <- function(patient_data){...}
# time_selector <- function(patient_data){...}
# data_modifier <- function(patient_data){...}
# 
# regimen <- regimen(treatment_allocator, time_selector, data_modifier)

## ----eval=FALSE---------------------------------------------------------------
# trial <- trial(...)
# trial$add_regime(regimen)

## ----eval=FALSE---------------------------------------------------------------
# trial <- trial(...)
# trial$add_arms(sample_ratio, soc, low_dose, high_dose)
# trial$add_regime(regimen)

## ----echo=FALSE, error=TRUE---------------------------------------------------
try({
msg <- tryCatch(
  {
    stop(' Member function trial$add_regimen() must be called before trial$add_arms(). ', 
         'A good practice is to call trial$add_regimen() immediately after trial() is executed. ')
    NULL
  },
  error = function(e) {
    cat('Error in trial$add_regimen(regimen) :\n', 
        e$message)
  }
)
})

## ----eval=FALSE---------------------------------------------------------------
# treatment_allocator <- function(patient_data){
#   ## add break point to develop and debug
#   # browser()
#   switch_to <- sample(c('low', 'high', 'stay'), nrow(patient_data),
#                         replace = TRUE, prob = c(.3, .4, .3))
#   data.frame(
#     patient_id = patient_data$patient_id,
#     new_treatment =
#       dplyr::case_when(
#         # patient die before progression cannot switch
#         patient_data$os == patient_data$pfs ~ NA_character_,
#         patient_data$arm == 'placebo' & switch_to == 'low' ~ 'low dose',
#         patient_data$arm == 'placebo' & switch_to == 'high' ~ 'high dose',
#         TRUE ~ NA_character_
#       )
#   )
# 
# }

## ----eval=FALSE---------------------------------------------------------------
# time_selector <- function(patient_data){
#   ## add break point to develop and debug
#   # browser()
#   data.frame(
#     patient_id = patient_data$patient_id,
#     ## all patient in patient_data progress before die
#     ## thus pfs < os and can switch.
#     ## See treatment_allocator()
#     switch_time = patient_data$pfs
#   )
# 
# }

## ----eval=FALSE---------------------------------------------------------------
# time_selector <- function(patient_data){
#   ## add break point to develop and debug
#   # browser()
#   data.frame(
#     patient_id = patient_data$patient_id,
#     ## all patient in patient_data progress before die
#     ## thus pfs < os and can switch.
#     ## See treatment_allocator()
#     switch_time = runif(nrow(patient_data), min = patient_data$pfs, max = patient_data$os)
#   )
# 
# }

## ----eval=FALSE---------------------------------------------------------------
# data_modifier <- function(patient_data){
#   ## add break point to develop and debug
#   # browser()
#   f <- ifelse(patient_data$new_treatment == 'low dose', 1.1, 1.15)
#   data.frame(
#     patient_id = patient_data$patient_id,
#     ## other_endpoint = ...,
#     os = patient_data$switch_time + f * (patient_data$os - patient_data$switch_time)
#   )
# 
# }

## ----eval=FALSE---------------------------------------------------------------
# treatment_allocator <- function(patient_data){
#   ## add break point to develop and debug
#   # browser()
#   data.frame(
#     patient_id = patient_data$patient_id,
#     new_treatment =
#       dplyr::case_when(
#         patient_data$arm == 'low dose' & patient_data$response == 0 ~ 'high dose',
#         TRUE ~ NA_character_
#       )
#   )
# 
# }

## ----eval=FALSE---------------------------------------------------------------
# time_selector <- function(patient_data){
#   ## add break point to develop and debug
#   # browser()
#   data.frame(
#     patient_id = patient_data$patient_id,
#     ## all patient in patient_data progress before die
#     ## thus pfs < os and can switch
#     switch_time = patient_data$response_readout
#   )
# 
# }

## -----------------------------------------------------------------------------
treatment_allocator <- function(patient_data){
  ## add break point to develop and debug
  # browser()
  data.frame(
    patient_id = patient_data$patient_id,
    new_treatment =
      dplyr::case_when(
        patient_data$arm == 'placebo' ~ 'new treatment',
        TRUE ~ NA_character_
      )
  )

}

## ----eval=FALSE---------------------------------------------------------------
# time_selector <- function(patient_data){
#   ## add break point to develop and debug
#   # browser()
#   data.frame(
#     patient_id = patient_data$patient_id,
#     ## all patient in patient_data progress before die
#     ## thus pfs < os and can switch
#     switch_time = ifelse(patient_data$os <= 1, .9 * patient_data$os, patient_data$os - 1)
#   )
# 
# }

## ----eval=FALSE---------------------------------------------------------------
# what <- list(allocator1, allocator2, allocator3)
# when <- list(selector1, selector2, selector3)
# how <- list(modifier1, modifier2, modifier3)
# 
# regimen <- regimen(what, when, how)
# 
# trial <- trial(...)
# trial$add_regimen(regimen)

