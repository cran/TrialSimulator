## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  cache = TRUE,
  cache.path = 'cache/defineNonTimeToEventEndpoints/',
  comment = '#>',
  dpi = 300,
  out.width = '100%'
)

## ----setup, message = FALSE, echo=FALSE---------------------------------------
library(dplyr)
library(TrialSimulator)
set.seed(12345)

## ----lljalk-------------------------------------------------------------------
## endpoints in placebo arm
tumor_cfb_pbo <- endpoint(name = 'cfb', type = 'non-tte', 
                          readout = c(cfb = 6),
                          generator = rnorm, mean = .8, sd = 3.2)
orr_pbo <- endpoint(name = 'orr', type = 'non-tte', 
                    readout = c(orr = 2),
                    generator = rbinom, size = 1, prob = .1)

## define the placebo arm
pbo <- arm(name = 'placebo')
pbo$add_endpoints(tumor_cfb_pbo, orr_pbo)

## endpoints in treatment arm
tumor_cfb_trt <- endpoint(name = 'cfb', type = 'non-tte', 
                          readout = c(cfb = 6), 
                          generator = rnorm, mean = -2.3, sd = 1.5)
orr_trt <- endpoint(name = 'orr', type = 'non-tte', 
                    readout = c(orr = 2), 
                    generator = rbinom, size = 1, prob = .25)

## define the treatment arm
trt <- arm(name = 'treatment')
trt$add_endpoints(tumor_cfb_trt, orr_trt)

## ----oheahs-------------------------------------------------------------------
dropout_pars <- weibullDropout(c(12, 18), c(.15, .30))
dropout_pars

## ----lihoehf------------------------------------------------------------------
accrual_rate <- data.frame(end_time = c(6, Inf), 
                           piecewise_rate = c(10, 20))

trial <- trial(
  name = 'Trial-31415', description = 'Example Clinical Trial', 
  n_patients = 420, duration = 30, 
  enroller = StaggeredRecruiter, accrual_rate = accrual_rate, 
  dropout = rweibull, scale = 30.636, shape = 1.939
)

## add arms to the trial
trial$add_arms(sample_ratio = c(1, 1), trt, pbo)
trial

## ----llieh--------------------------------------------------------------------
interim <- milestone(name = 'interim', 
                     when = eventNumber(endpoint = 'orr', n = 60), 
                     action = doNothing)

random <- milestone(name = 'random', 
                    when = 
                      calendarTime(time = 10) & 
                      (eventNumber(endpoint = 'cfb', n = 100) | 
                         eventNumber(endpoint = 'orr', n = 180)
                       ), 
                    action = doNothing)

final <- milestone(name = 'final', 
                   when = calendarTime(time = 30), 
                   action = doNothing)

## ----ihhf---------------------------------------------------------------------
## register milestones to the listener
listener <- listener()
listener$add_milestones(interim, random, final)

## run the trial
controller <- controller(trial, listener)
controller$run()

## ----ihegea-------------------------------------------------------------------
interim_data <- trial$get_locked_data(milestone_name = 'interim')
random_data <- trial$get_locked_data(milestone_name = 'random')
final_data <- trial$get_locked_data(milestone_name = 'final')
head(interim_data)

## ----eqfee--------------------------------------------------------------------
not_ready_at_interim <- 
  interim_data %>% 
  dplyr::filter(is.na(cfb) & 
                  is.na(orr) & 
                  enroll_time + 6 < dropout_time) %>% 
  head() %>% 
  print()

random_data %>% 
  dplyr::filter(patient_id %in% not_ready_at_interim$patient_id) %>% 
  print()

