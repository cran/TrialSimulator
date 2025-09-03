## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  cache.path = 'cache/wrapper/',
  comment = '#>',
  dpi = 300,
  out.width = '100%'
)

## ----setup, echo = FALSE, message = FALSE-------------------------------------
library(dplyr)
library(TrialSimulator)

## ----ioeajfls-----------------------------------------------------------------
## time-to-event endpoint
pfs <- endpoint(name = 'pfs', type = 'tte', generator = rexp, rate = .07)
## continuous endpoint
cep <- endpoint(name = 'cep', type = 'non-tte', 
                readout = c(cep = 0), generator = rnorm)
## binary endpoint
bep <- endpoint(name = 'bep', type = 'non-tte', 
                readout = c(bep = 0), generator = rbinom, size = 1, prob = .1)

## biomarker
bm <- endpoint(name = 'biomarker', type = 'non-tte', 
               readout = c(biomarker = 0), generator = rbinom, 
               size = 1, prob = .7)

## covariate
covar <- endpoint(name = 'x', type = 'non-tte', 
                  readout = c(x = 0), generator = rnorm)

pbo <- arm(name = 'pbo')
pbo$add_endpoints(pfs, cep, bep, bm, covar)

## ----eualgha, echo=FALSE------------------------------------------------------
## time-to-event endpoint
pfs <- endpoint(name = 'pfs', type = 'tte', generator = rexp, rate = .06)
## continuous endpoint
cep <- endpoint(name = 'cep', type = 'non-tte', 
                readout = c(cep = 0), generator = rnorm, mean = 1.2)
## binary endpoint
bep <- endpoint(name = 'bep', type = 'non-tte', 
                readout = c(bep = 0), generator = rbinom, size = 1, prob = .2)

## biomarker
bm <- endpoint(name = 'biomarker', type = 'non-tte', 
               readout = c(biomarker = 0), generator = rbinom, 
               size = 1, prob = .7)

## covariate
covar <- endpoint(name = 'x', type = 'non-tte', 
                  readout = c(x = 0), generator = rnorm)

low <- arm(name = 'low')
low$add_endpoints(pfs, cep, bep, bm, covar)

## time-to-event endpoint
pfs <- endpoint(name = 'pfs', type = 'tte', generator = rexp, rate = .04)
## continuous endpoint
cep <- endpoint(name = 'cep', type = 'non-tte', 
                readout = c(cep = 0), generator = rnorm, mean = 1.3)
## binary endpoint
bep <- endpoint(name = 'bep', type = 'non-tte', 
                readout = c(bep = 0), generator = rbinom, size = 1, prob = .35)

## biomarker
bm <- endpoint(name = 'biomarker', type = 'non-tte', 
               readout = c(biomarker = 0), generator = rbinom, 
               size = 1, prob = .7)

## covariate
covar <- endpoint(name = 'x', type = 'non-tte', 
                  readout = c(x = 0), generator = rnorm)

high <- arm(name = 'high')
high$add_endpoints(pfs, cep, bep, bm, covar)

accrual_rate <- data.frame(end_time = c(10, Inf),
                           piecewise_rate = c(30, 50))
trial <- trial(
  name = 'Trial-3415', n_patients = 300,
  seed = 1727811904, duration = 1000,
  enroller = StaggeredRecruiter, accrual_rate = accrual_rate,
  dropout = rexp, rate = -log(1 - 0.1)/18, ## 10% by month 18
  silent = TRUE
)

trial$add_arms(sample_ratio = c(1, 1, 1), pbo, low, high)

final <- milestone(name = 'final', action = doNothing, when = calendarTime(1000))

listener <- listener()
listener$add_milestones(final)

controller <- controller(trial, listener)

## ----ieoajf-------------------------------------------------------------------
controller$run(n = 1, plot_event = FALSE, silent = TRUE)
locked_data <- trial$get_locked_data('final')
head(locked_data)

table(locked_data$arm)

## ----ghkaljf------------------------------------------------------------------
## adjust for covariate x
fitCoxph(Surv(pfs, pfs_event) ~ arm + x, placebo = 'pbo', 
         data = locked_data, alternative = 'less', 
         scale = 'hazard ratio')

fitLogrank(Surv(pfs, pfs_event) ~ arm, placebo = 'pbo', 
           data = locked_data, alternative = 'less')

## more details
fitLogrank(Surv(pfs, pfs_event) ~ arm, placebo = 'pbo', 
           data = locked_data, alternative = 'less', tidy = FALSE)

## with strata
fitLogrank(Surv(pfs, pfs_event) ~ arm + strata(biomarker), placebo = 'pbo', 
           data = locked_data, alternative = 'less')

## analyze a subset
fitCoxph(Surv(pfs, pfs_event) ~ arm + strata(biomarker), placebo = 'pbo', 
         data = locked_data, alternative = 'less', 
         scale = 'log hazard ratio', 
         x > -2 & x < 3) ## define a subset

## ----saljfh-------------------------------------------------------------------
## ATE accounting for covariate x
fitLinear(cep ~ arm * x, placebo = 'pbo', 
          data = locked_data, alternative = 'greater')

## marginal model
fitLinear(cep ~ arm, placebo = 'pbo', 
          data = locked_data, alternative = 'greater')

## analyze a sub-group
fitLinear(cep ~ arm, placebo = 'pbo', 
          data = locked_data, alternative = 'greater', 
          biomarker == 1) ## define the subgroup

## ----pqeir--------------------------------------------------------------------
## compute regression coefficient of arm
fitLogistic(bep ~ arm * x + biomarker, placebo = 'pbo', 
            data = locked_data, alternative = 'greater', 
            scale = 'coefficient')

## compute odds ratio (ATE)
fitLogistic(bep ~ arm + x*biomarker, placebo = 'pbo', 
            data = locked_data, alternative = 'greater', 
            scale = 'odds ratio')

## compute risk ratio (ATE)
fitLogistic(bep ~ arm + x + biomarker, placebo = 'pbo', 
            data = locked_data, alternative = 'greater', 
            scale = 'risk ratio')

## ----uyta---------------------------------------------------------------------
## compute risk difference (ATE)
fitLogistic(bep ~ arm + x * biomarker, placebo = 'pbo', 
            data = locked_data, alternative = 'greater', 
            scale = 'risk difference')

## compute risk difference without covariate
fitLogistic(bep ~ arm, placebo = 'pbo', 
            data = locked_data, alternative = 'greater', 
            scale = 'risk difference')

## analyze a sub-group
fitLogistic(bep ~ arm, placebo = 'pbo', 
            data = locked_data, alternative = 'greater', 
            scale = 'risk difference', 
            x < 2 & biomarker != 1) ## define a subgroup

## analyze the same sub-group using the FM test,
## same estimate but different p-values
fitFarringtonManning(endpoint = 'bep', placebo = 'pbo', 
                     data = locked_data, alternative = 'greater', 
                     x < 2 & biomarker != 1)

