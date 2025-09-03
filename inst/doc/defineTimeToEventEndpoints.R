## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  cache = TRUE,
  cache.path = 'cache/defineTimeToEventEndpoints/',
  comment = '#>',
  dpi = 300,
  out.width = '100%'
)

## ----setup, message = FALSE, echo = FALSE-------------------------------------
library(dplyr)
library(survival)
library(survminer)
library(simdata)
library(TrialSimulator)
set.seed(12345)

## ----lkioie-------------------------------------------------------------------
pfs_pbo <- endpoint(name = 'PFS', type = 'tte', 
                    generator = rexp, rate = log(2)/5.6)

## ----adfla, message=FALSE-----------------------------------------------------
test_set <- pfs_pbo$test_generator(n = 1e5)
head(test_set)
median(test_set$PFS) ## should be close to 5.6

## ----ieaofdsaf, results='asis'------------------------------------------------
pfs_pbo

## ----aldj---------------------------------------------------------------------
pfs_trt <- endpoint(name = 'PFS', type = 'tte', 
                    generator = rexp, rate = log(2)/6.4)
median(pfs_trt$test_generator(n = 1e5)$PFS) ## should be close to 6.4

## ----alfjd--------------------------------------------------------------------
pbo <- arm(name = 'placebo')
pbo$add_endpoints(pfs_pbo)
trt <- arm(name = 'treatment')
trt$add_endpoints(pfs_trt)

## ----aldjfba------------------------------------------------------------------
risk_pbo <- data.frame(
  end_time = c(2, 8, 10), 
  piecewise_risk = c(1, 0.48, 0.25) * exp(-1)
)

pfs_pbo <- endpoint(name = 'PFS', type = 'tte', 
                    generator = PiecewiseConstantExponentialRNG, 
                    risk = risk_pbo, 
                    endpoint_name = 'PFS')
risk_trt <- risk_pbo %>% 
  mutate(hazard_ratio = c(1, .6, .7))

pfs_trt <- endpoint(name = 'PFS', type = 'tte', 
                    generator = PiecewiseConstantExponentialRNG, 
                    risk = risk_trt, 
                    endpoint_name = 'PFS')

test_set <- rbind(pfs_pbo$test_generator(n = 1e4) %>% mutate(arm = 'pbo'), 
                  pfs_trt$test_generator(n = 1e4) %>% mutate(arm = 'trt'))

sfit <- survfit(Surv(time = PFS, event = PFS_event) ~ arm, test_set)
ggsurvplot(sfit, data = test_set, palette = c("blue", "red"))

## ----ldao---------------------------------------------------------------------
head(test_set %>% slice_sample(prop = 1))

## ----eegj---------------------------------------------------------------------
pbo <- arm(name = 'placebo')
pbo$add_endpoints(pfs_pbo)

trt <- arm(name = 'treatment')
trt$add_endpoints(pfs_trt)

## ----laiojb-------------------------------------------------------------------
os_pbo <- endpoint(name = 'OS', type = 'tte', 
                   generator = rexp, rate = log(2)/7.2)
os_trt <- endpoint(name = 'OS', type = 'tte', 
                   generator = rexp, rate = log(2)/8.5)

median(os_pbo$test_generator(n = 1e5)$OS) ## should be close to 7.2
median(os_trt$test_generator(n = 1e5)$OS) ## should be close to 8.5

## add endpoint to existing arms
pbo$add_endpoints(os_pbo)
trt$add_endpoints(os_trt)

## ----albaiola-----------------------------------------------------------------
custom_generator <- function(n, pfs_rate, os_rate, corr){
  
  dist <- list()
  dist[['PFS']] <- function(x) qexp(x, rate = pfs_rate)
  dist[['OS']] <- function(x) qexp(x, rate = os_rate)
  dsgn = simdata::simdesign_norta(cor_target_final = corr, 
                                dist = dist, 
                                transform_initial = data.frame,
                                names_final = names(dist), 
                                seed_initial = 1)
  
  simdata::simulate_data(dsgn, n_obs = n) %>% 
    mutate(PFS_event = 1, OS_event = 1) ## event indicators
}

## ----ojhonln------------------------------------------------------------------
corr <- matrix(c(1, .6, .6, 1), nrow = 2)
eps_pbo <- endpoint(name = c('PFS', 'OS'), type = c('tte', 'tte'), 
                    generator = custom_generator, 
                    pfs_rate = log(2)/5.6, os_rate = log(2)/7.2, 
                    corr = corr)

eps_trt <- endpoint(name = c('OS', 'PFS'), type = c('tte', 'tte'), 
                    generator = custom_generator, 
                    pfs_rate = log(2)/6.4, os_rate = log(2)/8.5, 
                    corr = corr)

test_set <- rbind(eps_pbo$test_generator(n = 1e5) %>% mutate(arm = 'pbo'), 
                  eps_trt$test_generator(n = 1e5) %>% mutate(arm = 'trt'))

with(test_set, cor(PFS, OS)) ## should be close to 0.6

## sample medians match to the parameters well
test_set %>% 
  group_by(arm) %>% 
  summarise(PFS = median(PFS), OS = median(OS))


## ----daadidal, results='asis'-------------------------------------------------
eps_pbo

## ----ieaofafa, results='asis'-------------------------------------------------
eps_trt

## ----lioe---------------------------------------------------------------------
pbo <- arm(name = 'placebo')
pbo$add_endpoints(eps_pbo)

trt <- arm(name = 'treatment')
trt$add_endpoints(eps_trt)

