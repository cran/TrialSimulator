## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  cache.path = 'cache/crossoverAtMilestone/',
  comment = '#>',
  dpi = 300,
  out.width = '100%'
)

## ----setup, echo = FALSE, message = FALSE-------------------------------------
library(TrialSimulator)

## ----eval=FALSE---------------------------------------------------------------
# action <- function(trial){
#   locked_data <- trial$get_locked_data('interim')
#   ## ... interim analysis / decision making ...
#   trial$crossover(what = what_fn, how = how_fn)
# }
# 
# interim <- milestone(name = 'interim',
#                      when = eventNumber(endpoint = 'pfs', n = 300),
#                      action = action)

## ----eval=FALSE---------------------------------------------------------------
# crossover(trial, what = what_fn, how = how_fn)

## ----eval=FALSE---------------------------------------------------------------
# what(patient_data)   # -> patient_id, new_treatment
# when(patient_data)   # -> patient_id, switch_time
# how(patient_data)    # -> patient_id, <modified endpoints>

## ----funnel, echo=FALSE, fig.width=7.2, fig.height=4.2, out.width="95%", fig.alt="Funnel: trial_data is filtered to eligible patients (passed to what()), then to switchers (passed to when()), then to switchers with their switch_time (passed to how())."----
op <- par(mar = c(0, 0, 0, 0))
plot.new(); plot.window(xlim = c(0, 10), ylim = c(0, 10))
cx <- 5
yb <- c(9.3, 7.1, 4.9, 2.7, 0.6)            # band boundaries (top -> bottom)
wt <- c(3.4, 2.6, 1.8, 1.0)                 # band top widths
wd <- c(2.6, 1.8, 1.0, 0.5)                 # band bottom widths
for (i in 1:4) {
  polygon(cx + c(-wt[i], wt[i], wd[i], -wd[i]) / 2,
          c(yb[i], yb[i], yb[i + 1], yb[i + 1]),
          col = "#e8eef5", border = "#33618f", lwd = 1.4)
}
ymid <- (head(yb, -1) + tail(yb, -1)) / 2
left <- c("all patients\n(trial_data)",
          "eligible\n(open endpoint at T)",
          "switchers\n(what() selects)",
          "+ switch_time")
text(cx - 2.0, ymid, left, adj = 1, cex = 0.76)
fn <- c("what()", "when()", "how()")
for (i in 1:3) {
  arrows(cx + wt[i + 1] / 2 + 0.1, ymid[i + 1], cx + 2.0, ymid[i + 1],
         length = 0.07, col = "#33618f")
  text(cx + 2.1, ymid[i + 1], paste0("passed to ", fn[i]), adj = 0, cex = 0.82)
}
par(op)

## ----eval=FALSE---------------------------------------------------------------
# time_selector <- function(patient_data){
#   data.frame(
#     patient_id  = patient_data$patient_id,
#     switch_time = pmax(patient_data$pfs,
#                        patient_data$earliest_crossover_time_from_enrollment)
#   )
# }

## ----eval=FALSE---------------------------------------------------------------
# data_modifier <- function(patient_data){
#   data.frame(
#     patient_id = patient_data$patient_id,
#     ## extend only the residual (post-switch) survival; leave os unchanged
#     ## for patients whose event is at or before the switch
#     os = ifelse(patient_data$os > patient_data$switch_time,
#                 patient_data$switch_time +
#                   1.2 * (patient_data$os - patient_data$switch_time),
#                 patient_data$os)
#   )
# }

## ----eval=FALSE---------------------------------------------------------------
# what <- function(patient_data){
#   browser()                       # pause here with the eligible pool in scope
#   switchers <- patient_data[patient_data$arm == 'control', ]
#   data.frame(patient_id = switchers$patient_id, new_treatment = 'experimental')
# }

## -----------------------------------------------------------------------------
what <- function(patient_data){
  ## return only the patients who switch (here, everyone on control)
  switchers <- patient_data[patient_data$arm == 'control', ]
  data.frame(
    patient_id    = switchers$patient_id,
    new_treatment = 'experimental'
  )
}

when <- function(patient_data){
  data.frame(
    patient_id  = patient_data$patient_id,
    switch_time = pmax(patient_data$pfs,
                       patient_data$earliest_crossover_time_from_enrollment)
  )
}

how <- function(patient_data){
  data.frame(
    patient_id = patient_data$patient_id,
    os = ifelse(patient_data$os > patient_data$switch_time,
                patient_data$switch_time +
                  1.3 * (patient_data$os - patient_data$switch_time),
                patient_data$os)
  )
}

crossover_action <- function(trial){
  trial$get_locked_data('interim')        # interim decision making can go here
  ## delay = 0 opens crossover at the interim; use delay > 0 for a wash-out
  trial$crossover(what = what, when = when, how = how, delay = 0)
}

## ----echo=FALSE---------------------------------------------------------------
os_e  <- endpoint(name = 'os',  type = 'tte', generator = rexp, rate = log(2) / 12)
pfs_e <- endpoint(name = 'pfs', type = 'tte', generator = rexp, rate = log(2) / 6)
control <- arm(name = 'control', pfs <= os); control$add_endpoints(pfs_e, os_e)
trt     <- arm(name = 'trt',     pfs <= os); trt$add_endpoints(pfs_e, os_e)

trial <- trial(name = 'demo', n_patients = 200, seed = 123, duration = 36,
               enroller = StaggeredRecruiter,
               accrual_rate = data.frame(end_time = Inf, piecewise_rate = 15),
               silent = TRUE)
trial$add_arms(sample_ratio = c(1, 1), control, trt)

lst <- listener(silent = TRUE)
lst$add_milestones(
  milestone(name = 'interim', when = calendarTime(time = 12), action = crossover_action),
  milestone(name = 'final',   when = calendarTime(time = 36))
)
invisible(controller(trial, lst)$run(n = 1, silent = TRUE, plot_event = FALSE))

## -----------------------------------------------------------------------------
final <- trial$get_locked_data('final')
head(final[final$arm == 'control',
           c('patient_id', 'arm', 'regimen_trajectory', 'n_switches')])

## ----eval=FALSE---------------------------------------------------------------
# trial$crossover(what = what, how = how, delay = 2) # opens 2 time units later

## -----------------------------------------------------------------------------
head(expandRegimen(final))

