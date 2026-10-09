## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = '#>',
  fig.width = 7,
  fig.height = 4.5,
  dpi = 120,
  out.width = '100%'
)

## ----setup, echo = FALSE, message = FALSE-------------------------------------
library(kableExtra)
library(TrialSimulator)

## ----multiple-crossings-------------------------------------------------------
cp_formula <- function(D, z, d, alpha, hazard_ratio, omega = 0.25){
  t <- d / D
  theta <- log(hazard_ratio)
  pnorm((qnorm(alpha) - sqrt(t) * z -
           theta * sqrt(omega * D) * (1 - t)) / sqrt(1 - t))
}

z_example <- -1.8
d_example <- 100
target_example <- 0.2
D_grid <- 101:10000
cp_grid <- cp_formula(D_grid, z = z_example, d = d_example,
                      alpha = 0.025, hazard_ratio = 0.98)

crossing_index <- which(diff(cp_grid >= target_example) != 0) + 1
D_grid[crossing_index]

## ----multiple-crossings-table-------------------------------------------------
selected_D <- c(106, 107, 182, 183, 5000, 8674, 8675)
data.frame(
  D = selected_D,
  conditional_power = round(
    cp_formula(selected_D, z = z_example, d = d_example,
               alpha = 0.025, hazard_ratio = 0.98),
    4
  ),
  reaches_target = cp_formula(selected_D, z = z_example, d = d_example,
                              alpha = 0.025,
                              hazard_ratio = 0.98) >= target_example
)

## ----multiple-crossings-plot--------------------------------------------------
plot(D_grid, cp_grid, type = 'l', log = 'x',
     xlab = 'Final event number D (log scale)',
     ylab = 'Conditional power', ylim = c(0, 0.3))
abline(h = target_example, lty = 2, col = 'firebrick')
abline(v = c(107, 183, 8675), lty = 3, col = 'grey40')

## ----define-trial-------------------------------------------------------------
pfs_pbo <- endpoint(name = 'pfs', type = 'tte', generator = rexp,
                    rate = log(2) / 10)
pbo <- arm(name = 'pbo')
pbo$add_endpoints(pfs_pbo)

pfs_trt <- endpoint(name = 'pfs', type = 'tte', generator = rexp,
                    rate = log(2) / 14)
trt <- arm(name = 'trt')
trt$add_endpoints(pfs_trt)

accrual_rate <- data.frame(end_time = Inf, piecewise_rate = 30)
trial <- trial(name = 'cp-example', n_patients = 400,
               seed = 31416, enroller = StaggeredRecruiter,
               accrual_rate = accrual_rate, silent = TRUE)
trial$add_arms(sample_ratio = c(1, 1), pbo, trt)

interim <- milestone(name = 'interim', when = calendarTime(time = 15))
final <- milestone(name = 'final', when = calendarTime(time = 40))

listener <- listener(silent = TRUE)
listener$add_milestones(interim, final)

controller <- controller(trial, listener)
controller$run(n = 1, silent = TRUE, plot_event = FALSE)

## ----interim-fit--------------------------------------------------------------
interim_fit <- fitLogrank(
  Surv(pfs, pfs_event) ~ arm,
  placebo = 'pbo', data = trial$get_locked_data('interim'),
  alternative = 'less', tidy = FALSE
)
interim_fit[, c('arm', 'placebo', 'z', 'info')]

## ----compare-effects----------------------------------------------------------
planned_D <- 300
final_alpha <- 0.022

cp_at <- function(effect){
  trial$conditionalPower(
    'interim', Surv(pfs, pfs_event) ~ arm,
    placebo = 'pbo', alternative = 'less',
    alpha = final_alpha, D = planned_D, effect = effect
  )$cp
}

data.frame(
  assumption = c('interim trend', 'null future drift', 'hazard ratio 0.75'),
  conditional_power = c(cp_at('trend'), cp_at('null'), cp_at(0.75))
)

## ----compare-directions-------------------------------------------------------
cp_less <- trial$conditionalPower(
  'interim', Surv(pfs, pfs_event) ~ arm,
  placebo = 'pbo', alternative = 'less',
  alpha = final_alpha, D = planned_D, effect = 0.75
)

cp_greater <- trial$conditionalPower(
  'interim', Surv(pfs, pfs_event) ~ arm,
  placebo = 'trt', alternative = 'greater',
  alpha = final_alpha, D = planned_D, effect = 1 / 0.75
)

data.frame(
  formulation = c('trt vs pbo, less', 'pbo vs trt, greater'),
  z = c(cp_less$z, cp_greater$z),
  conditional_power = c(cp_less$cp, cp_greater$cp)
)

## ----reassess-events----------------------------------------------------------
reassess_500 <- trial$eventNumberReestimationFromConditionalPower(
  'interim', Surv(pfs, pfs_event) ~ arm,
  placebo = 'pbo', alternative = 'less',
  alpha = final_alpha, target_cp = 0.9,
  effect = 0.75, D_cap = 500
)

reassess_600 <- trial$eventNumberReestimationFromConditionalPower(
  'interim', Surv(pfs, pfs_event) ~ arm,
  placebo = 'pbo', alternative = 'less',
  alpha = final_alpha, target_cp = 0.9,
  effect = 0.75, D_cap = 600
)

rbind(cap_500 = reassess_500, cap_600 = reassess_600)[,
  c('d', 'D', 'D_cap', 'target_cp', 'achieved_cp', 'target_reached')]

## ----reassess-unbounded-------------------------------------------------------
trial$eventNumberReestimationFromConditionalPower(
  'interim', Surv(pfs, pfs_event) ~ arm,
  placebo = 'pbo', alternative = 'less',
  alpha = final_alpha, target_cp = 0.9, effect = 0.75
)[, c('d', 'D', 'D_cap', 'achieved_cp', 'target_reached')]

## ----promising-zone-----------------------------------------------------------
lower_cp <- 0.30
target_cp <- 0.90
cp_planned <- cp_at(0.75)

zone <- if(cp_planned < lower_cp){
  'unfavorable'
}else if(cp_planned < target_cp){
  'promising'
}else{
  'favorable'
}

data.frame(planned_D, cp_planned, zone,
           reassessed_D = reassess_600$D,
           target_reached = reassess_600$target_reached)

## ----promising-zone-action----------------------------------------------------
promising_zone_action <- function(trial){
  ## Pre-specified design parameters. The promising zone is defined using
  ## conditional power at the originally planned final analysis.
  planned_D <- 300
  lower_cp <- 0.30
  target_cp <- 0.90
  D_cap <- 600

  ## conditionalPower() obtains the data locked at the interim milestone
  ## and fits the treatment-vs-placebo logrank comparison internally.
  cp <- trial$conditionalPower(
    'interim', Surv(pfs, pfs_event) ~ arm,
    placebo = 'pbo', alternative = 'less',
    alpha = 0.022, D = planned_D, effect = 0.75
  )$cp

  if(cp < lower_cp){
    zone <- 'unfavorable'
    decision <- 'unfavorable: retain planned event number'
    selected_D <- planned_D
  }else if(cp < target_cp){
    zone <- 'promising'

    ## Only a promising result invokes event-number reassessment. The solver
    ## searches for the earliest solution after the observed interim events.
    reassessment <- trial$eventNumberReestimationFromConditionalPower(
      'interim', Surv(pfs, pfs_event) ~ arm,
      placebo = 'pbo', alternative = 'less',
      alpha = 0.022, target_cp = target_cp,
      effect = 0.75, D_cap = D_cap
    )

    ## This policy permits increases only. Use the exact solution when it is
    ## above planned_D; otherwise use the pre-specified cap as the fallback.
    if(reassessment$target_reached && reassessment$D > planned_D){
      selected_D <- reassessment$D
      decision <- 'promising: reassess to target conditional power'
    }else{
      selected_D <- reassessment$D_cap
      decision <- if(reassessment$target_reached){
        'promising: use cap because the earliest solution is not an increase'
      }else{
        'promising: use cap because no solution was found through the cap'
      }
    }

  }else{
    zone <- 'favorable'
    decision <- 'favorable: retain the planned design'
    selected_D <- planned_D
  }

  ## update_milestone() must be called from a milestone action. It queues a
  ## replacement for the registered, not-yet-triggered final milestone. The
  ## new event target takes effect after this action returns and applies only
  ## to the current replicate; the original trigger is restored before the
  ## next replicate. No update is needed when planned_D is retained.
  if(selected_D > planned_D){
    trial$update_milestone(
      'final',
      when = eventNumber(endpoint = 'pfs', n = selected_D,
                         arms = c('pbo', 'trt'))
    )
  }

  ## Re-evaluate conditional power at the event number that will actually be
  ## used as the final target, including a cap selected as a fallback.
  selected_cp <- trial$conditionalPower(
    'interim', Surv(pfs, pfs_event) ~ arm,
    placebo = 'pbo', alternative = 'less',
    alpha = 0.022, D = selected_D, effect = 0.75
  )$cp

  ## Values saved in an action function become simulation-output columns.
  trial$save(zone, 'promising_zone')
  trial$save(cp, 'conditional_power_at_planned_D')
  trial$save(selected_D, 'selected_final_event_number')
  trial$save(selected_cp, 'conditional_power_at_final_D')
  trial$save(selected_cp >= target_cp, 'target_reached_at_selected_D')
  trial$save(decision, 'promising_zone_decision')
}

## ----promising-zone-final-action----------------------------------------------
promising_zone_final_action <- function(trial){
  fit <- fitLogrank(
    Surv(pfs, pfs_event) ~ arm,
    placebo = 'pbo', data = trial$get_locked_data('final'),
    alternative = 'less', tidy = FALSE
  )

  trial$save(fit$info, 'observed_final_event_number')
  trial$save(fit$z, 'final_z')
  trial$save(fit$p, 'final_p')
  trial$save(fit$p <= 0.022, 'reject_final')
}

## ----promising-zone-program---------------------------------------------------
pfs_pbo <- endpoint(name = 'pfs', type = 'tte', generator = rexp,
                    rate = log(2) / 10)
pbo <- arm(name = 'pbo')
pbo$add_endpoints(pfs_pbo)

pfs_trt <- endpoint(name = 'pfs', type = 'tte', generator = rexp,
                    rate = 0.75 * log(2) / 10)
trt <- arm(name = 'trt')
trt$add_endpoints(pfs_trt)

accrual_rate <- data.frame(end_time = Inf, piecewise_rate = 30)
trial <- trial(
  name = 'promising-zone', n_patients = 700,
  seed = 20260831, enroller = StaggeredRecruiter,
  accrual_rate = accrual_rate, silent = TRUE
)
trial$add_arms(sample_ratio = c(1, 1), pbo, trt)

interim <- milestone(
  name = 'interim', action = promising_zone_action,
  when = eventNumber(endpoint = 'pfs', n = 150,
                     arms = c('pbo', 'trt'))
)
final <- milestone(
  name = 'final', action = promising_zone_final_action,
  when = eventNumber(endpoint = 'pfs', n = 300,
                     arms = c('pbo', 'trt'))
)

listener <- listener(silent = TRUE)
listener$add_milestones(interim, final)

controller <- controller(trial, listener)
controller$run(n = 1, silent = TRUE, plot_event = FALSE)

## ----promising-zone-single-output, echo = FALSE-------------------------------
single_output <- controller$get_output(tidy = TRUE)
single_output_table <- single_output[, c(
  'promising_zone', 'conditional_power_at_planned_D',
  'selected_final_event_number', 'conditional_power_at_final_D',
  'target_reached_at_selected_D', 'observed_final_event_number',
  'final_p', 'reject_final'
)]
single_output_table$final_p <- format.pval(
  single_output_table$final_p, digits = 3, eps = 0.001
)

single_output_table |>
  kable(
    col.names = c(
      'Zone', 'CP at planned D', 'Selected final D', 'CP at final D',
      'Target reached', 'Observed final D', 'Final p-value', 'Reject'
    ),
    digits = 3, align = c('l', 'r', 'r', 'r', 'c', 'r', 'r', 'c'),
    caption = 'Decision path and final analysis for one replicate'
  ) |>
  kable_styling(
    bootstrap_options = 'striped', full_width = FALSE, position = 'left'
  ) |>
  scroll_box(width = '100%')

## ----promising-zone-run, eval = FALSE-----------------------------------------
# controller$reset()
# controller$run(n = 1000, silent = TRUE, plot_event = FALSE)
# conditional_power_output <- controller$get_output()

## ----promising-zone-load-output, echo = FALSE---------------------------------
conditional_power_output <-
  TrialSimulator:::getConditionalPowerOutput()

## ----promising-zone-output-head, echo = FALSE---------------------------------
output_columns <- c(
  'promising_zone', 'conditional_power_at_planned_D',
  'selected_final_event_number', 'conditional_power_at_final_D',
  'target_reached_at_selected_D', 'observed_final_event_number',
  'final_p', 'reject_final'
)

output_preview <- conditional_power_output[seq_len(5), output_columns]
output_preview$final_p <- format.pval(
  output_preview$final_p, digits = 3, eps = 0.001
)

output_preview |>
  kable(
    col.names = c(
      'Zone', 'CP at planned D', 'Selected final D', 'CP at final D',
      'Target reached', 'Observed final D', 'Final p-value', 'Reject'
    ),
    digits = 3, align = c('l', 'r', 'r', 'r', 'c', 'r', 'r', 'c'),
    caption = 'Decision and final-analysis results for the first five replicates'
  ) |>
  kable_styling(
    bootstrap_options = 'striped', full_width = FALSE, position = 'left'
  ) |>
  scroll_box(width = '100%')

## ----promising-zone-zone-summary, echo = FALSE--------------------------------
zone_levels <- c('unfavorable', 'promising', 'favorable')
zone <- factor(conditional_power_output$promising_zone,
               levels = zone_levels)
zone_count <- table(zone)
mean_D_by_zone <- tapply(
  conditional_power_output$selected_final_event_number,
  zone, mean
)

zone_summary <- data.frame(
  zone = zone_levels,
  replicates = as.integer(zone_count),
  percent = round(100 * as.numeric(prop.table(zone_count)), 1),
  mean_selected_D = round(as.numeric(mean_D_by_zone), 1)
)

zone_summary |>
  kable(
    col.names = c('Zone', 'Replicates', 'Percent (%)', 'Mean selected D'),
    digits = 1,
    caption = 'Promising-zone frequency and selected event number'
  ) |>
  kable_styling(
    bootstrap_options = 'striped', full_width = FALSE, position = 'left'
  )

## ----promising-zone-overall-summary, echo = FALSE-----------------------------
overall_summary <- data.frame(
  metric = c(
    'final event target increased',
    'target CP reached at selected D',
    'final null hypothesis rejected'
  ),
  percent = round(100 * c(
    mean(conditional_power_output$selected_final_event_number > 300),
    mean(conditional_power_output$target_reached_at_selected_D),
    mean(conditional_power_output$reject_final)
  ), 1)
)

overall_summary |>
  kable(
    col.names = c('Operating characteristic', 'Percent (%)'), digits = 1,
    caption = 'Operating characteristics of the promising-zone design'
  ) |>
  kable_styling(
    bootstrap_options = 'striped', full_width = FALSE, position = 'left'
  )

## ----multiple-arms-setup, include = FALSE-------------------------------------
trial <- local({
  pfs_pbo <- endpoint(name = 'pfs', type = 'tte', generator = rexp,
                      rate = log(2) / 10)
  pbo <- arm(name = 'pbo')
  pbo$add_endpoints(pfs_pbo)

  pfs_low <- endpoint(name = 'pfs', type = 'tte', generator = rexp,
                      rate = 0.80 * log(2) / 10)
  low <- arm(name = 'low')
  low$add_endpoints(pfs_low)

  pfs_high <- endpoint(name = 'pfs', type = 'tte', generator = rexp,
                       rate = 0.70 * log(2) / 10)
  high <- arm(name = 'high')
  high$add_endpoints(pfs_high)

  accrual_rate <- data.frame(end_time = Inf, piecewise_rate = 40)
  trial <- trial(
    name = 'multiple-arms', n_patients = 1200,
    seed = 27183, enroller = StaggeredRecruiter,
    accrual_rate = accrual_rate, silent = TRUE
  )
  trial$add_arms(sample_ratio = c(1, 1, 1), pbo, low, high)

  interim <- milestone(
    name = 'interim',
    when = eventNumber(endpoint = 'pfs', n = 240,
                       arms = c('pbo', 'low', 'high'))
  )
  listener <- listener(silent = TRUE)
  listener$add_milestones(interim)

  controller <- controller(trial, listener)
  controller$run(n = 1, silent = TRUE, plot_event = FALSE)
  trial
})

## ----multiple-arms-output-----------------------------------------------------
trial$conditionalPower(
  'interim', Surv(pfs, pfs_event) ~ arm,
  placebo = 'pbo', alternative = 'less',
  D = c(low = 350, high = 400),
  alpha = c(high = 0.01, low = 0.015),
  effect = 0.75
)

trial$eventNumberReestimationFromConditionalPower(
  'interim', Surv(pfs, pfs_event) ~ arm,
  placebo = 'pbo', alternative = 'less',
  alpha = c(low = 0.015, high = 0.01),
  target_cp = c(high = 0.90, low = 0.85),
  effect = 'trend',
  D_cap = c(low = 600, high = 700)
)

