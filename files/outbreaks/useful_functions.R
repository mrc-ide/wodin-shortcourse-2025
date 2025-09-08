get_counts_from_cumulative <- function(t, cum, dt)
{
  tail_t <- t[t > max(t, na.rm = TRUE) - dt]
  n_tail <- length(tail_t)
  c(rep(NA, n_tail), cum[-seq_len(n_tail)] - head(cum, length(cum)-n_tail))
}

sse <- function(mod, dat)
{
  (mod - dat) ^ 2
}

pois_ll <- function(mod, dat)
{
  # using minus because optim does minimisation and we want to maximise the likelihood
  - dpois(round(dat), lambda = mod, log = TRUE)
}

error <- function(mod_output, output_to_fit, data_to_fit, fun_to_optimise = sse)
{
  res <- 0
  for(i in seq_len(length(output_to_fit)))
  {
    # find exact time step in output and throw error if doesn't exist
    times <- match(data_to_fit[[output_to_fit[i]]]$t, mod_output$t)
    if(any(is.na(times))) stop("Some data time steps don't have a corresponding time step in model output.")
    res_i <- fun_to_optimise(mod_output[[output_to_fit[i]]][times], data_to_fit[[output_to_fit[i]]]$x)
    res_i <- sum(res_i, na.rm = TRUE)
    res <- res + res_i
  }
  res
}

error_model <- function(params, 
                       params_names = c("I0", "mu_h_before", "cfr", "R0"),
                       data_to_fit,
                       output_to_fit = c("weekly_onset", "weekly_death_h"),
                       gen, 
                       tmax,
                       dt,
                      time_output, 
                      fun_to_optimise = sse)
{
  #print(params)
  ## TO DO: add check as params should be same length as params_names
  
  params_list <- setNames(as.list(params), params_names)
  mod <- do.call(gen, params_list) 
  
  mod_output <- post_process(mod, time_output)
  
  res <- error(mod_output, output_to_fit, data_to_fit, fun_to_optimise)
  if(res == Inf) res <- 1e6 # because sometimes 0 likelihood and 
  print(c(params, res))
  res
}

post_process <- function(mod, time_output)
{
  
  ### run the model and transform model output
  tmp_mod_output <- mod$run(time_output)
  mod_output <- mod$transform_variables(tmp_mod_output)
  
  # this is commented out as has been incorporated within the odin code now 
  # (in a better way)
  
  ### compute weekly number of incident cases by date of onset
  #mod_output$weekly_onset_old <- get_counts_from_cumulative(mod_output$t,
  #                                                          mod_output$cumul_onset, 7) 
  
  ### compute weekly number of deaths
  #mod_output$weekly_death_h_old <- get_counts_from_cumulative(mod_output$t,
  #                                                            mod_output$cumul_death_h, 7) 
  
  mod_output
}

plot_summary_model <- function(mod_output, dat = NULL, with_cumulative = FALSE)
{
  if(with_cumulative) par(mfrow = c(2, 2)) else par(mfrow = c(2, 1)) 
  ### plot epi curve by date of onset
  # cumulative
  if(with_cumulative) 
    plot(mod_output$t, mod_output$cumul_onset, type = "l",
         xlab = "Time (days)", ylab = "Cumulative number of cases (by onset date)")
  # weekly
  if(!is.null(dat))
  {
    ymax <- max(c(dat$weekly_onset$incidence, mod_output$weekly_onset), 
                na.rm = TRUE)
  } else
  {
    ymax <- max(mod_output$weekly_onset, 
                na.rm = TRUE)
  }
  plot(mod_output$t, mod_output$weekly_onset, type = "l", ylim = c(0, ymax),
       xlab = "Time (days)", ylab = "Weekly number of cases (by onset date)") 
  if(!is.null(dat))
    points(dat$weekly_onset$date, dat$weekly_onset$incidence, col = "red")
  
  ### plot epi curve for the dead by date of death
  # cumulative
  if(with_cumulative) 
    plot(mod_output$t, mod_output$cumul_death_h, type = "l",
         xlab = "Time (days)", ylab = "Cumulative number of deaths")
  # weekly
  if(!is.null(dat))
  {
    ymax <- max(c(dat$weekly_death$incidence, mod_output$weekly_death_h), 
                na.rm = TRUE)
  } else
  {
    ymax <- max(mod_output$weekly_death_h,  
                na.rm = TRUE)
  }
  plot(mod_output$t, mod_output$weekly_death_h, type = "l", ylim = c(0, ymax),
       xlab = "Time (days)", ylab = "Weekly number of deaths") 
  if(!is.null(dat)) points(dat$weekly_death$date, dat$weekly_death$incidence, col = "red")
}


post_process_dat <- function(dat_raw, reporting_prob = 1)
{
  dat <- dat_raw
  dat$weekly_onset$date <- dat$weekly_onset$date + 1 # to have one step of model run before starting to fit to incidence 
  # (otherwise initial incidence in model is zero which we never manage to fit to non zero observed incidence)
  dat$weekly_death$date <- dat$weekly_death$date + 1 # same as above
  dat$weekly_onset$incidence <- dat_raw$weekly_onset$incidence / reporting_prob
  dat$weekly_death$incidence <- dat_raw$weekly_death$incidence / reporting_prob
  dat
}
