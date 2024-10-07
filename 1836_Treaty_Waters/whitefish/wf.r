rm(list = ls())
gc()
library(RTMB)

# source data script
source("1836_Treaty_Waters/whitefish/data_script.r")

# source functions
source("functions/ddirmultinom.r")
source("functions/selectivity.r")


##################
## Parameters ####
##################
par = list()
par$log_sel = matrix(0, nrow = 2, ncol = data$n_fleet) # selectivity parameters
par$log_sel[1] = rep(-3, data$n_fleet) # selectivity p1
par$log_sel[2] = rep(6, data$n_fleet) # selectivity p2
par$log_q = rep(0.1, data$n_fleet) # initial value of catchability
par$log_theta = rep(log(0.5), data$n_fleet) # Dirichlet multinomial parameter
par$log_M = log(0.2) # natural mortality
par$log_n_init = log(exp(12) * exp(-exp(par$log_M) * (1:3))) # initial abundance at age
par$log_r = rep(0, data$n_year) # recruitment

# catchability deviations
# standard deviations for catchability deviations
if(data$qt_ctl == "time-varying")
{ 
  par$log_qt_devs = matrix(0, nrow = data$n_year - 1, ncol = data$n_fleet)
  par$log_qt_sd = rep(log(0.05), data$n_fleet)
  par$log_qt_devs = matrix(0, nrow = data$n_year - 1, ncol = data$n_fleet)
  par$log_qt_sd = rep(log(0.05), data$n_fleet)
}else 
{
  par$log_qt_devs = numeric(0)
  par$log_qt_sd = numeric(0)
  par$log_qt_devs = numeric(0)
  par$log_qt_sd = numeric(0)
}
par$log_ct_sd = rep(log(0.05), data$n_fleet) # standard deviations for catch 
par$log_r_sd = 0 # standard deviation for recruitment deviations
par$log_ct_sd = rep(log(0.05), data$n_fleet) # standard deviations for catch 
par$log_r_sd = 0 # standard deviation for recruitment deviations
# autocorrelation parameter
if(data$rec_ctl == "AR1") 
{
  par$t_phi = 0
  par$t_phi = 0
}else 
{
  par$t_phi = numeric(0)
  par$t_phi = numeric(0)
}


########################
## Assessment model ####
########################
f = function(par) 
f = function(par) 
{
  # data list from global to local environment
  # parameters (par) from argument to local environment
  getAll(data, par)
  
  # back transform parameters
  sel_pars = exp(log_sel)
  q = exp(log_q)
  theta = exp(log_theta)
  M = exp(log_M)
  n_init = exp(log_n_init)
  if(qt_ctl == "time-varying") qt_sd = exp(log_qt_sd)
  ct_sd = exp(log_ct_sd)
  r_sd = exp(log_r_sd)
  if(rec_ctl == "AR1") phi = 2 * plogis(t_phi) - 1

  # objective function:
  jnll = 0 # joint negative log likelihood
  nll = numeric(4) # individual negative log likelihoods - catchability, recruitment, catch, proportions at age
  jnll = 0 # joint negative log likelihood
  nll = numeric(4) # individual negative log likelihoods - catchability, recruitment, catch, proportions at age

  # selectivity
  sel = array(0, dim = c(n_year, n_age, n_fleet))
  for(f in 1:n_fleet) 
  {
    for(t in 1:n_year) 
    {
      sel[t,,f] = selectivity(sel_ctl = sel_ctl[f],
              input = la[t,],
              sel_pars = sel_pars[,f])
    }
  }
  

  # catchability
  log_qt = matrix(0, nrow = n_year, ncol = n_fleet)
  log_qt = matrix(0, nrow = n_year, ncol = n_fleet)
  for(f in 1:n_fleet)
  {
    if(qt_ctl == "time-varying") # time-varying
    {
      # initialize q
      log_qt[1,f] = log_q[f]
      log_qt[1,f] = log_q[f]
      # walk through remaining timesteps for q
      for(t in 2:n_year) 
      {
        log_qt[t,f] = log_qt[t - 1, f] + log_qt_devs[t - 1, f]
        log_qt[t,f] = log_qt[t - 1, f] + log_qt_devs[t - 1, f]
      }
      # likelihood for catchability deviations
      nll[1] = nll[1] - sum(dnorm(log_qt[,f], 0, qt_sd[f], log = TRUE))
    }else if(qt_ctl == "single") # constant q
    {
      for(t in 2:n_year) 
      {
        log_qt[t,f] = log_q[f] 
      }
    }
  } 

  # calculate mortalities
  F = array(0, dim = c(n_year, n_age, n_fleet))
  F = array(0, dim = c(n_year, n_age, n_fleet))
  for(f in 1:n_fleet) 
  {
    F[,,f] = exp(log_qt[,f]) * obs_eff[,f] * sel[,,f]
    F[,,f] = exp(log_qt[,f]) * obs_eff[,f] * sel[,,f]
  }
  Z = apply(F, MARGIN = c(1,2), FUN = sum) + M
  Z = apply(F, MARGIN = c(1,2), FUN = sum) + M

  # recruitment
  if (rec_ctl == "WN") # WN
  { 
    nll[2] = nll[2] - sum(dnorm(log_r, 0, r_sd, TRUE))
  }else if (rec_ctl == "RW") # RW
  { 
    for (t in 2:n_year) 
    {
      nll[2] = nll[2] - dnorm(log_r[t], log_r[t - 1], r_sd, TRUE)
      nll[2] = nll[2] - dnorm(log_r[t], log_r[t - 1], r_sd, TRUE)
    }
  }else if (rec_ctl == "AR1") # AR-1
  { 
    stationary_sd = sqrt(r_sd * r_sd / (1 - phi * phi))
    nll[2] = nll[2] - dautoreg(log_r, phi = phi, scale = stationary_sd, log = TRUE)
    stationary_sd = sqrt(r_sd * r_sd / (1 - phi * phi))
    nll[2] = nll[2] - dautoreg(log_r, phi = phi, scale = stationary_sd, log = TRUE)
  }

  # population model, including plus group
  # initialize log_n
  log_n = matrix(-20, nrow = n_year, ncol = n_age)
  log_n[, 1] = log_r
  log_n[1, 2:(length(log_n_init) + 1)] = log_n_init
  log_n = matrix(-20, nrow = n_year, ncol = n_age)
  log_n[, 1] = log_r
  log_n[1, 2:(length(log_n_init) + 1)] = log_n_init

  # population model, including plus group
  for (t in 2:n_year) 
  {
    for (a in 2:n_age) 
    {
      log_n[t, a] = log_n[t - 1, a - 1] - Z[t - 1, a - 1]
      log_n[t, a] = log_n[t - 1, a - 1] - Z[t - 1, a - 1]
    }
    #advancing + ( * ) there already
    log_n[t, n_age] = log(exp(log_n[t, n_age]) +     # advancing
    log_n[t, n_age] = log(exp(log_n[t, n_age]) +     # advancing
                           exp(log_n[t - 1, n_age]) * exp(-Z[t - 1, n_age]))
                                # there already
  }

  # Baranov's equation to calculate catch and proportions
  ct = array(0, dim = c(n_year, n_age, n_fleet))
  ct_total = matrix(0, nrow = n_year, ncol = n_fleet)
  pa = array(0, dim = c(n_year, n_age, n_fleet))
  ct = array(0, dim = c(n_year, n_age, n_fleet))
  ct_total = matrix(0, nrow = n_year, ncol = n_fleet)
  pa = array(0, dim = c(n_year, n_age, n_fleet))
  for(f in 1:n_fleet) 
  {
    ct[,,f] = F[,,f] / Z * (1 - exp(-Z)) * exp(log_n)
    ct_total[,f] = rowSums(ct[,,f])
    pa[,,f] = ct[,,f] / ct_total[,f]
    ct[,,f] = F[,,f] / Z * (1 - exp(-Z)) * exp(log_n)
    ct_total[,f] = rowSums(ct[,,f])
    pa[,,f] = ct[,,f] / ct_total[,f]
  }

  # likelihoods
  neff_dm = matrix(0, nrow = n_year, ncol = n_fleet)
  neff_dm = matrix(0, nrow = n_year, ncol = n_fleet)
  for (f in unique(fleet)) 
  {
    for (y in unique(year)) 
    {
      t = which(years == y)
      t = which(years == y)
      # catch likelihood
      idx = which(fleet == f & year == y & type == "ct")
      nll[3] = nll[3] - dnorm(obs[idx], ct_total[t,f], ct_sd[f], log = TRUE)
      idx = which(fleet == f & year == y & type == "ct")
      nll[3] = nll[3] - dnorm(obs[idx], ct_total[t,f], ct_sd[f], log = TRUE)
      # proportions at age likelihood
      idx = which(fleet == f & year == y & type == "pa")
      idx = which(fleet == f & year == y & type == "pa")
      if (length(idx) != 0 & sum(obs[idx]) > 0) 
      {
        if (bio_samp[t,f] > 100) 
        {
          nll[4] = nll[4] - ddirmultinom(obs[idx], pa[t,,f] + 1e-8, bio_samp[t,f], theta[f])
          neff_dm[t,f] = 1 / (1 + theta[f]) + ess[t,f] * (theta[f] / (1 + theta[f]))
          nll[4] = nll[4] - ddirmultinom(obs[idx], pa[t,,f] + 1e-8, bio_samp[t,f], theta[f])
          neff_dm[t,f] = 1 / (1 + theta[f]) + ess[t,f] * (theta[f] / (1 + theta[f]))
        }
      }
    }
  }

  jnll = sum(nll)
  jnll = sum(nll)

  # reporting
  REPORT(log_n)
  REPORT(sel)
  REPORT(F)
  REPORT(Z)
  REPORT(ct_total)
  REPORT(pa)

  return(jnll)
}

map = list()
map$log_M = factor(NA)
# map$log_qt_sd = factor(c(NA,NA))
map$log_ct_sd = factor(c(NA,NA))
map$log_r_sd = factor(NA)

obj = MakeADFun(f, par, map = map)
obj$fn()
obj$gr()

opt = nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e3, iter.max = 1e3))
sdr = sdreport(obj)
opt = nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e3, iter.max = 1e3))
sdr = sdreport(obj)
print(opt)
print(sdr)

# pl = obj$report(opt$par)
# matplot(exp(pl$log_n), type = "b")

# plot(data$obs_ct[,1], pch = 19)
# lines(pl$ct_total[,1])
