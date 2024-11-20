rm(list = ls())
gc()
library(RTMB)

# source data script
source("1836_Treaty_Waters/whitefish/data_script.r")

# source functions
source("functions/ddirmultinom.r")
source("functions/selectivity.r")
source("functions/recruitment.r")
source("functions/catchability.r")


##################
## Parameters ####
##################
par = list()
par$log_sel = matrix(0, nrow = 2, ncol = data$n_fleet) # selectivity parameters
par$log_sel[1] = rep(-3, data$n_fleet) # selectivity p1
par$log_sel[2] = rep(6, data$n_fleet) # selectivity p2
par$log_theta = rep(log(0.5), data$n_fleet) # Dirichlet multinomial parameter
par$log_M = log(0.2) # natural mortality
par$log_n_init = log(exp(12) * exp(-exp(par$log_M) * (1:8))) # initial abundance at age
par$log_r = rep(0, data$n_year) # recruitment

# catchability deviations
# standard deviations for catchability deviations
if(data$qt_ctl == "time-varying")
{ 
  par$log_q = numeric(0)
  par$log_qt = matrix(0, nrow = data$n_year, ncol = data$n_fleet)
  par$log_qt_sd = rep(log(0.05), data$n_fleet)
}else 
{
  par$log_q = rep(0, data$n_fleet) 
  par$log_qt = numeric(0)
  par$log_qt_sd = numeric(0)
}
par$log_ct_sd = rep(log(0.05), data$n_fleet) # standard deviations for catch 
par$log_r_sd = 0 # standard deviation for recruitment deviations
# autocorrelation parameter
if(data$rec_ctl == "AR1") 
{
  par$t_phi = 0
}else 
{
  par$t_phi = numeric(0)
}


########################
## Assessment model ####
########################
f = function(par) 
{
  # data list from global to local environment
  # parameters (par) from argument to local environment
  getAll(data, par)
  nobs = length(obs)
  
  # back transform parameters
  sel_pars = exp(log_sel)
  q = exp(log_q)
  theta = exp(log_theta)
  M = exp(log_M)
  n_init = exp(log_n_init)
  if(qt_ctl == "time-varying") qt_sd = exp(log_qt_sd)
  ct_sd = exp(log_ct_sd)
  r_sd = exp(log_r_sd)
  if(rec_ctl == "AR1")
  {
    phi = 2 * plogis(t_phi) - 1
  }else 
  {
    phi = NULL
  }

  # objective function:
  jnll = 0 # joint negative log likelihood
  nll = numeric(4) # individual negative log likelihoods - catchability, recruitment, catch, proportions at age

  # selectivity
  sel = array(0, dim = c(n_year, n_age, n_fleet))
  for(f in 1:n_fleet) 
  {
    for(t in 1:n_year) 
    {
      sel[t,,f] = selectivity(sel_ctl = sel_ctl[f],
                              input = ages,
                              sel_pars = sel_pars[ ,f])
    }
  }
  
  # catchability
  for(f in 1:n_fleet)
  {
    if(qt_ctl == "time-varying") # time-varying
    {
      nll_qt = catchability(qt_ctl = qt_ctl,
                            log_q = NULL,
                            log_qt = log_qt[,f],
                            qt_sd = qt_sd[f],
                            n_year = n_year)
      nll[1] = nll[1] + nll_qt
    
    }else if(qt_ctl == "single")
    {
      log_qt = catchability(qt_ctl = qt_ctl,
                            log_q = log_q[f],
                            log_qt = NULL,
                            qt_sd = NULL,
                            n_year = n_year)
    }
  }

  # calculate mortalities
  F = array(0, dim = c(n_year, n_age, n_fleet))
  for(f in 1:n_fleet) 
  {
    F[,,f] = exp(log_qt[,f]) * obs_eff[,f] * sel[,,f]
  }
  Z = apply(F, MARGIN = c(1,2), FUN = sum) + M

  # recruitment
  nll[2] = recruitment(rec_ctl = rec_ctl,
                        log_r = log_r,
                        r_sd = r_sd,
                        phi = phi)

  # population model, including plus group
  # initialize log_n
  log_n = matrix(-20, nrow = n_year, ncol = n_age)
  log_n[, 1] = log_r
  log_n[1, 2:(length(log_n_init) + 1)] = log_n_init

  # population model, including plus group
  for (t in 2:n_year) 
  {
    for (a in 2:n_age) 
    {
      log_n[t, a] = log_n[t - 1, a - 1] - Z[t - 1, a - 1]
    }
    #advancing + ( * ) there already
    log_n[t, n_age] = log(exp(log_n[t, n_age]) +     
                           exp(log_n[t - 1, n_age]) * exp(-Z[t - 1, n_age]))
  }

  # Baranov's equation to calculate catch and proportions
  log_pred = numeric(length(obs))
  for(i in 1:length(obs))
  {
    y = year_v[i] - min(year_v) + 1
    f = fleet_v[i]
    a = age_v[i] - min(age_v, na.rm = TRUE) + 1

    # first - if the catch data is summarized across ages
    # second - if there is catch at age data
    if(is.na(a)) 
    {
      n_est = exp(log_n[y, ] - log(Z[y, ]) + log(1 - exp(-Z[y, ])) + log(F[y, , f]))
      log_pred[i] = log(sum(n_est))
    }else 
    {
      log_pred[i] = log_n[y, a] - log(Z[y, a]) + log(1 - exp(-Z[y, a])) + log(F[y, a, f])
    }
  }

  # likelihoods
  neff_dm = matrix(0, nrow = n_year, ncol = n_fleet)
  for (f in unique(fleet_v)) 
  {
    for (y in unique(year_v)) 
    {
      t = which(years == y)
      # catch likelihood
      idx = which(fleet_v == f & year_v == y & type_v == "ct")
      nll[3] = nll[3] - dnorm(obs[idx], log_pred[idx], ct_sd[f], log = TRUE)
      
      # proportions at age likelihood
      idx = which(fleet_v == f & year_v == y & type_v == "pa")
      if (length(idx) != 0 & sum(obs[idx]) > 0) 
      {
        if (bio_samp[t,f] > 100) 
        {
          log_est = exp(log_pred[idx]) / sum(exp(log_pred[idx])) # convert to proportions at age
          nll[4] = nll[4] - ddirmultinom(obs[idx], log_est, bio_samp[t,f], theta[f])
          neff_dm[t,f] = 1 / (1 + theta[f]) + ess[t,f] * (theta[f] / (1 + theta[f]))
        }
      }
    }
  }

  # objective function
  jnll = sum(nll)

  # for reporting and plotting
  ct_total = matrix(0, nrow = n_year, ncol = n_fleet)
  pa = array(0, c(n_year, n_age, n_fleet))
  for(f in 1:n_fleet) 
  {
    for (y in unique(year_v)) 
    {
      t = which(years == y)
      idx = which(year_v == y & fleet_v == f & type_v == "ct")
      ct_total[t, f] = exp(log_pred[idx])
      idx = which(year_v == y & fleet_v == f & type_v == "pa")
      pa[t,,f] = exp(log_pred[idx])
    }
  }

  # reporting
  REPORT(log_n)
  REPORT(sel)
  REPORT(F)
  REPORT(Z)
  REPORT(ct_total)
  REPORT(pa)

  return(jnll)
}


#################
## Run model ####
#################
map = list()
map$log_M = factor(NA)
map$log_ct_sd = factor(c(NA))
map$log_r_sd = factor(NA)

obj = MakeADFun(f, par, map = map)
obj$fn()
obj$gr()

opt = nlminb(obj$par, obj$fn, obj$gr, 
            control = list(eval.max = 1e4, iter.max = 1e4))
sdr = sdreport(obj)
print(opt)
print(sdr)

