rm(list = ls())
gc()
library(RTMB)

############
## Data ####
############
load("1836_Treaty_Waters/whitefish/sim_data.Rdata")
data <- list()
data$n_age <- dat$aux$n_age # number of ages
data$n_year <- dat$aux$n_year # number of years
data$n_fleet <- dat$aux$n_fleet # number of fleets
data$ages <- dat$aux$ages # vector of ages
data$years <- dat$aux$years # vector of years
data$fleets <- dat$aux$fleets # vector of fleets
data$la <- dat$la # length at age
data$wa <- dat$wa # weight at age
data$mat <- dat$mat # maturity at age
data$obs_eff <- dat$obs_eff # observed effort
data$obs_ct <- dat$obs_ct # observed catch
data$bio_samp <- dat$samp # initial effective sample size
data$ess <- dat$samp # initial effective sample size
data$obs_pa <- dat$obs_pa # observed proportions at age (age composition data)

# vectorize all the observed data
# catch - gillnet
obs1 <- matrix(data$obs_ct[, 1])
obs1 <- data.frame(obs1)
names(obs1) <- "obs"
obs1$year <- data$years
obs1$age <- NA
obs1$type <- "ct"
obs1$fleet <- 1
obs1$obs <- log(obs1$obs)
# proportions at age - gillnet
obs2 <- data.frame(data$obs_pa[, , 1], check.names = FALSE)
obs2 <- reshape(
  data = obs2,
  direction = "long",
  varying = c(1:data$n_age),
  v.name = "obs",
  times = data$ages,
  timevar = "age",
  ids = data$years,
  idvar = "year"
)
rownames(obs2) <- NULL
obs2$fleet <- 1
obs2$type <- "pa"
obs2 <- obs2[, names(obs1)]
# catch - trapnet
obs3 <- matrix(data$obs_ct[, 2])
obs3 <- data.frame(obs3)
names(obs3) <- "obs"
obs3$year <- data$years
obs3$age <- NA
obs3$type <- "ct"
obs3$fleet <- 2
obs3$obs <- log(obs3$obs)
# proportions at age - trapnet
obs4 <- data.frame(data$obs_pa[, , 2], check.names = FALSE)
obs4 <- reshape(
  data = obs4,
  direction = "long",
  varying = c(1:data$n_age),
  v.name = "obs",
  times = data$ages,
  timevar = "age",
  ids = data$years,
  idvar = "year"
)
rownames(obs4) <- NULL
obs4$fleet <- 2
obs4$type <- "pa"
obs4 <- obs4[, names(obs1)]
obs_all <- rbind(obs1, obs2, obs3, obs4)

# incorporate into data list
data$obs <- obs_all$obs
data$year <- obs_all$year
data$age <- obs_all$age
data$type <- obs_all$type
data$fleet <- obs_all$fleet

# controls for alternative functions
# recruitment
# 0 = WN, 1 = RW, 2 = AR1
data$rec_ctl <- 1 
# selectivity
# 0 = logistic, 1 = gamma
data$sel_ctl <- c(0,0) 
# time-varying catchability
# 0 - single, 1 - time-varying
data$qt_ctl <- 0


##################
## Parameters ####
##################
par <- list()
par$log_sel <- matrix(0, nrow = 2, ncol = data$n_fleet) # selectivity parameters
par$log_sel[, 1] <- rep(-3, data$n_fleet) # selectivity p1
par$log_sel[, 2] <- rep(6, data$n_fleet) # selectivity p2
par$log_q <- rep(0.1, data$n_fleet) # initial value of catchability
par$log_theta <- rep(log(0.5), data$n_fleet) # Dirichlet multinomial parameter
par$log_M <- log(0.2) # natural mortality
par$log_n_init <- log(exp(12) * exp(-exp(par$log_M) * (1:3))) # initial abundance at age
par$log_r <- rep(0, data$n_year) # recruitment
# catchability deviations
# standard deviations for catchability deviations
if(data$qt_ctl == 1)
{ 
  par$log_qt_devs <- matrix(0, nrow = data$n_year - 1, ncol = data$n_fleet)
  par$log_qt_sd <- rep(log(0.05), data$n_fleet)
} else 
{
  par$log_qt_devs <- numeric(0)
  par$log_qt_sd <- numeric(0)
}
par$log_ct_sd <- rep(log(0.05), data$n_fleet) # standard deviations for catch 
par$log_r_sd <- 0 # standard deviation for recruitment deviations
# autocorrelation parameter
if(data$rec_ctl == 2) 
{
  par$t_phi <- 0
} else 
{
  par$t_phi <- numeric(0)
}


############################
## Additional functions ####
############################
# Dirichelt multinomial likelihood
ddirmultinom <- function(obs, pred, input_n, theta) 
{
  dir_param <- theta * input_n
  nll <- lgamma(input_n + 1) + lgamma(dir_param) - lgamma(input_n + dir_param)
  nll2 <- 0
  for (a in 1:length(pred)) 
  {
    nll2 <- nll2 - lgamma(obs[a] * input_n + 1) -
      lgamma(obs[a] * input_n + dir_param * pred[a]) +
      lgamma(dir_param * pred[a])
  }
  nll <- nll - nll2
  nll
}


########################
## Assessment model ####
########################
f <- function(par) 
{
  # get data and pars
  getAll(data, par)

  # back transform parameters
  sel_pars <- exp(log_sel)
  q <- exp(log_q)
  theta <- exp(log_theta)
  M <- exp(log_M)
  n_init <- exp(log_n_init)
  if(qt_ctl == 1) qt_sd <- exp(log_qt_sd)
  ct_sd <- exp(log_ct_sd)
  r_sd <- exp(log_r_sd)
  if(rec_ctl == 2) phi <- 2 * plogis(t_phi) - 1

  # objective function:
  jnll <- 0 # joint negative log likelihood
  nll <- numeric(4) # individual negative log likelihoods - catchability, recruitment, catch, proportions at age

  # selectivity
  sel <- array(0, dim = c(n_year, n_age, n_fleet))
  for(f in 1:n_fleet) 
  {
    # logistic
    if(sel_ctl[f] == 0) 
    {
      sel[,,f] <- 1 / (1 + exp(-sel_pars[f,1] * (la - sel_pars[f,2])))
    } else if(sel_ctl[f] == 1)  # gamma
    {
      p <- 0.5 * (sqrt(sel_pars[f,2]^2 + 4 * sel_pars[f,1]^2) - sel_pars[f,2])
      sel[,,f] <- (la / sel_pars[f,2])^(sel_pars[f,2] / p) * exp((sel_pars[f,2] - la) / p)
    }
  }

  # catchability
  log_qt <- matrix(0, nrow = n_year, ncol = n_fleet)
  for(f in 1:n_fleet)
  {
    if(qt_ctl == 1) # time-varying
    {
      # initialize q
      log_qt[1,f] <- log_q[f]
      # walk through remaining timesteps for q
      for(t in 2:n_year) 
      {
        log_qt[t,f] <- log_qt[t - 1, f] + log_qt_devs[t - 1, f]
      }
      # likelihood for catchability deviations
      nll[1] <- nll[1] - sum(dnorm(log_qt[,f], 0, qt_sd[f], log = TRUE))
    } else if(qt_ctl == 0) # constant q
    {
      for(t in 2:n_year) 
      {
        log_qt[t,f] <- log_q[f] 
      }
    }
  } 

  # calculate mortalities
  F <- array(0, dim = c(n_year, n_age, n_fleet))
  for(f in 1:n_fleet) 
  {
    F[,,f] <- exp(log_qt[,f]) * obs_eff[,f] * sel[,,f]
  }
  Z <- apply(F, MARGIN = c(1,2), FUN = sum) + M

  # recruitment
  if (rec_ctl == 0) # WN
  { 
    nll[2] <- nll[2] - sum(dnorm(log_r, 0, r_sd, TRUE))
  } else if (rec_ctl == 1) # RW
  { 
    for (t in 2:n_year) {
      nll[2] <- nll[2] - dnorm(log_r[t], log_r[t - 1], r_sd, TRUE)
    }
  } else if (rec_ctl == 2) # AR-1
  { 
    stationary_sd <- sqrt(r_sd * r_sd / (1 - phi * phi))
    nll[2] <- nll[2] - dautoreg(log_r, phi = phi, scale = stationary_sd, log = TRUE)
  }

  # population model, including plus group
  # initialize log_n
  log_n <- matrix(-20, nrow = n_year, ncol = n_age)
  log_n[, 1] <- log_r
  log_n[1, 2:(length(log_n_init) + 1)] <- log_n_init

  # population model, including plus group
  for (t in 2:n_year) 
  {
    for (a in 2:n_age) 
    {
      log_n[t, a] <- log_n[t - 1, a - 1] - Z[t - 1, a - 1]
    }
    log_n[t, n_age] <- log(
      exp(log_n[t, n_age]) + # advancing
        exp(log_n[t - 1, n_age]) * exp(-Z[t - 1, n_age]) # there already
    )
  }

  # Baranov's equation to calculate catch and proportions
  ct <- array(0, dim = c(n_year, n_age, n_fleet))
  ct_total <- matrix(0, nrow = n_year, ncol = n_fleet)
  pa <- array(0, dim = c(n_year, n_age, n_fleet))
  for(f in 1:n_fleet) 
  {
    ct[,,f] <- F[,,f] / Z * (1 - exp(-Z)) * exp(log_n)
    ct_total[,f] <- rowSums(ct[,,f])
    pa[,,f] <- ct[,,f] / ct_total[,f]
  }

  # likelihoods
  neff_dm <- matrix(0, nrow = n_year, ncol = n_fleet)
  for (f in unique(fleet)) 
  {
    for (y in unique(year)) 
    {
      t <- which(years == y)
      # catch likelihood
      idx <- which(fleet == f & year == y & type == "ct")
      nll[3] <- nll[3] - dnorm(obs[idx], ct_total[t,f], ct_sd[f], log = TRUE)
      # proportions at age likelihood
      idx <- which(fleet == f & year == y & type == "pa")
      if (length(idx) != 0 & sum(obs[idx]) > 0) 
      {
        if (bio_samp[t,f] > 100) 
        {
          nll[4] <- nll[4] - ddirmultinom(obs[idx], pa[t,,f] + 1e-8, bio_samp[t,f], theta[f])
          neff_dm[t,f] <- 1 / (1 + theta[f]) + ess[t,f] * (theta[f] / (1 + theta[f]))
        }
      }
    }
  }

  jnll <- sum(nll)

  # reporting
  REPORT(log_n)
  REPORT(sel)
  REPORT(F)
  REPORT(Z)
  REPORT(ct_total)
  REPORT(pa)

  jnll
}

map <- list()
map$log_M <- factor(NA)
# map$log_qt_sd <- factor(c(NA,NA))
map$log_ct_sd <- factor(c(NA,NA))
map$log_r_sd <- factor(NA)
obj <- MakeADFun(f, par, map = map)
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e3, iter.max = 1e3))
sdr <- sdreport(obj)
opt
sdr

# pl <- obj$report(opt$par)
# matplot(exp(pl$log_n), type = "b")

# plot(data$obs_ct[,1], pch = 19)
# lines(pl$ct_total[,1])
