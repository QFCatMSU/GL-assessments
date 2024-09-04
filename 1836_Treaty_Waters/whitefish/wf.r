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
data$ess <- dat$samp # initial effective sample size
data$obs_pa <- dat$obs_pa # observed proportions at age (age composition data)

# vectorize all the observed data
obs1 <- matrix(data$obs_ct[,1])
obs1 <- data.frame(obs1)
names(obs1) <- "obs"
obs1$year <- data$years
obs1$age <- NA
obs1$type <- "ct"
obs1$fleet <- 1
obs1$obs <- log(obs1$obs)

obs2 <- data.frame(data$obs_pa[,,1], check.names = FALSE)
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

obs3 <- matrix(data$obs_ct[,2])
obs3 <- data.frame(obs3)
names(obs3) <- "obs"
obs3$year <- data$years
obs3$age <- NA
obs3$type <- "ct"
obs3$fleet <- 2
obs3$obs <- log(obs3$obs)

obs4 <- data.frame(data$obs_pa[,,2], check.names = FALSE)
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

data$obs <- obs_all$obs
data$year <- obs_all$year
data$age <- obs_all$age
data$type <- obs_all$type
data$fleet <- obs_all$fleet


##################
## Parameters ####
##################
par <- list()
par$log_sel <- matrix(0, nrow = 2, ncol = data$n_fleet)
par$log_sel[,1] <- rep(-3, data$n_fleet) # p1
par$log_sel[,2] <- rep(6, data$n_fleet) # p2
par$log_q <- rep(0, data$n_fleet)
par$log_q_sd <- rep(log(0.05), data$n_year)
par$log_ct_sd <- rep(log(0.05), data$n_fleet)
par$log_theta <- rep(log(0.5), data$n_fleet)
par$log_qt_devs <- matrix(0, nrow = (data$n_year - 1), ncol = data$n_fleet)
par$log_M <- log(0.2)
par$log_r_sd <- 0
par$log_r_devs <- rep(0, data$n_year - 1)
par$log_r_init <- 10
par$log_n_init <- log(exp(par$log_r_init) * exp(-exp(par$log_M) * (1:5)))


############################
## Additional functions ####
############################
# Dirichelt multinomial likelihood
ddirmultinom <- function(obs, pred, input_n, theta) {
  dir_param <- theta * input_n
  nll <- lgamma(input_n + 1) + lgamma(dir_param) - lgamma(input_n + dir_param)
  nll2 <- 0
  for (a in 1:length(pred)) {
    nll2 <- nll2 - lgamma(obs[a] * input_n + 1) -
      lgamma(obs[a] * input_n + dir_param * pred[a]) +
      lgamma(dir_param * pred[a])
  }
  nll <- nll - nll2
  nll
}


#---------------
# assessment model 
f <- function(par) {
  # get data and pars
  getAll(data, par)

  # back transform, set up vectors and matrices for use later on
  sel_gill_p1 <- exp(log_sel_gill_p1)
  sel_gill_p2 <- exp(log_sel_gill_p2)
  sel_trap_p1 <- exp(log_sel_trap_p1)
  sel_trap_p2 <- exp(log_sel_trap_p2)
  sig <- exp(log_sig)
  phi <- 2 * plogis(t_phi) - 1
  M <- exp(log_M)
  sel_gill <- sel_trap <- matrix(0, nrow = n_year, ncol = n_age)
  log_n <- matrix(-10, nrow = n_year, ncol = n_age)

  # objective function:
  # ! - test case with jnll (no NaNs or NAs)
  jnll <- 0
  
  # catchability 
  for (t in 2:n_year) { # RW 
    jnll <- jnll - dnorm(log_qt_gill[t], log_qt_gill[t - 1], 
                         sig * prior_q_gill_sd, TRUE)
    jnll <- jnll - dnorm(log_qt_trap[t], log_qt_trap[t - 1],
                         sig * prior_q_trap_sd, TRUE)
  }
  # selectivity 
  # ! - test case here - 
  # values to test for sel_ctl == 0
  # sel_gill_p2/trap_p2 <- mean(la[15,])
  # sel_gill_p1/trap_p1 <- 0.02
  # should look like a logistic curve for each year
  # values to test for sel_ctl == 1
  # sel_gill_p2/trap_p2 <- mean(la[15,])
  # sel_gill_p1/trap_p1 <- 90
  if (sel_ctl == 0) { # logistic 
    sel_gill <- 1 / (1 + exp(-sel_gill_p1 * (la - sel_gill_p2)))
    sel_trap <- 1 / (1 + exp(-sel_trap_p1 * (la - sel_trap_p2)))
  } else if (sel_ctl == 1) { # gamma 
    p <- 0.5 * (sqrt(sel_gill_p2^2 + 4 * sel_gill_p1^2) - sel_gill_p2)
    sel_gill <- (la / sel_gill_p2)^(sel_gill_p2 / p) * exp((sel_gill_p2 - la) / p)
    
    p <- 0.5 * (sqrt(sel_trap_p2^2 + 4 * sel_trap_p1^2) - sel_trap_p2)
    sel_trap <- (la / sel_trap_p2)^(sel_trap_p2 / p) * exp((sel_trap_p2 - la) / p)
  }
  
  # recruitment 
  sdr <- sig * prior_sdr
  if (rec_ctl == 0) { # WN 
    jnll <- jnll - sum(dnorm(log_r, 0, sdr, TRUE))
  } else if (rec_ctl == 1) { # RW 
    for (t in 2:n_year) {
      jnll <- jnll - dnorm(log_r[t], log_r[t - 1], sdr, TRUE)
    }
  } else if (rec_ctl == 2) { # AR-1 
    stationary_sd <- sqrt(sdr * sdr / (1 - phi * phi))
    jnll <- jnll - dautoreg(log_r, phi = phi, scale = stationary_sd, log = TRUE)
  }

  # calculate gear specific F and overall Z
  F_gill <- exp(log_qt_gill) * obs_eff_gill * sel_gill
  F_trap <- exp(log_qt_trap) * obs_eff_trap * sel_trap
  F_total <- F_gill + F_trap
  Z <- F_total + M

  # initialize log_n
  log_n[, 1] <- log_r
  log_n[1, 2:(length(log_ninit) + 1)] <- log_ninit

  # population model, including plus group
  for (t in 2:n_year) {
    for (a in 2:n_age) {
      log_n[t, a] <- log_n[t - 1, a - 1] - Z[t - 1, a - 1]
    }
    log_n[t, n_age] <- log(
      exp(log_n[t, n_age]) + # advancing
        exp(log_n[t - 1, n_age]) * exp(-Z[t - 1, n_age]) # there already
    )
  }

  # Baranov's equation to calculate catch and proportions
  ct_gill <- (F_gill / Z) * (1 - exp(-Z)) * exp(log_n)
  ct_gill_total <- rowSums(ct_gill)
  pa_gill <- ct_gill / ct_gill_total
  biomass_gill <- mn_wt_gill * ct_gill
  ct_trap <- (F_trap / Z) * (1 - exp(-Z)) * exp(log_n)
  ct_trap_total <- rowSums(ct_trap)
  pa_trap <- ct_trap / ct_trap_total
  biomass_trap <- mn_wt_trap * ct_trap

  # log catches
  idx <- which(type == "obs_ct_gill")
  jnll <- jnll - sum(dnorm(obs[idx], log(ct_gill_total),
                           sig * prior_ct_gill_sd, TRUE))

  idx <- which(type == "obs_ct_trap")
  jnll <- jnll - sum(dnorm(obs[idx], log(ct_trap_total),
                           sig * prior_ct_trap_sd, TRUE))

  # age comps for each fleet (gear)
  for (f in unique(fleet)) {
    for (y in unique(year)) {
      t <- which(years == y)
      idx1 <- which(fleet == f & year == y & type == "obs_pa_gill")
      if (length(idx1) != 0) {
        jnll <- jnll - dmultinom(obs[idx1], prob = pa_gill[t, ], log = TRUE)
      }
      idx2 <- which(fleet == f & year == y & type == "obs_pa_trap")
      if (length(idx2) != 0) {
        jnll <- jnll - dmultinom(obs[idx2], prob = pa_trap[t, ], log = TRUE)
      }
    }
  }

  # reporting
  REPORT(log_n)
  REPORT(sel_gill)
  REPORT(sel_trap)
  REPORT(F_trap)
  REPORT(F_gill)
  REPORT(F_total)
  REPORT(Z)
  REPORT(ct_gill_total)
  REPORT(ct_trap_total)
  REPORT(pa_gill)
  REPORT(pa_trap)

  jnll
}

data$rec_ctl <- 1 # 0 = WN, 1 = RW, 2 = AR1
if (data$rec_ctl != 2) {
  map <- list(log_M = factor(NA), t_phi = factor(NA))
} else {
  map <- list(log_M = factor(NA))
}

data$sel_ctl <- 0 # 0 = logistic, 1 = gamma

obj <- MakeADFun(f, par, map = map)
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e3, iter.max = 1e3))
sdr <- sdreport(obj)
opt
sdr