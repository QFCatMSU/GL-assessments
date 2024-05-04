library(RTMB)

#---------------
# assessment model 
f <- function(par) {
  # get data and pars, specify observations, ADoverload trick
  getAll(data, par)
  obs <- OBS(obs)
  "[<-" <- ADoverload("[<-")

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
  jnll <- 0
  
  # catchability 
  for (t in 2:n_year) { # RW 
    jnll <- jnll - dnorm(log_qt_gill[t], log_qt_gill[t - 1], 
                         sig * prior_q_gill_sd, TRUE)
    jnll <- jnll - dnorm(log_qt_trap[t], log_qt_trap[t - 1],
                         sig * prior_q_trap_sd, TRUE)
  }
  # selectivity 
  if (sel_ctl == 0) { # logistic 
    sel_gill <- 1 / (1 + exp(-sel_gill_p1 * (la - sel_gill_p2)))
    sel_trap <- 1 / (1 + exp(-sel_trap_p1 * (la - sel_trap_p2)))
  }
  if (sel_ctl == 1) { # gamma 
    p <- 0.5 * (sqrt(sel_gill_p2 + 4 * sel_gill_p1^2) - sel_gill_p2)
    sel_gill <- (la / sel_gill_p2)^(sel_gill_p2 / p) * exp((sel_gill_p2 - la) / p)
    
    p <- 0.5 * (sqrt(sel_trap_p2 + 4 * sel_trap_p1^2) - sel_trap_p2)
    sel_trap <- (la / sel_trap_p2)^(sel_trap_p2 / p) * exp((sel_trap_p2 - la) / p)
  }
  
  # recruitment 
  sdr <- sig * prior_sdr
  if (rec_ctl == 0) { # WN 
    jnll <- jnll - sum(dnorm(log_r, 0, sdr, TRUE))
  }
  if (rec_ctl == 1) { # RW 
    for (t in 2:n_year) {
      jnll <- jnll - dnorm(log_r[t], log_r[t - 1], sdr, TRUE)
    }
  }
  if (rec_ctl == 2) { # AR-1 
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


#---------------
# retrospective analysis
nrun <- 5 # number of peels
output <- list()
sdout <- list()

for(peel in 1:nrun) {
  message("Minus ", peel, " years(s)")
  load("data/sim_data.RData")
  data <- list()
  data$n_age <- length(dat$aux$fage:dat$aux$lage)
  data$n_year <- length(dat$aux$fyear:(dat$aux$lyear - peel))
  peel_idx <- 1:data$n_year # peel index
  data$ages <- dat$aux$fage:dat$aux$lage
  data$years <- dat$aux$fyear:(dat$aux$lyear - peel)
  data$la <- dat$la[peel_idx, ]
  data$wa <- dat$wt[peel_idx, ]
  data$obs_eff_trap <- dat$effort_trap[peel_idx]
  data$obs_eff_gill <- dat$effort_gill[peel_idx] * dat$depth_adj[peel_idx]
  data$mn_wt_trap <- dat$mean_wt_kg_trap[peel_idx]
  data$mn_wt_gill <- dat$mean_wt_kg_gill[peel_idx]
  data$prior_sdr <- dat$aux$rhoSR
  data$prior_ct_trap_sd <- dat$aux$rhoCT
  data$prior_q_trap_sd <- dat$aux$rhoET
  data$prior_ct_gill_sd <- dat$aux$rhoCG
  data$prior_q_gill_sd <- dat$aux$rhoEG
  data$ess_trap <- dat$ESS_trap[peel_idx]
  data$ess_gill <- dat$ESS_gill[peel_idx]

  obs1 <- matrix((dat$biomass_gill[peel_idx] / dat$mean_wt_kg_gill[peel_idx]) / dat$gill_adj[peel_idx])
  obs1 <- data.frame(obs1)
  names(obs1) <- "obs"
  obs1$year <- data$years
  obs1$age <- NA
  obs1$type <- "obs_ct_gill"
  obs1$fleet <- 1
  obs1$obs <- log(obs1$obs)

  obs2 <- data.frame(dat$pa_gill[peel_idx, ] * dat$ESS_gill[peel_idx], check.names = FALSE)
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
  obs2$type <- "obs_pa_gill"
  obs2 <- obs2[, names(obs1)]

  obs3 <- matrix((dat$biomass_trap[peel_idx] / dat$mean_wt_kg_trap[peel_idx]) / dat$trap_adj[peel_idx])
  obs3 <- data.frame(obs3)
  names(obs3) <- "obs"
  obs3$year <- data$years
  obs3$age <- NA
  obs3$type <- "obs_ct_trap"
  obs3$fleet <- 2
  obs3$obs <- log(obs3$obs)

  obs4 <- data.frame(dat$pa_trap[peel_idx, ] * dat$ESS_trap[peel_idx], check.names = FALSE)
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
  obs4$type <- "obs_pa_trap"
  obs4 <- obs4[, names(obs1)]
  obs_all <- rbind(obs1, obs2, obs3, obs4)

  data$obs <- obs_all$obs
  data$year <- obs_all$year
  data$age <- obs_all$age
  data$type <- obs_all$type
  data$fleet <- obs_all$fleet

  par <- list()
  par$log_sel_trap_p1 <- -3
  par$log_sel_trap_p2 <- 5.7
  par$log_sel_gill_p1 <- -3
  par$log_sel_gill_p2 <- 6
  par$log_sig <- 0
  par$log_qt_trap <- rep(0, data$n_year)
  par$log_qt_gill <- rep(0, data$n_year)
  par$log_ninit <- rep(0, 3)
  par$log_r <- rep(0, data$n_year)
  par$t_phi <- 0
  par$log_M <- log(0.2)

  data$rec_ctl <- 1 # 0 = WN, 1 = RW, 2 = AR1
  if (data$rec_ctl != 2) {
    map <- list(log_M = factor(NA), t_phi = factor(NA))
  } else {
    map <- list(log_M = factor(NA))
  }
  data$sel_ctl <- 0 # 0 = logistic, 1 = gamma

  obj <- MakeADFun(f, par, map = map)
  opt <- nlminb(obj$par, obj$fn, obj$gr,
    control = list(eval.max = 1e3, iter.max = 1e3)
  )
  output[[peel]] <- obj$report(opt$par)
  sdout[[peel]] <- sdreport(obj)
}

# plots
# abundance
abund <- lapply(output, function(x) rowSums(exp(x$log_n)))
year_minus <- 1:(length(dat$aux$fyear:dat$aux$lyear) - 1)
plot(year_minus, seq(0, max(unlist(abund)) * 1.1, length.out = length(year_minus)),
  type = "n", ylab = "numbers at age", xlab = "years"
)
col_vec <- rainbow(length(abund))
for (i in 1:length(abund)) {
  lines(abund[[i]], col = col_vec[i])
}

# recruitment
rec <- lapply(output, function(x) exp(x$log_n[, 1]))
plot(year_minus, seq(0, max(unlist(rec)), length.out = length(year_minus)),
  type = "n", ylab = "recruitment", xlab = "years"
)
col_vec <- rainbow(length(rec))
for (i in 1:length(rec)) {
  lines(rec[[i]], col = col_vec[i])
}

# F (mean)
Ftot <- lapply(output, function(x) rowMeans(x$F_total))
plot(year_minus, seq(0, max(unlist(Ftot)), length.out = length(year_minus)),
  type = "n", ylab = "total fishing mortality (mean)", xlab = "years"
)
col_vec <- rainbow(length(Ftot))
for (i in 1:length(Ftot)) {
  lines(Ftot[[i]], col = col_vec[i])
}
