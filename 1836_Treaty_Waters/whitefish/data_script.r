############
## Data ####
############
#! - before comments, use space not tab
load("1836_Treaty_Waters/whitefish/sim_data.Rdata")
source("functions/format_data.r")
data = list()
data$n_age = dat$aux$n_age        # number of ages
data$n_year = dat$aux$n_year      # number of years
data$n_fleet = dat$aux$n_fleet    # number of fleets
data$ages = dat$aux$ages          # vector of ages
data$years = dat$aux$years # vector of years
data$fleets = dat$aux$fleets # vector of fleets
data$la = dat$la # length at age
data$wa = dat$wa # weight at age
data$mat = dat$mat # maturity at age
obs_eff = matrix(0, nrow = data$n_year, ncol = data$n_fleet)
obs_eff[,1] = dat$obs_eff
data$obs_eff = obs_eff # observed effort
data$obs_ct = dat$obs_ct # observed catch
bio_samp = matrix(0, nrow = data$n_year, ncol = data$n_fleet)
bio_samp[,1] = dat$samp
data$bio_samp = bio_samp # initial effective sample size
ess = matrix(0, nrow = data$n_year, ncol = data$n_fleet)
ess[,1] = dat$samp
data$ess = ess # initial effective sample size
data$obs_pa = dat$obs_pa # observed proportions at age (age composition data)

# vectorize all the observed data
obs_list = list()
type_fleet = data.frame(matrix(0, ncol = 2))
colnames(type_fleet) = c("type", "fleet")
# catch - gillnet
obs_list[[1]] = matrix(data$obs_ct[, 1])
obs_list[[1]] = data.frame(obs_list[[1]])
rownames(obs_list[[1]]) = data$years
type_fleet[1,] = c("ct", 1)

# proportions at age - gillnet
obs_list[[2]] = data.frame(data$obs_pa, check.names = FALSE)
rownames(obs_list[[2]]) = data$years
colnames(obs_list[[2]]) = data$ages
type_fleet[2,] = c("pa", 1)

# combine datasets
obs_all = format_data(d = obs_list, type_fleet = type_fleet)
# incorporate into data list
data$obs = obs_all$obs
data$year = obs_all$year
data$age = obs_all$age
data$type = obs_all$type
data$fleet = obs_all$fleet

# controls for alternative functions
#! define these as characters instead
# recruitment
# 0 = WN, 1 = RW, 2 = AR1
data$rec_ctl = "RW"
# selectivity
# 0 = logistic, 1 = gamma
data$sel_ctl = c("logistic") 
# time-varying catchability
# 0 - single, 1 - time-varying
data$qt_ctl = "time-varying"

# remove all except data list
rm(list = setdiff(ls(), c("data")))
