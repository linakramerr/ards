## This script simulates a data set to illustrate the joint modeling mediation 
## analysis approache for the paper titled:
## "Interleukin-6 is a mediator of therapeutic efficacy in acute respiratory 
## distress syndrome: an individual patient data meta-analysis of RCTs"

## Author: Lina Kramer

################################################################################
## Simulating data 
#######################################################+#########################

## Note ########################################################################
# This simulation code is an adapted version of code written by Dimitris Rizopoulos
# and Nina van Gerwen. The original code is available here: 
# https://github.com/drizopoulos/JMbayes2/blob/master/Development/MCMC/Surv_Model/simulate_fast.R
################################################################################

# Aim: simulate a data set that is similar to real data 
# We use parameters from the ARMA trial data set 


##  Variables to be generated: 
# record.id: unique patient identifier
# randomized_group: levels "Control" or "Intervention"
# day: measurement time for biomarker sample (range 0-3)
# conc_log10: log10 biomarker concentration
# time_mort28: censoring time (range 0-28)
# death_d28: event indicator (0 = alive or 1 = death)

## Two step approach:

# 1. Simulate the longitudinal biomarker trajectories.
# 2. Simulate the time-to-event data.


## 1. Simulating longitudinal data #############################################

set.seed(15) 

n <- 603 # number of subjects
K <- 2 # number of measurements per subject
t_max <- 28 # maximum follow-up time

# we construct a data frame with the design:
# everyone has a baseline measurement, and then measurements at random 
# follow-up times up to day 3
DF <- data.frame(record.id = rep(seq_len(n), each = K),
                 day = c(replicate(n, c(0, sort(sample(1:3))))),
                 randomized_group = rep(gl(2, n/2, 
                                           labels = c("Control", "Intervention")),
                                        each = K))

# create a new variable: "day_gt0" (day > 0), 
# for modeling effects only after day 0 (so groups are the same at baseline)
DF$day_gt0 <- ifelse(DF$day > 0, DF$day, 0)

# design matrices for the fixed and random effects
X <- model.matrix(~ randomized_group:day_gt0 + day_gt0, data = DF)
Z <- model.matrix(~ day, data = DF)

betas <- c(2.5324, -0.0382, -0.1153) # fixed effects coefficients
sigma <- 0.3971 # errors' standard deviation
D11 <- 0.3653 # variance of random intercepts
D22 <- 0.0063 # variance of random slopes

# we simulate random effects
b <- cbind(rnorm(n, sd = sqrt(D11)), rnorm(n, sd = sqrt(D22)))
# linear predictor
eta_conc_log10 <- as.vector(X %*% betas + rowSums(Z * b[DF$record.id, ]))
# we simulate normal longitudinal data
DF$conc_log10 <- rnorm(n * K, mean = eta_conc_log10, sd = sigma)


## 2. Simulating survival data #################################################
# We assume a Weibull function for the baseline hazard and use the Inverse
# Sampling Technique for survival data to obtain event times.
# The target event rate is 30 %.

# design matrix for the survival model
W <- model.matrix(~ randomized_group, data = DF[!duplicated(DF$record.id), ])
gammas <- c("Intercept" = -9.9, "randomized_group" = -0.0876) # gamma coeffs

eta_t <- as.vector(W %*% gammas)
alpha <- 1.5749
shape_wb <- 2
upp_Cens <- 28

## Function for Inverse Sampling Technique for Survival Data:
invS <- function (t, i, u) {
  # i denotes the subject
  randomized_group_i <- W[i, 2L]
  h <- function (s) { # hazard at time s
    X_at_s <- cbind(1, randomized_group_i, s)
    Z_at_s <- cbind(1, s)
    f <- as.vector(X_at_s %*% betas + rowSums(Z_at_s * b[rep(i, nrow(Z_at_s)), ]))
    exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] + f * alpha)
  }
  ## Then integrate this hazard out over time, and add the
  ## log of a random uniform variable to this
  integrate(h, lower = 0, upper = t)$value + log(u)
  ## When this equation is 0, you get the event time
}


u <- runif(n)
true_eventTimes <- numeric(n)

for(i in seq_len(n)){
  ## random upper limit
  Up <- 100
  Root <- try(uniroot(invS, interval = c(1e-05, Up), i = i, u = u[i])$root,
              TRUE)
  true_eventTimes[i] <- if (!inherits(Root, "try-error")) Root else 29
}

## The observed event time is the minimum of the true event time
## and the upper censoring time (for now)
obs_eventTimes <- pmin(true_eventTimes, upp_Cens)
status_ind <- as.numeric(true_eventTimes <= upp_Cens)

# check event rate (target ~ 30 %)
#summary(status_ind) 
# 27.36 %

DF$time_mort28 <- obs_eventTimes[DF$record.id]
DF$death_d28 <- status_ind[DF$record.id]
## remove longitudinal measurements past the event time 
DF <- DF[DF$day <= DF$time_mort28, ]

# save data
saveRDS(DF, "dataexample.rds")
