## This script uses a fictional data set to illustrate the joint modeling mediation 
## analysis approache for the paper titled:
## "Interleukin-6 is a mediator of therapeutic efficacy in acute respiratory 
## distress syndrome: an individual patient data meta-analysis of RCTs"


## Author: Lina Kramer

################################################################################
## Joint modeling 
################################################################################

data_long <- readRDS("exampledata.rds")

library(JMbayes2)

data_surv <- data_long[!duplicated(data_long$record.id), ]

# fit lme model for biomarker concentration over time
lmefit.example <- lme(conc_log10~ day:randomized_group + day,
                      random = ~ day | record.id, 
                      data = data_long, 
                      control = lmeControl(opt = "optim"),
                      na.action = na.omit)

# fit Cox-PH model
coxfit.example <- coxph(Surv(time_mort28, death_d28) ~ randomized_group, 
                        data = data_surv)

summary(coxfit.example) # results of the Cox-PH model


jointfit.example <- jm(coxfit.example, lmefit.example, 
                       time_var = "day",
                       n_iter = 50000L,
                       n_burnin = 5000L, n_chains = 2L, n_thin = 7L, cores= 2)

summary(jointfit.example) # results of the joint model


################################################################################
### obtain indirect, direct, and total effects 
################################################################################

# Function to compute direct, indirect and total effects #######################
# Inputs a joint model object, the used survival model, and the intervention 

get_effects <- function(jointmodel, coxmodel, intervention){
  
  require(dplyr)
  
  time_var = jointmodel$model_info$var_names$time_var
  
  ### for direct 
  # get gamma point estimate
  g_hat = jointmodel$statistics$Mean$gammas
  
  # mcmc samples for gamma
  g_samples = jointmodel$mcmc$gammas  %>% 
    lapply(function(mcmc_chain) mcmc_chain[, intervention]) %>% 
    unlist()
  
  # get upp and lower for gamma
  g_ci = quantile(g_samples, probs = c(0.025, 0.975))
  
  ### for indirect
  # get beta
  b_hat = jointmodel$statistics$Mean$betas1[paste0(time_var, ":", intervention)]  
  # get alpha
  a_hat = jointmodel$statistics$Mean$alphas  
  
  # indirect effect point estimate on log hazard scale
  ind_hat = b_hat*a_hat
  
  # obtain MCMC samples of ind
  b_vecs <- lapply(jointmodel$mcmc$betas1, function(mcmc_chain) { # samples beta
    mcmc_chain[, paste0(time_var, ":", intervention)]
  })
  
  b_samples <- unlist(b_vecs) %>% as.matrix()
  
  a_samples = jointmodel$mcmc$alphas %>% unlist() %>% as.matrix() #samples alpha
  
  # use MCMC samples to obtain CIs
  ind_samples = a_samples*b_samples
  ind_ci = quantile(ind_samples, probs = c(0.025, 0.975)) # get upp and lower
  
  ### for total
  total_hat = coxmodel$coefficients[intervention]
  total_ci = confint(coxmodel)[intervention, ]
  
  # combine
  res = data.frame(effect = c("direct", 
                              "indirect", 
                              "total (Cox-PH)"),
                   est= c(g_hat, 
                          ind_hat, 
                          total_hat),
                   CI_lower = c(g_ci[1], 
                                ind_ci[1], 
                                total_ci[1]),
                   CI_upper = c(g_ci[2], 
                                ind_ci[2], 
                                total_ci[2]))
  
  return(as_tibble(res))
}

# example results for direct, indirect, and total effects
res <- get_effects(jointmodel =jointfit.example, 
                   coxmodel =coxfit.example, 
                   intervention="randomized_groupIntervention")
res

