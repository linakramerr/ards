### This script includes all functions used for 



## get_effects() computes direct, indirect, and total effects ################## 

get_effects <- function(jointmodel, coxmodel, intervention){
  
  require(dplyr)
  
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
  b_hat = jointmodel$statistics$Mean$betas1[paste0("day:", intervention)]  
  # get alpha
  a_hat = jointmodel$statistics$Mean$alphas  
  
  # indirect effect point estimate on log hazard scale
  ind_hat = b_hat*a_hat
  
  # obtain MCMC samples of ind
  b_vecs <- lapply(jointmodel$mcmc$betas1, function(mcmc_chain) { # samples beta
    mcmc_chain[, paste0("day:", intervention)]
  })
  
  b_samples <- unlist(b_vecs) %>% as.matrix()
  
  a_samples = jointmodel$mcmc$alphas %>% unlist() %>% as.matrix() #samples alpha
  
  # use MCMC samples to obtain CIs
  ind_samples = a_samples*b_samples
  ind_ci = quantile(ind_samples, probs = c(0.025, 0.975)) # get upp and lower
  
  ### for total
  total_hat = coxmodel$coefficients[intervention]
  total_ci = confint(coxmodel)[intervention, ]
  
  # also get the total effect estimate in joint model
  total_hat_jm = ind_hat+g_hat
  
  # use MCMC samples to obtain CIs
  total_jm_samples = ind_samples + g_samples
  total_jm_ci = quantile(total_jm_samples, probs = c(0.025, 0.975)) # get upp and lower
  
  # combine
  res = data.frame(effect = c("direct", 
                              "indirect", 
                              "total (Cox-PH)",
                              "total (JM)"),
                   est= c(g_hat, 
                          ind_hat, 
                          total_hat,
                          total_hat_jm),
                   CI_lower = c(g_ci[1], 
                                ind_ci[1], 
                                total_ci[1],
                                total_jm_ci[1]),
                   CI_upper = c(g_ci[2], 
                                ind_ci[2], 
                                total_ci[2],
                                total_jm_ci[2]))
  
  return(as_tibble(res))
}


get_int <- function(lmemodel, intervention){
  
  # get beta for interaction
  a_hat = fixef(lmemodel)[paste0("day:", intervention)]
  # get ci
  a_ci = intervals(lmemodel, which = "fixed")
  a_lower = a_ci$fixed[paste0("day:", intervention), "lower"]
  a_upper = a_ci$fixed[paste0("day:", intervention), "upper"]
  
  # combine
  res_int = data.frame(effect = "time*treat",
                   est= a_hat,
                   CI_lower = a_lower,
                   CI_upper = a_upper)
  
  return(as_tibble(res_int))
}

get_alpha <- function(jointmodel, endpoint){
  
  alpha_est = data.frame(model =  "Association",
                         endpoint = paste(endpoint),
                         est= jointmodel$statistics$Mean$alphas,
                         se = jointmodel$statistics$SD$alphas,
                         CI_lower = quantile(unlist(jointmodel$mcmc$alphas),
                                               probs = c(0.025, 0.975))[[1]],
                         CI_upper = quantile(unlist(jointmodel$mcmc$alphas), 
                                             probs = c(0.025, 0.975))[[2]])
  return(alpha_est)
}
