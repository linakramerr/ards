## functions used for analyses of 5 trial data sets



### function to fit joint models ##############################################
fit_joint_model<- function(data_long, 
                           biom, 
                           event = "death_d28",
                           time_cens = "time_mort28",
                           covariates_0 = NULL,
                           long_out = "conc_log10",
                           fixed_eff = "day + day:randomized_group",
                           intervention = "randomized_group.y",
                           id_var= "record.id",
                           time_var = "day",
                           n_iter = 50000L,
                           n_burnin = 5000L,
                           n_thin = 5L, 
                           n_chains =2L){
  
  data_long[["randomized_group"]] <- ifelse(data_long[["randomized_group"]] == levels(data_long[["randomized_group"]])[1], 0, 1)
  
  
  ds_long <- data_long %>%
    filter(biomarker == biom) %>%
    filter(!is.na(.data[[event]])) %>%
    filter(!is.na(.data[[long_out]])) %>%
    filter(if_all(all_of(covariates_0), ~ !is.na(.))) %>%
    filter(if_all(all_of(covariates_0), ~ is.finite(.))) %>%
    filter(is.finite(.data[[long_out]]))
  
  ds_surv <- ds_long %>%
    distinct(.data[[id_var]], .keep_all = TRUE) %>%
    mutate(randomized_group.y = randomized_group)

  ### Longitudinal submodel ###
  
  # this part will generate the string for the lme formula, based on included covariates
  fixed_part <- ifelse(length(covariates_0) >0, 
                       paste(fixed_eff, "+", 
                             paste(covariates_0, collapse = " + ")),
                       fixed_eff)
  
  fixed_formula <- paste(long_out, "~", fixed_part)
  
  random_formula <- paste("~",time_var, "|", id_var)
     
  # fit mixed model
  lmefit <- lme(as.formula(fixed_formula),
                as.formula(random_formula), 
                data = ds_long, 
                control = lmeControl(opt = "optim", maxIter = 200),
                na.action = na.omit)
  
  ### Survival submodel ###
  
  # this part will generate the string for the surv formula, based on included covariates
  right_part <- ifelse(length(covariates_0) >0, 
                       paste(intervention, "+", 
                             paste(covariates_0, collapse = " + ")),
                       intervention)
  
  # create formula string 
  surv_formula <- paste("Surv(", time_cens, ",", event, ") ~", right_part)
  
  # fit Cox-PH model
  coxfit <- coxph(as.formula(surv_formula), 
                  data = ds_surv)

  
  # fit joint model
  jointfit<- jm(coxfit, lmefit, 
                time_var = time_var,
                n_iter = n_iter,
                n_burnin = n_burnin, 
                n_chains = n_chains, 
                n_thin = n_thin, 
                cores = 2)
  
  study <- data_long$study
  
  return(list(study = study,
              biomarker = biom,
              coxfit = coxfit,
              jointfit = jointfit,
              ds_long = ds_long,
              baseline_cov = covariates_0))
}

## function to show missing IL-6 data


create_missing_table <- function(df, 
                                 biomarker_prefix, 
                                 days = c(0, 2, 3, 5),
                                 time_var = "time_mort28",
                                 death_var = "death_d28") {
  require(gt)
  # Create variable names
  all_vars <- paste0(biomarker_prefix, "_", days)
  
  # Create the data with missingness indicators
  summary_data <- df %>%
    mutate(
      all_missing = if_all(all_of(all_vars), ~ is.na(.), 1, 0), 
      include_in_analysis = (!all_missing & !is.na(time_mort28))
    )

  # Create alive indicators for each timepoint
  for (d in days) {
    alive_var <- paste0("alive_day", d)
    summary_data <- summary_data %>%
      mutate(!!alive_var := is.na(.data[[time_var]]) | .data[[time_var]] >= d)
  }
  
  # Build the first row (missing all)
  first_row <- tibble(
    Measure = paste0("Missing all ", biomarker_prefix, " measurements"),
    N = nrow(summary_data),
    Available = sum(!summary_data$all_missing),
    Missing = sum(summary_data$all_missing),
    `% Missing` = round(sum(summary_data$all_missing) / nrow(summary_data) * 100, 1),
    Deaths = 0
  )
  
  # Build rows for each timepoint
  timepoint_rows <- map_dfr(days, function(d) {
    var_name <- paste0(biomarker_prefix, "_", d)
    alive_var <- paste0("alive_day", d)
    
    # N at this timepoint (alive AND have at least one measurement)
    n_at_risk <- sum(summary_data[[alive_var]] & summary_data$include_in_analysis, na.rm = TRUE)
    
    # Available measurements
    available <- sum(!is.na(summary_data[[var_name]]) & 
                       summary_data[[alive_var]] & 
                       summary_data$include_in_analysis, na.rm = TRUE)
    
    # Missing measurements
    missing <- sum(is.na(summary_data[[var_name]]) & 
                     summary_data[[alive_var]] & 
                     summary_data$include_in_analysis, na.rm = TRUE)
    
    # Deaths by this timepoint (among those with at least one measurement)
    deaths <- sum(summary_data$include_in_analysis) - n_at_risk
    
    tibble(
      Measure = paste0("Day ", d, if_else(d == min(days), " (Baseline)", "")),
      N = n_at_risk,
      Available = available,
      Missing = missing,
      `% Missing` = round(missing / n_at_risk * 100, 1),
      Deaths = deaths
    )
  })
  
  # Combine rows
  table_data <- bind_rows(first_row, timepoint_rows)
  
  # Create gt table
  table_data %>%
    gt() %>%
    tab_header(
      title = paste0(biomarker_prefix, " missing data summary"), 
      subtitle = paste0(sum(summary_data$all_missing), " patients with no measurements and ", sum(is.na(summary_data$death_d28)),  " with no mortality data excluded from day 0 onwards")
    ) %>%
    cols_label(
      Measure = "Time Point",
      N = "N at Risk",
      Deaths = "Cumulative Deaths"
    ) %>%
    tab_style(
      style = cell_fill(color = "lightgray"),
      locations = cells_body(rows = 1)
    ) %>%
    tab_footnote(
      footnote = paste0("N excludes patients missing all ", biomarker_prefix, 
                        " measurements and decreases over time as patients die"),
      locations = cells_column_labels(columns = N)
    ) %>%
    fmt_number(
      columns = c(N, Available, Missing, Deaths),
      decimals = 0
    )
}


## get_effects() computes direct, indirect, and total effects ################## 

get_effects <- function(jointmodel, coxmodel, intervention.y, intervention){
  
  require(dplyr)
  
  ### for direct 
  # get gamma point estimate
  g_hat = jointmodel$statistics$Mean$gammas["randomized_group.y"]
  
  # mcmc samples for gamma
  g_samples = jointmodel$mcmc$gammas %>% 
    lapply(function(mcmc_chain) mcmc_chain[, intervention.y]) %>% 
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
  total_hat = coxmodel$coefficients[intervention.y]
  total_ci = confint(coxmodel)[intervention.y, ]
  
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



library(dplyr)
library(JMbayes2)

pool_joint_models <- function(jointfit_list) {
  # jointfit_list: list of fitted JMbayes2 joint models (from multiple imputations)
  m <- length(jointfit_list)
  # function to extract coefficients and SEs
  extract_coefs <- function(fit, type = c("Longitudinal", "Survival", "Association")) {
    type <- match.arg(type)
  if(type == "Longitudinal") {
    beta <- fit[["statistics"]][["Mean"]][["betas1"]]
    se <- fit[["statistics"]][["SE"]][["betas1"]]
  }else if(type == "Association") {
    beta <- fit[["statistics"]][["Mean"]][["alphas"]]
    se  <- fit[["statistics"]][["SE"]][["alphas"]]
  }else if(type == "Survival"){
    beta <- fit[["statistics"]][["Mean"]][["gammas"]]
    se <- fit[["statistics"]][["SE"]][["gammas"]]
  }
  list(beta = beta, se = se)
}

# Pooling function for a list of coefficient vectors
rubins_rules <- function(coef_list) {
  betas <- do.call(rbind, lapply(coef_list, `[[`, "beta"))
  ses <- do.call(rbind, lapply(coef_list, `[[`, "se"))
  p <- nrow(betas)
  
  # 1) Pooled estimate
  beta_bar <- colMeans(betas)
  
  # 2) Within-imputation variance
  U_bar <- diag(colMeans(ses^2))
  
  # 3) Between-imputation variance
  B <- cov(betas)
  
  # 4) Total variance
  T_mat <- U_bar + (1 + 1/m) * B
  
  # 5) Pooled SEs
  pooled_se <- sqrt(diag(T_mat))
  
  # 6) 95% CI
  ci_lower <- beta_bar - 1.96 * pooled_se
  ci_upper <- beta_bar + 1.96 * pooled_se
  
  # 7) Z-values and p-values
  z <- beta_bar / pooled_se
  pvals <- 2 * (1 - pnorm(abs(z)))
  
  data.frame(
    Estimate = beta_bar,
    SE = pooled_se,
    CI_lower = ci_lower,
    CI_upper = ci_upper,
    z = z,
    p = pvals,
    row.names = names(beta_bar)
  )
}
  
  # Extract coefficients for longitudinal submodel
  long_list <- lapply(jointfit_list, extract_coefs, type = "Longitudinal")
  long_pooled <- rubins_rules(long_list)
  
  # Extract coefficients for survival submodel
  surv_list <- lapply(jointfit_list, extract_coefs, type = "Survival")
  surv_pooled <- rubins_rules(surv_list)
  
  # coefficients for association estimate
  asso_list <- lapply(jointfit_list, extract_coefs, type = "Association")
  asso_pooled <- rubins_rules(asso_list)
  
  list(
    Longitudinal = long_pooled,
    Survival = surv_pooled,
    Association = asso_pooled
  )
}



add_name_variables <- function(df_list) {
  
  parsed <- lapply(names(df_list), function(nm) {
    parts <- strsplit(nm, "_")[[1]]
    
    study <- parts[1]
    
    if (length(parts) == 2) {
      class <- "all"
      biomarker <- parts[2]
    } else if (length(parts) == 3) {
      class <- parts[2]
      biomarker <- parts[3]
    } else {
      stop(paste("Unexpected name format:", nm))
    }
    
    data.frame(name = nm, study = study, class = class, biomarker = biomarker,
               stringsAsFactors = FALSE)
  })
  
  bind_rows(parsed)
}

get_alpha <- function(jointmodel){#, endpoint){
  
  alpha_est = data.frame(model =  "Association",
                      #   endpoint = paste(endpoint),
                         est= jointmodel$statistics$Mean$alphas,
                         se = jointmodel$statistics$SD$alphas,
                         CI_lower = quantile(unlist(jointmodel$mcmc$alphas),
                                             probs = c(0.025, 0.975))[[1]],
                         CI_upper = quantile(unlist(jointmodel$mcmc$alphas), 
                                             probs = c(0.025, 0.975))[[2]])
  return(alpha_est)
}

