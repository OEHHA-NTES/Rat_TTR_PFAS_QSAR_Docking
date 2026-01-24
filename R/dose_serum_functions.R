#### DOSE TO SERUM FUNCTIONS ####


## Half life from elimination constant
t_1_2fnx <- function(kelim){
  t_1_2 = log(2)/kelim
  return(t_1_2)
}


### Steady state serum concentration function from dose

CssFnx <- function(dose_mg_kg, F_frac, ClTot_L_kg_day, tau_day){
  Css_mg_L = (dose_mg_kg * F_frac)/(ClTot_L_kg_day * tau_day)
  return(Css_mg_L)
}




##### AUC infinity estimation ###
AUC_t_inf_fnx <- function(C_t, k_elim){
  AUC_t_inf = C_t / k_elim 
  
  return(AUC_t_inf)
}

AUC_inf_fnx <- function(AUC_0_t, AUC_t_inf){
  AUC_inf = AUC_0_t + AUC_t_inf 
  
  return(AUC_inf)
}

#### EXAMPLE ###
# k_elim_male_gavage = 0.18 #hr-1
# Ct_male_gavage = 3.2e-4 #umol/L @ 24hr in serum
# AUC_t_inf_male_gavage <- AUC_t_inf_fnx(Ct_male_gavage, k_elim_male_gavage)
# AUC_inf_male_gavage <- AUC_inf_fnx(AUC_IG_male, AUC_t_inf_male_gavage)
# # iv
# k_elim_male_iv = 0.20 #hr-1
# Ct_male_iv = 3.e-5 #umol/L @ 24hr in serum
# AUC_t_inf_male_iv <- AUC_t_inf_fnx(Ct_male_iv, k_elim_male_iv)
# AUC_inf_male_iv <- AUC_inf_fnx(AUC_IV_male, AUC_t_inf_male_iv)

#### Bioavailability Function from AUC
Ffnx <- function(AUC_IG, AUC_IV, Dose_IV, Dose_IG){
  F_IG = (AUC_IG / AUC_IV) * (Dose_IV / Dose_IG)
  return(F_IG)
}

#logNormal Fnx
logNormalfnx <- function(m,std){
  y = 1+std^2/m^2
  mu = log(m/sqrt(y))
  sigma = sqrt(log(y))
  
  logNormal <- list(mu,sigma)
  
  return(logNormal)
}


##### MONTE CARLO #####

## generate distribution of values based on mean and SD
# Monte Carlo function that includes bioavailable fraction
monte_carlo_F <- function(df, n_sim) {
  simulations <- map_dfr(1:n_sim, function(i) {
    F_df <- df %>% 
      distinct(sex, chem, 
               value_AUC_t_24hr_IG, sd_AUC_t_24hr_IG,
               value_AUC_t_24hr_IV, sd_AUC_t_24hr_IV) %>% 
      mutate(
        # Calculate meanlog and sdlog
        mean_AUC_IG = mu_logNormal_function(value_AUC_t_24hr_IG, sd_AUC_t_24hr_IG),
        sigma_AUC_IG = sigma_logNormal_function(value_AUC_t_24hr_IG, sd_AUC_t_24hr_IG),
        mean_AUC_IV = mu_logNormal_function(value_AUC_t_24hr_IV, sd_AUC_t_24hr_IV),
        sigma_AUC_IV = sigma_logNormal_function(value_AUC_t_24hr_IV, sd_AUC_t_24hr_IV)
      ) %>% 
      mutate(
        IG_sim = rlnorm(n(), mean_AUC_IG, sigma_AUC_IG),
        IV_sim = rlnorm(n(),  mean_AUC_IV, sigma_AUC_IV),
        F_sim = IG_sim / (IV_sim * 10)) 
    
    df %>%
      left_join(F_df, by = c("sex", "chem")) %>% 
      mutate(
        meanlog_clearance_renal = mu_logNormal_function(value_clearance_renal, sd_clearance_renal),
        sdlog_clearance_renal = sigma_logNormal_function(value_clearance_renal, sd_clearance_renal),
        #simulate renal clearance using log-normal distribution
        Clrenal_L_kg_day_sim = rlnorm(n(),meanlog_clearance_renal,sdlog_clearance_renal),
        Css_mg_L = CssFnx(dose_mg_kg = dose_mg_kg,
                          ClTot_L_kg_day = Clrenal_L_kg_day_sim,
                          F_frac = F_sim,
                          tau_day = tau_day)
      ) %>%
      #  filter(F_sim <= 1.000) %>% 
      select(chem, sex, dose_mg_kg, Css_mg_L, 
             Clrenal_L_kg_day_sim, F_sim, IV_sim, IG_sim) %>%
      mutate(simulation = i) %>% 
      filter(Css_mg_L > 0)
  })
  
  simulations
}


# Monte Carlo function that excludes bioavailable fraction
monte_carlo_noF <- function(df, n_sim) {
  simulations <- map_dfr(1:n_sim, function(i) {
    F_df <- df %>% 
      # distinct(sex, chem, 
      #          value_Intravenous, sd_Intravenous,
      #          value_Intragastric, sd_Intragastric) %>% 
      mutate(
        # IV_sim = rnorm(n(), value_Intravenous, sd_Intravenous),
        #   IG_sim = rnorm(n(), value_Intragastric, sd_Intragastric),
        F_sim = runif(n(), 1, 1)) %>%  #uniform distribution of 1s
      distinct(sex, chem, F_sim)
    
    
    df %>%
      left_join(F_df, by = c("sex", "chem")) %>% 
      mutate(
        # Calculate meanlog and sdlog
        meanlog_clearance_renal = mu_logNormal_function(value_clearance_renal, sd_clearance_renal),
        sdlog_clearance_renal = sigma_logNormal_function(value_clearance_renal, sd_clearance_renal),
        #simulate renal clearance using log-normal distribution
        Clrenal_L_kg_day_sim = rlnorm(n(),meanlog_clearance_renal,sdlog_clearance_renal),
        #estimate Css
        Css_mg_L = CssFnx(dose_mg_kg = dose_mg_kg,
                          ClTot_L_kg_day = Clrenal_L_kg_day_sim,
                          F_frac = F_sim,
                          tau_day = tau_day)
      ) %>%
      #  filter(F_sim <= 1.000) %>% 
      select(chem, sex, dose_mg_kg, Css_mg_L, 
             Clrenal_L_kg_day_sim, F_sim) %>%
      mutate(simulation = i) #%>% 
    #  filter(Css_mg_L > 0)
  })
  
  simulations
}


##### EXAMPLE ####
# Run the Monte Carlo simulation
# monte_carlo_F_results <- monte_carlo_F(Fujii_anon_wide, n_sim)
# 
# saveRDS(monte_carlo_F_results , 
#         "../../output/data/dose_to_serum/monte_carlo_F_results.rds")