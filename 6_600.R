library(tidyverse)
library(lme4)
library(survival)
library(nlme)
library(openxlsx)


#Initial N, sim_repetition, sim_year
sample_N = 600
repetition = 20
sim_year=5
measure_time = seq(0.1, sim_year, by=0.1)


#Program Output file
dementia_group_coef = NULL
dementia_group_pvalue = NULL
gcp_group_coef = NULL
gcp_group_pvalue = NULL
simN_each = NULL
result_summary = NULL


#Load base NHANES Data
basenhanes_data = read.csv("selected_nhanespop.csv")
basenhanes_data = basenhanes_data %>% rename(pred_gcp = predict_1_fix) %>%
  mutate(gender=tolower(gender))

#Load Model parameters: stroke, gcp, dementia, death
## 2nd stroke model parameters, use stroke_lambda_avg, stroke_nu_avg 
stk_fit_para = read.xlsx("ModelParameters_2024Nov.xlsx",sheet=1, rowNames = TRUE)
stk_weibull = read.xlsx("ModelParameters_2024Nov.xlsx", sheet=2, rows = 1:3,rowNames = TRUE)

stroke_sbpcenter_beta = stk_fit_para["sbp_center", "coef"]
stroke_fgcenter_beta = stk_fit_para["fg_center", "coef"]
stroke_ldlcenter_beta = stk_fit_para["ldl_center", "coef"]

stroke_lambda = stk_weibull["strokeM_lambda", "Baseline_in_the_program"]
stroke_rou = stk_weibull["strokeM_rou", "Baseline_in_the_program"]

## GCP model
gcp_mod_coef = read.xlsx("ModelParameters_2024Nov.xlsx", sheet=3, rowNames = TRUE)
gcp_mod_random = read.xlsx("ModelParameters_2024Nov.xlsx", sheet=4, rows = 1:4, rowNames = TRUE)

gcpM_intercept = gcp_mod_coef["(Intercept)", "Value"]
gcpM_time_beta = gcp_mod_coef["time", "Value"]
gcpM_2ndstroke_beta = gcp_mod_coef["stroke_2", "Value"]
gcpM_2ndstroketime_beta = gcp_mod_coef["stroke_2_time", "Value"]
gcpM_agecenter_beta = gcp_mod_coef["age_Center", "Value"]
gcpM_edu2_beta = gcp_mod_coef["educ_fpsc2", "Value"]
gcpM_edu3_beta = gcp_mod_coef["educ_fpsc3", "Value"]
gcpM_edu4_beta = gcp_mod_coef["educ_fpsc4", "Value"]
gcpM_female_beta = gcp_mod_coef["gender2female", "Value"]
gcpM_sbpcenter_beta = gcp_mod_coef["sbp_center", "Value"]
gcpM_fgcenter_beta = gcp_mod_coef["fg_center", "Value"]
gcpM_ldlcenter_beta = gcp_mod_coef["ldl_center", "Value"]
gcpM_CHS_beta = gcp_mod_coef["studynamechs", "Value"]
gcpM_FOS_beta = gcp_mod_coef["studynamefos", "Value"]
gcpM_REGARDS_beta = gcp_mod_coef["studynameregards", "Value"]
gcpM_timeagec_beta = gcp_mod_coef["time:age_Center", "Value"]
gcpM_timeldlc_beta = gcp_mod_coef["time:ldl_center", "Value"]
gcpM_timesbpc_beta = gcp_mod_coef["time:sbp_center", "Value"]
gcpM_timefgc_beta = gcp_mod_coef["time:fg_center", "Value"]

gcpM_studyname_beta = (gcpM_CHS_beta+gcpM_FOS_beta+gcpM_REGARDS_beta)/4

R_std = gcp_mod_random["residual_forprogram", "StdDev"]
std_inter = gcp_mod_random["(Intercept)", "StdDev"]
std_time = gcp_mod_random["time", "StdDev"]
corr = as.numeric(gcp_mod_random["time", "Corr"])

## Dementia model
dementia_fit_para = read.xlsx("ModelParameters_2024Nov.xlsx", sheet=5, rowNames = TRUE)
dementia_weibull = read.xlsx("ModelParameters_2024Nov.xlsx", sheet=6, rows = 1:3, rowNames = TRUE)

dementia_agecenter_beta = dementia_fit_para["age_Center", "coef"]
dementia_female_beta = dementia_fit_para["gender2female", "coef"]
dementia_gcppred_beta = dementia_fit_para["gcp_pred", "coef"]
dementia_istroke_beta = dementia_fit_para["recurrent_stk", "coef"]

dementia_lambda = dementia_weibull["dementiaM_lambda", "Baseline_in_the_program"]
dementia_rou = dementia_weibull["dementiaM_rou", "Baseline_in_the_program"]

## Death model
death_fit_para = read.xlsx("ModelParameters_2024Nov.xlsx", sheet=7, rowNames = TRUE)
death_weibull = read.xlsx("ModelParameters_2024Nov.xlsx", sheet=8, rows = 1:3, rowNames = TRUE)

death_fgcenter_beta = death_fit_para["fg_center", "coef"]
death_female_beta = death_fit_para["gender2female", "coef"]
death_agecenter_beta = death_fit_para["age_Center", "coef"]
death_istroke_beta = death_fit_para["event_re_stk", "coef"]
death_idementia_beta = death_fit_para["event_dem", "coef"]
death_gcppred_beta = death_fit_para["gcp_pred", "coef"]

death_lambda = death_weibull["deathM_lambda", "Baseline_in_the_program"]
death_rou = death_weibull["deathM_rou", "Baseline_in_the_program"]


#Load Simulation functions
get_sim = function(data, n, error_matrix){
  miu_df = data.frame(mean_loga1c = mean(data$log_a1c),
                      #mean_logfg = mean(data$log_fasting_glucose),
                      mean_logsbp = mean(data$log_sbp),
                      mean_logldl = mean(data$log_ldl),
                      mean_stkage = mean(data$stroke_age))
  #print(miu_df)
  miu_df =  bind_rows(replicate(n, miu_df, simplify = FALSE))
  miu_df = as.matrix(miu_df, nrow=n, ncol=4)
  
  error = error_matrix %*% matrix(rnorm(4*n), nrow=4, ncol=n)
  error = t(error)
  #print(error)
  
  sim = error+miu_df
  sim = as.data.frame(sim) %>% mutate(a1c = exp(log_a1c),
                                      #fasting_glucose = exp(log_fasting_glucose),
                                      sbp = exp(log_sbp),
                                      ldlcholesterol = exp(log_ldl))
  
  return(sim)
}

sim_nhanes = function(nhanes_data, sample_size, treatment){
  data = nhanes_data
  data = data %>% mutate(log_ldl = log(ldlcholesterol),
                         log_sbp = log(sbp),
                         log_fasting_glucose = log(fasting_glucose),
                         log_a1c = log(a1c))
  
  edu1_female = data %>% filter(gender=="female", educ_fpsc==1)
  edu2_female = data %>% filter(gender=="female", educ_fpsc==2)
  edu3_female = data %>% filter(gender=="female", educ_fpsc==3)
  edu4_female = data %>% filter(gender=="female", educ_fpsc==4)
  edu1_male = data %>% filter(gender=="male", educ_fpsc==1)
  edu2_male = data %>% filter(gender=="male", educ_fpsc==2)
  edu3_male = data %>% filter(gender=="male", educ_fpsc==3)
  edu4_male = data %>% filter(gender=="male", educ_fpsc==4)
  
  min_a1c = min(data$a1c)
  max_a1c = max(data$a1c)
  min_sbp = min(data$sbp)
  max_sbp = max(data$sbp)
  min_ldl = min(data$ldlcholesterol)
  max_ldl = max(data$ldlcholesterol)
  min_fg = min(data$fasting_glucose)
  max_fg = max(data$fasting_glucose)
  min_stkage = min(data$stroke_age)
  max_stkage = max(data$stroke_age)
  
  
  #MVN model, use t_chol_decompo for later simulation
  mvn_mod = lm(cbind(log_a1c, log_sbp, log_ldl, stroke_age)~educ_fpsc+gender, data = data)
  residual = residuals(mvn_mod)
  covariance_matrix = as.matrix(var(residual))
  #print(covariance_matrix)
  
  chol_decompo = chol(covariance_matrix)
  nvars = dim(chol_decompo)[1]
  
  t_chol_decompo = t(chol_decompo)
  t_chol_decompo[1,1] = 0.03
  #print(t_chol_decompo)
  
  #get n by gender*education category
  n_bygroup = data %>% group_by(gender,educ_fpsc) %>% summarise(n =n_distinct(id)) %>%
    ungroup() %>%
    mutate(total = sum(n), prop = n/total, sim_n = round(sample_size*prop))
  n_female1 = n_bygroup$sim_n[1]
  n_female2 = n_bygroup$sim_n[2]
  n_female3 = n_bygroup$sim_n[3]
  n_female4 = n_bygroup$sim_n[4]
  n_male1 = n_bygroup$sim_n[5]
  n_male2 = n_bygroup$sim_n[6]
  n_male3 = n_bygroup$sim_n[7]
  n_male4 = n_bygroup$sim_n[8]
  
  #simulate by gender*education, use n above and get_sim() function
  sim_edu1_female = get_sim(edu1_female, n_female1, t_chol_decompo) %>% mutate(edu_fpsc=1, gender="female")
  sim_edu1_male = get_sim(edu1_male,n_male1, t_chol_decompo) %>% mutate(edu_fpsc=1, gender="male")
  sim_edu2_female = get_sim(edu2_female, n_female2, t_chol_decompo) %>% mutate(edu_fpsc=2, gender="female")
  sim_edu2_male = get_sim(edu2_male, n_male2, t_chol_decompo) %>% mutate(edu_fpsc=2, gender="male")
  sim_edu3_female = get_sim(edu3_female, n_female3, t_chol_decompo) %>% mutate(edu_fpsc=3, gender="female")
  sim_edu3_male = get_sim(edu3_male, n_male3, t_chol_decompo) %>% mutate(edu_fpsc=3, gender="male")
  sim_edu4_female = get_sim(edu1_female, n_female4, t_chol_decompo) %>% mutate(edu_fpsc=4, gender="female")
  sim_edu4_male = get_sim(edu1_male,n_male4, t_chol_decompo) %>% mutate(edu_fpsc=4, gender="male")
  
  sim_data = rbind(sim_edu1_female, sim_edu1_male)
  sim_data = rbind(sim_data, sim_edu2_female)
  sim_data = rbind(sim_data, sim_edu2_male)
  sim_data = rbind(sim_data, sim_edu3_female)
  sim_data = rbind(sim_data, sim_edu3_male)
  sim_data = rbind(sim_data, sim_edu4_female)
  sim_data = rbind(sim_data, sim_edu4_male)
  
  deduct = ifelse(treatment=="control", mean(sim_data$a1c)-7.5, ifelse(treatment=="intensive", mean(sim_data$a1c)-6.25, 0))
  deduct_min = ifelse(treatment=="control", 7, ifelse(treatment=="intensive", 5.5, 0))
  deduct_max = ifelse(treatment=="control", 8, ifelse(treatment=="intensive", 7, 0))
  
  sim_data = sim_data %>% mutate(sbp = ifelse(sbp<=min_sbp, min_sbp, ifelse(sbp>=max_sbp, max_sbp, sbp)),
                                 ldlcholesterol = ifelse(ldlcholesterol<=min_ldl, min_ldl, ifelse(ldlcholesterol>=max_ldl, max_ldl, ldlcholesterol)),
                                 #fasting_glucose = ifelse(fasting_glucose<=min_fg, min_fg, ifelse(fasting_glucose>=max_fg, max_fg, fasting_glucose)),
                                 a1c = ifelse(a1c<=min_a1c, min_a1c, ifelse(a1c>=max_a1c, max_a1c, a1c)),
                                 a1c_cut = a1c-deduct,
                                 a1c_cut = ifelse(a1c_cut<=deduct_min, deduct_min, ifelse(a1c_cut>=deduct_max, deduct_max, a1c_cut)),
                                 fasting_glucose_cut = a1c_cut*28.7-46.7,
                                 fasting_glucose = a1c*28.7-46.7,
                                 stroke_age = ifelse(stroke_age<=min_stkage, min_stkage, ifelse(stroke_age>=max_stkage, max_stkage, stroke_age)),
                                 group = treatment,
                                 id = paste0(row_number(), "_",group))  
  
  return(sim_data)
}

get_long = function(base_data, treatment){
  seq = data.frame(time = seq(0, 5, by=0.1)) 
  id_df = base_data %>% select(id) %>% slice(rep(1:n(), each=(5*10+1)))
  id_rep = cbind(id_df, seq) 
  
  a1clong = base_data %>% 
    select(id, a1c, a1c_cut) %>% 
    #mutate(a1c_cut = ifelse(a1c_cut<6.5, 6.5, ifelse(a1c_cut>7.5, 7.5, a1c_cut))) %>%
    pivot_longer(a1c:a1c_cut, names_to = "a1c_type", values_to = "a1c") %>%
    mutate(time = ifelse(a1c_type=="a1c", 0, 0.1)) %>%
    select(id, time, a1c) %>%
    arrange(id, time) %>%
    full_join(id_rep, c("id", "time")) %>%
    arrange(id, time) %>%
    
    group_by(id) %>%
    fill(c("a1c"), .direction = "down") %>%
    ungroup() %>%
    
    group_by(id) %>%
    mutate(cumean_a1c = cummean(a1c)) %>%
    ungroup() %>%
    
    mutate(fasting_glucose=a1c*28.7-46.7) %>%
    group_by(id) %>%
    mutate(cumean_fg = cummean(fasting_glucose))
  
  sbpldllong = base_data %>% 
    select(id, sbp,ldlcholesterol) %>%
    mutate(time=0) %>%
    full_join(id_rep, c("id", "time")) %>%
    arrange(id, time) %>%
    group_by(id) %>%
    fill(c("sbp", "ldlcholesterol"), .direction = "down") %>%
    ungroup() %>%
    group_by(id) %>%
    mutate(cumean_sbp = cummean(sbp),
           cumean_ldl = cummean(ldlcholesterol)) %>%
    ungroup() %>%
    select(id, time, ldlcholesterol, sbp, cumean_ldl, cumean_sbp) 
  
  
  rf_long = a1clong %>% full_join(sbpldllong, c("id", "time")) %>% mutate(group = treatment)
  return(rf_long)
}


#Simulation Begin
for(m in 1:repetition){
  print(m)
  sim_data = NULL
  
  #Simulate data based on sample_N and basenhanes_data
  sim_base_control = sim_nhanes(basenhanes_data, sample_N, "control")
  sim_base_intensive = sim_nhanes(basenhanes_data, sample_N, "intensive")
  longrf_control = get_long(sim_base_control, 'control')
  longrf_intensive = get_long(sim_base_intensive, "intensive")
  
  simN_each_m = n_distinct(sim_base_control$id)
  simN_each = c(simN_each, simN_each_m)
  
  sim_base_control=sim_base_control %>% mutate(dropout_time = runif(n(), min=0, max=1)*50,
                                               dropout_time = round(dropout_time, 1)) %>%
    mutate(i_dropout_time = ifelse(dropout_time>5, 0, 1))
  
  mean(sim_base_control$i_dropout_time)
  ## attrition rate does not depend on the event rate? yes or no ? 
  ## the first file, change the attrition rate to 1%
  ## write code to store the actual drop out rate
  
  sim_base_intensive = sim_base_intensive %>% mutate(dropout_time = runif(n(), min=0, max=1)*50,
                                                     dropout_time = round(dropout_time, 1)) %>%
    mutate(i_dropout_time = ifelse(dropout_time>5, 0, 1))
    
  nrow_control = nrow(sim_base_control)
  nrow_intensive = nrow(sim_base_intensive)
  random_mean = c(0, 0)
  random_cv = matrix(c(std_inter**2, std_inter*std_time*corr, std_inter*std_time*corr, std_time**2), ncol=2)
  gcp_random_control = MASS::mvrnorm(nrow_control, random_mean, random_cv)
  gcp_random_intensive = MASS::mvrnorm(nrow_intensive, random_mean, random_cv)
  ### the correlation of random intercept and random slope may influence the power
 ## ( when corrlation is big, the sd for the beta coefficient is small?)
  ## use this MASS package to consider the corr MASS::mvrnorm(nrow_control, random_mean, random_cv)
  
  
  sim_base_control = sim_base_control %>% mutate(random_intercept = gcp_random_control[,1],
                                                 random_slope = gcp_random_control[,2])
  sim_base_intensive = sim_base_intensive %>% mutate(random_intercept = gcp_random_intensive[,1],
                                                     random_slope = gcp_random_intensive[,2])
  
  baseline_data= rbind(sim_base_control, sim_base_intensive) %>% 
    rename(baseline_fg = fasting_glucose,
           baseline_a1c = a1c,
           baseline_ldl = ldlcholesterol,
           baseline_sbp = sbp,
           educ_fpsc = edu_fpsc) %>%
    mutate(enroll_time=0.5, 
           time=0)
  
  long_data = rbind(longrf_control, longrf_intensive) %>% select(-group, -a1c, -fasting_glucose, -ldlcholesterol, -sbp)
  longdata_time = unique(long_data$time)
  
  #Prepare simulation
  sim_next = baseline_data %>% 
    mutate(stroke_entered = FALSE, 
           stroke_time=0, 
           death=FALSE, 
           death_status=FALSE, 
           dementia_entered = FALSE, 
           dementia_time=0, 
           death_time=0) %>%
    mutate(gcp_fix = gcpM_intercept+gcpM_time_beta*(0.5)+
             gcpM_2ndstroke_beta*0+gcpM_2ndstroketime_beta*(0)+
             gcpM_agecenter_beta*(stroke_age-74)+
             gcpM_edu2_beta*(educ_fpsc==2)+gcpM_edu3_beta*(educ_fpsc==3)+gcpM_edu4_beta*(educ_fpsc==4) + 
             gcpM_female_beta*(gender=="female")+
             gcpM_sbpcenter_beta*(baseline_sbp-132) + gcpM_fgcenter_beta*(baseline_fg-105)+gcpM_ldlcenter_beta*(baseline_ldl-103) + 
             gcpM_studyname_beta+
             gcpM_timeagec_beta*(0.5)*(stroke_age-74)+gcpM_timeldlc_beta*(0.5)*(baseline_ldl-103)+
             gcpM_timesbpc_beta*(0.5)*(baseline_sbp-132)+
             gcpM_timefgc_beta*(0.5)*(baseline_fg-105),
           
           
           random_error = rnorm(n(), mean=0 , sd = R_std),
           gcp_noerror = gcp_fix+random_intercept+0.5*random_slope,
           gcp_werror = gcp_noerror+random_error,
           
           gcp_noerror = ifelse(gcp_noerror<4.041921, 4.041921, ifelse(gcp_noerror>73.03935, 73.03935, gcp_noerror)),
           gcp_werror = ifelse(gcp_werror<4.041921, 4.041921, ifelse(gcp_werror>73.03935, 73.03935, gcp_werror)),
           
           last_gcp_noerror = gcp_noerror,
           last_gcp_werror = gcp_werror,
           repetition=m)
  temp = sim_next %>% select(id, time, enroll_time, repetition, group, stroke_age, educ_fpsc,
                             stroke_entered, stroke_time, gcp_fix,
                             gcp_noerror, gcp_werror, last_gcp_werror, last_gcp_noerror,
                             dementia_entered, dementia_time, death, death_status, death_time, dropout_time, i_dropout_time)
  
  sim_data = rbind(sim_data, temp)
  k=0
  
  for(i in measure_time){
    k = k+1
    
    #get time i risk factors
    long_i = long_data %>% filter(time==longdata_time[k]) %>% rename(rf_time = time) 
    sim_next = sim_next %>% left_join(long_i, "id") %>% mutate(time = rf_time+0.1)
    print(paste("time:", unique(sim_next$time)))
    print(paste0("rf_time:", unique(sim_next$rf_time)))
    
    
    #2ND Stroke Module
    sim_next = sim_next %>% mutate(xbeta_stroke= stroke_sbpcenter_beta*(cumean_sbp-132)+stroke_fgcenter_beta*(cumean_fg-105)+stroke_ldlcenter_beta*(cumean_ldl-103),
                                   prob_stroke = 1-exp(((rf_time+enroll_time)^stroke_rou)*exp(xbeta_stroke)*stroke_lambda-((rf_time+enroll_time+0.1)^stroke_rou)*exp(xbeta_stroke)*stroke_lambda),
                                   stroke = rbinom(n(), 1, prob_stroke)&(!death_status))
    
    #GCP Module
    sim_next = sim_next %>% mutate(last_gcp_noerror = gcp_noerror,
                                   last_gcp_werror = gcp_werror,
                                   random_error = rnorm(n(), 0 , sd = R_std)) 
    sim_next = sim_next %>% mutate(gcp_fix = gcpM_intercept+gcpM_time_beta*(time+enroll_time)+
                                     gcpM_2ndstroke_beta*stroke_entered+
                                     stroke_entered*gcpM_2ndstroketime_beta*(time-stroke_time)+
                                     gcpM_agecenter_beta*(stroke_age-74)+
                                     gcpM_edu2_beta*(educ_fpsc==2)+gcpM_edu3_beta*(educ_fpsc==3)+gcpM_edu4_beta*(educ_fpsc==4) + 
                                     gcpM_female_beta*(gender=="female")+
                                     gcpM_sbpcenter_beta*(baseline_sbp-132) + gcpM_fgcenter_beta*(baseline_fg-105)+gcpM_ldlcenter_beta*(baseline_ldl-103) + 
                                     gcpM_studyname_beta+
                                     gcpM_timeagec_beta*(time+enroll_time)*(stroke_age-74)+gcpM_timeldlc_beta*(time+enroll_time)*(cumean_ldl-103)+
                                     gcpM_timesbpc_beta*(time+enroll_time)*(cumean_sbp-132)+
                                     gcpM_timefgc_beta*(time+enroll_time)*(cumean_fg-105))
    
    sim_next = sim_next %>% mutate(gcp_noerror = gcp_fix+random_intercept+(time+enroll_time)*random_slope,
                                   gcp_werror = gcp_noerror+random_error,
                                   
                                   gcp_noerror = ifelse(gcp_noerror<4.041921, 4.041921, ifelse(gcp_noerror>73.03935, 73.03935, gcp_noerror)),
                                   gcp_werror = ifelse(gcp_werror<4.041921, 4.041921, ifelse(gcp_werror>73.03935, 73.03935, gcp_werror)))
    
    #Dementia module
    sim_next = sim_next %>% mutate(xbeta_dementia = dementia_agecenter_beta*(stroke_age-74)+dementia_female_beta*(gender=="female")+
                                     dementia_gcppred_beta*last_gcp_noerror+dementia_istroke_beta*stroke_entered,
                                   prob_dementia = 1-exp(((rf_time+enroll_time)^dementia_rou)*exp(xbeta_dementia)*dementia_lambda-((rf_time+enroll_time+0.1)^dementia_rou)*exp(xbeta_dementia)*dementia_lambda),
                                   dementia = rbinom(n(),1, prob_dementia)&(!death_status))
    
    
    #Death module
    sim_next = sim_next %>% mutate(xbeta_death = death_fgcenter_beta*(cumean_fg-105)+death_female_beta*(gender=="female")+death_agecenter_beta*(stroke_age-74)+
                                     death_istroke_beta*stroke_entered+death_idementia_beta*dementia_entered+death_gcppred_beta*last_gcp_noerror,
                                   prob_death = 1-exp(((rf_time+enroll_time)^death_rou)*exp(xbeta_death)*death_lambda-((rf_time+enroll_time+0.1)^death_rou)*exp(xbeta_death)*death_lambda),
                                   death = rbinom(n(), 1, prob_death)&(!death_status),
                                   death_status = ifelse(death_status==TRUE, death_status, ifelse(death==TRUE, death, FALSE)))
    
    #Update final status
    sim_next = sim_next %>% mutate(death_time = ifelse(death_time==0, ifelse(death==1, time, 0), death_time),
                                   stroke_entered = ifelse(stroke_entered==TRUE, stroke_entered, ifelse(stroke==TRUE, stroke, FALSE)),
                                   stroke_time = ifelse(stroke_time==0, ifelse(stroke==1, time, 0), stroke_time),
                                   dementia_entered = ifelse(dementia_entered==TRUE, dementia_entered, ifelse(dementia==TRUE, dementia, FALSE)),
                                   dementia_time = ifelse(dementia_time==0, ifelse(dementia==1, time, 0), dementia_time))
    
    
    temp = sim_next %>% select(id, time, enroll_time, repetition, group, stroke_age, educ_fpsc,
                               stroke_entered, stroke_time, gcp_fix,
                               gcp_noerror, gcp_werror, last_gcp_werror, last_gcp_noerror,
                               dementia_entered, dementia_time, death, death_status, death_time, dropout_time, i_dropout_time)
    
    sim_next = sim_next %>% select(-gcp_fix, -cumean_fg, -cumean_sbp, -cumean_ldl, -rf_time, -cumean_a1c)
    sim_data = rbind(sim_data, temp)
  }
  
  #GCP analysis model
  gcpmodel_data = sim_data %>% filter(!death_status) %>%
    filter(time<=dropout_time) %>%
    mutate(analysis_time = enroll_time+time) %>%
    mutate(time = as.character(time)) %>%
    filter(time %in% c("0", "0.5", "1", "1.5", "2", "2.5", "3", "3.5", "4", "4.5", "5")) %>%
    mutate(stroke_age_gcpM = stroke_age+0.5,
           stroke_age_gcpM_c = stroke_age_gcpM-74)
  gcp_bygroup = lme(gcp_werror~analysis_time+ group+analysis_time:group + stroke_age_gcpM +as.factor(educ_fpsc),
                    
                    random=~1+analysis_time|id,
                    control = lmeControl(opt='optim'),
                    data = gcpmodel_data)
  gcp_bygroup_beta = as.numeric(summary(gcp_bygroup)$coefficients$fixed[8])
  vcov_gcp_bygroup = vcov(gcp_bygroup)
  gcp_bygroup_se = sqrt(vcov_gcp_bygroup[4,4])
  gcp_pvalue_m = summary(gcp_bygroup)$tTable[8, 5]
  gcp_group_coef = c(gcp_group_coef, gcp_bygroup_beta)
  gcp_group_pvalue = c(gcp_group_pvalue, gcp_pvalue_m)
  
  #dementia model
  dementia_data = sim_data %>%
    group_by(id) %>%
    mutate(i = row_number(), n = n()) %>%
    filter(i==n) %>%
    mutate(death_time_a = ifelse(death_status, death_time+enroll_time, ifelse(i_dropout_time==1, dropout_time+enroll_time, 5+enroll_time)),
           dementia_time_a = ifelse(dementia_entered==1, dementia_time+enroll_time, death_time_a)) %>%
    mutate(dementia_time_a1 = ifelse(dementia_time_a<=0.5, 0.5,
                                     ifelse(dementia_time_a<=1, 1,
                                            ifelse(dementia_time_a<=1.5, 1.5,
                                                   ifelse(dementia_time_a<=2, 2,
                                                          ifelse(dementia_time_a<=2.5, 2.5,
                                                                 ifelse(dementia_time_a<=3, 3,
                                                                        ifelse(dementia_time_a<=3.5, 3.5, ifelse(dementia_time_a<=4, 4,
                                                                                                                 ifelse(dementia_time_a<=4.5, 4.5,
                                                                                                                        ifelse(dementia_time_a<=5, 5, 5.5)))))))))),
           death_time_a1 = ifelse(death_time_a<=0.5, 0.5,
                                  ifelse(death_time_a<=1, 1,
                                         ifelse(death_time_a<=1.5, 1.5,
                                                ifelse(death_time_a<=2, 2,
                                                       ifelse(death_time_a<=2.5, 2.5,
                                                              ifelse(death_time_a<=3, 3,
                                                                     ifelse(death_time_a<=3.5, 3.5, ifelse(death_time_a<=4, 4,
                                                                                                           ifelse(death_time_a<=4.5, 4.5,
                                                                                                                  ifelse(death_time_a<=5, 5, 5.5))))))))))) %>%
    select(id, dementia_entered, death_status, dementia_time_a,death_time_a, dementia_time_a1, death_time_a1, group, stroke_age, educ_fpsc) %>%
    mutate(i_dementia = as.numeric(dementia_entered),
           i_dementia = ifelse(i_dementia==1&death_status==TRUE, ifelse(death_time_a1==dementia_time_a1, 0, i_dementia), i_dementia),
           t_dementia = ifelse(i_dementia==1, dementia_time_a1, death_time_a1),
           stroke_ageM = stroke_age+0.5) 
  
  #adjust for age and education level
  dementia_bygroup = coxph(Surv(t_dementia, i_dementia)~as.factor(group)+stroke_ageM+as.factor(educ_fpsc), data = dementia_data)
  
  dementia_bygroup_beta = as.numeric(summary(dementia_bygroup)$coefficients[1])
  dementia_bygroup_se = as.numeric(summary(dementia_bygroup)$coefficients[3])
  dementia_pvalue_m = summary(dementia_bygroup)$coefficients[1,5]
  dementia_group_coef = c(dementia_group_coef, dementia_bygroup_beta)
  dementia_group_pvalue = c(dementia_group_pvalue, dementia_pvalue_m)
  
  #Incidence rate
  re_stroke = sim_data %>% 
    filter(time<=dropout_time) %>%
    group_by(group, time, repetition) %>%
    summarise(rep_mean = mean(stroke_entered)) %>%
    ungroup() %>%
    group_by(group, time) %>%
    summarise(stroke_rate = mean(rep_mean))%>%
    ungroup()
  re_dementia = sim_data %>% 
    filter(time<=dropout_time) %>%
    group_by(group, time, repetition) %>%
    summarise(rep_mean = mean(dementia_entered)) %>%
    ungroup() %>%
    group_by(group, time) %>%
    summarise(dementia_rate = mean(rep_mean)) %>%
    ungroup()
  re_death = sim_data %>% 
    filter(time<=dropout_time) %>%
    group_by(group, time, repetition) %>%
    summarise(rep_mean = mean(death_status)) %>%
    ungroup() %>%
    group_by(group, time) %>%
    summarise(death_rate = mean(rep_mean)) %>%
    ungroup()
  compare_gcp = sim_data %>% arrange(id, time) %>%
    filter(time<=dropout_time) %>%
    filter(!death_status) %>%
    #filter(time<=dropout_time) %>%
    #mutate(time = as.character(time)) %>%
    #filter(time %in% c("0", "0.5"," 1", "1.5", "2", "2.5", "3", "3.5", "4", "4.5", "5", "5.5")) %>%
    select(id, time, gcp_werror,gcp_noerror, gcp_fix, group) %>%
    group_by(group, time) %>%
    summarise(avg_gcp_werror = mean(gcp_werror),
              avg_gcp_noerror = mean(gcp_noerror),
              avg_gcp_fix = mean(gcp_fix))%>%
    ungroup() #%>%
  #mutate(time = as.numeric(time)) 
  result_m = re_stroke %>% full_join(re_dementia,c("group","time")) %>%
    full_join(re_death, c("group", "time")) %>%
    full_join(compare_gcp, c("group", "time")) %>%
    mutate(repetition=m)
  
  result_summary = rbind(result_summary, result_m)
  
  
  
  rm(sim_next)
  rm(sim_data)
  rm(sim_base_control)
  rm(sim_base_intensive)
  rm(longrf_control)
  rm(longrf_intensive)
  rm(baseline_data)
  rm(long_data)
}


array_index = Sys.getenv("SLURM_ARRAY_TASK_ID")

result_df = data.frame(
  dementia_coef = dementia_group_coef,
  dementia_pvalue = dementia_group_pvalue,
  gcp_coef = gcp_group_coef,
  gcp_pvalue = gcp_group_pvalue,
  sim_N = simN_each,
  case = paste0(sample_N)
)
result_df = result_df %>% mutate(iteration = row_number(),
                                 array_id = array_index)



write.csv(result_df, file=paste0("6_600_summary/6_600_" , array_index , ".csv"))
write.csv(result_summary, file=paste0("6_600_eventsummary/6_600_event_" , array_index , ".csv"))






