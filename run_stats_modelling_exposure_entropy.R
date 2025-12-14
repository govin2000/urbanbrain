# Load required libraries
library(lmerTest)  # for lmer()
library(dplyr)     # for mutate, select, group_by, ungroup, bind_rows, across
library(readr)     # for read_csv, write_csv
library(tidyr)     # for pivot_longer, pivot_wider, unite
library(tibble)    # 

smart_round <- function(x) {
  ifelse(abs(x) < 0.01, round(x, 5), round(x, 3))
}

run_lmer_models <- function(data, formula_text,region,exposure_var) {
  model <- lmer(as.formula(formula_text), data = data, weights = weights)
  
  model_summary <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE, conf.level = 0.95)
  
  model_summary <- model_summary %>%
    dplyr::mutate(Region = region, Exposure = exposure_var) %>%
    dplyr::select(Region, Exposure, term, estimate, conf.low, conf.high, statistic, p.value)
  return(list(model = model, estimates = model_summary))
}

reshape_table<-function(estimates_table){
  
  estimates_wide <- estimates_table %>% mutate(CI = paste0("(", smart_round(conf.low), ", ", smart_round(conf.high), ")"),
                                               p_value_fmt = ifelse(p.value < 0.05, paste0(smart_round(p.value)), as.character(round(p.value, 3))),
                                               Estimate_CI = paste0(smart_round(estimate), " ", CI),
                                               t_value_fmt = as.character(round(statistic,3), 2) ) %>%
    dplyr::select(Region, Exposure, Term = term, Estimate_CI, t_value = t_value_fmt, p_value_fmt)
  
  estimates_wide_pivot <- estimates_wide %>%
    pivot_longer(cols = c(Estimate_CI, t_value, p_value_fmt), names_to = "Stat", values_to = "Value") %>%
    unite("Term_Stat", Term, Stat) %>%
    pivot_wider(names_from = Term_Stat, values_from = Value) 
  
  estimates_wide_pivot <- estimates_wide_pivot %>%
    mutate(across(where(is.list), ~ sapply(., `[`, 1)))
  
   
  return(estimates_wide_pivot)
}

setwd("~/Documents/Personal/Research/EPOCH/SydneyMAS/LongitudinalAnalysis/Script")

# Load data
df <- read_csv("../FinalData/AllData_ent.csv")

df <- df %>%
  mutate(
    Sex = factor(Sex),
    NESB_Status = factor(NESB_Status),
    living = factor(living),
    hearaid = factor(hearaid),
    Scanner_Type=factor(Scanner_Type),
    agec = age - 70 + tisyear,
    agec2 = agec*agec,
    z_eTIV = as.numeric(scale(EstimatedTotalIntraCranialVol))
  )



direct_effect_covar <- c("Sex", "education", "NESB_Status",  "living","moves","pp_dn_n10",
                      "ses_av_n10","hearaid","z_eTIV","Scanner_Type","luco_pr_n10","lupa_pr_n10",
                      "blu_pr_c10","tree_pr_n10","tp_cn_n10", "ap_pm25_2", "ap_no2_2")

exposure_var<-"entropy"

# Create average and centered exposure variables to capture both between- and within-person effects

df <-df %>%
  group_by(MAS_ID) %>%
  mutate(temp_var = get(exposure_var),
         obs_n = sum(!is.na(temp_var)),
         avg_val = mean(temp_var, na.rm = TRUE),
         centered_val = ifelse(obs_n >= 2, temp_var - avg_val, NA_real_),
         !!paste0(exposure_var, "_avg") := as.numeric(avg_val),
         !!paste0(exposure_var, "_centered") := as.numeric(centered_val)
  )%>%
  select(-temp_var, -obs_n, -avg_val, -centered_val)%>%
  ungroup()

cat("\nRunning model for exposure:", exposure_var, "\n")


exposure_centered_var <- paste0(exposure_var, "_centered")
exposure_avg_var <- paste0(exposure_var, "_avg")




  estimates_table_direct <- tibble()
  df <-df %>%mutate(
    tail=(rh_Hippocampal_tail+lh_Hippocampal_tail),
    body=(rh_Whole_hippocampal_body+lh_Whole_hippocampal_body),
    head=(rh_Whole_hippocampal_head+lh_Whole_hippocampal_head)
  )
  

  hippocampal_regions<-c(
    "head", "body", "tail"
  )
 
  region_type='head_body_tail'
  
  for (region in hippocampal_regions) {
          
     
      
      formula_text2 <- paste0(region, "~ agec+ agec2+", exposure_avg_var,"+",exposure_centered_var,
                              "+agec:",exposure_avg_var,"+agec:",exposure_centered_var,
                              "+agec2:",exposure_avg_var,"+agec2:",exposure_centered_var,"+",
                              paste(direct_effect_covar, collapse = " + "),
                              " + (1+agec|MAS_ID)")
      lmer_output2<-run_lmer_models(df,  formula_text2,region, exposure_var) 
      estimates_table_direct<-bind_rows(estimates_table_direct,lmer_output2$estimates)
      
       
    }
 

  
  

  estimates_wide_pivot5<-reshape_table(estimates_table_direct)
  selected_colnames <- names(estimates_wide_pivot5)[grepl(exposure_var, names(estimates_wide_pivot5))]
  selected_estimates_wide_pivot5 <- estimates_wide_pivot5[, c("Region",selected_colnames)]
  selected_estimates_wide_pivot5$entropy_avg_p_value_adj<-p.adjust(selected_estimates_wide_pivot5[["entropy_avg_p_value_fmt"]],method = "fdr")
  selected_estimates_wide_pivot5$agec_p_value_adj<-p.adjust(selected_estimates_wide_pivot5[["agec:entropy_avg_p_value_fmt"]],method = "fdr")
  selected_estimates_wide_pivot5$agec2_p_value_adj<-p.adjust(selected_estimates_wide_pivot5[["agec2:entropy_avg_p_value_fmt"]],method = "fdr")
  
  selected_estimates_wide_pivot5$entropy_centered_p_value_adj<-p.adjust(selected_estimates_wide_pivot5[["entropy_centered_p_value_fmt"]],method = "fdr")
  selected_estimates_wide_pivot5$agec_centered_p_value_adj<-p.adjust(selected_estimates_wide_pivot5[["agec:entropy_centered_p_value_fmt"]],method = "fdr")
  selected_estimates_wide_pivot5$agec2_centered_p_value_adj<-p.adjust(selected_estimates_wide_pivot5[["agec2:entropy_centered_p_value_fmt"]],method = "fdr")
  
  
  df2_sel <- selected_estimates_wide_pivot5 %>%
    select(Region,
           main_effect_CI = entropy_avg_Estimate_CI,
           main_effect_t_value=entropy_avg_t_value,
           main_effect_p_value = entropy_avg_p_value_adj,
           linear_CI = `agec:entropy_avg_Estimate_CI`,
           linear_t_value=`agec:entropy_avg_t_value`,
           linear_p_value = agec_p_value_adj,
           nonlinear_CI = `agec2:entropy_avg_Estimate_CI`,
           nonlinear_t_value=`agec2:entropy_avg_t_value`,
           nonlinear_p_value = agec2_p_value_adj
    )
  write_csv(df2_sel, paste0("../FinalData/Table3_lmer_model_direct",region_type,'_', exposure_var, ".csv"))
  
  
  df3_sel <- selected_estimates_wide_pivot5 %>%
    select(Region,
           main_effect_CI = entropy_centered_Estimate_CI,
           main_effect_t_value=entropy_centered_t_value,
           main_effect_p_value = entropy_centered_p_value_adj,
           linear_CI = `agec:entropy_centered_Estimate_CI`,
           linear_t_value=`agec:entropy_centered_t_value`,
           linear_p_value = agec_centered_p_value_adj,
           nonlinear_CI = `agec2:entropy_centered_Estimate_CI`,
           nonlinear_t_value=`agec2:entropy_centered_t_value`,
           nonlinear_p_value = agec2_centered_p_value_adj
    )
  write_csv(df3_sel, paste0("../FinalData/Table4_lmer_model_direct",region_type,'_', exposure_var, ".csv"))
  
  
  
  