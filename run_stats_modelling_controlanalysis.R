# Load required libraries
library(lmerTest)  # for lmer()
library(dplyr)     # for mutate, select, group_by, ungroup, bind_rows, across
library(readr)     # for read_csv, write_csv
library(tidyr)     # for pivot_longer, pivot_wider, unite
library(tibble)    # 

setwd("~/Documents/Personal/Research/EPOCH/SydneyMAS/LongitudinalAnalysis/Script")

# Load data
df_raw <- read_csv("../FinalData/AllData_Modelled.csv")
df_raw <- df_raw %>%
  mutate(
    Sex = factor(Sex),
    NESB_Status = factor(NESB_Status),
    living = factor(living),
    hearaid = factor(hearaid),
    Scanner_Type=factor(Scanner_Type),
    total_hippo=head+body+tail)
    

#1. Compute IQR-based bounds for tail within each wave
hippo_bounds <- df_raw %>%
  group_by(wave) %>%
  summarise(
    Q1    = quantile(total_hippo, 0.25, na.rm = TRUE),
    Q3    = quantile(total_hippo, 0.75, na.rm = TRUE),
    IQR   = Q3 - Q1,
    lower = Q1 - 1.5 * IQR,
    upper = Q3 + 1.5 * IQR,
    .groups = "drop"
  )

# 1. Compute IQR-based bounds for tail within each wave
# tail_bounds <- df_raw %>%
#   group_by(wave) %>%
#   summarise(
#     Q1    = quantile(tail, 0.25, na.rm = TRUE),
#     Q3    = quantile(tail, 0.75, na.rm = TRUE),
#     IQR   = Q3 - Q1,
#     lower = Q1 - 1.5 * IQR,
#     upper = Q3 + 1.5 * IQR,
#     .groups = "drop"
#   )

# # 2. Join bounds back to the main data and flag outliers
# df_flagged <- df_raw %>%
#   left_join(tail_bounds, by = "wave") %>%
#   mutate(
#     tail_outlier = !is.na(tail) & (tail < lower | tail > upper)
#   )

# 2. Join bounds back to the main data and flag outliers
df_flagged <- df_raw %>%
  left_join(hippo_bounds, by = "wave") %>%
  mutate(
    hippo_outlier = !is.na(total_hippo) & (total_hippo < lower | total_hippo > upper)
  )

# 3. Quick check: how many outliers per wave?
df_flagged %>%
  count(wave, hippo_outlier)

df_clean <- df_flagged %>%
  filter(!hippo_outlier) %>%           # drop outlier rows
  select(-Q1, -Q3, -IQR, -lower, -upper, -hippo_outlier)

region<-c( "tail")


direct_effect_covar <- c("Sex", "education", "NESB_Status",  "living","moves","pp_dn_n10",
                         "ses_av_n10","hearaid","z_eTIV","Scanner_Type","luco_pr_n10","lupa_pr_n10",
                         "blu_pr_c10","tree_pr_n10","tp_cn_n10", "ap_pm25_2", "ap_no2_2")

exposure_var<-"si_dn_n10"
exposure_centered_var <- paste0(exposure_var, "_centered")
exposure_avg_var <- paste0(exposure_var, "_avg")

estimates_table_basic <- tibble()
formula_text2 <- paste0(region, "~ agec+ agec2+", exposure_avg_var,"+",exposure_centered_var,
                        "+agec:",exposure_avg_var,"+agec:",exposure_centered_var,
                        "+agec2:",exposure_avg_var,"+agec2:",exposure_centered_var,"+",
                        paste(direct_effect_covar, collapse = " + "),
                        " + (1+agec|MAS_ID)")
lmer_output2<-run_lmer_models(df_clean,  formula_text2,region, exposure_var) 

estimates_wide_pivot5<-reshape_table(lmer_output2$estimates)
selected_colnames <- names(estimates_wide_pivot5)[grepl(exposure_var, names(estimates_wide_pivot5))]
selected_estimates_wide_pivot5 <- estimates_wide_pivot5[, c("Region",selected_colnames)]

write_csv(selected_estimates_wide_pivot5,"../FinalData/ControlAnalysis_OutlierRemoved_Tail.csv")


#control analysis for subject which are in all waves
# Find how many waves exist in the dataset
required_waves <- 3

# Filter subjects who have all waves
df <- df_raw %>%
  group_by(MAS_ID) %>%
  filter(n_distinct(wave) == required_waves) %>%
  arrange(MAS_ID, wave) %>%
  ungroup()
   
region<-c( "tail")


direct_effect_covar <- c("Sex", "education", "NESB_Status",  "living","moves","pp_dn_n10",
                      "ses_av_n10","hearaid","z_eTIV","Scanner_Type","luco_pr_n10","lupa_pr_n10",
                      "blu_pr_c10","tree_pr_n10","tp_cn_n10", "ap_pm25_2", "ap_no2_2")

exposure_var<-"si_dn_n10"
exposure_centered_var <- paste0(exposure_var, "_centered")
exposure_avg_var <- paste0(exposure_var, "_avg")

estimates_table_basic <- tibble()
formula_text2 <- paste0(region, "~ agec+ agec2+", exposure_avg_var,"+",exposure_centered_var,
                              "+agec:",exposure_avg_var,"+agec:",exposure_centered_var,
                              "+agec2:",exposure_avg_var,"+agec2:",exposure_centered_var,"+",
                              paste(direct_effect_covar, collapse = " + "),
                              " + (1+agec|MAS_ID)")
lmer_output2<-run_lmer_models(df,  formula_text2,region, exposure_var) 
esti<-reshape_table(lmer_output2$estimates)

esti$si_dn_n10_avg_Estimate_CI
esti$`agec2:si_dn_n10_avg_Estimate_CI`
esti$`agec:si_dn_n10_avg_Estimate_CI`


## Keep only Wave 1 and Wave 2
df_w1 <- df_raw %>% 
  filter(wave %in% c(1)) %>% 
  arrange(MAS_ID, wave)

df_w1 <- df_w1 %>%
  mutate(
    dage  = age - mean(age, na.rm = TRUE),
    dage2 = dage^2
  )

# Identify outliers using IQR rule
Q1  <- quantile(df_w1$tail, 0.25, na.rm = TRUE)
Q3  <- quantile(df_w1$tail, 0.75, na.rm = TRUE)
IQR <- Q3 - Q1

lower_bound <- Q1 - 1.5 * IQR
upper_bound <- Q3 + 1.5 * IQR

# Filter out outliers
df_w1_no_out <- df_w1 %>%
  filter(tail >= lower_bound & tail <= upper_bound)

# Check number removed
nrow(df_w1) - nrow(df_w1_no_out)
model_cross <- lm(
  tail ~ si_dn_n10 + Sex + dage+dage2+education + NESB_Status + living +
    pp_dn_n10 + ses_av_n10 + hearaid + z_eTIV + Scanner_Type +
    luco_pr_n10 + lupa_pr_n10 + blu_pr_c10 + tree_pr_n10 +
    tp_cn_n10 + ap_pm25_2 + ap_no2_2,
  data = df_w1
)

s <- summary(model_cross)

estimate <- s$coefficients["si_dn_n10", "Estimate"]
t_value  <- s$coefficients["si_dn_n10", "t value"]
p_value  <- s$coefficients["si_dn_n10", "Pr(>|t|)"]

df_raw %>%
  ggplot(aes(x = tail)) +
  geom_histogram(bins = 30, color = "black", fill = "skyblue", alpha = 0.7) +
  labs(
    title = "Distribution of Hippocampal Tail Volume (Wave 1)",
    x = "Hippocampal Tail Volume",
    y = "Count"
  ) +
  theme_minimal(base_size = 14)


## Keep only Wave 1 and Wave 2
df_w12 <- df_raw %>% 
  filter(wave %in% c(1, 2)) %>% 
  arrange(MAS_ID, wave)

df_w1 

lmer_output3<-run_lmer_models(df_w12,  formula_text2,region, exposure_var) 





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

smart_round <- function(x) {
  ifelse(abs(x) < 0.01, round(x, 5), round(x, 3))
}