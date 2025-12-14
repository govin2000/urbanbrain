library(lme4)
library(ggeffects)
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(readr)
library(lmerTest)


setwd("~/Documents/Personal/Research/EPOCH/SydneyMAS/LongitudinalAnalysis/Script")

# Load data
df <- read_csv("../FinalData/AllData_Modelled.csv")



regions<-c("head", "body", "tail")


# Function to fit model and extract prediction with CIs for each region
get_region_predictions <- function(region_var) {
  
  # Fit model
  formula_str <- paste0(region_var, " ~ agec + I(agec^2) + Sex + education + living + ",
                        "moves + z_eTIV + Scanner_Type+hearaid + ses_av_n10 + pp_dn_n10 + (1 + agec | MAS_ID)")
  model <- lmer(as.formula(formula_str), data = df, weights = weights)
  
  # Get population-level predictions using ggpredict
  preds <- ggpredict(
    model,
    terms = "agec [0:25 by=0.5]",  # corresponds to age 70–95
    condition = list(
      Sex = "Female",
      living = "community - with spouse",
      hearaid = "no",
      education = mean(df$education, na.rm = TRUE),
      moves = mean(df$moves, na.rm = TRUE),
      z_eTIV = mean(df$z_eTIV, na.rm = TRUE),
      ses_av_n10 = mean(df$ses_av_n10, na.rm = TRUE),
      pp_dn_n10 = mean(df$pp_dn_n10, na.rm = TRUE)
    )
  )
  
  # Calculate percent of baseline volume
  preds$age <- preds$x + 70
  baseline <- preds$predicted[preds$age == 70]
  preds$percent <- 100 * (preds$predicted / baseline)
  preds$ci.low_percent <- 100 * (preds$conf.low / baseline)
  preds$ci.high_percent <- 100 * (preds$conf.high / baseline)
  preds$Region <- region_var
  
  return(preds)
}

# Get predictions for all regions
all_preds <- map_dfr(regions, get_region_predictions)

ages_of_interest <- c(75, 80, 85, 90, 95)

ci_points <- all_preds %>%
  filter(age %in% ages_of_interest)

pd <- position_dodge(width = 0.8)

ggplot(all_preds, aes(x = age, y = percent, color = Region)) +
  geom_line(size = 1.3) +
  geom_errorbar(
    data = ci_points,
    aes(ymin = ci.low_percent, ymax = ci.high_percent),
    width = 0.6,
    size  = 0.8,
    position = pd  )+
  labs(
    title = "Predicted Change in Tail Volume",
    x = "Age (years)",
    y = "Volume Change (% of Baseline)"
  ) +
  scale_y_continuous(limits = c(NA, 100)) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid = element_blank()
  )
ggsave(
  filename = "../FinalData/hippocampal_volume_retention.tiff",
  plot = last_plot(),  # or replace with your plot object, e.g., plot = my_plot
  width = 7,            # inches
  height = 5,           # inches
  dpi = 600,            # dots per inch
  units = "in",
  compression = "lzw"   # recommended for TIFF
)




#Generate Interaction Plots
# Population-level predictions at -1.5 SD, Mean, +1.5 SD
library(multcomp)

exposure_var<-"si_dn_n10"

exposure_centered_var <- paste0(exposure_var, "_centered")
exposure_avg_var <- paste0(exposure_var, "_avg")

formula_text<-"tail~ agec+ agec2+si_dn_n10_avg+si_dn_n10_centered+
              agec:si_dn_n10_avg+agec:si_dn_n10_centered+agec2:si_dn_n10_avg+agec2:si_dn_n10_centered+
              Sex + education + NESB_Status + living + moves + pp_dn_n10 + ses_av_n10 + hearaid + z_eTIV+Scanner_Type+
              luco_pr_n10 +lupa_pr_n10 + blu_pr_c10 + tree_pr_n10 + tp_cn_n10 + ap_pm25_2 + ap_no2_2 + (1+agec|MAS_ID)"


model <- lmer(as.formula(formula_text), data = df, weights = weights)

# partial residuals
X     <- model.matrix(model)
beta  <- fixef(model)
yhat_full  <- as.vector(X %*% beta)  # fixed-effects prediction
resid_full <- residuals(model)

beta_x <- beta["si_dn_n10_avg"]

df$partial_resid <- resid_full + beta_x * df$si_dn_n10_avg


resid_plot <- ggplot(df, aes(x = si_dn_n10_avg, y = partial_resid)) +
  geom_point(alpha = 0.4, size = 1, colour = "grey40") +
  geom_smooth(method = "lm", se = TRUE, colour = "#08519C", fill = "#87CEEB") +
  labs(
    title = "Association Between SID and Tail Volume",
    x = "Street intersection density (person-level, /km²)",
    y = "Total tail volume (Partial Residuals) (mm³)"
  ) +
  theme_classic(base_size = 16)
ggsave(
  filename = "../FinalData/tail_volume_vs_si.tiff",
  plot = last_plot(),  # or replace with your plot object, e.g., plot = my_plot
  width = 7,            # inches
  height = 5,           # inches
  dpi = 600,            # dots per inch
  units = "in",
  compression = "lzw"   # recommended for TIFF
)

# Get model matrix and fixed effects only
X <- model.matrix(model)
beta <- fixef(model)

# Predicted values from fixed effects only
y_hat_fixed <- as.vector(X %*% beta)

y <- df$tail

resid_full <- residuals(model)

partial_resid <- resid_full + (beta["si_dn_n10_avg"] * df$si_dn_n10_avg)

df$partial_resid <- partial_resid
df$partial_x <- df$si_dn_n10_avg

# Get mean and SD for SID
sid_mean <- mean(df[[exposure_avg_var]], na.rm = TRUE)
sid_sd <- sd(df[[exposure_avg_var]], na.rm = TRUE)



# SID levels
sid_levels <- c("Low" = sid_mean - 1.5 * sid_sd,
                "Mean" = sid_mean,
                "High" = sid_mean + 1.5 * sid_sd)

# Ages (agec, since model uses centered age)
ages <- c( 5, 15,  25)  # corresponds to age 75, 80, 85

# Container for all contrast results
contrast_df <- data.frame()

for (sid_label in names(sid_levels)) {
  sid_val <- sid_levels[[sid_label]]
  
  for (a in ages) {
    a2 <- a^2
    
    # Construct contrast expression
    contrast_expr <- paste0(
      "agec*", a,
      " + agec2*", a2,
      " + agec:si_dn_n10_avg*", sid_val, "*", a,
      " + agec2:si_dn_n10_avg*", sid_val, "*", a2,
      " = 0"
    )
    
    # Run contrast test
    glht_out <- glht(model, linfct = c(contrast_expr))
    ci <- confint(glht_out)$confint
    est <- ci[1, "Estimate"]
    lower <- ci[1, "lwr"]
    upper <- ci[1, "upr"]
    pval <- summary(glht_out)$test$pvalues[1]
    
    # Store in dataframe
    contrast_df <- bind_rows(contrast_df, data.frame(
      SID = sid_label,
      Age = a + 70,
      Estimate = est,
      CI_low = lower,
      CI_high = upper,
      p_value = pval
    ))
  }
}

# Plot
p<-ggplot(contrast_df, aes(x = Age, y = Estimate, color = SID)) +
  geom_line(size = 1.2) +
  geom_point(position = position_dodge(width = 0.3), size = 3) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2,
                position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(
    title = "Marginal Effect of SID on Tail Volume",
    x = "Age (years)",
    y = "Change from Age 70 (mm³)",
    color = "SID Level"
  ) +
  theme_classic(base_size = 14)

ggsave("../FinalData/InteractionPlot.png", plot = p, width = 6, height = 5, dpi = 600)


