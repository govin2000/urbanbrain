#!/usr/bin/env Rscript

# ===============================
# Publication-quality MRIQC plots
# ===============================
# Usage (from terminal):
#   Rscript make_mriqc_plots.R group_T1w.tsv output_plots/


suppressPackageStartupMessages({
  library(tidyverse)
})

# -------- Helper: nice theme --------
theme_pub <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid       = element_blank(),
      axis.title       = element_text(face = "bold"),
      axis.text        = element_text(color = "black"),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text       = element_text(face = "bold"),
      plot.title       = element_text(face = "bold", hjust = 0),
      legend.position  = "none"
    )
}

setwd("/Users/govindapoudel/Documents/Personal/Research/EPOCH/SydneyMAS/LongitudinalAnalysis/Script")

mriqc_csv<-"../FinalData/to_upload/bias_field.csv" 

out_dir = "../FinalData/mriqc_plots"
metrics = c("Min_Bias",	"Max_Bias",	"Mean_Bias",	"SD_Bias")
id_col = "subid"
  
df <- read_csv(mriqc_csv, show_col_types = FALSE)
  
df<-df%>%mutate(
  wave = case_when(
    session == "ses-1" ~ 1L,
    session == "ses-2" ~ 2L,
    session == "ses-3" ~ 4L,
    TRUE ~ NA_integer_
  )
)
# Identify grouping variable if present
group_var <- "wave"
  
  
df[[group_var]] <- as.factor(df[[group_var]])
  
  
message("Metrics to plot: ", paste(metrics, collapse = ", "))
  
  for (metric in metrics) {
    message("Processing metric: ", metric)
    
    # Skip if metric is all NA or non-numeric
    if (!is.numeric(df[[metric]])) {
      warning("Metric '", metric, "' is not numeric. Skipping.")
      next
    }
    
    metric_df <- df %>%
      dplyr::select(all_of(c(id_col, group_var, metric))) %>%
      dplyr::rename(value = !!metric) %>%
      dplyr::filter(!is.na(value))
    
    if (nrow(metric_df) == 0) {
      warning("No non-NA values for metric '", metric, "'. Skipping.")
      next
    }
    
    # Compute global z-score for outlier highlighting
    metric_df <- metric_df %>%
      mutate(
        z_value  = as.numeric(scale(value)),
        outlier  = abs(z_value) > 3
      )
    
    # Base plot
    if (!is.na(group_var)) {
      p <- ggplot(metric_df, aes_string(x = group_var, y = "value")) +
        geom_violin(trim = FALSE, fill = "grey90", color = "black") +
        geom_boxplot(width = 0.15, outlier.shape = NA) +
        geom_jitter(
          
          width = 0.1,
          alpha = 0.8,
          size = 1.5
        ) +
        scale_color_manual(
          values = c(`FALSE` = "black", `TRUE` = "red")
        ) +
        labs(
          x = "Wave",
          y = metric,
          title = paste0("Distribution of ", metric, " by ", group_var)
        )
    } else {
      # Single-group plot
      p <- ggplot(metric_df, aes(x = factor(1), y = value)) +
        geom_violin(trim = FALSE, fill = "grey90", color = "black") +
        geom_boxplot(width = 0.15, outlier.shape = NA) +
        geom_jitter(
          aes(color = outlier),
          width = 0.05,
          alpha = 0.8,
          size = 1.5
        ) +
        scale_color_manual(
          values = c(`FALSE` = "black", `TRUE` = "black")
        ) +
        labs(
          x = NULL,
          y = metric,
          title = paste0("Distribution of ", metric)
        ) +
        scale_x_discrete(labels = NULL)
    }
    
    p <- p + theme_pub()
    
    # Save plot
    out_file <- file.path(out_dir, paste0("biasfield_", metric, ".png"))
    ggsave(
      filename = out_file,
      plot     = p,
      width    = 7,
      height   = 5,
      units    = "in",
      dpi      = 300
    )
    
    message("  → Saved: ", out_file)
  }
  

