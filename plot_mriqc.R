#!/usr/bin/env Rscript

# ===============================
# Publication-quality MRIQC plots
# ===============================
# Now: all metrics in one figure, one facet per metric
# Usage:
#   Rscript make_mriqc_plots.R

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
      legend.position  = "bottom"
    )
}

setwd("/Users/govindapoudel/Documents/Personal/Research/EPOCH/SydneyMAS/LongitudinalAnalysis/Script")

mriqc_csv <- "../FinalData/mriqc_qc_ex.csv"
out_dir   <- "../FinalData/mriqc_plots"

# Create output dir if needed
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

metrics <- c("cjv", "cnr", "wm2max", "fber", "qi_2")
id_col  <- "subid"
group_var <- "wave"

# -------------------------
# Load + reshape the data
# -------------------------
df <- read_csv(mriqc_csv, show_col_types = FALSE) %>%
  mutate(
    wave = case_when(
      session == "ses-1" ~ 1L,
      session == "ses-2" ~ 2L,
      session == "ses-3" ~ 4L,
      TRUE ~ NA_integer_
    ),
    wave = factor(wave)
  )

message("Metrics to plot: ", paste(metrics, collapse = ", "))

# Keep only the relevant columns and pivot to long format
long_df <- df %>%
  select(all_of(c(id_col, group_var, metrics))) %>%
  pivot_longer(
    cols      = all_of(metrics),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  group_by(metric) %>%
  # z-score per metric to mark outliers
  mutate(
    z_value = as.numeric(scale(value)),
    outlier = abs(z_value) > 3
  ) %>%
  ungroup()

# ---------------
# Build the plot
# ---------------
p <- ggplot(long_df, aes(x = !!sym(group_var), y = value)) +
  # Violin per wave, coloured by wave
  geom_violin(
    aes(fill = !!sym(group_var)),
    trim  = FALSE,
    color = "black",
    alpha = 0.7
  ) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  # Points, coloured by outlier status
  geom_jitter(
    aes(color = outlier),
    width = 0.1,
    alpha = 0.8,
    size  = 1.5
  ) +
  scale_fill_brewer(palette = "Set2", name = "Wave") +
  scale_color_manual(
    values = c(`FALSE` = "black", `TRUE` = "red"),
    name   = "Outlier (|z| > 3)"
  ) +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  labs(
    x     = "Wave",
    y     = "MRIQC metric value",
    title = "MRIQC metrics by wave"
  ) +
  theme_pub()


# ---------------
# Save the figure (vector + optional PNG)
# ---------------

# Vector PDF (recommended for Nature)
out_file_pdf <- file.path(out_dir, "mriqc_all_metrics.pdf")
ggsave(
  filename = out_file_pdf,
  plot     = p,
  device   = cairo_pdf,   # needs Cairo; gives nice embedded fonts
  width    = 10,
  height   = 6,
  units    = "in"
)

# Optional: EPS version (also vector, sometimes requested)
out_file_eps <- file.path(out_dir, "mriqc_all_metrics.eps")
ggsave(
  filename = out_file_eps,
  plot     = p,
  device   = cairo_ps,
  width    = 10,
  height   = 6,
  units    = "in"
)

# Optional: keep a PNG for quick viewing
out_file_png <- file.path(out_dir, "mriqc_all_metrics.png")
ggsave(
  filename = out_file_png,
  plot     = p,
  width    = 10,
  height   = 6,
  units    = "in",
  dpi      = 300
)

message("→ Saved PDF: ", out_file_pdf)
message("→ Saved EPS: ", out_file_eps)
message("→ Saved PNG: ", out_file_png)
