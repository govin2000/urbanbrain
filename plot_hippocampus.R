#!/usr/bin/env Rscript

library(tidyverse)

# ----------------------------------------------------------
# 1. Load data
# ----------------------------------------------------------
data <- read_csv("../FinalData/AllData_ent.csv")

data <- data %>%
  mutate(
    wave = factor(wave, levels = c(1, 2, 4),
                  labels = c("Wave 1", "Wave 2", "Wave 4"))
  )

# ----------------------------------------------------------
# 2. Create total head, body, tail volumes
#    (edit column names here if needed)
# ----------------------------------------------------------
data <- data %>%
  mutate(
    total_head = (rh_Whole_hippocampal_head + lh_Whole_hippocampal_head),
    total_body = (rh_Whole_hippocampal_body + lh_Whole_hippocampal_body),
    total_tail = (rh_Hippocampal_tail        + lh_Hippocampal_tail)
  )

# ----------------------------------------------------------
# 3. Long format
# ----------------------------------------------------------
long_dat <- data %>%
  select(wave, total_head, total_body, total_tail) %>%
  pivot_longer(
    cols      = c(total_head, total_body, total_tail),
    names_to  = "region",
    values_to = "volume"
  ) %>%
  mutate(
    region = factor(
      region,
      levels = c("total_head", "total_body", "total_tail"),
      labels = c("Total Head Volume",
                 "Total Body Volume",
                 "Total Tail Volume")
    )
  )

# ----------------------------------------------------------
# 4. Outlier detection (within wave × region)
# ----------------------------------------------------------
is_outlier_iqr <- function(x) {
  q1  <- quantile(x, 0.25, na.rm = TRUE)
  q3  <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  (x < (q1 - 1.5 * iqr)) | (x > (q3 + 1.5 * iqr))
}

long_dat <- long_dat %>%
  group_by(region, wave) %>%
  mutate(outlier = is_outlier_iqr(volume)) %>%
  ungroup()



# ----------------------------------------------------------
# 5. Publication-quality violin plots (colourful)
# ----------------------------------------------------------
# Custom palette for waves (soft but distinct)
wave_cols <- c(
  "Wave 1" = "#4E79A7",  # blue
  "Wave 2" = "#59A14F",  # green
  "Wave 4" = "#F28E2B"   # orange
)

p <- ggplot(long_dat, aes(x = wave, y = volume)) +
  geom_violin(aes(fill = wave), trim = FALSE, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white",
               color = "black", linewidth = 0.4) +
  geom_jitter(aes(color = outlier), width = 0.06, alpha = 0.8, size = 1.3) +
  scale_fill_manual(values = wave_cols, name = "Wave") +
  scale_color_manual(
    values = c(`FALSE` = "black", `TRUE` = "#D62728"),
    breaks = c(FALSE, TRUE),
    labels = c("Normal", "Outlier"),
    name = "Outlier status"
  ) +
  facet_wrap(~ region, scales = "free_y", nrow = 1) +
  labs(x = "Wave", y = "Total volume (L+R) (mm³)") +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 14),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.4)
  )



print(p)

# ----------------------------------------------------------
# 6. Save high-resolution versions
# ----------------------------------------------------------
ggsave(
  filename = "../FinalData/hipp_total_head_body_tail_violin.png",
  plot     = p,
  width    = 11,
  height   = 4,
  dpi      = 600
)

# (Optional) TIFF for journals that prefer it
ggsave(
  filename = "../FinalData/hipp_total_head_body_tail_violin.tiff",
  plot     = p,
  width    = 11,
  height   = 4,
  dpi      = 600,
  compression = "lzw"
)
