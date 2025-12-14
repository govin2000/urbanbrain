suppressPackageStartupMessages(library(tidyverse))

# ---- Load + tidy ----
dat <- read_csv("../FinalData/AllData.csv") %>%
  mutate(
    wave = factor(wave, levels = c(1, 2, 4), labels = c("Wave 1", "Wave 2", "Wave 4")),
    sid  = si_dn_n10
  ) %>%
  select(wave, sid) %>%
  filter(!is.na(sid))

# ---- Outliers per wave (z > 3) ----
dat <- dat %>%
  group_by(wave) %>%
  mutate(outlier = abs(as.numeric(scale(sid))) > 3) %>%
  ungroup()

ggplot(dat, aes(wave, sid)) +
  geom_violin(aes(fill = wave), trim = FALSE, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.4) +
  geom_jitter(aes(color = outlier), width = 0.06, alpha = 0.8, size = 1.3) +
  scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "#D62728"), guide = "none") +
  labs(x = "Wave", y = "Street intersection density (/km²)") +
  theme_classic(base_size = 14) +
  theme(
    axis.title   = element_text(face = "bold"),
    axis.text    = element_text(color = "black"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
  )