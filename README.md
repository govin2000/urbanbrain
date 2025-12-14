#Key R scripts for modelling and generating images used in the manuscript detailed below:
"Neighbourhood street connectivity and hippocampus volume in older adults"

run_stats_modelling.R - Linear mixed effects modelling script
run_stats_modelling_exposure_entropy.R - LME for second exposure


## Requirements

### R version
- R >= 4.2 recommended

### R packages
Most scripts rely on:
- `tidyverse`
- `lme4`, `lmerTest` (mixed-effects models)
- `broom.mixed` (tidy model outputs)
- `emmeans` (marginal effects / contrasts; optional)
- `performance`, `car` (diagnostics; optional)
- `patchwork` or `cowplot` (plot layouts; optional)
- `readr`, `stringr`, `glue` (usually via tidyverse)

Install packages:
```r
install.packages(c(
  "tidyverse","lme4","lmerTest","broom.mixed","emmeans","performance","car",
  "patchwork","cowplot"
))
