
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
