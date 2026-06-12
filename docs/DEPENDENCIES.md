# Software dependencies

Author: Wang Linglong

The public release contains R and Python scripts. Package use was preserved from the original analysis code.

## R

Core packages used across scripts include:

- `ggplot2`, `dplyr`, `tidyr`, `readr`, `tidyverse`, `forcats`
- `vegan`, `ape`, `ggrepel`, `ggsci`, `showtext`, `grid`
- `ComplexHeatmap`, `circlize`, `cols4all`, `randomcoloR`
- `sf`, `ggforce`, `ggalluvial`, `igraph`, `ggraph`
- `MuMIn`, `broom`, `mgcv`, `pheatmap`

See `docs/r_dependency_inventory.csv` for script-level package calls extracted from the code.

## Python

The Mantel analysis script uses:

- `pandas`
- `numpy`
- `scipy`
- `matplotlib`

## Reproducibility note

The code release standardizes script names, encoding, authorship metadata, and public directory structure. The original statistical and plotting logic is preserved. Some figure scripts require input tables or map files that must be deposited in `data/public_figures/` or supplied through `JEV_FIGURE_INPUT_DIR`.
