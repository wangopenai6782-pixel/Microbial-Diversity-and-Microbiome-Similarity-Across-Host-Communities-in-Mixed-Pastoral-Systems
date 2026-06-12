# Public code release

Author: Wang Linglong

This repository contains the analysis and figure-generation code required to reproduce the analyses in the associated manuscript.
The public scripts were organized from the author's original `.txt` code files and converted to UTF-8 `.R`/`.py` files with stable English file names.

## Code availability statement

All code required to reproduce the analyses in this paper is available on GitHub at [add GitHub repository URL here].

## Repository layout

- `scripts/public_figures/`: figure and statistical-analysis scripts converted from the original manuscript code.
- `scripts/run_all_public_figures.R`: helper script to run all R figure scripts in lexical order.
- `data/public_figures/`: expected location for public input tables and map files.
- `results/public_figures/`: default output location for regenerated figures and tables.
- `docs/script_manifest.csv`: mapping from public file names to original local `.txt` source files.
- `docs/path_inventory.csv`: inventory of portable input/output paths referenced by the scripts.
- `docs/DEPENDENCIES.md`: R/Python dependency summary.

## Quick start

Install R package dependencies, place the required input files under `data/public_figures/`, and run:

```r
Sys.setenv(JEV_FIGURE_INPUT_DIR = "data/public_figures")
Sys.setenv(JEV_FIGURE_OUTPUT_DIR = "results/public_figures")
source("scripts/run_all_public_figures.R", encoding = "UTF-8")
```

Individual scripts can also be run directly, for example:

```bash
Rscript scripts/public_figures/figure1_species_density_violin.R
```

The Mantel distance script is Python-based:

```bash
python scripts/public_figures/supplementary_mantel_distance_relationship.py
```

## Notes for release

Before publishing on GitHub, replace `[add GitHub repository URL here]` with the final repository link.
If large or restricted datasets cannot be committed to GitHub, add accession numbers or repository links in this README and keep the code paths documented in `docs/path_inventory.csv`.

