# Data files for public figure scripts

Place manuscript figure input tables and geospatial files in this directory before running the scripts.
The scripts were converted from the author's original local `.txt` files and still document the original
file names in `docs/path_inventory.csv`.

Recommended layout:

- `data/public_figures/*.csv` for tabular inputs.
- `data/public_figures/China_shapefile/` for China boundary shapefiles, if using the original map script.
- `data/public_figures/nmg_shapefile/` for Inner Mongolia boundary shapefiles, if using the original map script.

Large raw sequencing data and privacy-sensitive field metadata should be deposited separately or documented
with accession numbers before GitHub release.
