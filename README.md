# ENDURE shipwreck deterioration — reproducible publication package

This directory is designed to reproduce the analyses and figures reported in the
ENDURE shipwreck deterioration manuscript without access to the original development
filesystem.

## Craig: one-time preparation on Ubuntu

1. Copy/unzip this repository to:

   `/media/syms/ExtremeSSD/ShipwreckENDURE_Publication`

2. Open `analysis/00_prepare_publication_data.ipynb`.
3. Run it top-to-bottom.

It reads the original development files from:

- `/media/syms/ExtremeSSD/ShipwreckENDURE/ABMsubset.nc`
- `/media/syms/ExtremeSSD/ShipwreckENDURE/wreck_agent_decay_for_gam.parquet`

and writes only:

- `data/model_inputs/publication_model_inputs.nc`
- `data/calibration/wreck_calibration.csv`
- `data_manifest.json`

The calibration CSV contains only the four analytical variables used in the GAM.
Wreck names, identifiers and coordinates are not released.

## Reproducing the paper

After data preparation, the remaining notebooks use only repository-relative paths
and should run on Linux, macOS or Windows:

1. `01_calibration.ipynb`
2. `02_generate_manuscript_outputs.ipynb`
3. `03_sensitivity_analysis.ipynb`
4. `04_make_figures.ipynb`

Notebook 01 regenerates `data/gam_age_curve.csv`.

Notebook 02 regenerates threshold maps, model-driver maps, representative trajectories
and the compact variance dataset used in the manuscript and supplement.

Notebook 03 regenerates the final sensitivity CSV products.

Notebook 04 regenerates the manuscript and supplementary figures from those outputs.

## Platform notes

The public analysis notebooks contain no `/media/syms/...` paths. Only Notebook 00,
which is a one-time data-extraction step for the author, knows about the original SSD.

A collaborator on macOS should be able to clone/copy the prepared repository and run
Notebooks 01–04 without editing filesystem paths.

## Important caveat

This package reconstructs the publication analysis from the frozen derived environmental
input. It does not attempt to reproduce the much larger upstream Copernicus / EMODnet
download and preprocessing archive. Those source products and preprocessing provenance
should be documented separately in the manuscript/README.
