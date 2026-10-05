# Environmental forcing and deterioration trajectories in submerged shipwreck heritage

Code and derived data supporting the manuscript:

**Syms, C., McGonigle, C., Quinn, R. & Gregory, D.**\
*Environmental forcing and deterioration trajectories in submerged
shipwreck heritage*.\
*Heritage Science* (under review).

This repository is the publication-scoped reproducibility package for
the manuscript. It contains the custom model code, the minimal derived
environmental input used by the model, an anonymised wreck calibration
dataset, and notebooks that reproduce the analyses and manuscript
figures.

It is not the complete ENDURE development environment. Ongoing
development of the shipwreck deterioration framework, including
interactive and practitioner-facing applications, is maintained
separately.

## Repository structure

``` text
src/
    decay_model.py
    model_config.py
    response_functions.py
    sensitivity_helpers.py
    wave_physics.py

analysis/
    01_calibration.ipynb
    02_generate_manuscript_outputs.ipynb
    03_sensitivity_analysis.ipynb
    04_make_figures.ipynb

data/
    model_inputs/
        publication_model_inputs.nc
    calibration/
        wreck_calibration.csv
    gam_age_curve.csv

outputs/
    figure_data/
    figure_maps/
    driver_maps/
    manuscript_figures/
```

## Reproducibility boundary

The repository starts from a frozen, derived regional environmental
dataset rather than the full upstream Copernicus and EMODnet download
archive.

`data/model_inputs/publication_model_inputs.nc` contains the
environmental fields required to reproduce the published model analyses
over the study domain. The much larger source datasets and preprocessing
archive are not duplicated here; their provenance and processing are
described in the manuscript.

The archaeological calibration dataset has been deliberately reduced to
the four variables used in the fitted GAM:

-   `Obs`
-   `chronological_age`
-   `log_CR_physics`
-   `log_CR_bio`

Wreck names, site identifiers and geographic coordinates have been
removed because they derive from restricted heritage records and are not
required to reproduce the reported analysis.

## Running the analyses

Run the notebooks in this order:

``` text
01_calibration.ipynb
02_generate_manuscript_outputs.ipynb
03_sensitivity_analysis.ipynb
04_make_figures.ipynb
```

### 01 --- Calibration

Fits the empirical age-condition GAM using the anonymised wreck
calibration dataset and writes:

``` text
data/gam_age_curve.csv
```

### 02 --- Manuscript model outputs

Runs the deterioration model over the regional environmental input and
regenerates the spatial threshold maps, environmental-driver rasters,
representative trajectories and associated figure data.

GeoTIFF outputs are written as north-up WGS84 rasters.

### 03 --- Sensitivity analysis

Runs the one-at-a-time parameter sensitivity analyses and writes the CSV
products used by the supplementary sensitivity figures.

### 04 --- Figures

Generates the manuscript and supplementary figures from the calibration
and model outputs produced above.

## One-time author data preparation

`00_prepare_publication_data.ipynb` is retained to document how the
public input datasets were derived from the development files. It is
primarily for the authors and requires access to the original ENDURE
development archive.

Readers reproducing the published analysis do **not** need to run
Notebook 00 because the resulting publication input files are included
in the repository.

## Software

The deterioration model is implemented in Python. Calibration and
publication graphics use R.

Python dependencies are listed in `requirements.txt`.

Principal R packages include:

``` text
mgcv
ggplot2
terra
tidyterra
sf
rnaturalearth
rnaturalearthdata
data.table
cowplot
```

The analysis notebooks use repository-relative paths and are intended to
run without modification on Linux, macOS or Windows once the required
Python and R environments are available.

## Data availability

The repository includes the minimal derived environmental dataset and anonymised archaeological calibration data needed to reproduce the analyses reported in the manuscript.

Exact archaeological site coordinates are not included because they derive from restricted heritage records and are not required for reproducibility of the reported analyses.

Data provenance and reuse conditions are described in `DATA_LICENSE.md`.

## Archived release

The publication reproducibility release (`v1.0.0`) is permanently archived in Zenodo:

**DOI:** <https://doi.org/10.5281/zenodo.23159324>

The GitHub repository may continue to be updated. The Zenodo `v1.0.0` archive provides the fixed version corresponding to the manuscript.

## Licence and citation

Original software code in this repository is released under the MIT License. Data included in the repository retain the provenance and reuse conditions described in `DATA_LICENSE.md`.

Citation metadata are provided in `CITATION.cff`. Please cite the archived Zenodo release when using the repository:

**Syms, C., McGonigle, C., Quinn, R. & Gregory, D. (2026).**
_Environmental forcing and deterioration trajectories in submerged shipwreck heritage: code and data_. Version 1.0.0. Zenodo.
<https://doi.org/10.5281/zenodo.23159324>

The article DOI will be added when available.
