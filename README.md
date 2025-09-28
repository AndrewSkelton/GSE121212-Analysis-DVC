# Analysis of GSE121212 RNA-seq Data

This repository contains exploratory analyses of the public RNA sequencing dataset GSE121212. It includes processed count matrices, reference metadata, notebooks, and workflow scripts used to reproduce figures and summary statistics.

## Repository Structure

- `Data/` – Alignment outputs, quantification results, and log files.
- `Notebooks/` – Reproducible analysis notebooks detailing QC and differential expression steps.
- `Pheno/` – Sample- and patient-level phenotype metadata.
- `Reference/` – Reference sequences, annotations, and indices used for alignment and quantification.
- `Workstreams/` – Workflow definitions and scripts for preprocessing and downstream analysis.

## Large Files

Git LFS is enabled to store large binary artefacts such as the RDS snapshot of the project state. If cloning this repository, ensure Git LFS is installed (`git lfs install`) so tracked data are fetched automatically.

## Environment Setup

The analysis is managed with [`renv`](https://rstudio.github.io/renv/). To create the project library, run:

```r
renv::restore()
```

This installs the package versions recorded in `renv.lock`. When adding new package dependencies, install them with `renv::install()` and update the lockfile via `renv::snapshot()` so collaborators can stay in sync.

## Getting Started

1. Clone the repository: `git clone <repo-url>`
2. Install Git LFS (if not already): `git lfs install`
3. Run `renv::restore()` to populate the project library.
4. Open the notebooks in RStudio or Jupyter to reproduce the analysis.

## Reproducibility

Environment management and additional workflow details are documented within the notebooks. Consider using renv or conda to capture library versions when expanding the analysis.
