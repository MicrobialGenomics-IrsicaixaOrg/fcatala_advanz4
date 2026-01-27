# Advanz4 

[![DOI](https://zenodo.org/badge/955869511.svg)](https://doi.org/10.5281/zenodo.18387903)

## Introduction 🦠

This repository contains the code and documentation for the analyses presented in the paper:

**Dolutegravir Restores Gut Microbiota in Late-stage HIV-1 Unlike Darunavir: An Open-Label, Randomized Clinical Trial**

(*ClicialTrials.gov; NCT02337322*)

## Study Summary 📋

Late HIV-1 presentation is associated with impaired immune reconstitution and gut microbiota dysbiosis. In this substudy of a multicentre, open-label and randomized clinical trial, 88 antiretroviral-naïve individuals with very advanced HIV were randomized to receive lamivudine/abacavir combined with either dolutegravir (DTG) or ritonavir-boosted darunavir (DRV/r) and followed for 2 years.

This repository provides all the code used to analyze the clinical and microbiome data generated in the study:

- Alpha diversity analyses

- Beta diversity analyses

- Differential abundance testing

- Correlation analyses with inflammation markers

- Functional abundance profiling

- Network analyses

Please refer to the article manuscript for the complete description of the study design, methods, results and conclusions.

## Repository Structure 📁

* `analysis/`: contains Quarto (`.qmd`) files with the full analysis worfkflows.

* `docs/`: contains rendered HTML reports from the Quarto files.

* `data/`: placeholder for datasets required by the analysis (not included due to privacy reasons).

* `.github/`: CD/CI workflows and repository configuration.

## Installation & Environment Setup 🔧

Analyses were conducted in **R(4.5.1)** with dependencies managed via `renv`. To reproduce the analyses, please follow the instructions below.

1. Clone this repository:

```bash
git clone https://github.com/MicrobialGenomics-IrsicaixaOrg/fcatala_advanz4.git
cd fcatala_advanz4
```

2. Open an R session and install `renv` (if not already installed).

```R
install.packages("renv")
```

3. Restore the project environment:

```R
renv::restore()
``` 

This will install all required packages in the versions used for the analyses at a local project library.

## Data Availability 🔐

Due to privacy reasons, raw sequencing data and clinical metadata required as input files are not included in this repository. Please contact the corresponding author for data access requests.

Once obtained, place the required files in the `data/raw/` directory following the structure indicated above. Additional intermediate files are generated during the analysis execution and stored in `data/processed/`.

## Reproducing the Analyses ⚙️

After setting up the environment and downloading the required data:

1. Open an R session in the project directory.

2. Render all Quarto analysis files. For example:

    ```bash
    quarto render analysis/01_alpha_diversity.qmd
    ```
    or, from R:

   ```R
   quarto::quarto_render("analysis/01_alpha_diversity.qmd")
   ```
   
> [!IMPORTANT]
> Analyses files should be rendered separately and in numerical order (e.g. 01, 02, 03...) to avoid data dependency issues.

3. Output files
    
    * Rendered HTML reports are saved in the `docs/`.

    * Intermediate data files are saved in `data/processed/`.

    * Generated figures are stored in `data/` under the corresponding analysis subdirectory (e.g., results from `02_alpha_diversity.qmd` will be saved in `data/02_alpha_diversity/`).


## Citation 📜

If you use this code or analyses, please cite the original paper:

> F. Català et al. *Dolutegravir Restores Gut Microbiota in Late-stage HIV-1 Unlike Darunavir: An Open-Label, Randomized Clinical Trial*. [Journal, Year, DOI - to be updated]

## Contact 📧

For inquiries about the analyses or data requests, please contact:

- **Analysis Lead**: Francesc Català-Moll

- **Email**: fcatala@irsicaixa.es

- **Affiliation**: Microbial Genomics Group -  IrsiCaixa AIDS Research Institute, Badalona, Spain.
