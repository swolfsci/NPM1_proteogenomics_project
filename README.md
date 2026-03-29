# NPM1 Proteogenomics Project

Analysis code for: **"A deregulated mitochondrial proteome defines a high-risk subtype of NPM1-mutant AML"**

## Data

Download the data folder from [figshare DOI: XXXX] and place its contents into `data/`.

The `data/` directory should contain:

| File | Description |
|------|-------------|
| `protein_intensities_raw.csv` | Raw protein-level intensities (10,507 proteins x 330 samples) |
| `protein_matrix_vsn_imputed.csv` | Analysis-ready matrix after QC, VSN, median centering, DreamAI imputation (7,223 x 284) |
| `clinical_data.csv` | De-identified clinical and NGS data |
| `sample_metadata.csv` | Sample-to-BioID mapping and MS run metadata |
| `cluster_mapping.csv` | Locked proteomic cluster assignments and diffusion map coordinates |
| `MitovsAll.csv` | Differential expression results from Jayavelu, Wolf, Buettner et al. Cancer Cell 2022 |
| `Human.MitoCarta3.0.xls` | MitoCarta 3.0 database |
| `AMLCellType_Genesets.gmt` | AML cell type signatures (Zeng et al.) |
| `Furtwaengler_signatures.csv` | Single-cell proteomics signatures (Furtwaengler et al.) |

Raw mass spectrometry data are deposited at PRIDE (accession: XXXX).

## Reproducing the analysis

```r
# 1. Install dependencies
renv::restore()

# 2. Run scripts in order
source("R/00_preprocessing.R")       # QC, normalization (or start from provided matrix)
source("R/01_differentiation.R")     # Figures 1, S1
source("R/02_mutations.R")           # Figures 2, S2
source("R/03_mito_score.R")          # Figures 3, S3
source("R/04_robustness.R")          # Table S2, robustness analyses
```

Scripts read from `data/` and write to `output/figures/` and `output/tables/`.

Each script sources `R/utils.R` for shared helper functions and explicitly loads its required data at the top. Scripts 01-03 save intermediate R objects to `output/` that downstream scripts read. Script 04 is self-contained.

## Software

Analysis was performed in R 4.4.2. All package versions are recorded in `renv.lock`.
