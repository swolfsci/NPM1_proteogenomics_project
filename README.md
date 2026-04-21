# NPM1 Proteogenomics Project

Analysis code for: **"A deregulated mitochondrial proteome defines a high-risk subtype of NPM1-mutant AML"**

## Data availability and scope of reproducibility

The AMLSG 09-09 trial cohort comprises individual-level clinical and genetic data
that are not publicly shareable under the trial's data-use agreement. This scopes
what can be reproduced from this repository alone:

**Reproducible from public data** (processed proteomics in Supplementary Dataset 1
of the manuscript, plus the reference resources listed below):

- Differentiation-state analyses (diffusion maps, metaclustering, GSVA scores)
- Mito-score derivation and ssGSEA scoring
- Comparative GSEA and pathway analyses
- All figure panels that do not depend on clinical outcomes

**Requires access to clinical and genetic data from the corresponding authors:**

- Survival analyses (Kaplan-Meier, Cox regression, cumulative incidence)
- Mutation-cluster association tests
- Multivariable models adjusting for clinical covariates

To support readers adapting this pipeline to other cohorts, the repository includes
`data/clinical_data.csv.example` — a schema file with column definitions and dummy
rows indicating the expected input format for all outcome and genetic analyses.
Running the survival scripts on this file will execute end-to-end without errors
but produce results that are not biologically meaningful.

## Data

Download the data folder from [figshare DOI: XXXX] and place its contents into `data/`.

The `data/` directory should contain:

| File | Description | Availability |
|------|-------------|--------------|
| `protein_intensities_raw.csv` | Raw protein-level intensities (10,507 proteins × 330 samples) | Public |
| `protein_matrix_vsn_imputed.csv` | Analysis-ready matrix after QC, VSN, median centering, DreamAI imputation (7,223 × 284) | Public (Supp. Dataset 1) |
| `cluster_mapping.csv` | Locked proteomic cluster assignments and diffusion map coordinates | Public |
| `MitovsAll.csv` | Differential expression results from Jayavelu, Wolf, Buettner et al. *Cancer Cell* 2022 | Public |
| `Human.MitoCarta3.0.xls` | MitoCarta 3.0 database | Public (external resource) |
| `AMLCellType_Genesets.gmt` | AML cell type signatures (Zeng et al.) | Public (external resource) |
| `Furtwaengler_signatures.csv` | Single-cell proteomics signatures (Furtwängler et al.) | Public (external resource) |
| `sample_metadata.csv` | Sample-to-BioID mapping and MS run metadata | Public |
| `clinical_data.csv` | De-identified clinical and NGS data | **Upon request** (see above) |
| `clinical_data.csv.example` | Schema file with dummy rows for testing | Public (included in repo) |

Raw mass spectrometry data are deposited at PRIDE (accession: XXXX).

## Reproducing the analysis

```r
# 1. Install dependencies
renv::restore()

# 2. Run scripts in order
source("R/00_preprocessing.R")       # QC, normalization (or start from provided matrix)
source("R/01_differentiation.R")     # Figures 1, S1
source("R/02_mutations.R")           # Figures 2, S2 (requires clinical_data.csv)
source("R/03_mito_score.R")          # Figures 3, S3 (survival analyses require clinical_data.csv)
source("R/04_robustness.R")          # Table S2, robustness analyses
```

Scripts read from `data/` and write to `output/figures/` and `output/tables/`.

Each script sources `R/utils.R` for shared helper functions and explicitly loads its
required data at the top. Scripts 01–03 save intermediate R objects to `output/`
that downstream scripts read. Script 04 is self-contained.

If `clinical_data.csv` is absent, scripts 02 and 03 fall back to `clinical_data.csv.example`
and print a warning; figures that depend on outcome data will be generated with placeholder
results.

## Software

Analysis was performed in R 4.4.2. All package versions are recorded in `renv.lock`.

## Contact

For access to individual-level clinical and genetic data, contact the corresponding
authors of the manuscript.