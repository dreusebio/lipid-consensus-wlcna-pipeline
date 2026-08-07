# Lipid Consensus WLCNA Pipeline

Consensus Weighted Lipid Correlation Network Analysis (WLCNA) of gestational lipidomics data, relating lipid co-expression modules to adverse pregnancy outcomes (APOs) and postpartum weight retention (PPWR).

---

## Project summary

**Data type:** Targeted/untargeted lipidomics, dried blood spot (451 identified lipids)

**Study / cohort:** GROWELL (NCT02904473)

**Timepoints:** Baseline (~10–16 weeks), 36–38 weeks gestation, postpartum

**Sample size:** n = 49 participants

**Analytic framework:** Consensus  WLCNA (module preservation across timepoints), differential lipid analysis, PLS-DA (`ropls` for statistics, MetaboAnalystR for biplots)

**Primary outputs:** Consensus lipid correlation networks/modules, module–trait correlations, differential lipid results, PLS-DA models, Table 1, publication-ready figures

**Associated publication:** Kuodza, G.E., Keeton, V.F., Williams, L.A. *et al.* Exploring cardiometabolic markers in adverse pregnancy outcomes: insights from the GROWell study. *npj Women's Health* (2026). [doi.org/10.1038/s44294-026-00158-3](https://doi.org/10.1038/s44294-026-00158-3)

---

## Project overview

This repository implements the full downstream analysis for the GROWELL Lipidomics Study, a secondary analysis nested within a randomized clinical trial (NCT02904473). Plasma lipid profiles from dried blood spots were measured at three gestational timepoints and related to two families of outcomes:

- **Adverse pregnancy outcomes (APOs):** composite APO, hypertensive disorders of pregnancy (HDP), gestational diabetes (GDM), preterm birth, and other APOs
- **Postpartum weight retention (PPWR)**

The pipeline builds a **consensus WLCNA network** across the three timepoints, evaluates module preservation between timepoints, tests module eigengene–trait associations, runs a differential lipid analysis, fits PLS-DA models, and generates the manuscript's figures and Table 1.

**Central hypothesis:** The purpose of this study was to build on existing evidence linking lipidomic profiles to APOs and PPWR at 6 months, and to identify groups of lipids early in pregnancy that could be associated with APOs in a remote study of women with pre-pregnancy overweight or obesity.

Using WLCNA, the analysis identified triglyceride (TG)-rich lipid signatures associated with APOs and PPWR: some in early-pregnancy TG networks and several individual TGs consistently associated with APOs, and some TG networks at 36–38 weeks and 3 months postpartum associated with PPWR at 6 months.

This repository is built primarily around the GROWELL analysis and manuscript; the config-driven design (`config/config.yml`) makes reuse with other trait sets possible, but reproducing the GROWELL results is the primary use case.

---

## Study design

- **Design:** Longitudinal, nested within a randomized clinical trial (NCT02904473)
- **Sample size:** 49 participants
- **Biological sample:** Dried blood spot (whole blood, fingerstick)
- **Timepoints:** Baseline (~10–16 wk), 36–38 wk (Wk36_38), postpartum
- **Features:** 451 identified lipids

**Outcomes:**

- APO traits: composite, HDP, GDM, preterm, other
- PPWR (continuous and/or dichotomized: excessive vs. not)

**Statistical framework:**

- Consensus WGCNA/WLCNA for module detection; module preservation (Zsummary, median rank) assessed pairwise across timepoints
- Spearman/Pearson/bicor module–trait correlation
- Differential lipid analysis with nominal (p < 0.05) and FDR (BH q = 0.05, q = 0.10) significance tiers
- PLS-DA via `ropls` (Q², R²X/R²Y, permutation testing); MetaboAnalystR used for biplot/score-plot figures only
- Sensitivity analyses: treatment-arm-adjusted and covariate-adjusted (age, race, BMI group) module–trait associations

---

## Analysis workflow

Scripts in `R/` run in numbered order; `scripts/run_all.R` (or `run_all.sh`) executes the full pipeline end-to-end.

### 0. Setup

- `00_load_packages.R` — loads/attaches required packages
- `00_figure_theme.R` — shared ggplot2 theme for manuscript figures

### 1–4. Data preparation

- `01_prepare_raw_data.R` — imports and QCs raw lipidomics data (identified/unidentified, by ESI mode and timepoint)
- `02_prepare_traits.R` — prepares demographic and APO/PPWR trait data
- `03_normalize_lipids.R` — sample and row normalization per timepoint (QC plots: density, PCA scores)
- `04_prepare_lipids.R` — builds the multi-timepoint lipid expression object (`multiExpr_ID.rds`) for consensus WGCNA

### 5. Consensus network construction

- `05_run_consensus_wgcna.R` — soft-thresholding power selection and consensus WLCNA network/module construction across timepoints
- `05b_module_preservation.R` — pairwise module preservation statistics across timepoints (Zsummary, median rank; Figure S6)

### 6–9. Module characterization and association

- `06_module_membership.R` — module membership and hub lipid identification per timepoint
- `07_annotate_lipids.R` — lipid class annotation and ESI mode summaries
- `08_module_trait_correlations.R` — module eigengene–trait correlation testing (APO and PPWR trait families, BH-corrected separately)
- `09_module_eigennode_plots.R` — eigennode dendrograms and network plots per timepoint

### 10–12. Differential lipids and classification

- `10_differential_lipids.R` — differential lipid analysis (raw p < 0.05, BH q = 0.05/0.10)
- `11_lipid_boxplots.R` — boxplots for individual significant lipids
- `12_plsda.R` — PLS-DA statistics via `ropls` (Q², permutation p-values)
- `12_plsda_metaboanalyst.R` — MetaboAnalystR biplots/score plots for PLS-DA figures

### 13–14. Sensitivity analyses

- `13_sensitivity_analyses.R` — treatment-arm-adjusted vs. unadjusted module–trait associations
- `14_sensitivity_analyses_covariates.R` — covariate-adjusted (age, race, BMI group) module–trait associations

### 15. Tables and figures

- `15_table1.R` — generates manuscript Table 1 (Word + Excel)
- `99_assemble_figures.R` — assembles final manuscript figures from individual panels
- `diagnose_mp.R` — diagnostic utility for troubleshooting module preservation runs

---

## Repository structure

```
lipid-consensus-wlcna-pipeline/
├── R/
│   ├── 00_load_packages.R
│   ├── 00_figure_theme.R
│   ├── 01_prepare_raw_data.R
│   ├── 02_prepare_traits.R
│   ├── 03_normalize_lipids.R
│   ├── 04_prepare_lipids.R
│   ├── 05_run_consensus_wgcna.R
│   ├── 05b_module_preservation.R
│   ├── 06_module_membership.R
│   ├── 07_annotate_lipids.R
│   ├── 08_module_trait_correlations.R
│   ├── 09_module_eigennode_plots.R
│   ├── 10_differential_lipids.R
│   ├── 11_lipid_boxplots.R
│   ├── 12_plsda.R
│   ├── 12_plsda_metaboanalyst.R
│   ├── 13_sensitivity_analyses.R
│   ├── 14_sensitivity_analyses_covariates.R
│   ├── 15_table1.R
│   ├── 99_assemble_figures.R
│   └── diagnose_mp.R
├── scripts/
│   ├── run_all.R
│   └── run_all.sh
├── config/
│   └── config.yml                     # study parameters, timepoints, trait definitions, WGCNA settings
├── reproducibility/
│   ├── pixi.toml
│   ├── pixi.lock
│   ├── install_r_packages.R
│   ├── session_info.txt
│   └── README.md
├── data/
│   ├── raw/                            # raw lipidomics + trait data (also archived at Metabolomics Workbench, ST005032)
│   └── processed/                      # intermediate .rds objects generated by the pipeline
├── results/                            # full analysis outputs (all timepoints, n = 49)
│   ├── 01_raw_data/
│   ├── 02_traits/
│   ├── 03_normalization/
│   ├── 04_lipids/
│   ├── 05_consensus_wgcna/
│   │   └── power14_split4_merge010/    # module preservation outputs (Figure S6)
│   ├── 06_module_membership/
│   ├── 07_annotations/
│   ├── 08_module_trait_cor/
│   ├── 09_module_eigennode_plots/
│   ├── 10_differential/
│   ├── 11_lipid_boxplots/
│   ├── 12_plsda/
│   ├── 12_plsda_metaboanalyst/
│   ├── 13_sensitivity/
│   ├── 14_sensitivity_covariates/
│   ├── 15_table1/
│   └── figures/
│       ├── final/                      # assembled manuscript figures (PDF/PNG/TIFF)
│       ├── panels/                     # individual figure panels 
│       ├── supplement/
│       └── Main_Figures_and_Tables/, Supplement/  # figures/tables as submitted to the journal
├── results_48_sample/                  # parallel run on the earlier 48-sample subset as an additional sensitivity analysis 
├── LICENSE
└── lipid-consensus-wlcna-pipeline.Rproj
```

---

## Quick start (reproducible environment)

This project uses **[Pixi](https://pixi.sh/latest/installation/)** for reproducible environments. If you don't have Pixi installed, follow the [official installation instructions](https://pixi.sh/latest/installation/) before proceeding.

### 1. Clone repository

```
git clone https://github.com/dreusebio/lipid-consensus-wlcna-pipeline.git
cd lipid-consensus-wlcna-pipeline
```

### 2. Configure environment

```
export PIXI_HOME=/path/to/your/pixi_home
export PIXI_CACHE_DIR=/scratch/$USER/pixi-cache
unset RATTLER_CACHE_DIR
unset XDG_CACHE_HOME
```

#### (cluster users, e.g. UC Davis Hive)

```
export PIXI_HOME=/quobyte/lasallegrp/George/.pixi
export PIXI_CACHE_DIR=/tmp/$USER/pixi-cache
unset RATTLER_CACHE_DIR
unset XDG_CACHE_HOME
```

### 3. Install environment

```
cd reproducibility
pixi install --concurrent-downloads 1 --concurrent-solves 1
```

### 4. Install R packages

```
pixi run Rscript reproducibility/install_r_packages.R
```

### 5. Verify installation

```
pixi run Rscript -e 'sessionInfo()'
```

Compare against `reproducibility/session_info.txt` for the expected package versions.

---

## Running the analysis

Run the full pipeline end-to-end:

```
pixi run Rscript scripts/run_all.R
```

or

```
bash scripts/run_all.sh
```

Or run individual steps interactively from `R/` in numbered order (`00` → `99`), configuring study parameters, timepoints, and trait definitions in `config/config.yml`.

---

## Reproducibility

This project is reproducible using:

```
reproducibility/
  pixi.toml
  pixi.lock
  install_r_packages.R
  session_info.txt
```

**Key features:**

- `pixi.lock` → guarantees identical environments across machines/clusters
- `install_r_packages.R` → handles CRAN + Bioconductor + GitHub dependencies (WGCNA, ropls, MetaboAnalystR, etc.)
- Config-driven design (`config/config.yml`): trait definitions, timepoint labels, WGCNA parameters, and file paths are specified in one place, not hard-coded
- BH correction applied separately within each trait family (APO, PPWR) to avoid inflating significance across unrelated hypothesis sets
- Soft-thresholding power selected as the minimum consensus power at which all timepoints achieve R² ≥ 0.80, not tuned post hoc to improve results
- PLS-DA statistics (Q², permutation p-values) computed via `ropls`, not `mixOmics`/MetaboAnalystR, for valid cross-validated performance metrics; MetaboAnalystR used only for biplot figures

---

## Platform support

| Platform         | Support      |
| ---------------- | ------------ |
| Linux (cluster)  | ✅ full       |
| macOS            | ✅ full       |
| Windows          | ⚠️ untested  |

---

## Troubleshooting

### Pixi install fails (network/cache issues)

```
export PIXI_CACHE_DIR=/tmp/$USER/pixi-cache
```

### Reinstall environment

```
rm -rf /tmp/$USER/pixi-cache
cd reproducibility
pixi install --concurrent-downloads 1 --concurrent-solves 1
```

### Missing R/Bioconductor packages

```
pixi run Rscript reproducibility/install_r_packages.R
```

---

## Data availability

The processed lipidomics data for the GROWELL study are also deposited at the **Metabolomics Workbench** (Study ID **ST005032**), the canonical, citable archive for this dataset:

[https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST005032](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST005032)

For convenience and full end-to-end reproducibility, this repository also directly includes `data/raw/`, `data/processed/`, and the complete `results/` tree (all pipeline outputs, at both the full n = 49 sample and the earlier 48-sample subset). Cloning the repository is sufficient to inspect or rerun the analysis without a separate data download.

Manuscript figures and tables as submitted to *npj Women's Health* are in `results/figures/Main_Figures_and_Tables/` and `results/figures/Supplement/`.

---

## Citation

If using this pipeline, please cite:

> Kuodza, G.E., Keeton, V.F., Williams, L.A. *et al.* Exploring cardiometabolic markers in adverse pregnancy outcomes: insights from the GROWell study. *npj Women's Health* (2026). https://doi.org/10.1038/s44294-026-00158-3

---

## Contact

George Eusebio Kuodza
LaSalle Lab, UC Davis
<gekuodza@health.ucdavis.edu>

---


---

## License

See `LICENSE`
