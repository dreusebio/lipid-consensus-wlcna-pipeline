# Lipid Consensus Weight Lipid Correlation Network Analysis (WLCNA) Pipeline: Module–Trait (Weight / PPWR / APO) Pipeline

This pipeline builds a **consensus WLCNA network** from lipidomics data across
multiple timepoints and relates lipid modules to sample-level traits
(e.g., **weight**, **PPWR**, **APO**).  
It is designed to be **reusable**: other users can plug in *their own* traits via
a simple YAML configuration file.

---

## Folder structure

```text
lipid-consensus-wlcna-pipeline/
│
├── config/
│   └── config.yml                    ← all paths & parameters
│
├── R/                                ← all analysis scripts
│   ├── 00_load_packages.R
│   ├── 01_prepare_traits.R
│   ├── 02a_normalize_lipids.R
│   ├── 02_prepare_lipids.R
│   ├── 03_run_consensus_wgcna.R
│   ├── 04_module_membership.R
│   ├── 05_module_trait_correlations.R
│   ├── 06_plot_module_trait_timepoints.R
│   ├── 07_qc_normalization.R
│   ├── 08_fdr_correction.R
│   ├── 09_plsda_permutation.R
│   └── 10_sensitivity_analyses.R
│
├── scripts/
│   └── run_all.R                     ← master orchestration only
│
├── reproducibility/
│   ├── pixi.toml
│   ├── pixi.lock
│   ├── install_r_packages.R
│   ├── session_info.txt
│   └── README.md
│
├── data_raw/                         ← in .gitignore (patient data)
│   ├── raw/                          ← unnormalized lipid CSVs
│   ├── normalized/                   ← MetaboAnalyst output CSVs
│   └── growell_curated_deid_george.csv
│
├── data_processed/                   ← in .gitignore (generated .rds files)
│
├── results/                          ← in .gitignore (generated outputs)
│   ├── 01_traits/
│   ├── 02a_normalization/
│   ├── 02_lipids/
│   ├── 03_consensus_wgcna/
│   ├── 04_module_membership/
│   ├── 05_module_trait_cor/
│   ├── 06_module_trait_plots/
│   ├── 07_differential/
│   ├── 08_fdr_correction/
│   ├── 09_plsda/
│   └── 10_sensitivity/
│
├── .gitignore
└── README.md
```
For the GROWELL lipidomics paper, please check this folder for the script that was used here **scripts_for_growell_lipidomics_analysis**
