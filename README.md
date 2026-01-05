# Lipid Consensus WLCNA Pipeline: Module–Trait (Weight / PPWR / APO) Pipeline

This pipeline builds a **consensus WGCNA network** from lipidomics data across
multiple timepoints and relates lipid modules to sample-level traits
(e.g., **weight**, **PPWR**, **APO**).  
It is designed to be **reusable**: other users can plug in *their own* traits via
a simple YAML configuration file.

---

## Folder structure

```text
project_root/
  R/
    00_setup_packages.R
    01_prepare_traits_demographic_data_bmi.R
    02_prepare_lipid_expression_multiExpr_ID.R
    03_run_consensus_WGCNA_ID.R
    04_setup_pheno_ID_pre_traits.R
    05_module_trait_correlations_traits.R
    03_utils_module_trait_timepoints_plots.R
  scripts/
    06_plot_module_trait_timepoints.R
    run_all.R
  config/
    trait_config.yml
  data_raw/
    ID_Updated_Combined_Traits_apo_analysis.csv
    identified_lipids_baseline_data_normalized.csv
    identified_lipids_TP36_data_normalized.csv
    identified_postpartum_data_normalized.csv
  data_processed/
    (generated intermediate files)
  results/
    tables/
    plots/

For the GROWELL lipidomics paper, please check this folder for the script that was used here **scripts_for_growell_lipidomics_analysis**
