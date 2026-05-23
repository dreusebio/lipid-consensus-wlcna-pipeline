cd /Users/gekuodza/Documents/GitHub/lipid-consensus-wlcna-pipeline

pixi run --manifest-path reproducibility/pixi.toml Rscript R/00_load_packages.R

pixi run --manifest-path reproducibility/pixi.toml Rscript R/01_prepare_traits.R

pixi run --manifest-path reproducibility/pixi.toml Rscript R/00b_prepare_raw_data.R

pixi run --manifest-path reproducibility/pixi.toml Rscript R/02_normalize_lipids.R

pixi run --manifest-path reproducibility/pixi.toml Rscript R/02b_prepare_lipids.R

pixi run --manifest-path reproducibility/pixi.toml Rscript R/03_run_consensus_wgcna.R

pixi run --manifest-path reproducibility/pixi.toml Rscript R/04_module_membership.R

pixi run --manifest-path reproducibility/pixi.toml Rscript R/04b_annotate_lipids.R

pixi run --manifest-path reproducibility/pixi.toml Rscript R/05_module_trait_correlations.R

pixi run --manifest-path reproducibility/pixi.toml Rscript R/06_plot_module_trait_timepoints.R

pixi run --manifest-path reproducibility/pixi.toml Rscript R/07_differential_lipids.R

pixi run --manifest-path reproducibility/pixi.toml Rscript R/08_fdr_correction.R

pixi run --manifest-path reproducibility/pixi.toml Rscript R/09_plsda_permutation.R

pixi run --manifest-path reproducibility/pixi.toml Rscript R/10_sensitivity_analyses.R