#!/bin/bash
# -------------------------------------------------------
# run_pipeline.sh
# Full GROWELL lipidomics pipeline execution
# Run from repo root:  bash run_pipeline.sh
# -------------------------------------------------------

set -e  # stop on first error
cd /Users/gekuodza/Documents/GitHub/lipid-consensus-wlcna-pipeline

PIXI="pixi run --manifest-path reproducibility/pixi.toml Rscript"

echo "======================================================"
echo " GROWELL Lipidomics Pipeline"
echo "======================================================"

# ── Infrastructure (sourced by other scripts, not run directly) ──
# 00_load_packages.R
# 00_figure_theme.R

# # ── Data preparation ──────────────────────────────────────
# echo "[01/13] Preparing raw lipidomics data..."
# $PIXI R/01_prepare_raw_data.R

# echo "[02/13] Preparing trait data..."
# $PIXI R/02_prepare_traits.R

# echo "[03/13] Normalizing lipids (MetaboAnalystR)..."
# $PIXI R/03_normalize_lipids.R

# echo "[04/13] Building multiExpr object for WGCNA..."
# $PIXI R/04_prepare_lipids.R

# # # ── Network analysis ──────────────────────────────────────
# echo "[05/13] Running consensus WLCNA..."
# $PIXI R/05_run_consensus_wgcna.R

# echo "[06/13] Computing module membership (kME)..."
# $PIXI R/06_module_membership.R

# ── Annotation (must run before 08 to generate module order) ──
echo "[07/13] Annotating lipids (RefMet/PubChem)..."
$PIXI R/07_annotate_lipids.R

# ── Statistical analyses ──────────────────────────────────
# echo "[08/13] Module-trait correlations..."
# $PIXI R/08_module_trait_correlations.R

echo "[09/13] Module eigennode plots..."
$PIXI R/09_module_eigennode_plots.R

# echo "[10/13] Differential lipid analysis (APO + PPWR)..."
# $PIXI R/10_differential_lipids.R

# echo "[11/13] Lipid boxplots (selected lipids of interest)..."
# $PIXI R/11_lipid_boxplots.R

# echo "[12/13a] PLS-DA permutation testing..."
# $PIXI R/12_plsda.R

# echo "[12/13b] PLS-DA permutation testing..."
# $PIXI R/12_plsda_metaboanalyst.R

# echo "[13/13] Sensitivity analyses..."
# $PIXI R/13_sensitivity_analyses.R

# ── Figure assembly ───────────────────────────────────────
echo "[Final] Assembling publication figures..."
$PIXI R/99_assemble_figures.R

echo "======================================================"
echo " Pipeline complete!"
echo " Results: results/"
echo " Figures: results/figures/final/"
echo "======================================================"