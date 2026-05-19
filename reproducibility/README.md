# Reproducibility

This folder contains everything needed to recreate the exact R environment used for this analysis.

## Setup

**Requirements:** pixi, Mac arm64 (M1/M2/M3), Xcode command line tools.

**One-time C++ compiler fix** (required for Bioconductor packages):
```bash
mkdir -p ~/.R
cat > ~/.R/Makevars << 'MAKEVARS'
CXX = clang++
CXX11 = clang++
CXX14 = clang++
CXX17 = clang++
CXX11STD = -std=c++14
CXX14STD = -std=c++14
CXX17STD = -std=c++17
CXXFLAGS = -std=c++14 -O2 -fPIC
MAKEVARS
```

**Install the environment:**
```bash
cd reproducibility
pixi install
pixi run install-r-pkgs
pixi run session-info
```

## Files

- `pixi.toml` — R version and all conda-installable packages
- `pixi.lock` — exact locked versions (auto-generated, do not edit)
- `install_r_packages.R` — installs Bioconductor + MetaboAnalystR
- `session_info.txt` — exact package versions used at analysis time

## Running the analysis

```bash
pixi run run-all
```
