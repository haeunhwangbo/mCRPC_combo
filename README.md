# mCRPC_combo

Author: Haeun (Hannah) Hwangbo

This repository contains source codes for collaborations between the Palmer Lab and Huang Lab.

## Getting Started

Installing all required packages. You would need conda to access the virtual environment.

```bash
git clone https://github.com/haeunhwangbo/mCRPC_combo.git
conda env create -f env/environment_short.yml

conda activate surv
```

Dependencies: numpy, scipy, pandas, scikit-learn, lifelines, matplotlib, seaborn, snakemake

## Reconstructing the analysis

Download all data and place them into a directory of your choice. The original data directory was organized as

```bash
data
├── PFS
├── experimental
├── rPFS
├── raw
│   ├── experimental
│   └── trials
└── waterfall

```

All data files should be placed in their approprite directories as described in `config.yaml` file.

You can reconstruct all tables and figures in the article by running the following code. We recommend using at least 4 cores because the code utilizes parallel computing. If you have limitied computation power, reducing the `NRUN` in the `src/predictive_power.py` to 100 or 1000 will significantly reduce the run time. 

```bash
# specify number of cores {N} you want to use
snakemake --cores {N} all
```

