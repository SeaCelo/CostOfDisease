This branch is the `OG-Core 0.15.5 / OG-ZAF 0.0.11` version of the Cost of Disease paper run.

Relative to the committed `OG-Core 0.14.3 / OG-ZAF 0.0.8` GBD-based branch, it differs in three ways:

- the environment is pinned to newer OG-Core and OG-ZAF GitHub commits
- the branch is intended to rerun the same GBD-based mortality mapping under those versions
- `main.py` uses `p.RC_TPI = 0.0075` so the reform scenarios complete under this version

Pinned package sources in `code/environment.yml`:

- OG-Core: `bfe4c5bef220dee8538508b5175c2759e497e903`
- OG-ZAF: `819fad9349d9c4643675fb0450e0cc653a7d331a`

Checked-in mortality input:

- `source/JDE/hiv-data/IHME-GBD_2023_DATA-ddf37f70-1/IHME-GBD_2023_DATA-ddf37f70-1.csv`
- `code/demographic_data/hiv_mortality_profile_gbd_sa_2023.csv`
- `code/build_hiv_mortality_profile.py` rebuilds the frozen profile from the raw GBD CSV if needed

Environment build:

```bash
conda env create -n cod-paper -f code/environment-og15-zaf11.lock.yml
conda activate cod-paper
```

Then run:

```bash
cd code
python main.py
```

Current manuscript/output status:

- `source/JDE/CostOfDisease_JDE_bins.tex` and the manuscript-facing GDP tables in `source/JDE/tables_figures_bins/` now reflect the current `OG-Core 0.15.5 / OG-ZAF 0.0.11` outputs
- the remaining figure assets in `source/JDE/tables_figures_bins/` are still carried forward from the committed `OG-Core 0.14.3 / OG-ZAF 0.0.8` GBD-based branch until they are regenerated on this branch

Caveat:

- baseline passes the older `RC_TPI = 0.001` standard
- the current stack still needs a looser `RC_TPI` than the `OG-Core 0.14.3 / OG-ZAF 0.0.8` branch
- the current rerun solved the High scenario at `p.RC_TPI = 0.0075`
