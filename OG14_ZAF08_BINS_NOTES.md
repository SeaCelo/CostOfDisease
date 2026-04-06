This branch is the `OG-Core 0.14.3 / OG-ZAF 0.0.8` version of the Cost of Disease paper run with a South Africa GBD-based HIV age-mortality profile.

Relative to the pre-binning `OG-Core 0.14.3 / OG-ZAF 0.0.8` version, it adds:

- South Africa 2023 HIV/AIDS mortality-rate data from the GBD Results Tool
- a frozen smooth age-specific HIV mortality profile constructed from the non-overlapping GBD age groups
- one calibrated mortality-shock scalar per scenario so realized year-5 excess deaths match the scenario totals
- the same environment pinning, low-scenario fix, reform demographic-path fix, and generated table cleanup already in the pre-binning branch
- the exact GBD CSV used to build the profile is checked into `source/JDE/hiv-data/IHME-GBD_2023_DATA-ddf37f70-1/`
- the runtime profile used by the paper run is checked into `code/demographic_data/hiv_mortality_profile_gbd_sa_2023.csv`
- `code/build_hiv_mortality_profile.py` rebuilds the frozen profile from the raw GBD CSV if we ever need to refresh it

Current manuscript/output status:

- `source/JDE/tables_figures_bins/mortality_rates.png` reflects the GBD-based method
- the refreshed results from the latest rerun are:
  - 20-year average GDP losses: `-10.71`, `-17.62`, `-25.49`
  - 2\% NPV losses: `-2,976.15`, `-4,755.79`, `-6,767.19`
- a point to flag for discussion is that, relative to the cleaned-up proportional-scaling reference version, this GBD-based mapping produces larger average GDP losses in the first 20 years and larger discounted GDP losses over the full 100-year horizon
- `source/JDE/CostOfDisease_JDE_bins.tex` describes the GBD-based age mapping rather than the old proportional-scaling rule

Environment build:

```bash
conda env create -n cod-paper -f code/environment-og14-zaf08.lock.yml
conda activate cod-paper
```

Then run:

```bash
cd code
python main.py
```
