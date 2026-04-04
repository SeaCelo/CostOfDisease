# Simulating the costs of disease

The scripts in this directory are set up to simulate the cost of disease for South Africa.  This is determined in `main.py` which calibrates OG-Core to South African data and then simulates the cost of disease for South Africa.

Descriptions of files in this directory:
* `main.py`: This is the script to execute. It runs the simulations and saves the results to disk.
* `get_pop_data.py`: This script has utilities to get the population data for the baseline and to compute the counterfactual population under increased mortality.
* `create_plots_tables.py`: This script has utilities to create plots and tables for the results.

## Reproducible environment

Use `environment-osx-arm64.lock.yml` to rebuild the exact tested package set on Apple Silicon Macs:

```bash
conda env create -p /tmp/cod-disease-paper-env -f code/environment-osx-arm64.lock.yml
```

Then verify the package stack directly:

```bash
MPLCONFIGDIR=/tmp/cod-disease-paper-env-mpl \
XDG_CACHE_HOME=/tmp/cod-disease-paper-env-cache \
/tmp/cod-disease-paper-env/bin/python -c "import ogcore, ogzaf, get_pop_data, create_plots_tables, main; import numpy; print(numpy.__version__)"
```

If you need a fresh solve from `environment.yml`, keep `numpy<2` and `pandas<3` in the specification. The pinned `ogcore==0.14.3` stack fails under NumPy 2 due scalar assignment behavior in the steady-state solver, and `ogzaf==0.0.8` is not compatible with `pandas` 3.x through the `pandas-datareader` import path.
