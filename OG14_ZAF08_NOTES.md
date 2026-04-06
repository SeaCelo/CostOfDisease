This branch is the `OG-Core 0.14.3 / OG-ZAF 0.0.8` version of the Cost of Disease paper run.

It includes the following code changes relative to `main`:

- reproducible environment pinning
- corrected low excess-deaths scenario
- reform demographic-path fix, including holding infant mortality on its baseline path
- generated table cleanup

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
