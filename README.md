# ion-functions

Python library of oceanographic data processing functions for the
[Ocean Observatories Initiative (OOI)](https://oceanobservatories.org/). It converts raw instrument output (L0) into
calibrated (L1) and derived (L2) data products, and applies corrections such as magnetic declination to products
computed onboard the instruments. Functions are organized by instrument family and take and return NumPy arrays.

**Documentation:** https://oceanobservatories.github.io/ion-functions/

**Repository:** https://github.com/oceanobservatories/ion-functions

## Installation

Install [Miniforge](https://github.com/conda-forge/miniforge), which uses the conda-forge channel by default, then
create the environment from `conda_env.yml` and install the package:

```bash
git clone https://github.com/oceanobservatories/ion-functions.git
cd ion-functions
conda env create -f conda_env.yml
conda activate ion
pip install --no-build-isolation --no-deps -e .
```

Runtime dependencies are listed in `conda_env.yml` (pinned to the OOI production versions) and `pyproject.toml`; the
test and documentation tools are optional extras in `pyproject.toml`. A C compiler is optional; see the
[installation guide](https://oceanobservatories.github.io/ion-functions/installation/) for details.

## Running the tests

```bash
pip install --no-build-isolation -e ".[test]"         # test tools (pytest, matplotlib)
pytest                                                # full suite
pytest ion_functions/data/test/test_ctd_functions.py  # one instrument family
```

## Project layout

```
ion_functions/
    data/        instrument family modules (ctd_functions.py, do2_functions.py, ...)
    data/test/   unit tests, one file per module
    qc/          quality control functions and Cython extension
    qc/test/     unit tests for the QC functions
extensions/      C source for the QC and polynomial calibration extensions
docs/            MkDocs documentation source
```

## Data product levels

| Level | Description | Example |
|-------|-------------|---------|
| L0 | Raw instrument output (counts, voltages, instrument-specific units) | Thermistor counts |
| L1 | Computed from L0 inputs from one instrument with calibration coefficients and metadata, or an instrument-computed value rescaled or corrected | Seawater temperature (&deg;C) |
| L2 | Derived value computed from at least one L1 or L2 product | Practical salinity (PSS-78) |

Algorithms come from the OOI Data Product Specifications (DPS), vendor documentation, and peer-reviewed literature.
Where sources disagree, the implemented code is authoritative and the discrepancy is noted in the documentation.

## Contributing

See [CONTRIBUTING.md](docs/CONTRIBUTING.md) for reporting bugs, proposing changes, and opening pull requests.

## License

Copyright 2013 UC Regents. Licensed under the [Apache License, Version 2.0](LICENSE).
