# ion-functions

Python library of oceanographic data processing functions for the
[Ocean Observatories Initiative (OOI)](https://oceanobservatories.org/).

ion-functions transforms raw instrument data into calibrated scientific
data products and applies additional corrections (e.g., magnetic declination) 
to already calculated data products. It covers a broad range of instrument
families, including CTD, ADCP, velocity, dissolved oxygen, CO<sub>2</sub>, fluorometry,
pH, pressure, optical, meteorology, and more.

**Documentation:** https://oceanobservatories.github.io/ion-functions/

**Repository:** https://github.com/oceanobservatories/ion-functions

## Overview

OOI instruments output raw data (L0) in instrument-specific units — counts,
voltages, or proprietary formats. ion-functions applies vendor calibration
coefficients and documented algorithms to produce:

- **L1** — calibrated physical quantities (e.g., temperature in &deg;C, pressure
  in dbar)

Instruments also output data in already calculated physical quantities. These
may require additional corrections (e.g., unit conversions, magnetic declination).
ion-functions applies these corrections either as part of the conversion from L0 
to L1, or to the already calculated physical quantities. Additional functions
are included which are used in the calculations of derived products, where 2 or
more L1 data products are combined to create a new L2 data product.

- **L2** — derived scientific products (e.g., practical salinity, dissolved
  oxygen concentration, partial pressure of CO<sub>2</sub>)

The functions are organized according to instrument family and are 
implemented as standalone Python modules under `ion_functions/data/`, with 
pure functions that take NumPy arrays as input and return NumPy arrays as 
output.

Algorithms are derived from OOI Data Product Specifications (DPS), vendor
documentation, and peer-reviewed literature. Where discrepancies exist
between and within sources, the implemented code is authoritative. Discrepancies
are noted in the documentation.

## Requirements

- Python = 3.12
- Cython >= 3.0 (required to build C extensions)
- NumPy >= 2.4
- SciPy >= 1.17
- [gsw](https://teos-10.github.io/GSW-Python/) >= 3.6 (TEOS-10 seawater
  thermodynamics)
- [ppigrf](https://github.com/IAGA-VMOD/ppigrf) >= 2.0 (geomagnetic field
  model)
- mkdocstrings (for the creation and editing of the documentation pages)
- mkdocstrings-python
- mkdocs-material

## Installation

The recommended approach uses conda to create an isolated environment with
all dependencies:

```bash
git clone https://github.com/oceanobservatories/ion-functions.git
cd ion-functions
conda env create -f conda_env.yml
conda activate ion
pip install -e .
```

### Building the C extensions

Two Cython extensions must be compiled before use:

- `ion_functions/qc/qc_extensions.pyx` — quality control algorithms (stuck,
  spike, gradient), wrapping C source in `extensions/`
- `ion_functions/data/polycals.pyx` — polynomial calibration calculations

`pip install -e .` will build them automatically. To rebuild explicitly:

```bash
python setup.py build_ext --inplace
```

## Running the tests

```bash
pytest
```

To run tests for a specific instrument family:

```bash
pytest ion_functions/data/test/test_ctd_functions.py
```

### C extension tests

A separate test suite exercises the C extensions directly:

```bash
make
extensions/test
```

## Project structure

```
ion_functions/
    data/           # Core data processing functions, one module per
                    # instrument family (ctd_functions.py, do2_functions.py, …)
    data/test/      # Unit tests, one file per data module
    data/perf/      # Performance benchmarks (not part of normal test suite)
    qc/             # Quality control functions and Cython extension
    qc/test/        # Unit tests for QC functions
extensions/         # C source for performance-critical QC algorithms
docs/               # MkDocs documentation source
```

## Data product levels

| Level | Description | Example |
|-------|-------------|---------|
| L0 | Raw instrument output (counts, voltages, proprietary units) | SBE thermistor counts |
| L1 | Calibrated physical quantity | Temperature (°C) |
| L2 | Derived scientific product | Practical salinity (PSS-78) |

## Contributing

Contributions are welcome. Please see [CONTRIBUTING.md](CONTRIBUTING.md) for
guidelines on submitting bug reports, proposing changes, and opening pull
requests.

## License

Copyright 2013 UC Regents. Licensed under the
[Apache License, Version 2.0](LICENSE.txt).
