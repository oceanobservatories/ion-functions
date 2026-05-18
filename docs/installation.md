# Installation

## Requirements

- Python 3.12
- Cython >= 3.0
- NumPy >= 2.4
- SciPy >= 1.17
- [gsw](https://teos-10.github.io/GSW-Python/) >= 3.6 (TEOS-10 seawater thermodynamics)
- [ppigrf](https://github.com/IAGA-VMOD/ppigrf) >= 2.0 (geomagnetic field model)

## Installing

The recommended approach uses conda to create an isolated environment with
all dependencies:

```bash
git clone https://github.com/oceanobservatories/ion-functions.git
cd ion-functions
conda env create -f conda_env.yml
conda activate ion
pip install -e .
```

`pip install -e .` will build the required Cython extensions automatically.

## Running the tests

```bash
pytest
```
