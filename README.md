# ion-functions

Python library of oceanographic data processing functions for the
[Ocean Observatories Initiative (OOI)](https://oceanobservatories.org/).

ion-functions transforms raw instrument data (L0) into calibrated scientific
data products at L1 and L2 levels, covering CTD, ADCP, velocity, dissolved
oxygen, CO2, fluorometry, pH, pressure, meteorology, and more.

**Documentation:** https://oceanobservatories.github.io/ion-functions/

## Installation

```bash
conda env create -f conda_env.yml
conda activate ion
python setup.py develop
```

## Running tests

```bash
pytest
```

## Running C-extension tests

```bash
make
extensions/test
```

## License

See [LICENSE.txt](LICENSE.txt).
