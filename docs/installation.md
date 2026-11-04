# Installation

## Requirements

ion-functions is installed into a conda environment built from conda-forge packages. Install
[Miniforge](https://github.com/conda-forge/miniforge), which uses the conda-forge channel by default. The runtime
dependencies are listed in `conda_env.yml` at the repository root, pinned to the versions used in OOI production, and,
for pip, in `pyproject.toml`. The test and documentation tools are optional extras in `pyproject.toml`.

## Installing

```bash
git clone https://github.com/oceanobservatories/ion-functions.git
cd ion-functions
conda env create -f conda_env.yml
conda activate ion
pip install --no-build-isolation --no-deps -e .
```

The `--no-build-isolation --no-deps` options make pip build and install the package with the packages already in the
conda environment, rather than downloading its build tools and dependencies from PyPI.

## Optional C extensions

Two compiled extensions, the QC functions (`ion_functions.qc`) and `ion_functions.data.polycals`, are built when a C
compiler is available (gcc on Linux and macOS, the Microsoft C++ Build Tools on Windows). Without one, the installation
still succeeds and every function in the instrument family modules works; only those two modules are unavailable.
Deployments that need them should confirm both built:

```bash
python -c "import ion_functions.qc.qc_extensions, ion_functions.data.polycals"
```

To build them, install a C compiler -- for example conda-forge's compiler package (`conda install c-compiler`, Linux
and macOS) or the system gcc (`sudo apt install gcc` on Debian and Ubuntu) -- and then reinstall the package. pip
reports a skipped extension only as a warning, so use the import check above, or `pip install -v ...`, to confirm the
extensions were built.

## Running the tests

Install the test tools (pytest and matplotlib), then run the suite:

```bash
pip install --no-build-isolation -e ".[test]"
pytest
```

## Contributing

To report a bug, propose a change, or set up a development fork, see [Contributing](CONTRIBUTING.md).
