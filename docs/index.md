# ION Functions

Python library of oceanographic data processing functions for the
[Ocean Observatories Initiative (OOI)](https://oceanobservatories.org/).

ION Functions transforms raw instrument data (L0) into calibrated scientific
data products at L1 and L2 levels, covering CTD, ADCP, velocity, dissolved
oxygen, CO2, fluorometry, pH, pressure, meteorology, and more.

## Instrument families

Documentation is being added progressively as modules are converted to the
NumPy docstring format. Linked entries have full API documentation; unlinked
entries are pending conversion.

| Module                              | Instruments                                  | Data products |
|-------------------------------------|----------------------------------------------|---------------|
| [CTD](api/ctd_functions.md)         | SBE 16Plus, SBE 37IM, SBE 52MP, glider CTDs | TEMPWAT, PRESWAT, CONDWAT, PRACSAL, DENSITY |
| [CO2](api/co2_functions.md)         | SAMI2-CO2, CO2-Pro                           | PCO2WAT, PCO2ATM, PCO2SSW, CO2FLUX |
| [Fluorometer](api/flo_functions.md) | ECO FLORD, FLORT, FLNTU                      | CHLAFLO, CDOMFLO, FLUBSCT |
| [Hydrophone](api/hyd_functions.md)  | HYDBB, HYDLF                                 | HYDAPBB, HYDAPLF |
| [pH](api/ph_functions.md)           | SAMI-II pH, Sea-Bird Deep SeapHOx V2         | PHWATER |
| Dissolved oxygen                    | SBE 43, Aanderaa Optode                      | DOXYGEN |
| Water velocity (single-point)       | Aquadopp, FSI, Nobska MAVS                   | VELPTTU |
| Ocean bottom seismometer            | OBSBB, OBSBK                                 | GRNDVEL, GRNDACC, SGRDVEL |
| Nitrate                             | SUNA                                         | NITROPT |
| Optical                             | SPECTIR, SPKIR, PARAD, FLORT                 | OPTATTN, OPTABSN, PAR |
| Pressure                            | SBE 26Plus, Nano                             | PRESWAT |
| Seafloor                            | TRHPH, RASFL                                 | THSPHTE, SULFIDE |
| ADCP                                | Workhorse, Pinnacle, VADCP                   | VELPROF, ECHOINT |
| Flux direct covariance              | FDCHP                                        | FDCHP products |
| Mass spectrometer                   | MASSP                                        | MASSP products |
| Meteorology                         | METBK                                        | WINDAVG, TEMPAIR, BARPRES, and others |
| Generic                             | —                                            | Shared utilities |
| QC functions                        | —                                            | Quality control flags |

## Data product levels

- **L0** — Raw instrument output (counts, voltages)
- **L1** — Converted/calibrated engineering units
- **L2** — Derived scientific products (e.g., practical salinity, dissolved oxygen concentration)

## Quick start

See the repository [README](https://github.com/oceanobservatories/ion-functions) for installation and usage.
