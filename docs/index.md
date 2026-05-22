# ION Functions

ION Functions is a Python library of oceanographic data processing algorithms
for the [Ocean Observatories Initiative (OOI)](https://oceanobservatories.org/),
a National Science Foundation-funded ocean observing network that operates
fixed and mobile platforms across the full depth range of the ocean, from
surface buoys to the seafloor.

OOI instruments output raw data (L0) in instrument-native units — counts,
voltages, or proprietary binary formats. ION Functions applies vendor
calibration coefficients and documented algorithms to convert those raw
outputs into calibrated physical quantities (L1) and derived scientific
products (L2), following the data product specifications developed by the
OOI project science team. Algorithms are derived from OOI Data Product
Specifications (DPS), vendor documentation, and peer-reviewed literature.
Where discrepancies exist between sources, they are noted in the documentation.

The library covers a broad range of instrument families — CTD,
dissolved oxygen, CO<sub>2</sub>, fluorometry, pH, nitrate, velocity,
passive acoustics, seismology, optics, pressure, meteorology, and more —
organized as standalone Python modules under `ion_functions/data/`. Each
module contains pure functions that accept NumPy arrays as input and return
NumPy arrays as output, with no internal state.

For installation and usage, see the [Installation](installation.md) page.

## Data products

The primary output of ION Functions is a set of standardized oceanographic
data products (associated with individual variables within a data set), each 
assigned a processing level that reflects how far the data has been transformed 
from the raw instrument output. L0 is raw, L1 is a calibrated physical quantity 
derived from a single sensor channel, and L2 is a derived scientific product 
computed from two or more L1 or L2 inputs. See [Data Product Levels](data_products.md) for 
full definitions.

## Instrument families

Documentation is being added progressively as modules are converted to the
NumPy docstring format. Linked entries have full API documentation; unlinked
entries are pending conversion. The instrument families were defined at the
beginning of the project were used to organize the modules. In two cases, PH
and ADCP, a different name was used for the module, with the original project
defined family name noted in parentheses.

| Family                          | Class and Description | Data Products                                                                            |
|---|---|---|
| [CO2](api/co2_functions.md)     | PCO2A — partial pressure of pCO<sub>2</sub> in air and surface water<br>PCO2W — partial pressure of CO<sub>2</sub> in water | PCO2WAT, PCO2ATM, PCO2SSW, CO2FLUX                                                       |
| [PH (CO2)](api/ph_functions.md) | PHSEN — seawater pH | PHWATER                                                                                  |
| [CTD](api/ctd_functions.md)     | CTDAV — CTD for AUV<br>CTDBP — CTD with pump<br>CTDGV — CTD for gliders<br>CTDMO — CTD for moorings<br>CTDPF — CTD for profilers | TEMPWAT, PRESWAT, CONDWAT, PRACSAL, DENSITY                                              |
| [DO2](api/do2_functions.md)     | DOSTA — dissolved oxygen, stable response<br>DOFST — dissolved oxygen, fast response | DOCONCS, DOCONCF, DOXYGEN                                                                |
| FDC                             | FDCHP — direct covariance flux system | FLUXHOT, FLUXMOM, FLUXWET, WINDTUR, TMPATUR, RELHUMI, MOISTUR                            |
| [FLO](api/flo_functions.md)     | FLORD — 2-wavelength fluorometer (Chl-a, backscatter)<br>FLORT — 3-wavelength fluorometer (Chl-a, CDOM, backscatter) | CHLAFLO, CDOMFLO, FLUBSCT                                                                |
| [HYD](api/hyd_functions.md)     | HYDBB — broadband passive acoustic receiver, water column<br>HYDLF — broadband passive acoustic receiver, seafloor | HYDAPBB, HYDAPLF                                                                         |
| MET                             | METBK — bulk meteorology | TEMPAIR, BARPRES, WINDAVG, RELHUMI, SHRTIRR, LONGIRR, PRECIPM, TEMPSRF, CONDSRF, SPECHUM |
| MSP                             | MASSP — mass spectrometer for dissolved gases in hydrothermal and seep fluids | MASSPEC                                                                                  |
| [NIT](api/nit_functions.md)     | NUTNR — optical nitrate sensor | NITRTSC                                                                                  |
| [OBS](api/obs_functions.md)     | OBSBB — broadband seismometer<br>OBSSP — short-period seismometer | GRNDVEL, GRNDACC, SGRDVEL                                                                |
| [OPT](api/opt_functions.md)     | OPTAA — spectrophotometer (attenuation and absorption)<br>PARAD — photosynthetically active radiation<br>SPKIR — downwelling spectral irradiance | OPTATTN, OPTABSN, OPTPARW, SPECTIR                                                       |
| PRS                             | BOTPT — bottom pressure and tilt<br>PRESF — seafloor pressure<br>PREST — seafloor pressure | BOTPRES, BOTTILT, BOTSFLU, SFLPRES                                                       |
| SFL                             | FLOBN — benthic flow meter<br>OSMOI — osmotic fluid sampler<br>THSPH — hydrothermal vent fluid chemistry (temperature, H<sub>2</sub>, H<sub>2</sub>S, pH)<br>TMPSF — seafloor diffuse flow thermistor array<br>TRHPH — vent fluid temperature and resistivity | BENTHFL, TEMPSFL, THSPHTE, THSPHHC, THSPHHS, THSPHPH, TRHPHTE, TRHPHEH, TRHPHCC          |
| [VEL](api/vel_functions.md)     | VELPT — single-point horizontal velocity<br>VEL3D — single-point 3-D velocity | VELPTMN, VELPTTU, VELTURB                                                                |
| ADCP (VEL)                      | ADCPA — velocity profile, 50 m range<br>ADCPS — velocity profile, 600–700 m range<br>ADCPT — velocity profile, 300 m range<br>VADCP — 3-D velocity profile for turbulence | VELPROF, ECHOINT                                                                         |
| WAV                             | WAVSS — directional and non-directional wave statistics | WAVSTAT                                                                                  |
| Generic                         | Shared utility functions | —                                                                                        |

