# ION Functions

ION Functions is a Python library of oceanographic data processing algorithms
for the [Ocean Observatories Initiative (OOI)](https://oceanobservatories.org/),
a National Science Foundation-funded ocean observing network that operates
fixed and mobile platforms across the full depth range of the ocean, from
surface buoys to the seafloor.

OOI records raw data in the instrument's native format, or wrapped in a format specific to the data logger system used
by the Marine Implementing Organization (MIO) operating the array. The OOI data system parses the raw data into named
parameters in datasets held in its internal data store. Depending on the instrument, the parsed parameters may be
defined as L0, L1, or L2 data products, or as engineering or other values that are not defined data products. When a
user requests a dataset, the OOI data system uses ION Functions to compute the L1 and/or L2 data products defined for
the dataset, including products the instrument may also report directly, applying scaling factors or vendor calibration
coefficients as needed and using the algorithms documented here. Algorithms are derived from OOI Data Product
Specifications (DPS), vendor documentation, and peer-reviewed literature. Where discrepancies exist between sources,
they are noted in the documentation.

The library covers a broad range of instrument families -- CTD,
dissolved oxygen, CO<sub>2</sub>, fluorometry, pH, nitrate, velocity,
passive acoustics, seismology, optics, pressure, meteorology, and more --
organized as standalone Python modules under `ion_functions/data/`. Each
module contains pure functions that accept NumPy arrays as input and return
NumPy arrays as output, with no internal state.

For installation and usage, see the [Installation](installation.md) page. To report a bug or
contribute changes, see [Contributing](CONTRIBUTING.md).

## Data products

The primary output of ION Functions is a set of standardized oceanographic data products (associated with individual 
variables within a data set), each assigned a processing level that reflects how far the data has been transformed 
from the instrument raw data. L0 is unprocessed data in instrument or sensor units; L1 is computed from one 
instrument's L0 inputs, together with calibration coefficients and metadata (if applicable), or is an 
instrument-computed value that has been rescaled or corrected; and L2 is a derived value computed with at least one L1 
or L2 product as an input. See [Data Product Levels](data_products.md) for full definitions.

## Instrument families

Each instrument family below links to its documentation page. The families follow the instrument groupings defined at
the start of the OOI project; in two cases, pH and ADCP, the module uses a different name, with the original family
name given in parentheses.

| Family                              | Class and Description | Data Products                                                                            |
|-------------------------------------|---|---|
| [CO2](api/co2_functions.md)         | PCO2A -- partial pressure of pCO<sub>2</sub> in air and surface water<br>PCO2W -- partial pressure of CO<sub>2</sub> in water | PCO2WAT, PCO2ATM, PCO2SSW, CO2FLUX                                                       |
| [PH (CO2)](api/ph_functions.md)     | PHSEN -- seawater pH | PHWATER                                                                                  |
| [CTD](api/ctd_functions.md)         | CTDAV -- CTD for AUV<br>CTDBP -- CTD with pump<br>CTDGV -- CTD for gliders<br>CTDMO -- CTD for moorings<br>CTDPF -- CTD for profilers | TEMPWAT, PRESWAT, CONDWAT, PRACSAL, DENSITY                                              |
| [DO2](api/do2_functions.md)         | DOSTA -- dissolved oxygen, stable response<br>DOFST -- dissolved oxygen, fast response | DOCONCS, DOCONCF, DOXYGEN                                                                |
| [FDC](api/fdc_functions.md)         | FDCHP -- direct covariance flux system | FLUXHOT, FLUXMOM, WINDTUR, TMPATUR                                                       |
| [FLO](api/flo_functions.md)         | FLORD -- 2-wavelength fluorometer (Chl-a, backscatter)<br>FLORT -- 3-wavelength fluorometer (Chl-a, CDOM, backscatter) | CHLAFLO, CDOMFLO, FLUBSCT                                                                |
| [HYD](api/hyd_functions.md)         | HYDBB -- broadband passive acoustic receiver, water column<br>HYDLF -- broadband passive acoustic receiver, seafloor | HYDAPBB, HYDAPLF                                                                         |
| [MET](api/met_functions.md)         | METBK -- bulk meteorology | TEMPAIR, BARPRES, WINDAVG, RELHUMI, SHRTIRR, LONGIRR, PRECIPM, TEMPSRF, CONDSRF, SPECHUM |
| [MSP](api/msp_functions.md)         | MASSP -- mass spectrometer for dissolved gases in hydrothermal and seep fluids | DISSGAS, CALRANG, TSTAMP, MSINLET, MSWATER, NAFEFF, GASMODE, SMPMODE, TOTLGAS            |
| [NIT](api/nit_functions.md)         | NUTNR -- optical nitrate sensor | NITRTSC                                                                                  |
| [OBS](api/obs_functions.md)         | OBSBB -- broadband seismometer<br>OBSSP -- short-period seismometer | GRNDVEL, GRNDACC, SGRDVEL                                                                |
| [OPT](api/opt_functions.md)         | OPTAA -- spectrophotometer (attenuation and absorption)<br>PARAD -- photosynthetically active radiation<br>SPKIR -- downwelling spectral irradiance | OPTATTN, OPTABSN, OPTPARW, SPECTIR                                                       |
| [PLIMS](api/plims_functions.md)     | PLIMS -- plankton imaging system (Imaging FlowCytobot) | FLOWRT |
| [PRS](api/prs_functions.md)         | BOTPT -- bottom pressure and tilt | BOTPRES, BOTTILT, BOTSFLU                                                                |
| [SFL](api/sfl_functions.md)         | FLOBN -- benthic flow meter<br>OSMOI -- osmotic fluid sampler<br>PRESF -- seafloor pressure (SBE 26plus)<br>PREST -- seafloor pressure (SBE 54)<br>THSPH -- hydrothermal vent fluid chemistry (temperature, H<sub>2</sub>, H<sub>2</sub>S, pH)<br>TMPSF -- seafloor diffuse flow thermistor array<br>TRHPH -- vent fluid temperature and resistivity | BENTHFL, TEMPSFL, THSPHTE, THSPHHC, THSPHHS, THSPHPH, TRHPHTE, TRHPHEH, TRHPHCC, SFLPRES |
| [VEL](api/vel_functions.md)         | VELPT -- single-point horizontal velocity<br>VEL3D -- single-point 3-D velocity | VELPTMN, VELPTTU                                                                         |
| [ADCP (VEL)](api/adcp_functions.md) | ADCPA -- velocity profile, 50 m range<br>ADCPS -- velocity profile, 600-700 m range<br>ADCPT -- velocity profile, 300 m range<br>VADCP -- 3-D velocity profile for turbulence | VELPROF, VELTURB, ECHOINT                                                                |
| [WAV](api/wav_functions.md)         | WAVSS -- directional and non-directional wave statistics | WAVSTAT                                                                                  |
| [Generic](api/generic_functions.md) | Shared utility functions | --                                                                                        |

## Quality control

Automated quality control tests (global and local range, spike, trend, stuck value, gradient, and combined
flags) that are applied to the data products above are documented on the [Quality Control](api/qc_functions.md)
page.
