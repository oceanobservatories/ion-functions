# Ocean Observatories Initiative Ion Functions Documentation

This is documentation for the Ocean Observatories Initiative (OOI) Ion Functions package, which contains code for calculating L1 and L2 properties from OOI L0 observational data. For information on the OOI, visit [the website](https://oceanobservatories.org/). To explore OOI data including fields calculated using these ion functions, visit the [OOI Data Explorer](https://dataexplorer.oceanobservatories.org/). The source code for the Ion Functions module is available on [GitHub](https://github.com/oceanobservatories/ion-functions), and is linked directly throughout the documentation.

The Ion Functions package houses functions to generate various calculated parameters from OOI data parameters. These calculated parameters are calculated at request time by the OOI infrastructure and are delivered to users alongside the rest of the measured parameters for each instrument.


## Project layout
 * `ion_functions`: top level package
   - `data`: contains many modules for each of the different instruments requiring calculated parameters/
   - `qc`: contains various quality controls tests and datasets
   - `test`: contains unit tests for the modules in `data`

## Project documentation

### Reference documentation for each of the ion functions modules.
  * [ADCP](adcp_functions.md)
  * [CTD](ctd_functions.md)
  * [CO2](co2_functions.md)
  * [DOSTA](do2_functions.md)
  * [Fluorometer](flo_functions.md)
  * [Hydrophone](hyd_functions.md)
  * [Optical Backscatter](obs_functions.md)
  * [Met](met_functions.md)
  * [Disolved Gas](msp_functions.md)
  * [Disolved Nitrogen](nit_functions.md)
  * [Seismometer](obs_functions.md)
  * [OPTAA](opt_functions.md)
  * [pH](ph_functions.md)
  * [Seafloor Pressure](prs_functions.md)
  * [Seafloor Hydrothermal Vent](sfl_functions.md)
  * [Velocity](vel_functions.md)
  * [WAVSS](wav_functions.md)
### Quality Control and Testing
  * [QC information](qc-tests.md)
