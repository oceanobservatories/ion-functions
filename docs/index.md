# Ocean Observatories Initiative Ion Functions Documentation

For information on the Ocean Observatories Initiative (OOI), visit [the website](https://oceanobservatories.org/). To explore OOI data including fields calculated using these ion functions, visit the [OOI Data Explorer](https://dataexplorer.oceanobservatories.org/). The source code is available on [GitHub](https://github.com/oceanobservatories/ion-functions).

The Ion Functions package houses functions to generate various calculated parameters from OOI data parameters. These calculated parameters are calculated on the fly by the OOI infrastructure and are delivered to users alongside the rest of the measured parameters for each instrument.


## Project layout
 * `ion_functions`: top level package
   - `data`: contains many modules for each of the different instruments requiring calculated parameters/
   - `qc`: contains various quality controls tests and datasets
   - `test`: contains unit tests for the modules in `data`

## Project documentation

### Reference documentation for each of the ion functions modules.
  * [CTD](ctd_functions.md)
  * [CO2](co2_functions.md)
### Quality Control and Testing
  * [QC information](qc-tests.md)
