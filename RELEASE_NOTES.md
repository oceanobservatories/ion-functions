# Release 2.4.9 2023-03-15

Issue #15558 - Add PCO2W/PHSEN voltage calculations 
Issue #15072 - Address ADCP in the Sky

# Release 2.4.8 2023-01-18

Issue #15346 - Change in PCO2W/PHSEN (SAMI) data output due to new Rev K board

# Release 2.4.7 2022-02-23

Issue #15167 - Glider Oxygen data fix

# Release 2.4.6 2020-06-09

Issue #14537 - Updates to use numpy=1.16

# Release 2.4.5 2020-03-31

Issue #13402 - fix assorted bugs in FDCHP data product algorithms
- Update FDCHP functions for changes to numpy/scipy libraries
- Update FDCHP and PRS tests for library changes
- Reduce PRS function/test footprint so tests can run without MemoryError

# Release 2.4.4 2019-11-07

Issue #14302 - ADCPTE101 - failure to generate data products
- Adds ability to compute a 3-beam solution when converting from beam to Earth coordinates

# Release 2.4.3 2019-08-19

Issue #11399 - Expose,use depth from pressure function

# Release 2.4.2 2019-02-14

Issue #13461 - corrected VEL3D-K data product algorithm.

# Release 2.4.1 2018-12-11

PCO2W corrections

# Release 2.4.0 2018-05-18

Issue #13276 - make stuck test handle invalid arguments
- Use the absolute value of num parameter if it is negative
- Use the default value of 10 for num if it is 0

# Version 0.0.2

* Issue 12615 Merged updates from ooici/ion-functions repository

# Version 0.0.1

* Initial version - see also ion_functions/version.py
