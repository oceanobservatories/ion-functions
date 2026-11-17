# Data Products

OOI data product levels are defined in the OOI Data Management Plan (DMP), Section 4.1.1 (OOI, 2014):

- **Instrument raw data** -- The data as received from the instrument, in an instrument-specific format. Depending on
  the instrument, they may interleave data from several sensors, be in native sensor units (e.g., counts, volts), or
  already include processing performed within the instrument (e.g., primary calibration). Instrument raw data are
  always archived by OOI, on the [OOI Raw Data Server](https://rawdata.oceanobservatories.org/files/), and are not
  formal data products. Example: a format 0 binary file from an SBE 37IM on a Global Flanking Mooring.

- **L0 (Unprocessed Data)** -- Data in instrument/sensor units and at instrument/sensor resolution, separated sensor
  by sensor (unpacked and de-interleaved), and available in OOI-supported formats (e.g., NetCDF). L0 data are
  sometimes mistakenly called "raw" data. They are always archived by OOI, which enables reprocessing and lets
  researchers apply their own quality control and calibration.

- **L1 (Basic Data)** -- L0 data that have been calibrated and are in well-defined scientific units; they may have
  some quality control (QC) applied.

- **L2 (Derived Data Products)** -- Data that have been manipulated and/or combined, and that include new values not
  present in the original data.

## Levels in ION Functions

OOI assigns a level to each individual variable, not to a dataset: a single dataset typically holds L0, L1, and L2
variables side by side, and the level of a variable cannot be inferred from the dataset or file name.

In practice, OOI records raw data in the instrument's native format, or wrapped in a format specific to the data logger
system used by the Marine Implementing Organization (MIO) operating the array.

When a user requests a dataset, the OOI data system uses ION Functions to compute the L1 and/or L2 data products defined
for that dataset from the parsed L0, L1, or L2 parameters. Some instruments also report L1 or L2 values computed
onboard; OOI may recompute those products independently, and a disagreement between the two can expose errors in the
calibration coefficients. ION Functions applies the levels as follows:

- **L1** -- Computed from L0 inputs from one instrument, together with calibration coefficients and metadata such
  as position and time (if applicable); or from an instrument-computed L1 value that has been rescaled or corrected 
- (e.g., for magnetic declination). L1 variables may be reported directly by the instrument (onboard calculation), 
- calculated by ION Functions, or both.

- **L2** -- A derived value computed with at least one L1 or L2 product as an input. The inputs may come from the
  same instrument (e.g., practical salinity, PRACSAL_L2, from the CTD's own conductivity, temperature, and pressure)
  or from another instrument (e.g., dissolved oxygen, DOXYGEN_L2, corrected using salinity from a co-located CTD).

## References

OOI (2014). Data Management Plan, Version 2-00. Document Control Number 1102-00000. (Not publicly posted.)
