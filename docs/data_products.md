# Data Products

OOI Data products, and their different levels, are defined for individual variables (or parameters) within a dataset,
not to the dataset as a whole. A single dataset may contain 0-N variables defined as an OOI data product. These may
exist at multiple processing levels simultaneously reflecting the different transformations applied to each measured
or derived quantity. The data product levels defined by OOI are:

- **L0** -- Raw instrument output exactly as recorded by the sensor, in whatever units the instrument natively 
  produces (counts, volts, or other unprocessed units). L0 values have no calibration or correction applied and
  are usually not directly interpretable as physical quantities.

- **L1** -- Computed from L0 inputs from one instrument, together with calibration coefficients and metadata such
  as position and time; or an instrument-computed value that has been rescaled or corrected (e.g., for magnetic
  declination). L1 variables may be reported directly by the instrument (onboard calculation), calculated by ION
  Functions, or both.

- **L2** -- A derived value, computed with at least one or more L1 or L2 product as an input. The inputs may come from 
  the same instrument (e.g., practical salinity, PRACSAL_L2, from the CTD's own conductivity, temperature, and
  pressure) or from another instrument (e.g., dissolved oxygen, DOXYGEN_L2, corrected using salinity from a co-located 
  CTD).
