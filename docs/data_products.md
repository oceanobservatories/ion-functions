# Data Product Levels

Data product levels apply to individual variables within a dataset, not to the
dataset as a whole. A single dataset will typically contain variables at
multiple processing levels simultaneously, reflecting the different
transformations applied to each measured or derived quantity.

- **L0** — Raw instrument output exactly as recorded by the sensor, in
  whatever units the instrument natively produces (counts, volts, or other
  unprocessed units). L0 values have no calibration or correction applied and
  are not directly interpretable as physical quantities.

- **L1** — A single variable transformed to produce a physically meaningful or
  standardized result. The transformation may be as simple as scaling units
  (e.g., converting velocity from mm/s to m/s), or may involve applying
  factory calibration coefficients and corrections (e.g., correcting for
  magnetic declination). The defining characteristic of an L1 product is that
  it is derived from a single input variable, though that input may itself be
  an L0 or another L1 product.

- **L2** — A variable derived by combining two or more L1 or L2 inputs through
  a more complex calculation. L2 products represent quantities that cannot be
  measured directly by a single sensor channel but are computed from multiple
  calibrated inputs. Practical salinity, for example, is an L2 product derived
  from conductivity, temperature, and pressure measurements.
