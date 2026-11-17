# PLIMS Functions

## Background

The Ocean Observatories Initiative deploys a McLane Labs Imaging FlowCytobot (IFCB) as its Plankton Imaging System
(PLIMS). The instrument is mounted at 7 m depth on the Near Surface Instrument Frame (NSIF) of the Central Surface
Mooring at the Pioneer MAB Array.

| Class | Hardware | Platform | Designator Meaning |
|---|---|---|---|
| PLIMS-A | McLane Labs Imaging FlowCytobot (IFCB) | Pioneer MAB Central Surface Mooring, NSIF (7 m) | Plankton Imaging System |

The IFCB creates three files for each sample, sharing the name `D[yyyymmdd]T[hhmmss]_IFCB[instrument_number]`:

- a header file (`.hdr`) with the instrument settings and a key to the format of the `.adc` file;
- an analog-to-digital converter file (`.adc`) with the data from the instrument sensors; and
- a region-of-interest file (`.roi`) with all the image data from the sample, stored as a binary stream.

Because of satellite telemetry bandwidth constraints, only the `.hdr` and `.adc` files are telemetered; the `.roi`
image files are available only after the instrument is recovered. The telemetered files provide diagnostic data and
metadata during the deployment.

`plims_functions.py` computes one data product, FLOWRT_L1. The other PLIMS data products are either parsed directly
from the header file or computed by the OOI data system from simple expressions defined in the OOI preload database,
not by ION Functions:

| Data Product | Units | Description | Source |
|---|---|---|---|
| HUMIDI_L1 | percent | Humidity within the instrument pressure housing | Parsed from the header file |
| TRGSEL_L1 | | Trigger selection: 2 to trigger on chlorophyll only (PMT B), 3 to trigger on both scattering and chlorophyll (PMT A and PMT B) | Parsed from the header file |
| SAMTIM_L1 | seconds since 1900-01-01 | Instrument time at the start of the sample run | Parsed from the header file |
| LOOKTI_L1 | s | Cumulative time spent capturing images during the sample run | Preload expression |
| FLOWRT_L1 | mL min-1 | Sample flow rate | `plims_flow_rate` |
| VOLANA_L2 | mL | Volume of sample analyzed during the run time | Preload expression |
| ROIPML_L2 | counts mL-1 | Number of regions of interest (ROIs) recorded per mL of sample analyzed | Preload expression |

This page documents the algorithm used to compute FLOWRT_L1 and how it feeds the other derived PLIMS products. For
locating the L0 inputs, see [Using These Functions with OOI Data](../using-ooi-data.md).

### Primary Sources

| Document | Title |
|---|---|
| 1341-00850 | OOI (2026). Data Product Specification: PLIMS, Version 0-06. Document Control Number 1341-00850. (Not publicly posted.) |

---

### FLOWRT_L1 -- Sample Flow Rate

FLOWRT_L1 is the rate at which the IFCB injects sample (mL min$^{-1}$). The inputs are instrument settings reported
in the header file, not calibration coefficients.

#### FLOWRT_L1 Algorithm

The standard injection rate is 0.25 mL min$^{-1}$. When the instrument is set to run fast (`runSampleFast`), the
standard rate is multiplied by the run fast factor (`RunFastFactor`):

$$\text{FLOWRT} = 0.25 \times f$$

where $f$ is the run fast factor when run fast is set, and 1 otherwise. Per DPS 1341-00850, run fast is always set
when sampling.

#### Data Product Chain

FLOWRT_L1 is one step in the chain that produces the volume analyzed and the ROI concentration. The other steps are
computed by the OOI data system, as specified in DPS 1341-00850:

$$\text{LOOKTI} = t_{run} - t_{inhibit}$$

$$\text{VOLANA} = \frac{\text{LOOKTI} \times \text{FLOWRT}}{60}$$

$$\text{ROIPML} = \frac{N_{roi}}{\text{VOLANA}}$$

where $t_{run}$ is the run time (s), $t_{inhibit}$ is the inhibit time (s), and $N_{roi}$ is the number of ROIs
recorded in the sample, all reported in the IFCB header and `.adc` files. The volume analyzed is usually less than
the volume sampled (typically 5 mL).

---

## Core Functions

::: ion_functions.data.plims_functions.plims_flow_rate

#### History
| Date | Author | Change |
|---|---|---|
| 2025-02-04 | Joffrey Peters | Initial code |
| 2026-10-30 | Samuel Dahlberg | Removed sample speed and sample volume from the calculation (Redmine #16237) |
| 2026-11-16 | Christopher Wingard | Converted to NumPy docstring format; added documentation page |

---

## References

McLane Labs. *Imaging FlowCytobot User Manual*, Rev 20.E.05.

OOI (2026). Data Product Specification: PLIMS, Version 0-06. Document Control Number 1341-00850. (Not publicly
posted.)
