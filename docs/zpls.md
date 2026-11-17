# Bioacoustic Sonar (ZPLSC, ZPLSG)

OOI measures water column acoustic backscatter with bioacoustic sonars at most of its Coastal and Global array sites:
ZPLSC on the coastal moorings and cabled platforms, and ZPLSG on the Global Apex Profiler Moorings. ION Functions
contains no functions for these instruments. Their data are processed separately, with the
[ooi-zpls-echograms](https://github.com/oceanobservatories/ooi-zpls-echograms) package, which is built on the
open-source [echopype](https://echopype.readthedocs.io/en/latest/) library (Lee et al., 2024). This page summarizes how
that processing works; see the package [README](https://github.com/oceanobservatories/ooi-zpls-echograms) for more information.

## Instruments

The sensor depends on whether the site is uncabled or cabled: the ASL Environmental Sciences Acoustic Zooplankton Fish 
Profiler (AZFP) at uncabled sites, and a Kongsberg Simrad echo sounder at the cabled sites.

| Sensor | Sites | Instrument raw data |
|---|---|---|
| [AZFP](https://www.aslenv.com/product-details.html#AZFP) | Uncabled Coastal and Global sites | `.01A` files, plus an instrument `.XML` configuration file |
| [EK60](https://www.kongsberg.com/discovery/products/discontinued/simrad-ek60/) | Cabled sites (CE02SHBP, CE04OSPS), before summer 2023 | `.raw` files (with `.bot` and `.idx`) |
| [EK80](https://www.kongsberg.com/what-we-do/ocean-space/ocean-science/ek-scientific-echo-sounder/) | Cabled sites (CE02SHBP, CE04OSPS), since summer 2023 | `.raw` files |

The Global Apex Profiler Moorings carry two AZFPs on the 150 m sphere, one looking up toward the surface and one
looking down.

Differences that affect processing:

- **AZFP** processing requires the instrument's `.XML` file, which holds the instrument-specific settings and
  calibration coefficients.
- **EK80** data are processed for the narrowband (CW) channels only, recorded with either power or complex sample
  encoding (detected per file). Broadband (LFM) channels are dropped; files with only LFM channels, and files recorded
  at a full broadband duty cycle, are skipped.
- The converted files are NetCDF for the AZFP and EK60, and Zarr for the EK80.

## Processing

The instrument raw data are available on the [OOI Raw Data Server](https://rawdata.oceanobservatories.org/files/).
ooi-zpls-echograms processes them in three steps:

1. **Converted** -- Each instrument raw data file is converted, one to one, with the echopype convert function to
   echopype's adaptation of the SONAR-netCDF4 convention (Macaulay and Pena, 2018).
2. **Calibrated** -- The echopype `compute_Sv` function computes the volume backscattering strength (Sv), the range
   from the sensor face, and the nominal depth of each range bin. Results are saved in daily NetCDF files.
3. **Averaged** -- Sv is averaged in time (15-minute means for the coastal ZPLSC, 60-minute means for the global
   ZPLSG), with gaps of up to three averaging windows filled by linear interpolation. Echograms are drawn from the
   averaged data.

Weekly echograms let users quickly assess the quality and usefulness of the data for their research. The processed
NetCDF files are preliminary and are provided in a standard format for further processing and analysis.

## References

Lee, W.-J., Setiawan, L., Tuguinay, C., Mayorga, E., and Staneva, V. (2024). Interoperable and scalable echosounder
data processing with Echopype. ICES Journal of Marine Science, 81(10), 1941-1951.
[doi:10.1093/icesjms/fsae133](https://doi.org/10.1093/icesjms/fsae133).

Macaulay, G. and Pena, H. (Eds.) (2018). The SONAR-netCDF4 convention for sonar data, Version 1.0. ICES Cooperative
Research Report No. 341, 33 pp. [doi:10.17895/ices.pub.4392](https://doi.org/10.17895/ices.pub.4392).

[ooi-zpls-echograms: OOI Bioacoustic Sonar Processing and Echogram
Generation.](https://github.com/oceanobservatories/ooi-zpls-echograms)
