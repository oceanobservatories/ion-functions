# Using These Functions with OOI Data

The instrument family pages document the algorithms used to compute each OOI data product. This page explains how to
connect those functions to the data you download from OOI: what a download contains, where the calibration
coefficients come from, how the variables in an OOI data file map to the function arguments, and how to reproduce a
data product from its inputs.

## What an OOI Data Download Contains

OOI data are delivered as NetCDF files, one per deployment and data stream. Each file contains the instrument's raw
(L0) inputs alongside the derived (L1 and L2) data products computed from them, so in most cases everything a function
needs is already in the file. For example, a beam-coordinate ADCP file contains the four beam velocities
(`velocity_beam1` through `velocity_beam4`), the percent-good values, `heading`, `pitch`, and `roll` together with the
`eastward_seawater_velocity`, `northward_seawater_velocity`, and `upward_seawater_velocity` products computed from
them.

The variable names in the file are the OOI parameter names, not the function argument names. The mapping between the
two is described [below](#mapping-file-variables-to-function-arguments).

Files produced by the OOI data system use `obs` as the record dimension. With xarray, use
`ds.swap_dims({'obs': 'time'})` to index the data by time.

## Timestamps

OOI stores all times as seconds since 1900-01-01 00:00:00 UTC (the NTP epoch), and the functions in this library
expect time arguments in that form. The one exception is `dataqc_solarelevation`
([QC Functions](api/qc_functions.md)), which takes seconds since 1970-01-01 (Unix time).
`ntp_to_unix_time` on the [Generic Functions](api/generic_functions.md) page converts NTP to Unix time.

xarray decodes the `time` variable to `datetime64` values when a file is opened. To recover NTP seconds:

```python
ntp = (ds['time'].values - np.datetime64('1900-01-01')) / np.timedelta64(1, 's')
```

## Calibration Coefficients

Calibration coefficients appear in this library as function arguments, and in the OOI data system as named values
prefixed with `CC_` (for example `CC_wl` or `CC_scale_factor1`). Each family page with calibration coefficients has a
**Calibration Coefficients** table listing the function argument, the `CC_` name, and where the value comes from.
See [Calibration Coefficients](calibrations.md) for how coefficients are obtained and managed.

### Where to find the values

The coefficient values can be obtained in three ways.

**The asset-management repository.** All coefficient values are maintained as CSV files in the
[`oceanobservatories/asset-management`](https://github.com/oceanobservatories/asset-management/) GitHub repository,
under `calibration/<CLASS><SERIES>/` (for example `calibration/NUTNRB/`). Each file is named for the instrument and
the date from which the values apply (for example `CGINS-NUTNRB-00239__20120920.csv`) and has four columns:

| Column | Contents |
|---|---|
| `serial` | Instrument serial number |
| `name` | Coefficient name (`CC_...`) |
| `value` | Coefficient value; arrays are written as bracketed lists |
| `notes` | Source of the value, typically the vendor calibration file |

**The OOI M2M API.** The asset-management service that the OOI data system itself queries is available through the
M2M API (OOINet credentials required):
`https://ooinet.oceanobservatories.org/api/m2m/12587/asset/cal?refdes=<reference designator>`. It returns the
deployment events for the instrument, each with the `sensor` -> `calibration` values described below.

**The provenance file.** A data request can include a provenance file for each deployment, named
`deploymentNNNN_<stream>_provenance_<start>-<end>.json` (for requests made through the OOI M2M API, set
`include_provenance=true`). Under each deployment event, `sensor` -> `calibration` lists every coefficient applied to
that deployment by name, with its value, the date from which it applies (`eventStartTime`), and the asset-management
file it came from (`dataSource`).

When a request spans more than one calibration, the OOI data system applies each set of values to the records that
fall within its period. When you call a function yourself, supply the values that apply to the records you pass.

### Location values

`CC_lat` (or `CC_latitude`), `CC_lon` (or `CC_longitude`), and `CC_depth` look like calibration coefficients but are
not stored in the calibration files. The OOI data system takes them from the deployment record for the instrument, and
they are therefore not listed in the family-page Calibration Coefficients tables. They can be obtained in three ways:

- **The asset-management deployment files.** The `deployment/` directory of the
  [`oceanobservatories/asset-management`](https://github.com/oceanobservatories/asset-management/) repository has one
  CSV file per site (for example `CE02SHBP_Deploy.csv`), with one row per instrument deployment. The `lat`, `lon`, and
  `deployment_depth` columns supply `CC_lat`, `CC_lon`, and `CC_depth`, for the deployment identified by the
  `Reference Designator` and `deploymentNumber` columns.
- **The OOI M2M API.** The same `asset/cal?refdes=<reference designator>` request that returns the calibration
  coefficients (see [above](#where-to-find-the-values)) returns the `location` (latitude, longitude, and depth) of
  each deployment.
- **The data file.** The deployment latitude and longitude appear in the data file as the `lat` and `lon` global
  attributes.

### Names that do not describe the value

A few `CC_` names are historical and do not describe the quantity they hold. For example, the chi factor used for
the FLUBSCT_L2 backscatter calculation is stored as `CC_angular_resolution` (see [FLO Functions](api/flo_functions.md)).
Where this occurs it is noted on the family page.

## Mapping File Variables to Function Arguments

The OOI data system defines which inputs are passed to each function in the
[`oceanobservatories/preload-database`](https://github.com/oceanobservatories/preload-database) repository:

- `csv/ParameterFunctions.csv` lists each parameter function (`PFIDnnn`): the Python function, its module, and its
  argument names.
- `csv/ParameterDefs.csv` lists each parameter (`PDnnn`): its `name` (the variable name in the data file), units, and,
  for computed parameters, the function used (`parameterfunctionid`) and the argument map (`parameterfunctionmap`).

The argument map is a dictionary from function argument to input. Inputs take four forms:

| Form | Example | Meaning |
|---|---|---|
| `PDnnn` | `"b1": "PD3811"` | A parameter in the same data stream; look up its `name` in `ParameterDefs.csv` (here `velocity_beam1`) |
| `CC_...` | `"lat": "CC_lat"` | A calibration coefficient or location value (see [above](#calibration-coefficients)) |
| `dpi_<PRODUCT>` | `"ctd_sp": "dpi_PRACSAL_L2"` | A data product from a co-located instrument (see [below](#inputs-from-co-located-instruments)) |
| `<stream>.PDnnn` | `"botpres": "botpt_nano_sample.PD848"` | A parameter from another named data stream from the same instrument |

### Inputs from co-located instruments

Some functions need a data product measured by another instrument, most often temperature, salinity, or pressure from
a nearby CTD (for example nitrate, dissolved oxygen, and the backscatter calculation for the fluorometers). The
argument map names these by data product (`dpi_TEMPWAT_L1`, `dpi_PRACSAL_L2`), sometimes with a list of
alternatives in order of preference (for example `["dpi_TEMPWAT_L1", "dpi_TEMPSRF_L1"]` for subsurface and surface
platforms). The OOI data system finds a co-located instrument that provides the product and interpolates its values
onto the timestamps of the instrument being processed; the interpolated values are included in the download (for
example `int_ctd_pressure` in an ADCP file). When calling a function yourself, interpolate the co-located data onto
your instrument's timestamps in the same way.

### Array and spectral inputs

Some instruments report an array for each record, such as the 256-channel UV spectrum of the nitrate sensor or the
beam velocities in each ADCP depth bin. These are passed with one row per record, so that the arrays have shape
(records, channels) or (records, bins). Calibration coefficient arrays are passed with the same leading dimension: the
OOI data system tiles each coefficient array to one row per record, which allows the coefficients to change between
deployments within a single call.

### Inputs supplied by separate streams

A few products use inputs that the OOI data system supplies from a separate data stream rather than from the
instrument itself. For example, the BOTSFLU seafloor uplift products use predicted tides from the
`instrument_predicted_tide` stream (see [PRS Functions](api/prs_functions.md)).

## Worked Example: ADCP Velocity

This example reproduces the upward seawater velocity in a beam-coordinate ADCP file (a cabled ADCPT-B) from the
inputs in the same file.

The data product `upward_seawater_velocity` is `PD3821` in `ParameterDefs.csv`. Its parameter function is `PFID30`,
`adcp_beam_vertical`, with the argument map:

```
{"b1": "PD3811", "b2": "PD3812", "b3": "PD3813", "b4": "PD3814",
 "pg1": "PD873", "pg2": "PD874", "pg3": "PD875", "pg4": "PD876",
 "h": "PD672", "p": "PD673", "r": "PD674", "vf": "PD617"}
```

Looking up each `PDnnn` in `ParameterDefs.csv` gives the variable names in the file:

| Argument | Parameter | Variable in the file | Units |
|---|---|---|---|
| `b1`-`b4` | PD3811-PD3814 | `velocity_beam1`-`velocity_beam4` | mm s-1 |
| `pg1`-`pg4` | PD873-PD876 | `percent_good_beam1`-`percent_good_beam4` | percent |
| `h` | PD672 | `heading` | cdegrees |
| `p` | PD673 | `pitch` | cdegrees |
| `r` | PD674 | `roll` | cdegrees |
| `vf` | PD617 | `sysconfig_vertical_orientation` | 1 |

```python
import numpy as np
import xarray as xr

from ion_functions.data.adcp_functions import adcp_beam_vertical

ds = xr.open_dataset('deployment0013_CE02SHBP-LJ01D-05-ADCPTB104-streamed-adcp_velocity_beam_'
                     '20261031T031524_011850-20261101T031134_114704.nc')
ds = ds.swap_dims({'obs': 'time'})

w = adcp_beam_vertical(
    ds['velocity_beam1'].values, ds['velocity_beam2'].values,
    ds['velocity_beam3'].values, ds['velocity_beam4'].values,
    ds['percent_good_beam1'].values, ds['percent_good_beam2'].values,
    ds['percent_good_beam3'].values, ds['percent_good_beam4'].values,
    ds['heading'].values, ds['pitch'].values, ds['roll'].values,
    ds['sysconfig_vertical_orientation'].values)

print(np.nanmax(np.abs(w - ds['upward_seawater_velocity'].values)))   # ~1e-6 m/s
```

The horizontal components use the same inputs plus the deployment location and time, which `adcp_beam_eastward`
(`PFID28`) maps as `"lat": "CC_lat"`, `"lon": "CC_lon"`, and `"dt": "PD7"` (`time`). The location comes from the
deployment record and is in the file's `lat` and `lon` global attributes; the time must be converted back to NTP
seconds (see [Timestamps](#timestamps)):

```python
from ion_functions.data.adcp_functions import adcp_beam_eastward

n = ds.sizes['time']
lat = np.full(n, ds.attrs['lat'])
lon = np.full(n, ds.attrs['lon'])
ntp = (ds['time'].values - np.datetime64('1900-01-01')) / np.timedelta64(1, 's')

u = adcp_beam_eastward(
    ds['velocity_beam1'].values, ds['velocity_beam2'].values,
    ds['velocity_beam3'].values, ds['velocity_beam4'].values,
    ds['percent_good_beam1'].values, ds['percent_good_beam2'].values,
    ds['percent_good_beam3'].values, ds['percent_good_beam4'].values,
    ds['heading'].values, ds['pitch'].values, ds['roll'].values,
    ds['sysconfig_vertical_orientation'].values, lat, lon, ntp)
```
