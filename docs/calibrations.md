# Calibration Coefficients

ION functions may require calibration coefficients as inputs to convert raw
instrument data into calibrated scientific data products.

## Sources and Management

All calibration coefficients used by ion-functions are supplied with
individual instruments and updated during the service and recalibration
process. Coefficient sources vary by instrument type -- some come from
vendor-supplied calibration documents, some are downloaded directly from
the instrument over a serial connection, and some are obtained from both
sources (cross-checked to verify agreement). Entry into the OOI Asset
Management database (maintained as CSV files in the
[`oceanobservatories/asset-management`](https://github.com/oceanobservatories/asset-management/)
GitHub repository) is through a mix of manual and scripted processes depending
on the source. All entries are reviewed and confirmed before use.

Coefficients are updated on a recalibration cycle tied to each instrument's
service schedule -- typically after each recovery and refurbishment for moored
instruments, or on a time-based interval for assets such as gliders that may
be deployed and recovered multiple times between calibrations. For most
instrument types the calibration is performed by the manufacturer or vendor;
a small number of instrument types are calibrated in-house or receive a secondary 
pre-deployment calibration.

## Use in the Library

Within ion-functions, calibration coefficients are passed directly as
function arguments -- scalars or arrays depending on the function. The library
has no built-in mechanism for loading or retrieving calibration values; it
performs no database lookups and reads no configuration files at runtime.
Users calling functions directly are responsible for supplying correct
coefficient values. In the OOI data pipeline, the system retrieves the
appropriate coefficients from the Asset Management database and passes them
to the functions automatically.

## Naming Conventions

Calibration coefficient names in each function's parameter list generally 
reflect the naming conventions used in the vendor's calibration documentation 
or the relevant OOI Data Product Specification. In most cases, the names match
the vendor's convention directly; where they differ, the name is defined on the 
relevant instrument family page. The meaning and units of each coefficient are
given where it is used in the Background section of the relevant instrument family page, and its OOI
calibration name (`CC_...`) and source are listed in that page's Calibration Coefficients table.

## Finding the Coefficients for an OOI Dataset

The coefficient values applied to a dataset can be obtained from the calibration CSV files in the asset-management
repository, through the OOI M2M API, or from the provenance file that can accompany an OOI data download. The location values `CC_lat`,
`CC_lon`, and `CC_depth` are taken from the deployment record rather than the calibration files. See
[Using These Functions with OOI Data](using-ooi-data.md#calibration-coefficients) for details.