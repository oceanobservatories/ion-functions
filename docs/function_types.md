# Function Types

Instrument family functions are grouped according to their usage within the
OOI system.

## Core Functions

These are the primary functions used in the calculation of oceanographic 
parameters by OOI from the different instrument families. They were built from 
multiple sources (e.g., vendor code and/or documentation, published algorithms, 
pseudocode developed by OOI Project Scientists, etc.). 

## Helper Functions

Some of the core functions require additional functions that are used
within the core function as part of the processing workflow. These are grouped 
as **helper** functions within each instrument family. For example, 
the [`flo_bback_total`](api/flo_functions.md#core-functions) function calls a 
total of four [helper functions](api/flo_functions.md#helper-functions): 
`flo_zhang_scatter_coeffs`, `flo_reactive_index`, `flo_isotherm_compress` and 
`flo_density_seawater`.   

Helper functions are documented alongside their core counterparts in each
instrument family page (where applicable), under a collapsed "Helper 
Functions" section.

## Wrapper Functions

Some modules expose thin wrapper functions in addition to the core and helper
functions. These wrappers exist because the OOI data management system 
historically imposed constraints on function outputs. For each function, the
system could only accept one output, while the core and helper functions
could return multiple outputs. External users should use the **core** and 
**helper** functions directly as they will be more efficient.

Two wrapper patterns appear in the codebase:

**Single-output wrappers** extract one value from a core function that returns
multiple outputs. The OOI data management system historically required one
output per function call; the core function is the right entry point for
external use since it returns all outputs at once.

**Named-product wrappers** call a single shared implementation function under
instrument-specific or product-specific names. For example, `flo_chla`,
`flo_cdom`, and `flo_beta` are all thin wrappers around `flo_scale_and_offset`,
which is the actual computation. The named wrappers exist so that each OOI data
product has a uniquely named function; external users should call
`flo_scale_and_offset` directly with the appropriate calibration coefficients.

Wrapper functions are documented alongside their core counterparts in each
instrument family page (where applicable), under a collapsed "Wrapper 
Functions" section.
