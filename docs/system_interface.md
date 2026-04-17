# OOI System Interface

Some modules expose thin wrapper functions in addition to the core functions
documented in the API Reference. These wrappers exist because the OOI data
management system historically imposed constraints on function outputs.
External users should use the **core functions** directly.

Two wrapper patterns appear in the codebase:

**Single-output wrappers** extract one value from a core function that returns
multiple outputs. The OOI data management system historically required one
output per function call; the core function is the right entry point for
external use since it returns all outputs at once.

**Named-product wrappers** call a single shared implementation function under
instrument-specific or product-specific names. For example, `flo_chla`,
`flo_cdom`, and `flo_beta` are all thin wrappers around `flo_scale_and_offset`,
which is the actual computation. The named wrappers exist so that each OOI data
product has a uniquely named callable; external users should call
`flo_scale_and_offset` directly with the appropriate calibration coefficients.

Wrapper functions are documented alongside their core counterparts in each
instrument family page, under a collapsed "OOI System Interface" section.
