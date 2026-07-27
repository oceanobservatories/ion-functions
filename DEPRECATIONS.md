# Deprecations

Functions marked with the `@deprecated` decorator emit a
`DeprecationWarning` at call time. They are retained for reference but
should not be used in new code. A future cleanup PR will remove them.

| Module | Function | Reason |
|---|---|---|
| `fdc_functions.py` | `fdc_flux_and_wind` | Burst assembly incompatible with current data system |
| `fdc_functions.py` | `fdc_tmpatur` | Burst assembly incompatible with current data system |
| `fdc_functions.py` | `fdc_windtur_north` | Burst assembly incompatible with current data system |
| `fdc_functions.py` | `fdc_windtur_up` | Burst assembly incompatible with current data system |
| `fdc_functions.py` | `fdc_windtur_west` | Burst assembly incompatible with current data system |
| `fdc_functions.py` | `fdc_fluxhot` | Burst assembly incompatible with current data system |
| `fdc_functions.py` | `fdc_fluxmom_alongwind` | Burst assembly incompatible with current data system |
| `fdc_functions.py` | `fdc_fluxmom_crosswind` | Burst assembly incompatible with current data system |
| `fdc_functions.py` | `fdc_time_L1` | Burst assembly incompatible with current data system |
| `fdc_functions.py` | `fdc_time_L2` | Burst assembly incompatible with current data system |
| `hyd_functions.py` | `hyd_bb_acoustic_pwaves` | Pipeline never implemented |
| `hyd_functions.py` | `hyd_lf_acoustic_pwaves` | Pipeline never implemented |
| `obs_functions.py` | `obs_bb_ground_velocity` | Pipeline never implemented |
| `obs_functions.py` | `obs_bb_ground_acceleration` | Pipeline never implemented |
| `obs_functions.py` | `obs_sp_ground_velocity` | Pipeline never implemented |
