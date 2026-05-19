# Deprecations

Functions marked with the `@deprecated` decorator emit a
`DeprecationWarning` at call time. They are retained for reference but
should not be used in new code. A future cleanup PR will remove them.

| Module | Function | Reason |
|---|---|---|
| `hyd_functions.py` | `hyd_bb_acoustic_pwaves` | Pipeline never implemented |
| `hyd_functions.py` | `hyd_lf_acoustic_pwaves` | Pipeline never implemented |
| `obs_functions.py` | `obs_bb_ground_velocity` | Pipeline never implemented |
| `obs_functions.py` | `obs_bb_ground_acceleration` | Pipeline never implemented |
| `obs_functions.py` | `obs_sp_ground_velocity` | Pipeline never implemented |
