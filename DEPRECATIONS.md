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
| `msp_functions.py` | `calc_l2_totlgas_smph2scon` | Pipeline never implemented |
| `msp_functions.py` | `calc_l2_totlgas_smpco2con` | Pipeline never implemented |
| `msp_functions.py` | `calc_timestamp_totlgas_smph2scon` | Pipeline never implemented |
| `msp_functions.py` | `calc_timestamp_totlgas_smpco2con` | Pipeline never implemented |
| `msp_functions.py` | `calc_l2_totlgas_bkgh2scon` | Pipeline never implemented |
| `msp_functions.py` | `calc_l2_totlgas_bkgco2con` | Pipeline never implemented |
| `msp_functions.py` | `calc_timestamp_totlgas_bkgh2scon` | Pipeline never implemented |
| `msp_functions.py` | `calc_timestamp_totlgas_bkgco2con` | Pipeline never implemented |
| `msp_functions.py` | `calc_l2_mswater_smpphval` | Pipeline never implemented |
| `msp_functions.py` | `calc_l2_mswater_bkgphval` | Pipeline never implemented |
| `msp_functions.py` | `calc_msinlet_smpphint` | Pipeline never implemented |
| `msp_functions.py` | `calc_msinlet_smpphint_timestamp` | Pipeline never implemented |
| `msp_functions.py` | `calc_msinlet_bkgphint` | Pipeline never implemented |
| `msp_functions.py` | `calc_msinlet_bkgphint_timestamp` | Pipeline never implemented |
| `msp_functions.py` | `calc_msinlet_cal1phint` | Pipeline never implemented |
| `msp_functions.py` | `calc_msinlet_cal1phint_timestamp` | Pipeline never implemented |
| `msp_functions.py` | `calc_msinlet_cal2phint` | Pipeline never implemented |
| `msp_functions.py` | `calc_msinlet_cal2phint_timestamp` | Pipeline never implemented |
| `msp_functions.py` | `calc_smpnafeff` | Pipeline never implemented |
| `msp_functions.py` | `calc_smpnafeff_timestamp` | Pipeline never implemented |
| `msp_functions.py` | `calc_dissgas_smpmethcon` | Pipeline never implemented |
| `msp_functions.py` | `calc_dissgas_smpethcon` | Pipeline never implemented |
| `msp_functions.py` | `calc_dissgas_smph2con` | Pipeline never implemented |
| `msp_functions.py` | `calc_dissgas_smparcon` | Pipeline never implemented |
| `msp_functions.py` | `calc_dissgas_smph2scon` | Pipeline never implemented |
| `msp_functions.py` | `calc_dissgas_smpo2con` | Pipeline never implemented |
| `msp_functions.py` | `calc_dissgas_smpco2con` | Pipeline never implemented |
| `msp_functions.py` | `calc_dissgas_bkgmethcon` | Pipeline never implemented |
| `msp_functions.py` | `calc_dissgas_bkgethcon` | Pipeline never implemented |
| `msp_functions.py` | `calc_dissgas_bkgh2con` | Pipeline never implemented |
| `msp_functions.py` | `calc_dissgas_bkgarcon` | Pipeline never implemented |
| `msp_functions.py` | `calc_dissgas_bkgh2scon` | Pipeline never implemented |
| `msp_functions.py` | `calc_dissgas_bkgo2con` | Pipeline never implemented |
| `msp_functions.py` | `calc_dissgas_bkgco2con` | Pipeline never implemented |
| `msp_functions.py` | `calc_dissgas_cal1methcon` | Pipeline never implemented |
| `msp_functions.py` | `calc_dissgas_cal1co2con` | Pipeline never implemented |
| `msp_functions.py` | `calc_dissgas_cal2methcon` | Pipeline never implemented |
| `msp_functions.py` | `calc_dissgas_cal2co2con` | Pipeline never implemented |
| `msp_functions.py` | `calc_timestamp_smpmethcon` | Pipeline never implemented |
| `msp_functions.py` | `calc_timestamp_smpethcon` | Pipeline never implemented |
| `msp_functions.py` | `calc_timestamp_smph2con` | Pipeline never implemented |
| `msp_functions.py` | `calc_timestamp_smparcon` | Pipeline never implemented |
| `msp_functions.py` | `calc_timestamp_smph2scon` | Pipeline never implemented |
| `msp_functions.py` | `calc_timestamp_smpo2con` | Pipeline never implemented |
| `msp_functions.py` | `calc_timestamp_smpco2con` | Pipeline never implemented |
| `msp_functions.py` | `calc_timestamp_bkgmethcon` | Pipeline never implemented |
| `msp_functions.py` | `calc_timestamp_bkgethcon` | Pipeline never implemented |
| `msp_functions.py` | `calc_timestamp_bkgh2con` | Pipeline never implemented |
| `msp_functions.py` | `calc_timestamp_bkgarcon` | Pipeline never implemented |
| `msp_functions.py` | `calc_timestamp_bkgh2scon` | Pipeline never implemented |
| `msp_functions.py` | `calc_timestamp_bkgo2con` | Pipeline never implemented |
| `msp_functions.py` | `calc_timestamp_bkgco2con` | Pipeline never implemented |
| `msp_functions.py` | `calc_timestamp_cal1methcon` | Pipeline never implemented |
| `msp_functions.py` | `calc_timestamp_cal1co2con` | Pipeline never implemented |
| `msp_functions.py` | `calc_timestamp_cal2methcon` | Pipeline never implemented |
| `msp_functions.py` | `calc_timestamp_cal2co2con` | Pipeline never implemented |
| `msp_functions.py` | `calc_calrang_smpmethcon` | Pipeline never implemented |
| `msp_functions.py` | `calc_calrang_smpethcon` | Pipeline never implemented |
| `msp_functions.py` | `calc_calrang_smph2con` | Pipeline never implemented |
| `msp_functions.py` | `calc_calrang_smparcon` | Pipeline never implemented |
| `msp_functions.py` | `calc_calrang_smph2scon` | Pipeline never implemented |
| `msp_functions.py` | `calc_calrang_smpo2con` | Pipeline never implemented |
| `msp_functions.py` | `calc_calrang_smpco2con` | Pipeline never implemented |
| `msp_functions.py` | `calc_calrang_bkgmethcon` | Pipeline never implemented |
| `msp_functions.py` | `calc_calrang_bkgethcon` | Pipeline never implemented |
| `msp_functions.py` | `calc_calrang_bkgh2con` | Pipeline never implemented |
| `msp_functions.py` | `calc_calrang_bkgarcon` | Pipeline never implemented |
| `msp_functions.py` | `calc_calrang_bkgh2scon` | Pipeline never implemented |
| `msp_functions.py` | `calc_calrang_bkgo2con` | Pipeline never implemented |
| `msp_functions.py` | `calc_calrang_bkgco2con` | Pipeline never implemented |
| `msp_functions.py` | `calc_calrang_cal1methcon` | Pipeline never implemented |
| `msp_functions.py` | `calc_calrang_cal1co2con` | Pipeline never implemented |
| `msp_functions.py` | `calc_calrang_cal2methcon` | Pipeline never implemented |
| `msp_functions.py` | `calc_calrang_cal2co2con` | Pipeline never implemented |
| `msp_functions.py` | `gas_concentration` | Pipeline never implemented |
| `msp_functions.py` | `GasModeDetermination` | Pipeline never implemented |
| `msp_functions.py` | `SmpModeDetermination` | Pipeline never implemented |
