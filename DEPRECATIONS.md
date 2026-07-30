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
| `fdc_functions.py` | `fdc_accelsclimode` | Burst assembly incompatible with current data system; not decorated, internal-only helper |
| `fdc_functions.py` | `fdc_alignwind` | Burst assembly incompatible with current data system; not decorated, internal-only helper |
| `fdc_functions.py` | `fdc_anglesclimodeyaw` | Burst assembly incompatible with current data system; not decorated, internal-only helper |
| `fdc_functions.py` | `fdc_despikesimple` | Burst assembly incompatible with current data system; not decorated, internal-only helper |
| `fdc_functions.py` | `fdc_detrend` | Burst assembly incompatible with current data system; not decorated, internal-only helper |
| `fdc_functions.py` | `fdc_filtcoef` | Burst assembly incompatible with current data system; not decorated, internal-only helper |
| `fdc_functions.py` | `fdc_grv` | Burst assembly incompatible with current data system; not decorated, internal-only helper |
| `fdc_functions.py` | `fdc_process_compass_data` | Burst assembly incompatible with current data system; not decorated, internal-only helper |
| `fdc_functions.py` | `fdc_quantize_data` | Burst assembly incompatible with current data system; not decorated, internal-only helper |
| `fdc_functions.py` | `fdc_sonic` | Burst assembly incompatible with current data system; not decorated, internal-only helper |
| `fdc_functions.py` | `fdc_trans` | Burst assembly incompatible with current data system; not decorated, internal-only helper |
| `fdc_functions.py` | `fdc_update` | Burst assembly incompatible with current data system; not decorated, internal-only helper |
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
| `msp_functions.py` | `average_mz` | Pipeline never implemented; not decorated, internal-only helper |
| `msp_functions.py` | `deconvolution_correction` | Pipeline never implemented; not decorated, internal-only helper |
| `msp_functions.py` | `rga_status_process` | Pipeline never implemented; not decorated, internal-only helper |
| `msp_functions.py` | `SamplePreProcess` | Pipeline never implemented; not decorated, internal-only helper |
| `msp_functions.py` | `BackgroundPreProcess` | Pipeline never implemented; not decorated, internal-only helper |
| `msp_functions.py` | `Cal1PreProcess` | Pipeline never implemented; not decorated, internal-only helper |
| `msp_functions.py` | `Cal2PreProcess` | Pipeline never implemented; not decorated, internal-only helper |
| `met_functions.py` | `seasurface_skintemp_correct` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm) |
| `met_functions.py` | `warmlayer` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm) |
| `met_functions.py` | `coare35vn` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm) |
| `met_functions.py` | `met_buoyfls` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm) |
| `met_functions.py` | `met_buoyflx` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm) |
| `met_functions.py` | `met_frshflx` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm) |
| `met_functions.py` | `met_heatflx` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm) |
| `met_functions.py` | `met_latnflx` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm) |
| `met_functions.py` | `met_mommflx` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm) |
| `met_functions.py` | `met_netlirr` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm) |
| `met_functions.py` | `met_rainflx` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm) |
| `met_functions.py` | `met_sensflx` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm) |
| `met_functions.py` | `met_sphum2m` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm) |
| `met_functions.py` | `met_stablty` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm) |
| `met_functions.py` | `met_tempa2m` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm) |
| `met_functions.py` | `met_tempskn` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm) |
| `met_functions.py` | `met_wind10m` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm) |
| `met_functions.py` | `met_heatflx_minute` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm) |
| `met_functions.py` | `met_latnflx_minute` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm) |
| `met_functions.py` | `met_netlirr_minute` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm) |
| `met_functions.py` | `met_sensflx_minute` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm) |
| `met_functions.py` | `met_timeflx` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm) |
| `met_functions.py` | `met_netsirr_hourly` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm) |
| `met_functions.py` | `met_rainrte` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm) |
| `met_functions.py` | `air_density` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm); not decorated, internal-only helper |
| `met_functions.py` | `airtemp_at_refheight` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm); not decorated, internal-only helper |
| `met_functions.py` | `calc_rain_rate` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm); not decorated, internal-only helper |
| `met_functions.py` | `gravity` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm); not decorated, internal-only helper |
| `met_functions.py` | `latent_heat_vaporization_pure_water` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm); not decorated, internal-only helper |
| `met_functions.py` | `net_longwave_up` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm); not decorated, internal-only helper |
| `met_functions.py` | `psit_26` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm); not decorated, internal-only helper |
| `met_functions.py` | `psiu_26` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm); not decorated, internal-only helper |
| `met_functions.py` | `rain_heat_flux` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm); not decorated, internal-only helper |
| `met_functions.py` | `sea_spechum` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm); not decorated, internal-only helper |
| `met_functions.py` | `spechum_at_refheight` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm); not decorated, internal-only helper |
| `met_functions.py` | `water_thermal_expansion` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm); not decorated, internal-only helper |
| `met_functions.py` | `windspeed_at_refheight` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm); not decorated, internal-only helper |
| `met_functions.py` | `charnock_wind` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm); not decorated, internal-only helper |
| `met_functions.py` | `coolskin_parameters` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm); not decorated, internal-only helper |
| `met_functions.py` | `effective_relwind` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm); not decorated, internal-only helper |
| `met_functions.py` | `obukhov_for_init` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm); not decorated, internal-only helper |
| `met_functions.py` | `obukhov_length_scale` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm); not decorated, internal-only helper |
| `met_functions.py` | `roughness_lengths` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm); not decorated, internal-only helper |
| `met_functions.py` | `roughness_lengths_for_init` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm); not decorated, internal-only helper |
| `met_functions.py` | `scaling_parameters` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm); not decorated, internal-only helper |
| `met_functions.py` | `condition_data` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm); not decorated, internal-only helper |
| `met_functions.py` | `make_hourly_data` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm); not decorated, internal-only helper |
| `met_functions.py` | `warmlayer_time_keys` | Superseded by the official TOGA-COARE implementation (github.com/NOAA-PSL/COARE-algorithm); not decorated, internal-only helper |
| `generic_functions.py` | `ntp_to_unix_time` | Superseded by the fixed 2208988800-second NTP/Unix epoch offset constant |
