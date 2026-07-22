# WAV Functions

## Background

The WAV instrument family consists of the TRIAXYS OEM Directional Wave Sensor
(WAVSS), manufactured by AXYS Technologies, Inc. The WAVSS measures
three-dimensional acceleration, rotation rate, and orientation in the Earth's
magnetic field. From these measurements the instrument reconstructs the motion
of the sea surface in the wave field. All wave statistics computations are
performed internally by the instrument's proprietary firmware; ion-functions
does not re-implement those algorithms.

The WAVSS acquires accelerations over a sampling interval (typically 20 minutes) 
and produces a single set of wave property values representing that interval, 
along with a time series of platform displacement over the course of the 
acquisition.

Wave statistics are also measured by Teledyne RDI acoustic Doppler current
profilers (ADCPs) mounted on the seafloor at the shallowest OOI Endurance
and Pioneer Mid-Atlantic Bight (MAB) locations. Those instruments collect raw 
ping-level data which RDI's WavesMon software processes into derived wave 
statistics. This measurement method is used at those sites because the moorings 
are shallow enough for a bottom-mounted ADCP to accurately measure wave parameters. 
The WavesMon outputs map to the same WAVSTAT product namespace as the WAVSS outputs; 
ion-functions is not involved in the WavesMon processing pipeline.

### Primary Sources

| DCN | Document |
|---|---|
| [1341-00450](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00450_Data_Product_SPEC_WAVSTAT_OOI.pdf) | OOI (2012). Data Product Specification for Wave Statistics. Document Control Number 1341-00450. |

### WAVSTAT Data Products

The WAVSS and WavesMon pipelines together produce a set of wave statistics
under the collective product name WAVSTAT. Most WAVSS products are extracted
directly from the instrument's NMEA output sentences without further
computation by ion-functions. WavesMon products are computed by RDI's
WavesMon software from raw ADCP data and are not processed by
ion-functions. The table below lists all WAVSTAT sub-products, their units,
and their source(s).

<table style="table-layout:fixed;width:100%;border-collapse:collapse;border-spacing:0">
<colgroup>
<col style="width:20%">
<col style="width:30%">
<col style="width:12%">
<col style="width:38%">
</colgroup>
<thead>
<tr>
<th style="white-space:normal;padding:6px 8px;border-bottom:2px solid #e1e4e5;text-align:left">Product ID</th>
<th style="white-space:normal;padding:6px 8px;border-bottom:2px solid #e1e4e5;text-align:left">Description</th>
<th style="white-space:normal;padding:6px 8px;border-bottom:2px solid #e1e4e5;text-align:left">Units</th>
<th style="white-space:normal;padding:6px 8px;border-bottom:2px solid #e1e4e5;text-align:left">Source</th>
</tr>
</thead>
<tbody>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-N0</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Number of zero crossings in displacement data (QC use)</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">--</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSS firmware; extracted from $TSPWA</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-HMAX</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Maximum wave height</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">m</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSS firmware; extracted from $TSPWA. WavesMon: Hmax</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-HAVG</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Mean significant wave height</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">m</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSS firmware; extracted from $TSPWA. WavesMon: Hm</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-TAVG</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Period associated with mean significant wave height</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">s</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSS firmware; extracted from $TSPWA. WavesMon: Tz</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-HSIG</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Significant wave height (average of highest 1/3)</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">m</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSS firmware; extracted from $TSPWA. WavesMon: H13</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-HM0</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Significant wave height from spectral moments</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">m</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSS firmware; extracted from $TSPWA. WavesMon: Hs</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-TSIG</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Significant wave period (average period of Hsig waves)</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">s</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSS firmware; extracted from $TSPWA. WavesMon: T13</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-H10</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Average height of highest 1/10 of waves</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">m</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSS firmware; extracted from $TSPWA. WavesMon: H1/10</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-T10</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Average period of H10 waves</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">s</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSS firmware; extracted from $TSPWA. WavesMon: T1/10</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-TP</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Peak wave period</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">s</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSS firmware; extracted from $TSPWA. WavesMon: Tp (see also WAVSTAT-TSIG, WAVSTAT-TP5)</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-TP5</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Peak wave period via Read method (alternative to TP)</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">s</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSS firmware; extracted from $TSPWA. WavesMon: Tp (see also WAVSTAT-TP, WAVSTAT-TSIG)</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-TMAX</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Maximum peak wave period</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">s</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WavesMon: Tmax</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-D_L0</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Mean wave direction (magnetic)</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">deg</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSS firmware; extracted from $TSPWA</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-D_L2</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Mean wave direction (true north)</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">deg</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5"><code>wav_triaxys_correct_mean_wave_direction</code>. WavesMon: Dp (reported to true north; no correction applied by ion-functions)</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-DMEAN</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Mean peak wave direction</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">deg</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WavesMon: Dmean</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-DS</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Mean directional spread of wave field</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">deg</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSS firmware; extracted from $TSPWA</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-DEPTH</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Water depth from pressure sensor</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">m</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WavesMon: Depth Water level</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-TP_SEA</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Peak period in sea region of power spectrum</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">s</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WavesMon: Tp_Sea</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-D_SEA</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Peak direction in sea region at peak period</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">deg</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WavesMon: Dp_Sea</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-HSIG_SEA</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Significant wave height in sea region of power spectrum</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">m</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WavesMon: Hs_Sea</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-TP_SWELL</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Peak period in swell region of power spectrum</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">s</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WavesMon: Tp_Swell</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-D_SWELL</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Peak swell direction at peak period in swell region</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">deg</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WavesMon: Dp_Swell</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-HSIG_SWELL</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Significant wave height in swell region of power spectrum</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">m</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WavesMon: Hs_Swell</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-PND</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Power spectral density, non-directional spectra</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">m<sup>2</sup> Hz<sup>-1</sup></td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSS firmware; extracted from $TSPNA</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-FND_L1</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Frequency values for non-directional spectral bins</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Hz</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5"><code>wav_triaxys_nondir_freq</code></td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-FDS_L1</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Frequency values for directional spectral bins</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Hz</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5"><code>wav_triaxys_dir_freq</code></td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-PDS</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Power spectral density, directional spectra</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">m<sup>2</sup> Hz<sup>-1</sup></td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSS firmware; extracted from $TSPMA. WavesMon: directional spectrum (mm<sup>2</sup> Hz<sup>-1</sup> per cycle, 90 directions x 128 frequencies)</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-DDS_L0</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Wave directions from directional spectra (magnetic)</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">deg</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSS firmware; extracted from $TSPMA</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-DDS_L2</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Wave directions from directional spectra (true north)</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">deg</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5"><code>wav_triaxys_correct_directional_wave_direction</code>. WavesMon: directional spectrum directions (reported to true north; no correction applied by ion-functions)</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-SDS</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Directional spread from directional spectra</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">deg</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSS firmware; extracted from $TSPMA</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-MOTX_L0</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Eastward buoy displacement (magnetic frame)</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">m</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSS firmware; extracted from $TSPHA</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-MOTX_L1</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Eastward buoy displacement (true east)</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">m</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5"><code>wav_triaxys_magcor_buoymotion_x</code></td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-MOTY_L0</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Northward buoy displacement (magnetic frame)</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">m</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSS firmware; extracted from $TSPHA</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-MOTY_L1</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Northward buoy displacement (true north)</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">m</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5"><code>wav_triaxys_magcor_buoymotion_y</code></td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-MOTZ</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Vertical buoy displacement</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">m</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSS firmware; extracted from $TSPHA</td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-MOTT_L1</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Time of each buoy displacement measurement</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">seconds since 1900-01-01</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5"><code>wav_triaxys_buoymotion_time</code></td></tr>
<tr><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WAVSTAT-VELPROF</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">Current magnitude at each depth level</td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">m s<sup>-1</sup></td><td style="white-space:normal;padding:6px 8px;border-bottom:1px solid #e1e4e5">WavesMon: depthlevel1..N magnitude</td></tr>
<tr><td style="white-space:normal;padding:6px 8px">WAVSTAT-DIRPROF</td><td style="white-space:normal;padding:6px 8px">Current direction at each depth level</td><td style="white-space:normal;padding:6px 8px">deg</td><td style="white-space:normal;padding:6px 8px">WavesMon: depthlevel1..N direction</td></tr>
</tbody>
</table>

### Frequency Vector Reconstruction

WAVSTAT-FND_L1 and WAVSTAT-FDS_L1 are the frequency axis vectors for the
non-directional and directional power spectral density arrays, respectively.
Rather than transmitting every frequency value, the instrument reports only
the number of bins, the initial frequency, and the uniform frequency spacing.
ion-functions reconstructs the full frequency vector from these three values.

For the non-directional spectrum, all data packets share the same number of
frequency bins, so the reconstruction is vectorized across packets. For the
directional spectrum, the number of active bins (nfreq_dir) can vary between
data packets as a function of ocean conditions -- it is always less than or
equal to the number of non-directional bins. The output array is therefore
sized to the (fixed) non-directional bin count, with fill values occupying
unused positions.

The reconstruction for both products follows:

$$f_i = f_0 + (i - 1) \times \Delta f, \quad i = 1, 2, \ldots, N$$

where $f_0$ is the initial frequency [Hz], $\Delta f$ is the frequency spacing
[Hz], and $N$ is the number of active bins for that data packet.

### Platform Motion Time Reconstruction

WAVSTAT-MOTT_L1 is the time axis corresponding to the WAVSTAT-MOTX,
WAVSTAT-MOTY, and WAVSTAT-MOTZ displacement time series. The instrument
reports the NTP timestamp of the data sentence, the elapsed time to the first
displacement sample, and the uniform sampling interval. ion-functions
reconstructs the absolute NTP time of each sample:

$$t_i = t_0 + t_{init} + (i - 1) \times \Delta t, \quad i = 1, 2, \ldots, N$$

where $t_0$ is the NTP timestamp from the $TSPHA sentence [s since 1900-01-01],
$t_{init}$ is the initial time offset to the first sample [s], $\Delta t$ is
the time spacing between samples [s], and $N$ is the number of displacement
measurements.

### Magnetic Declination Correction

The WAVSS contains an internal compass. Wave directions (WAVSTAT-D and
WAVSTAT-DDS) and buoy displacement components (WAVSTAT-MOTX and WAVSTAT-MOTY)
are reported in the magnetic reference frame. ion-functions corrects these to
true north using the IGRF-14 magnetic declination model via
`generic_functions.magnetic_declination`.

The direction correction applies a rotation modulo 360 deg:

$$D_{true} = (D_{mag} + \theta + 360) \bmod 360$$

where $\theta$ is the magnetic declination [deg] at the instrument's location
and time. The same formula applies element-wise to the WAVSTAT-DDS array;
fill-valued positions are preserved through the operation.

The buoy displacement correction rotates the (X, Y) coordinate pair by
$\theta$ using a 2-D rotation matrix:

$$\begin{align}
X_{true} &= X_{mag} \times \cos\theta + Y_{mag} \times \sin\theta \\
Y_{true} &= -X_{mag} \times \sin\theta + Y_{mag} \times \cos\theta
\end{align}$$

The WAVSS is a surface sensor; the depth parameter for the declination
calculation defaults to 0 (sea level).

Full algorithm derivations, processing references, and source documentation
are listed in the [References](#references) section.

---

## Core Functions

::: ion_functions.data.wav_functions.wav_triaxys_nondir_freq

#### History
| Date | Author | Change |
|---|---|---|
| 2014-04-03 | Russell Desiderio | Initial code |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.wav_functions.wav_triaxys_dir_freq

#### Additional Notes

The number of active directional frequency bins (nfreq_dir) can vary between
data packets as a function of measured ocean conditions at fixed instrument
settings. As a result, the size of the WAVSTAT-FDS_L1 output array, and
correspondingly of the WAVSTAT-PDS and WAVSTAT-SDS arrays and the
WAVSTAT-DDS_L2 product, is not fixed across packets. The output is sized to
the (fixed) non-directional bin count, with fill values in unused positions.

#### History
| Date | Author | Change |
|---|---|---|
| 2014-04-03 | Russell Desiderio | Initial code |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.wav_functions.wav_triaxys_buoymotion_time

#### History
| Date | Author | Change |
|---|---|---|
| 2014-04-07 | Russell Desiderio | Initial code |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.wav_functions.wav_triaxys_correct_mean_wave_direction

#### History
| Date | Author | Change |
|---|---|---|
| 2014-04-08 | Russell Desiderio | Initial code |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.wav_functions.wav_triaxys_correct_directional_wave_direction

#### Additional Notes

The variable-length behavior of nfreq_dir described in the Additional Notes
for `wav_triaxys_dir_freq` applies equally here: the WAVSTAT-DDS_L2 output
array is sized to the non-directional bin count, with fill values in unused
positions. The correction preserves fill-valued positions by converting them
to NaN before applying the rotation and restoring fills afterward.

#### History
| Date | Author | Change |
|---|---|---|
| 2014-04-09 | Russell Desiderio | Initial code |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Helper Functions

::: ion_functions.data.wav_functions.magnetic_correction_einsum

#### Additional Notes

The executable code in this function is identical to
`magnetic_correction_vctrzd` in `adcp_functions.py`. It was written
separately to handle the vectorized (i, j) case of one magnetic declination
per ensemble of (u, v) pairs without a for loop. The ADCP DPS citations in
the docstring reflect the origin of the rotation algorithm; the function
operates identically for (X, Y) displacement coordinates as for (u, v)
velocity components.

#### History
| Date | Author | Change |
|---|---|---|
| 2014-04-04 | Russell Desiderio | Initial code |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Wrapper Functions

::: ion_functions.data.wav_functions.wav_triaxys_magcor_buoymotion_x

#### History
| Date | Author | Change |
|---|---|---|
| 2014-04-10 | Russell Desiderio | Initial code |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.wav_functions.wav_triaxys_magcor_buoymotion_y

#### History
| Date | Author | Change |
|---|---|---|
| 2014-04-10 | Russell Desiderio | Initial code |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## References

[OOI (2012). Data Product Specification for Wave Statistics. Document Control
Number 1341-00450.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00450_Data_Product_SPEC_WAVSTAT_OOI.pdf)

[OOI (2020). Data Product Specification for Velocity Profile and Echo
Intensity. Document Control Number 1341-00750.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00750_Data_Product_SPEC_VELPROF_ECHOINT_OOI.pdf)

[OOI (2013). Data Product Specification for Turbulent Velocity Profile and
Echo Intensity. Document Control Number 1341-00760.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00760_Data_Product_SPEC_VELTURB_ECHOINT_OOI.pdf)

Alken, P., Thebault, E., Beggan, C.D., et al. (2021). International
Geomagnetic Reference Field: the thirteenth generation. *Earth Planets Space*,
73, 49. <https://doi.org/10.1186/s40623-020-01288-x>

Strom, K.M., and Reistad, H. (2024). ppigrf: Python package for computing
the International Geomagnetic Reference Field (IGRF).
<https://github.com/IAGA-VMOD/ppigrf>

AXYS Technologies Inc. *TRIAXYS OEM Directional Wave Sensor, User's Manual.*
Sidney, BC: AXYS Technologies Inc. May 2005.

Teledyne RD Instruments. *WavesMon v3.08 User's Guide.* P/N 957-6232-00.
Teledyne RD Instruments.

Teledyne RD Instruments. *Waves Primer: Wave Measurements and the RDI ADCP
Waves Array Technique.* Teledyne RD Instruments.
