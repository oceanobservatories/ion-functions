# HYD Functions

## Background

The Ocean Observatories Initiative deploys two classes of passive acoustic
hydrophone instruments as part of the cabled portion of the Coastal Endurance 
Array and elsewhere in the Regional Cabled Array (RCA) to record underwater 
sound pressure waves in the frequency ranges relevant to seismic, biological, 
and anthropogenic sources.

| Class | Hardware | Platform | Designator meaning |
|---|---|---|---|
| HYDBB | Ocean Sonics icListen HF | Moored (fixed depth) | Hydrophone, Broadband |
| HYDLF | Low-frequency hydrophone | Moored (fixed depth) | Hydrophone, Low Frequency |

`hyd_functions.py` converts raw L0 data from both instrument classes into
their respective L1 data products: Broadband Acoustic Pressure Waves
(HYDAPBB_L1) and Low Frequency Acoustic Pressure Waves (HYDAPLF_L1). All
calibration coefficients are from factory calibration sheets supplied with
individual instruments. Within the OOI data system, HYD functions fall
under the Water Column and Seafloor/Crust science regimes and the Passive
Acoustics category.

### Primary Sources

| DCN | Document |
|---|---|
| 1341-00820 | OOI (2013). Data Product Specification for Broadband Acoustic Pressure Waves. (Not released.) |
| 1341-00821 | [OOI (2013). Data Product Specification for Low Frequency Acoustic Pressure Waves.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00821_Data_Product_SPEC_HYDAPLF_OOI.pdf) |

---

### HYDAPBB_L1 — Broadband Acoustic Pressure Waves

The HYDBB instrument digitizes the hydrophone signal into a 24-bit WAV
data stream at up to 200 kHz. Data are streamed over a 10/100BaseT Ethernet
connection using the icListen RIFF/WAV format with a full-scale range of
$\pm 3\ \text{V}$.

The L0 input (`wav`) arrives as a normalized floating-point value in the
range $[-1, 1]$ as encoded in the WAV format. The L1 product is computed in
two steps.

**Step 1 — Recover voltage from WAV scaling:**

$$V(t) = \text{wav}(t) \times 3$$

which maps the normalized WAV values back to the instrument's $\pm 3\ \text{V}$
full-scale range.

**Step 2 — Remove external gain:**

The instrument header encodes the external gain $G$ in dBV. The gain is
converted from dBV to a linear scale factor and divided out:

$$G_\text{linear} = 10^{G/20}$$

$$S(t) = \frac{V(t)}{G_\text{linear}}$$

where $S(t)$ is the gain-compensated time-series voltage (HYDAPBB_L1) in
Volts.

The OCVR (Open Circuit Voltage Response, or Sensitivity) of the icListen HF
is a function of frequency and cannot be applied to a broadband time-series
without frequency-domain filtering. Therefore, the L1 product is
gain-compensated voltage rather than calibrated pressure. Signal levels are
accurate to within 1 dB re $1\ \mu\text{Pa/V}$.

HYDBB instruments are deployed on the RSN system at Hydrate Ridge (PN1A,
PN1B) and Axial (PN3A, PN3B) Hybrid Mooring subsites, and on the Endurance
Array Hybrid Mooring subsites PN1C and PN1D.

---

### HYDAPLF_L1 — Low Frequency Acoustic Pressure Waves

The HYDLF hydrophone is attached electrically and physically to the co-located
OBSBB or OBSBK broadband seismometer. Its analog output is digitized by the
Guralp DM24 digitizer at up to 1000 samples per second at 24-bit depth. The
digitized signal is transmitted via SEEDlink protocol as SEED blockettes,
routed through a US Navy data diversion switch, and made available to OOI via
an Antelope Orbserver export.

The L0 input is raw digital counts. The L1 product converts counts to Volts
using the DM24 bit weight (gain), which has units of $\mu\text{V/count}$:

$$S(t) = x(t) \times G \times 10^{-6}$$

where $x(t)$ is the raw count time-series (HYDAPLF_L0), $G$ is the Guralp
DM24 fixed bit weight in $\mu\text{V/count}$ (default 3.2), and the
$10^{-6}$ factor converts $\mu\text{V}$ to $\text{V}$. $S(t)$ is the L1
time-series (HYDAPLF_L1) in Volts.

As with HYDAPBB, the frequency-dependent OCVR of the HYDLF hydrophone cannot
be applied to a broadband time-series without prior frequency-domain
filtering, so the L1 product is gain-compensated voltage. Resolution is
better than 0.1 dB re $1\ \mu\text{Pa}$; accuracy is $\pm 1\ \text{dB}$
re $1\ \mu\text{Pa}$.

The acoustic signal is digitized at the same rate and with the same
time-stamp as the co-located seismic channels, facilitating separation of
the waterborne pressure signal from the solid-earth seismic signal.

HYDLF instruments are deployed on the RSN system at Hydrate Ridge (PN1A,
PN1B) and Axial (PN3A, PN3B) subsites.

Full algorithm derivations, calibration procedures, and source references
are listed in the [References](#references) section.

---

## Core functions

::: ion_functions.data.hyd_functions.hyd_bb_acoustic_pwaves

#### Additional Notes
The WAV format encodes audio samples as normalized floating-point values
in the range $[-1, 1]$. The multiply-by-3 step in the code recovers the
physical voltage by scaling to the icListen HF full-scale range of
$\pm 3\ \text{V}$. The external gain in dBV is taken from the instrument
event header; if no external gain is applied the gain value is 0 dB
(linear factor of 1).

#### History
| Date | Author | Change |
|---|---|---|
| 2014-05-16 | Christopher Wingard | Initial code |
| 2023-08-15 | Samuel Dahlberg | Removed use of numexpr; renamed local variables |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.hyd_functions.hyd_lf_acoustic_pwaves

#### History
| Date | Author | Change |
|---|---|---|
| 2014-07-09 | Christopher Wingard | Initial code |
| 2023-08-15 | Samuel Dahlberg | Removed use of numexpr |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## References

Urick, R. J. (1983). Principles of Underwater Sound, 3rd ed., pp. 44-53.
McGraw-Hill, New York.

OOI (2013). Data Product Specification for Broadband Acoustic Pressure
Waves. Document Control Number 1341-00820. (Not released.)

[OOI (2013). Data Product Specification for Low Frequency Acoustic
Pressure Waves. Document Control Number 1341-00821.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00821_Data_Product_SPEC_HYDAPLF_OOI.pdf)
