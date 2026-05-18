# OBS Functions

## Background

The Ocean Observatories Initiative deploys two classes of ocean bottom
seismometer instruments on the Regional Cabled Array (RCA) to record
seismic signals at the seafloor: the Broadband Ocean Bottom Seismometer
(OBSBB) and the Short Period Ocean Bottom Seismometer (OBSSP). Both
classes are cabled, streaming data in real time to shore and to the
Incorporated Research Institutions for Seismology (IRIS).

`obs_functions.py` converts raw L0 digital counts from both instrument
classes into their respective L1 data products: Broadband Ground Velocity
(GRNDVEL_L1), Broadband Ground Acceleration (GRNDACC_L1), and Short Period
Ground Velocity (SGRDVEL_L1). All calibration coefficients are from factory
calibration sheets supplied with individual instruments. Within the OOI
data system, OBS functions fall under the Seafloor/Crust science regime
and the Ground Motion category.

### Primary Sources

| DCN | Document |
|---|---|
| 1341-00090 | [OOI (2013). Data Product Specification for Broadband Ground Velocity.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00090_Data_Product_SPEC_GRNDVEL_OOI.pdf) |
| 1341-00100 | [OOI (2013). Data Product Specification for Broadband Ground Acceleration.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00100_Data_Product_SPEC_GRNDACC_OOI.pdf) |
| 1341-00110 | [OOI (2013). Data Product Specification for Short Period Ground Velocity.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00110_Data_Product_SPEC_SGRDVEL_OOI.pdf) |

---

### OBSBB — Broadband Ocean Bottom Seismometer

OBSBB instruments provide characterization of seismicity and earthquake
activity along tectonic plate boundaries. The OOI array of broadband
seismometers is primarily deployed to detect earthquakes along the
subduction zone of the Oregon margin and at Axial Volcano. Because the
instruments are cabled and connected directly to shore, they provide
real-time earthquake detection.

The OBSBB hardware is the Guralp CMG-1T (360 s to 50 Hz) paired with a
CMG-5T strong motion accelerometer, a DM24/7-EAM digitizer/interface, and
a co-located Low Frequency Hydrophone (HTI-90-U). The instrument provides
an Ethernet (10/100) interface and is synchronized via pulse-per-second
timing. In sedimented areas (Slope Base, Southern Hydrate Ridge, and the
base of Axial Seamount), instruments are buried beneath the seafloor in
60 cm deep by 60 cm diameter caissons filled with silica beads to improve
seismic coupling. At the summit of Axial Seamount, where sediment is
absent, the instrument is installed on basement rock and surrounded by
gravel-filled bags.

The CMG-1T seismometer measures ground velocity; its analog output is
digitized by the DM24 ADC at up to 1000 samples per second at 24-bit
resolution and passed through a cascade of FIR filters to produce output
data streams at 200, 40, and 1 samples per second. The CMG-5T strong
motion accelerometer measures ground acceleration and produces output data
streams at 200 and 1 samples per second.

The digitized signals are formatted for transmission via SEEDlink protocol
as SEED blockettes, routed through a US Navy data diversion switch, and
made available to OOI via an Antelope Orbserver export.

---

### GRNDVEL_L1 — Broadband Ground Velocity

GRNDVEL_L1 is the broadband ground velocity (m/s) time-series measured by
the Guralp CMG-1T seismometer channel of the OBSBB instrument. The L0
input is the raw digital count output of the Guralp DM24 ADC. The L1
product is computed by applying the DM24 bit weight (gain) and the CMG-1T
sensor sensitivity:

$$GRNDVEL = raw \times \frac{gain \times 10^{-6}}{2 \times sensitivity}$$

where $raw$ is the L0 count time-series, $gain$ is the DM24 fixed bit
weight ($\mu$V/count, default 3.2), $10^{-6}$ converts $\mu$V to V, and
$sensitivity$ is the CMG-1T sensor sensitivity (V/(m/s), default 1500.0).
The result is the L1 ground velocity time-series in m/s.

---

### GRNDACC_L1 — Broadband Ground Acceleration

GRNDACC_L1 is the broadband ground acceleration (m/s$^2$) time-series
measured by the Guralp CMG-5T strong motion accelerometer channel of the
OBSBB instrument. The L0 input is the raw digital count output of the
Guralp DM24 ADC. The L1 product is computed by applying the DM24 bit
weight (gain) and the CMG-5T sensor sensitivity:

$$GRNDACC = raw \times \frac{gain \times 10^{-6}}{2 \times sensitivity}$$

where $raw$ is the L0 count time-series, $gain$ is the DM24 fixed bit
weight ($\mu$V/count, default 3.2), $10^{-6}$ converts $\mu$V to V, and
$sensitivity$ is the CMG-5T sensor sensitivity (V/(m/s$^2$), default
0.508). The result is the L1 ground acceleration time-series in m/s$^2$.

The GRNDACC and GRNDVEL signals are digitized at the same rate and with
the same time-stamp as the co-located HYDAPLF acoustic signal, facilitating
correlation across data products.

---

### OBSSP — Short Period Ocean Bottom Seismometer

OBSSP instruments detect vibrations from small, locally generated
earthquakes in the frequency range 0.1 Hz to 100 Hz. These events arise
from local phenomena such as melt movement beneath volcanoes and upward
flow of hydrothermal fluids in conduits feeding hydrothermal vent systems.
OBSSP instruments enable imaging of seismic energy traveling through the
seafloor.

OBSSP instruments are deployed on the RCA at Hydrate Ridge (PN1B, Hydrate
Summit) and Axial Volcano (PN3B, Axial Caldera). Each deployment includes
three orthogonal sensors: one oriented vertically and two oriented
horizontally. All instruments stream data in real time to IRIS.

The OBSSP digitized signal is formatted for transmission via SEEDlink
protocol, routed through the same US Navy data diversion switch used by
the OBSBB instruments, and made available to OOI via an Antelope Orbserver
export.

---

### SGRDVEL_L1 — Short Period Ground Velocity

SGRDVEL_L1 is the short period ground velocity (m/s) time-series measured
by the OBSSP instrument. The L0 input is the raw digital count output of
the Guralp DM24 ADC. The L1 product is computed using the same formula
structure as GRNDVEL:

$$SGRDVEL = raw \times \frac{gain \times 10^{-6}}{2 \times sensitivity}$$

where $raw$ is the L0 count time-series, $gain$ is the DM24 fixed bit
weight ($\mu$V/count, default 2.84), $10^{-6}$ converts $\mu$V to V, and
$sensitivity$ is the sensor sensitivity (V/(m/s), default 1200.0). The
result is the L1 short period ground velocity time-series in m/s.

The DM24 ADC digitizes the signal at up to 1000 samples per second at
24-bit resolution; FIR filters produce output data streams at 200, 40, and
1 samples per second. Output accuracy is 5% in amplitude and 0.1% in
precision across the instrument passband; time-stamping accuracy is 100
microseconds UTC.

Full algorithm derivations, calibration procedures, and source references
are listed in the [References](#references) section.

---

## Core Functions

::: ion_functions.data.obs_functions.obs_bb_ground_velocity

#### History
| Date | Author | Change |
|---|---|---|
| 2014-07-09 | Christopher Wingard | Initial code |
| 2023-08-15 | Samuel Dahlberg | Removed use of Numexpr library |
| 2025-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.obs_functions.obs_bb_ground_acceleration

#### History
| Date | Author | Change |
|---|---|---|
| 2014-07-09 | Christopher Wingard | Initial code |
| 2023-08-15 | Samuel Dahlberg | Removed use of Numexpr library |
| 2025-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.obs_functions.obs_sp_ground_velocity

#### History
| Date | Author | Change |
|---|---|---|
| 2014-07-09 | Christopher Wingard | Initial code |
| 2023-08-15 | Samuel Dahlberg | Removed use of Numexpr library |
| 2025-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## References

[OOI (2013). Data Product Specification for Broadband Ground Velocity.
Document Control Number 1341-00090.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00090_Data_Product_SPEC_GRNDVEL_OOI.pdf)

[OOI (2013). Data Product Specification for Broadband Ground Acceleration.
Document Control Number 1341-00100.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00100_Data_Product_SPEC_GRNDACC_OOI.pdf)

[OOI (2013). Data Product Specification for Short Period Ground Velocity.
Document Control Number 1341-00110.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00110_Data_Product_SPEC_SGRDVEL_OOI.pdf)

Scherbaum, F. (2007). *Of Poles and Zeros: Fundamentals of Digital
Seismology*, 2nd ed. Springer, Dordrecht.
