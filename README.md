# Crystal orientation fabric anisotropy causes directional hardening of the Northeast Greenland Ice Stream
## Specfab COF evolution along EGRIP flow line

These scripts were used to simulate the COF evolution along a flowline starting at EGRIP in the Northeast Greenland Ice Stream (NEGIS) which are published 
in Gerber et al., (2023), Crystal orientation fabric anisotropy causes directional hardening of the Northeast Greenland Ice Stream, Nature Communications

## Script overview

1)  calculate strain-rate and spin tensor from satellite velocities:

`extract_strainrate_spin.m` 
- calculates strain-rate and spin components from surface velocities **downstream_flowline_rawdata.shp** and rotates it into the flowframe
- saves output in **strainrate_spin.csv**

2) Iteratetively update fabric tensor along flow line:

`EGRIP_fabricEvolution.py`
- read the strain and spin tensor data from **strainrate_spin.csv**
- initialize the fabric at EGRIP by loading **nlm_EGRIP_n.csv**
- calculate COF evolution with `sf.dndt_LATROT(nlm_prev, D, W)` and update fabric from previous step
- save output in **eigenvalues_downstream.csv** and COF plots in **images/fabric_stepXXXX.png**
