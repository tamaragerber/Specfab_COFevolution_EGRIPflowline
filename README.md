# Crystal orientation fabric anisotropy causes directional hardening of the Northeast Greenland Ice Stream
## Specfab COF evolution along EGRIP flow line

These scripts were used to infer horizontal anisotropy from airborne radar crosspoints in the Northeast Greenland Ice Stream (NEGIS) which are published 
in Gerber et al., (2023), Crystal orientation fabric anisotropy causes directional hardening of the Northeast Greenland Ice Stream, Nature Communications

## Script overview

1)  calculate strain-rate and spin tensor from satellite velocities

`extract_strainrate_spin.m` 
- calculates strain-rate and spin components from surface velocities **downstream_flowline_rawdata.shp** and rotates it into the flowframe
- saves output in **strainrate_spin.csv**

2) Iteratetively update fabric tensor along flow line.

EGRIP_fabricEvolution.py, using nlm_EGRIP_n.csv from NEGIS_fabric_aicorr.py/example_plot

output --> eigenvalues_downstream.csv
images/fabric_stepXXXX.png
png2gif.py --> EGRIP_fabricEvolution.gif
plot Eigenvalues and cumulative strain along flowline with Eigenvalues2D.m 
