# README #

Source code for Rieger et al. 2017.

### Latest Software Release
[![DOI](https://zenodo.org/badge/105279206.svg)](https://zenodo.org/badge/latestdoi/105279206)

### Required Software
This code has been tested in MATLAB 2016b (9.1). Execution requires MATLAB's Statistics, Global Optimization, Optimization, and SimBiology toolboxes.

### Significant Functions and Scripts
* **a_make_paper_figs.m** - function that reads the contents of txtout/ and makes figures similar to the manuscript figures.
* **a_run_vpop_fit.m** - main generation script. Will cycle through methods and iterates to generate data for manuscript.

### Included Packages
The code makes use of three packages from MATLAB Central:
* allcomb: https://www.mathworks.com/matlabcentral/fileexchange/10064-allcomb-varargin-
* distributionPlot: https://www.mathworks.com/matlabcentral/fileexchange/23661-violin-plots-for-plotting-multiple-distributions--distributionplot-m-
* error_ellipse: http://www.mathworks.com/matlabcentral/fileexchange/4705-error-ellipse

Since all of the packages' licenses permit redistribution with attribution, they are included as subdirectories with this repository, however, they have individual licenses and copyrights, which are independent from vpop-gen's license and copyright.

### Contact
* Ted Rieger: ted.rieger@pfizer.com
