# README #

Source code for Rieger et al. 2017.

### Latest Software Release
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1117870.svg)](https://doi.org/10.5281/zenodo.1117870)


### Required Software
This code has been tested in MATLAB 2016b (9.1). Execution requires MATLAB's Statistics, Global Optimization, Optimization, and SimBiology toolboxes. Presently, there appears to be a minor problem with plotting for versions of MATLAB before 2016b.

### Significant Functions and Scripts
* **a_make_paper_figs.m** - function that reads the contents of txtout/ and makes figures similar to the manuscript figures.
* **a_run_vpop_fit.m** - function for generation of simulation results. Based on inputs, will cycle through methods and iterates to generate data for manuscript.

### Included Packages
The code makes use of three packages from MATLAB Central:
* allcomb: https://www.mathworks.com/matlabcentral/fileexchange/10064-allcomb-varargin-
* distributionPlot: https://www.mathworks.com/matlabcentral/fileexchange/23661-violin-plots-for-plotting-multiple-distributions--distributionplot-m-
* error_ellipse: http://www.mathworks.com/matlabcentral/fileexchange/4705-error-ellipse

Since all of the packages' licenses permit redistribution with attribution, they are included as subdirectories with this repository. The user should be aware of the individual licenses, credits, and copyrights, which are independent from vpop-gen.

### Contact
* Ted Rieger: ted.rieger@pfizer.com

![alt text](https://github.com/openPfizer/DigitalHealthData/blob/master/img/osbypfizer.png)
