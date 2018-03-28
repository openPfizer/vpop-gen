# VPop-Gen #

VPop-Gen is the MATLAB source code for Rieger et al. 2018. With this code, the user should be able to reproduce all the figures and results from the manuscript. This code is also a natural starting point for extending the methods described in the manuscript to other models, however, it has only been tested against the Van de Pas model to date.  

### Latest Software Release
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1117870.svg)](https://doi.org/10.5281/zenodo.1117870)

### Required Software
This code has been tested in MATLAB 2016b (9.1). Required MathWorks toolboxes:
* Statistics and Machine Learning;
* Global Optimization;
* Optimization;
* SimBiology.

For versions of MATLAB before 2016b, there is a problem with generating the final plots.

### Significant Functions
* **a_make_paper_figs.m** - function that reads the contents of txtout/ and makes figures similar to the manuscript figures.
* **a_run_vpop_fit.m** - function for generation of simulation results. Based on inputs, will cycle through methods and iterates to generate data for manuscript.

### Included Packages
Vpop-Gen uses of three packages from MATLAB Central:
* allcomb: https://www.mathworks.com/matlabcentral/fileexchange/10064-allcomb-varargin-
* distributionPlot: https://www.mathworks.com/matlabcentral/fileexchange/23661-violin-plots-for-plotting-multiple-distributions--distributionplot-m-
* error_ellipse: http://www.mathworks.com/matlabcentral/fileexchange/4705-error-ellipse

All packages are included as subdirectories with this repository.

![alt text](https://github.com/openPfizer/DigitalHealthData/blob/master/img/osbypfizer.png)
