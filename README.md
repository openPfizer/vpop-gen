# README #

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
The code makes use of three packages from MATLAB Central:
* allcomb: https://www.mathworks.com/matlabcentral/fileexchange/10064-allcomb-varargin-
* distributionPlot: https://www.mathworks.com/matlabcentral/fileexchange/23661-violin-plots-for-plotting-multiple-distributions--distributionplot-m-
* error_ellipse: http://www.mathworks.com/matlabcentral/fileexchange/4705-error-ellipse

![alt text](https://github.com/openPfizer/DigitalHealthData/blob/master/img/osbypfizer.png)
