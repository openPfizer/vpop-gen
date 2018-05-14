# VPop-Gen #
VPop-Gen is the MATLAB source code for Rieger et al. 2018. With this code, the user should be able to reproduce all the figures and results from the manuscript. This code is also a natural starting point for extending the methods described in the manuscript to other models, however, it has only been tested against the Van de Pas model to date.  

### Latest Software Release
[![DOI](https://zenodo.org/badge/105279206.svg)](https://zenodo.org/badge/latestdoi/105279206)

### Required Software
[MATLAB](https://www.mathworks.com) and several of its toolboxes are required to execute Vpop-Gen. Vpop-Gen was written and tested using MATLAB 2016b (9.1). Execution of the code requires these additional toolboxes:
* [Global Optimization](https://www.mathworks.com/products/global-optimization.html)
* [Optimization](https://www.mathworks.com/products/optimization.html)
* [SimBiology](https://www.mathworks.com/products/simbiology.html)
* [Statistics and Machine Learning](https://www.mathworks.com/products/statistics.html)

For versions of MATLAB before 2016b there is a problem with generating the final plots.

### Functions for Running Vpop-Gen
* `a_make_paper_figs.m` - function that reads the contents of txtout/ and makes figures similar to the manuscript figures.
* `a_run_vpop_fit.m` - function for generation of simulation results. Based on inputs, will cycle through methods and iterates to generate data for manuscript.

### Included Packages
Vpop-Gen uses three packages from [MATLAB Central](https://www.mathworks.com/matlabcentral/):
* [allcomb](https://www.mathworks.com/matlabcentral/fileexchange/10064-allcomb-varargin-)
* [distributionPlot](https://www.mathworks.com/matlabcentral/fileexchange/23661-violin-plots-for-plotting-multiple-distributions--distributionplot-m-i)
* [error_ellipse](http://www.mathworks.com/matlabcentral/fileexchange/4705-error-ellipse)

All packages are included as subdirectories with this repository.

![alt text](https://github.com/openPfizer/DigitalHealthData/blob/master/img/osbypfizer.png)
