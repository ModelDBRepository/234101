Simulation software and model files used in 
Danner SM, Shevtsova NA, Frigon A, Rybak IA. Computational modeling of spinal circuits controlling limb coordination and gaits in quadrupeds. eLife 2017;6:e31050.



Source code needs to be compiled with c++11 extensions turned on. 
Boost header files need to be located in the header search paths.
libnetcdf_c++4 needs to be linked.
A makefile is supplied. Software was tested on macOS 10.13.1, Cygwin on Windows 10, and Ubuntu 16.04.

Files in directory 'models' specify the network structure and all parameters of the model.

usage: executable -f config_file [-o output_file] [-u name value] [-a alpha]
                                 [-U varname tstep value1 value2 [tstep2] [-V varname]]
 -f config_file: text file specifying the neural network model
 -o output_file: the path where the simulationresults should be written to. The file will
                  be written in CDF-4 file format.
 -u name value:  updates a variable (name) specified in the config_file to the value given
 -a alpha: sets alpha to a constant value (overrides the configuration file)
 -U varname tstep value1 value2 [tstep2]: sets variable varname to value1 and changes
                  it to value2 at time tstep, if tstep2 is specified, the variable will be
                  changed back to value1 at time tstep2
 -V varname: adds an additional variable to the update process of -U (has no effect if -U
                  is not specified.

Output will be written as a netcdf file that can be open natively with matlab or using the netcdf library in various languages. It contains a single matrix (under /data) first column is the time, all other columns are V of the neurons (same sequence as in model file) followed by the h parameter of the NaP current.

Three Matlab scripts are provided that perform simulations presented in the paper and plot the appropriate figures.

create_bfdiag.m  -  bifurcation diagramscreate_noise_diags.m - simulations with increased noisedynamicDriveChange.m - simulations with abrupt drive changes

Usage of the scripts is described in more detail in the comments within the files.