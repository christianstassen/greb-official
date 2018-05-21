# README #

### What is this repository for? ###

* Quick summary

This repository contains the GREB code for the paper 'A Hydrological Cycle Model for the Globally Resolved Energy Balance Model (GREB)' by Stassen et al. 2018.

### How do I get set up? ###

* Summary of set up

 1. To download the code via github click on the 'clone or download' tab on the top.
 2. The downloaded zip-archive contains the model source code and all input files.
 3. The main file to execute is 'run/greb.web-public.hydro.com'.
 4. To execute the main file try '.run/greb.web-public.hydro.com' or 'sh run/greb.web-public.hydro.com'
 5. The main file also contains several settings (i.e. what experiment, how many years to simulate, etc)

* Configuration

 1. log_exp  -> Experiment identifier (30=ENSO, 40=climate change)
 2. irain    -> Precipitation parameterisation (0=best GREB, -1=orig. GREB, 1=rq only, 2=omega only, 3=rq and omega)
 3. ieva     -> Evaporation parameterisation (0=best GREB, -1=orig. GREB, 1=parameter fitting, 2=parameter fitting and external wind speed)
 4. iadv     -> Advection parameterisation (0=best GREB, -1=orig. GREB)
 5. iqflux   -> q-flux corrections (0=orig. GREB, 1=annual mean spec. humidity, 2=annual mean spec. hum. and surface temperature)
 6. iclim    -> external boundary conditions (0=NCEP (orig. GREB), 1=ERA-Interim)

* Dependencies

 1. External files:
     1. topography, soil moisture, solar radiation, ocean mix-layer depth, omega, omega standard deviation, surface temperature, surface humidity,
     meridional wind, zonal wind, wind speed, glacier mask, cloud cover
     2. For the experiments (i.e. climate change) additional external files are needed.
 2. Libraries:
     1. NetCDF library (if wanted)
     2. netcdf-python (if wanted)

* Compiling options

As a standard only one compiler option is used (-Ofast), which speeds up the model simulations.

### Additional scripts and tests
You can use the script 'plot_control.py' to plot the binary output with python.
This is only a simple script and considered as test to see if the model is running properly.

If you prefer netcdf you can run the python script bin2netcdf.py to transform all files in /output to netcdf.
This requires additional libraries for python (i.e. netcdf-python).

### Important
The model has only been tested on Linux (Ubuntu) and Mac Os (High Sierra). It might not be running on other operating systems and/or versions.
Good luck with running the model.

### Contribution guidelines ###

### Who do I talk to? ###

* Repo owner or admin: Christian Stassen (christian.stassen@monash.edu)

* Other community or team contact: Dietmar Dommenget (dietmar.dommenget@monash.edu)
