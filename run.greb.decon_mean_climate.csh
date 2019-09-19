#!/bin/csh
#
# run script for deconstruct mean climate experiments with the Globally Resolved Energy Balance (GREB) Model
#
# author: Tobias Bayr and Dietmar Dommenget

# create work directory if does not already exist
if (! -d work ) mkdir work

# else clean up work directory
if (-d work ) rm -f work/*



#####################
# BEGIN USER INPUT! #
#####################

# switches to turn on (1) or off (0) the different processes
# see mscm.dkrz.de => deconstruct 2xco2 response for details
set LOG_ICE   = 1	# ice
set LOG_CLOUD = 1	# clouds
set LOG_OCEAN = 1	# ocean
set LOG_ATMOS = 1	# atmosphere
set LOG_HDIF  = 1	# heat diffusion
set LOG_HADV  = 1	# heat advection
set LOG_CO2   = 1	# CO2
set LOG_HYDRO = 1	# hydrology
set LOG_VDIF  = 1	# vapour diffusion
set LOG_VADV  = 1	# vapour advection
set LOG_QFLUX = 1	# model correction

### compile GREB model (uncomment one of these three options)
### gfortran compiler (Linux (e.g. Ubuntu), Unix or MacBook Air)
#gfortran -fopenmp -march=native -O3 -ffast-math -funroll-loops greb.model.mscm.f90 greb.shell.mscm.f90 -o greb.x
gfortran -Ofast -ffast-math -funroll-loops -fopenmp greb.model.mscm.f90 greb.shell.mscm.f90 -o greb.x
### ifortran compiler (Mac)
# ifort -assume byterecl -O3 -xhost -align all -fno-alias greb.model.mscm.f90 greb.shell.mscm.f90 -o greb.x
### g95 compiler (other Linux)
# g95 greb.model.mscm.f90 greb.shell.mscm.f90 -o greb.x

###################
# END USER INPUT! #
###################

setenv OMP_NUM_THREADS 2
setenv KMP_AFFINITY verbose,none

set SCENARIO='greb.mean.decon.exp-'
set NUMBER=${LOG_QFLUX}${LOG_ICE}${LOG_CLOUD}${LOG_VADV}${LOG_VDIF}${LOG_HYDRO}${LOG_OCEAN}${LOG_CO2}${LOG_HADV}${LOG_HDIF}${LOG_ATMOS}
set FILENAME=${SCENARIO}${NUMBER}
echo 'EXPERIMENT: '${FILENAME}

# move complied files to work directory
mv greb.x work/.
mv *.mod work/.

# change to work directory
cd work

#  generate namelist
cat >namelist <<EOF
&NUMERICS
time_flux = 3  		! length of flux corrections run [yrs]
time_ctrl = 50 		! length of control run [yrs]
/
&PHYSICS
log_exp = 1 		! deconstruct mean state
log_cloud_dmc   = $LOG_CLOUD
log_ocean_dmc   = $LOG_OCEAN
log_atmos_dmc   = $LOG_ATMOS
log_co2_dmc     = $LOG_CO2
log_hydro_dmc   = $LOG_HYDRO
log_qflux_dmc   = $LOG_QFLUX
log_ice    	= $LOG_ICE
log_hdif   	= $LOG_HDIF
log_hadv   	= $LOG_HADV
log_vdif   	= $LOG_VDIF
log_vadv   	= $LOG_VADV

log_conv    = -1
log_rain    = -1
log_eva     = -1
log_clim    = -1
/
EOF

# run model
./greb.x

# postprocessing
# create output directory if does not already exist
if (! -d ../output ) mkdir ../output
# rename control run output and move it to output folder
mv control.bin ../output/${FILENAME}.bin

# create description file
cat >../output/${FILENAME}.ctl <<EOF
dset ^${FILENAME}.bin
undef 9.e27
xdef  96 linear 0 3.75
ydef  48 linear -88.125 3.75
zdef   1 linear 1 1
tdef 12 linear 15jan0  1mo
vars 8
tsurf  1 0 tsurf
tatmos 1 0 tatmos
tocean 1 0 tocean
vapor  1 0 vapour
ice    1 0 ice
precip 1 0 precip
eva 1 0 eva
qcrcl 1 0 qcrcl
endvars
EOF

exit
