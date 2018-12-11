#!/bin/csh
#
# run script for deconstruct 2xCO2 experiments with the Globally Resolved Energy Balance (GREB) Model
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
set LOG_TOPO  = 1	# topography
set LOG_CLOUD = 1	# clouds
set LOG_HUMID = 1	# humidity
set LOG_HDIF  = 1	# heat diffusion
set LOG_HADV  = 1	# heat advection
set LOG_ICE   = 1	# ice albedo feedback
set LOG_OCEAN = 1	# deep ocean heat uptake
set LOG_HYDRO = 1	# hydrology
set LOG_VDIF  = 1	# vapour diffusion
set LOG_VADV  = 1	# vapour advection

# Hydro cycle decomp
set LOG_EVA        = 2 # 0=no-eva;   1=greb-eva; 2=forced-eva
set LOG_OMEGA_EXT  = 1 # 0=no-omega; 1=??? ;     2=forced-omega
set LOG_HWIND_EXT  = 1 # 0 and 1 would be set in LOG_DIFF and LOG_ADV; 2=forced-crcl
set LOG_TSURF_EXT  = 1

# length of sensitivity experiment in years
set YEARS=50


### compile GREB model (uncomment one of these three options)
### gfortran compiler (Linux (e.g. Ubuntu), Unix or MacBook Air)
#gfortran -fopenmp -march=native -O3 -ffast-math -funroll-loops greb.model.mscm.f90 greb.shell.mscm.f90 -o greb.x
gfortran -Ofast -ffast-math -funroll-loops -fopenmp greb.model.mscm.f90 greb.shell.mscm.f90 -o greb.x # Christian
### ifortran compiler (Mac)
# ifort -assume byterecl -O3 -xhost -align all -fno-alias greb.model.mscm.f90 greb.shell.mscm.f90 -o greb.x
### g95 compiler (other Linux)
# g95 greb.model.mscm.f90 greb.shell.mscm.f90 -o greb.x

###################
# END USER INPUT! #
###################

setenv OMP_NUM_THREADS 2
setenv KMP_AFFINITY verbose,none

set SCENARIO='climatechange'
set NUMBER=${LOG_OCEAN}${LOG_VADV}${LOG_VDIF}${LOG_HYDRO}${LOG_ICE}${LOG_HADV}${LOG_HDIF}${LOG_HUMID}${LOG_CLOUD}${LOG_TOPO}
set HYDRO_NUMBER=${LOG_TSURF_EXT}${LOG_EVA}${LOG_OMEGA_EXT}${LOG_HWIND_EXT}
set FILENAME=${NUMBER}HYDRO${HYDRO_NUMBER}.${SCENARIO}
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
time_ctrl = 3 		! length of control run [yrs]
time_scnr = $YEARS ! length of scenario run [yrs]
/
&PHYSICS
log_exp = 230 		! deconstruct climate change
log_topo_drsp   = $LOG_TOPO
log_cloud_drsp  = $LOG_CLOUD
log_humid_drsp  = $LOG_HUMID
log_ocean_drsp  = $LOG_OCEAN
log_hydro_drsp  = $LOG_HYDRO
log_ice    	= $LOG_ICE
log_hdif   	= $LOG_HDIF
log_hadv   	= $LOG_HADV
log_vdif   	= $LOG_VDIF
log_vadv   	= $LOG_VADV

log_eva     = $LOG_EVA
log_tsurf_ext = $LOG_TSURF_EXT
log_hwind_ext = $LOG_HWIND_EXT
log_omega_ext = $LOG_OMEGA_EXT
/
EOF

# run model
chmod +x greb.x
./greb.x

# postprocessing
# create output directory if does not already exist
if (! -d ../output ) mkdir ../output

# rename control run output and move it to output folder
mv control.bin ../output/control.exp-${FILENAME}.bin
mv scenario.bin ../output/response.exp-${FILENAME}.bin
mv scenario.gmean.bin ../output/response.gmean.exp-${FILENAME}.bin

# calculate months of scenario run for header file
@ MONTHS = $YEARS * 12

# control run
cat >../output/control.exp-${FILENAME}.ctl <<EOF
dset ^control.exp-${FILENAME}.bin
undef 9.e27
xdef  96 linear 0 3.75
ydef  48 linear -88.125 3.75
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

# scenario run
cat >../output/response.exp-${FILENAME}.ctl <<EOF
dset ^response.exp-${FILENAME}.bin
undef 9.e27
xdef  96 linear 0 3.75
ydef  48 linear -88.125 3.75
tdef $MONTHS linear 15jan0  1mo
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

cat >../output/response.gmean.exp-${FILENAME}.ctl <<EOF
dset ^response.gmean.exp-${FILENAME}.bin
undef 9.e27
xdef 12 linear 0 3.75
ydef  1 linear -88.125 3.75
zdef  $YEARS linear 1 1
tdef  1 linear 15jan0  1yr
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
