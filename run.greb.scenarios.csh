#!/bin/csh
#
# run script for scenario experiments with the Globally Resolved Energy Balance (GREB) Model
#
# author: Tobias Bayr and Dietmar Dommenget

# create work directory if does not already exist
if (! -d work ) mkdir work

# else clean up work directory
if (-d work ) rm -f work/*


# possible sensitivity experiments and suggested/maximum experiment length in years
#
#  EXP = 20  2xCO2					[ 50 years]
#  EXP = 21  4x   CO2					[ 50 years]
#  EXP = 22  10x  CO2					[ 50 years]
#  EXP = 23  0.5x CO2					[ 50 years]
#  EXP = 24  0x   CO2					[ 50 years]
#
#  EXP = 25  CO2-wave 30yrs-period			[100 years]
#  EXP = 26  2xCO2 30yrs followed by 70yrs CO2-ctrl	[100 years]
#  EXP = 27  solar constant +27W/m2 (~2xCO2 warming)	[ 50 years]
#  EXP = 28  11yrs solar cycle				[ 50 years]
#
#  EXP = 30  paleo solar 231 kyr BP & CO2=200ppm	[ 50 years]
#  EXP = 31  paleo solar 231 kyr BP 			[ 50 years]
#  EXP = 32  paleo CO2=200ppm 231 kyr BP 		[ 50 years]
#
#  EXP = 35  solar radiation obliquity changes		[ 50 years]
#  EXP = 36  solar radiation eccentricity changes	[ 50 years]
#  EXP = 37  solar radiation radius changes		[ 50 years]
#
#  EXP = 40  partial 2xCO2 Northern hemisphere 		[ 50 years]
#  EXP = 41  partial 2xCO2 Southern hemisphere 		[ 50 years]
#  EXP = 42  partial 2xCO2 Tropics 			[ 50 years]
#  EXP = 43  partial 2xCO2 Extratropics 		[ 50 years]
#  EXP = 44  partial 2xCO2 Ocean 			[ 50 years]
#  EXP = 45  partial 2xCO2 Land 			[ 50 years]
#  EXP = 46  partial 2xCO2 Boreal Winter 		[ 50 years]
#  EXP = 47  partial 2xCO2 Boreal Summer 		[ 50 years]
#
#  EXP = 95  IPCC A1B scenario				[150 years]
#  EXP = 96  IPCC RCP26 scenario			[550 years]
#  EXP = 97  IPCC RCP45	scenario			[550 years]
#  EXP = 98  IPCC RCP60	scenario			[550 years]
#  EXP = 99  IPCC RCP85	scenario			[550 years]
#
#  EXP = 100 run model with your own CO2 scenario
#
#
#  EXP = 230 run a climate change experiment with forced boundary conditions
#            (surface temperature, hodrizontal winds and omega) of the CMIP5
#            rcp85 ensemble mean response
#
#  EXP = 240 & 241 run a El Nino (La Nina) experiment with forced boundary conditions
#            (surface temperature, hodrizontal winds and omega) of the ERA-Interim
#            composite mean response
#
#
# some general remarks to the sensitivity experiments:
# - all scenarios will start in 1950
# - EXP 20-24 are abrupt climate change experiment, that will reach
#   the new equilibrium climate after ~50 years
# - EXP 26 is a climate change experiment with two abrupt changes
# - EXP 25,28 are climate change experiments where the boundary conditions
#   change contiously
# - EXP 27 the incoming solar rediation is increased
# - EXP 30-32 are paleo experiment with boundary conditions 231 kyr before present
# - EXP 35-37 are experiments with changed orbital parameters
# - EXP 40-47 are experiments where CO2 is doubled in parts of the globe or seasons
# - EXP 95-99 are IPCC scenarios, where the CO2 data is available
#   for the years given in brackets from 1950 onward
# - EXP 100: you can run the model with your own CO2 forcing,
#   which should be in the format [year co2] like for the IPCC scenarios
#   and variable "FILENAME" below should be the same name as your CO2 forcing file, but without '.txt'


#####################
# BEGIN USER INPUT! #
#####################

# settings for scenario
# scenario number from list above
set EXP=20
# length of sensitivity experiment in years
set YEARS=50

# for EXP = 35 choose here a value between -250 and 900 (with an increment of 25) for the obliquity:
# => possible range: [-250 (= -25deg),  900 (= +90deg)], todays value 225 (=22.5deg)
set OBL=0

# for EXP = 36 choose here a value between -30 and 30 (with an increment of 1) for the eccentricity:
# => possible range: [-30 (= -0.3), 30 (= +0.3)], todays value ~2 (=0.02)
set ECC=0

# for EXP = 37 give here the deviation of the earths radius around the sun in %
# suggested range [-20:+20], todays value 0
set DRAD=0

# for EXP='100', give here the name of input CO2 forcing data set without '.txt'
set CO2input=none

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

# move complied files to work directory
mv greb.x work/.
mv *.mod work/.

# change to work directory
cd work

# link solar forcing for paleo and orbital scenarios
set SOLSCEN='nosolfile'
touch nosolfile
set INDIR='../input/solar_forcing_scenarios/'
if ( $EXP == 30 ) set SOLSCEN=${INDIR}'greb.solar.231K_hybers.corrected.bin'
if ( $EXP == 31 ) set SOLSCEN=${INDIR}'greb.solar.231K_hybers.corrected.bin'
if ( $EXP == 35 ) set SOLSCEN=${INDIR}'greb.solar.obliquity.'${OBL}'.bin'
if ( $EXP == 36 ) set SOLSCEN=${INDIR}'greb.solar.eccentricity.'${ECC}'.bin'
# link solar forcing scenario
ln -s $SOLSCEN solar_scenario

# link CO2 forcing for IPCC RCP scenarios
set CO2='noco2file'
touch noco2file
if ( $EXP == 96 ) set CO2='../input/ipcc.scenario.rcp26.forcing.txt'
if ( $EXP == 97 ) set CO2='../input/ipcc.scenario.rcp45.forcing.txt'
if ( $EXP == 98 ) set CO2='../input/ipcc.scenario.rcp6.forcing.txt'
if ( $EXP == 99 ) set CO2='../input/ipcc.scenario.rcp85.forcing.txt'
if ( $EXP == 100 ) set CO2='../input/'${CO2input}'.txt'
# link CO2 forcing file
ln -s $CO2 co2forcing

#  generate namelist
cat >namelist <<EOF
&NUMERICS
time_flux = 3  		! length of flux corrections run [yrs]
time_ctrl = 3 		! length of control run [yrs]
time_scnr = $YEARS  	! length of scenario run [yrs]
/
&PHYSICS
 log_exp = $EXP 	! sensitivity run as set above
 dradius = $DRAD	! deviations from the earth radius around the sun in %
/
EOF

# run model
./greb.x

# postprocessing
# create output directory if does not already exist
if (! -d ../output ) mkdir ../output

# create filename
set FILENAME=exp-${EXP}
if ( $EXP == 20 ) set FILENAME=exp-${EXP}.2xCO2
if ( $EXP == 21 ) set FILENAME=exp-${EXP}.4xCO2
if ( $EXP == 22 ) set FILENAME=exp-${EXP}.10xCO2
if ( $EXP == 23 ) set FILENAME=exp-${EXP}.0.5xCO2
if ( $EXP == 24 ) set FILENAME=exp-${EXP}.0xCO2
if ( $EXP == 25 ) set FILENAME=exp-${EXP}.CO2-wave
if ( $EXP == 26 ) set FILENAME=exp-${EXP}.CO2-step-function
if ( $EXP == 27 ) set FILENAME=exp-${EXP}.dS0.27W
if ( $EXP == 28 ) set FILENAME=exp-${EXP}.solar-cycle
if ( $EXP == 30 ) set FILENAME=exp-${EXP}.solar.231K.CO2-200ppm
if ( $EXP == 31 ) set FILENAME=exp-${EXP}.solar.231K
if ( $EXP == 32 ) set FILENAME=exp-${EXP}.231K.CO2-200ppm
if ( $EXP == 35 ) set FILENAME=exp-${EXP}.obliquity.${OBL}
if ( $EXP == 36 ) set FILENAME=exp-${EXP}.eccentricity.${ECC}
if ( $EXP == 37 ) set FILENAME=exp-${EXP}.radius.${DRAD}
if ( $EXP == 40 ) set FILENAME=exp-${EXP}.partial.2xCO2.n-hemis
if ( $EXP == 41 ) set FILENAME=exp-${EXP}.partial.2xCO2.s-hemis
if ( $EXP == 42 ) set FILENAME=exp-${EXP}.partial.2xCO2.tropics
if ( $EXP == 43 ) set FILENAME=exp-${EXP}.partial.2xCO2.extrop
if ( $EXP == 44 ) set FILENAME=exp-${EXP}.partial.2xCO2.ocean
if ( $EXP == 45 ) set FILENAME=exp-${EXP}.partial.2xCO2.land
if ( $EXP == 46 ) set FILENAME=exp-${EXP}.partial.2xCO2.winter
if ( $EXP == 47 ) set FILENAME=exp-${EXP}.partial.2xCO2.summer
if ( $EXP == 95 ) set FILENAME=exp-${EXP}.IPCC.A1B
if ( $EXP == 96 ) set FILENAME=exp-${EXP}.IPCC.RCP26
if ( $EXP == 97 ) set FILENAME=exp-${EXP}.IPCC.RCP45
if ( $EXP == 98 ) set FILENAME=exp-${EXP}.IPCC.RCP60
if ( $EXP == 99 ) set FILENAME=exp-${EXP}.IPCC.RCP85
if ( $EXP == 100 ) set FILENAME=exp-${EXP}.${CO2input}

# rename scenario run output and move it to output folder
mv scenario.bin ../output/scenario.${FILENAME}.bin
mv scenario.gmean.bin ../output/scenario.gmean.${FILENAME}.bin

# calculate months of scenario run for header file
@ MONTHS = $YEARS * 12

# scenario run
cat >../output/scenario.${FILENAME}.ctl <<EOF
dset ^scenario.${FILENAME}.bin
undef 9.e27
xdef  96 linear 0 3.75
ydef  48 linear -88.125 3.75
zdef   1 linear 1 1
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

cat >../output/scenario.gmean.${FILENAME}.ctl <<EOF
dset ^scenario.gmean.${FILENAME}.bin
undef 9.e27
xdef 12 linear 0 3.75
ydef  1 linear -88.125 3.75
zdef  $YEARS linear 1 1
tdef  1 linear 15jan0  1mo
vars 8
tsurf  $YEARS 0 tsurf
tatmos $YEARS 0 tatmos
tocean $YEARS 0 tocean
vapor  $YEARS 0 vapour
ice    $YEARS 0 ice
precip $YEARS 0 precip
eva $YEARS 0 eva
qcrcl $YEARS 0 qcrcl
endvars
EOF

exit
