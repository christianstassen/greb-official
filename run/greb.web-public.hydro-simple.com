#
#----------------------------------------------------------
#   The Globally Resolved Energy Balance (GREB) Model
#----------------------------------------------------------
#
#   Authors; Dietmar Dommenget & Christian Stassen
#
#   Reference: - Conceptual Understanding of Climate Change with a Globally Resolved Energy Balance Model
#                by Dietmar Dommenget and Janine Flöter, submitted to Climate Dynamics 2010.
#              - A Hydrological Cycle Model for the Globally Resolved Energy Balance Model (GREB)
#                by Christian Stassen & Dietmar Dommenget
#
#   Describtion: This routine runs 50 year long A1b climate change scenario.
#
#

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# echo New GREB model (Climate Change A1B)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
log_exp=44
irain=0
ieva=0
iadv=0
iqflux=0
iclim=1


# compile GREB model
gfortran -Ofast ../src/greb.model.f90 ../src/greb.shell.web-public.f90 -o ../src/greb.x

#  namelist
cat >namelist <<EOF
&NUMERICS
time_flux = 3  ! length of flux corrections run [yrs]
time_ctrl = 5  ! length of control run  [yrs]
time_scnr = 5  ! length of scenario run [yrs]
/
&PHYSICS
 log_exp  = $log_exp ! complete GREB model; 2xCO2 forcing
 log_rain = $irain
 log_eva  = $ieva
 log_qflux = $iqflux
 log_adv_scheme = $iadv
 log_clim = $iclim
/
EOF

# run model
../src/greb.x
