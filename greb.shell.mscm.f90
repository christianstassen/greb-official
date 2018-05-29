program  greb_shell

! initialisation of the greb model

USE mo_numerics
USE mo_physics

! declare output fields
real, dimension(xdim,ydim,ndays_yr) :: Tc1, Ta1, q1, ap1
real, dimension(xdim,ydim,ndays_yr) :: Tc2, Ta2, q2, ap2

integer, dimension(ndays_yr)::  t = (/(i,i=1,ndays_yr)/) ! jday index

100 FORMAT('climate: ',F9.2, 5E12.4)

print*,'% start climate shell'

! open input files
open(10,file='namelist')
open(11,file='../input/ncep.tsurf.1948-2007.clim.bin',  	ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
open(12,file='../input/ncep.zonal_wind.850hpa.clim.bin',      ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
open(13,file='../input/ncep.meridional_wind.850hpa.clim.bin', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
open(14,file='../input/atmospheric_humidity.clim.bin',        ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
open(15,file='../input/isccp.cloud_cover.clim.bin',     	ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
open(16,file='../input/ncep.soil_moisture.clim.bin',   	ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
open(17,file='../input/Tocean.clim.bin',			ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
open(18,file='../input/woce.ocean_mixed_layer_depth.clim.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
open(19,file='../input/global.topography.bin',      		ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
open(20,file='../input/greb.glaciers.bin',   			ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
open(21,file='../input/solar_radiation.clim.bin', 		ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*ydim*nstep_yr)
open(22,file='../input/erainterim.windspeed.bin',      ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
open(23,file='../input/erainterim.omega.vertmean.bin',      ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
open(24,file='../input/erainterim.omega_std.vertmean.bin',      ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)

! read namelist
read(10,numerics)
read(10,physics)

! read climatologies
do n=1,nstep_yr
   read(11,rec=n) tclim(:,:,n)
   read(12,rec=n) uclim(:,:,n)
   read(13,rec=n) vclim(:,:,n)
   read(14,rec=n) qclim(:,:,n)
   read(15,rec=n) cldclim(:,:,n)
   read(16,rec=n) swetclim(:,:,n)
   read(17,rec=n) Toclim(:,:,n)
   read(18,rec=n) mldclim(:,:,n)
   read(22,rec=n) wsclim(:,:,n)
   read(23,rec=n) omegaclim(:,:,n)
   read(24,rec=n) omegastdclim(:,:,n)
end do

! read fix data
read(19,rec=1)  z_topo
read(20,rec=1)  glacier
read(21,rec=1)  sw_solar_ctrl

! read scenario solar forcing for paleo scenarios or oribital forcings
if ( log_exp .eq. 30 .or. log_exp .eq. 31 .or. log_exp .eq. 35 .or. log_exp .eq. 36 ) then
open(25,file='solar_scenario', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*ydim*nstep_yr)
read(25,rec=1)  sw_solar_scnr
end if

! open CO2 forcing file for IPCC RCP scenarios (CO2 is read in forcing subroutine)
if ( log_exp .ge. 96 .and. log_exp .le. 100 ) then
open(26,file='co2forcing')
end if

! open external forcing for climate change (ensemble mean) (it is read in forcing subroutine)
if ( log_exp .eq. 230 ) then
  open(31,file='../input/tsurf.response.cmip5.ensmean', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(32,file='../input/zonal.wind.response.cmip5.ensmean', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(33,file='../input/meridional.wind.response.cmip5.ensmean', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(34,file='../input/omega.response.cmip5.ensmean', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(35,file='../input/windspeed.response.cmip5.ensmean', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  do i=1,nstep_yr ! Read in the anomalies
    read(31,rec=i) Tclim_anom_cc(:,:,i)
    read(32,rec=i) uclim_anom_cc(:,:,i)
    read(33,rec=i) vclim_anom_cc(:,:,i)
    read(34,rec=i) omegaclim_anom_cc(:,:,i)
    read(35,rec=i) wsclim_anom_cc(:,:,i)
  end do
end if

! ENSO forcing
if ( log_exp .eq. 240 .or. log_exp .eq. 241 ) then
  ! open external forcing for El Nino (era-interim composite mean) (it is read in forcing subroutine)
  if ( log_exp .eq. 240 ) then
    open(41,file='../input/erainterim.tsurf.elnino.forcing.bin', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    open(42,file='../input/erainterim.zonal.wind.elnino.forcing.bin', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    open(43,file='../input/erainterim.meridional.wind.elnino.forcing.bin', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    open(44,file='../input/erainterim.omega.elnino.forcing.bin', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    open(45,file='../input/erainterim.windspeed.elnino.forcing.bin', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  ! open external forcing for La Nina (era-interim composite mean) (it is read in forcing subroutine)
  else if ( log_exp .eq. 241 ) then
    open(41,file='../input/erainterim.tsurf.lanina.forcing.bin', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    open(42,file='../input/erainterim.zonal.wind.lanina.forcing.bin', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    open(43,file='../input/erainterim.meridional.wind.lanina.forcing.bin', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    open(44,file='../input/erainterim.omega.lanina.forcing.bin', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    open(45,file='../input/erainterim.windspeed.lanina.forcing.bin', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  end if
  do i=1,nstep_yr ! Read in the anomalies
    read(41,rec=i) Tclim_anom_enso(:,:,i)
    read(42,rec=i) uclim_anom_enso(:,:,i)
    read(43,rec=i) vclim_anom_enso(:,:,i)
    read(44,rec=i) omegaclim_anom_enso(:,:,i)
    read(45,rec=i) wsclim_anom_enso(:,:,i)
  end do
end if ! ENSO forcing


! start greb_model run
print*,'% time flux/control/scenario: ', time_flux, time_ctrl, time_scnr
call greb_model

end program greb_shell
