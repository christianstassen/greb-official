!
!----------------------------------------------------------
!   The Globally Resolved Energy Balance (GREB) Model
!----------------------------------------------------------
!
!   Authors: Dietmar Dommenget, Janine Flöter, Tobias Bayr and Christian Stassen
!
!   References:	- Conceptual Understanding of Climate Change with a
!              	Globally Resolved Energy Balance Model
!               by Dietmar Dommenget and Janine Flöter, J. Clim Dyn (2011) 37: 2143.
!               doi:10.1007/s00382-011-1026-0

!               - A hydrological cycle model for the Globally Resolved Energy Balance (GREB) model v1.0
!               by Christian Stassen, Dietmar Dommenget & Nicholas Loveday.
!               Geosci. Model Dev., 12, 425-440, https://doi.org/10.5194/gmd-12-425-2019, 2019.
!
!		            - The Monash Simple Climate Model Experiments: An interactive database
!		            of the mean climate, climate change and scenarios simulations
!		            by Dietmar Dommenget, Kerry Nice, Tobias Bayr, Dieter Kasang, Christian Stassen
!		            and Mike Rezny, submitted to Geoscientific Model Development
!
!
!  input fields: The GREB model needs the following fields to be specified before
!                the main subroutine greb_model is called:
!
!    Tclim(xdim,ydim,nstep_yr):   mean Tsurf                       [K]
!    uclim(xdim,ydim,nstep_yr):   mean zonal wind speed            [m/s]
!    vclim(xdim,ydim,nstep_yr):   mean meridional wind speed       [m/s]
!    qclim(xdim,ydim,nstep_yr):   mean atmospheric humidity        [kg/kg]
!  cldclim(xdim,ydim,nstep_yr):   total cloud cover                [0-1]
! swetclim(xdim,ydim,nstep_yr):   soil wetness, fraction of total  [0-1]
!   Toclim(xdim,ydim,nstep_yr):   mean deep ocean temperature      [K]
!  mldclim(xdim,ydim,nstep_yr):   mean ocean mixed layer depth     [m]
!   z_topo(xdim,ydim):            topography (<0 are ocean points) [m]
!  glacier(xdim,ydim):            glacier mask ( >0.5 are glacier points )
! sw_solar(ydim,nstep_yr):        24hrs mean solar radiation       [W/m^2]
!
!
! possible experiments:
!
!  log_exp = 1  deconstruct mean climate
!  		you can switch on or off climate components as given in the
!		module physics in 'deconstruct mean state switches' section
!
!  log_exp = 10  deconstruct 2xCO2 response
!  		 you can switch on or off climate components as given in the
!		 module physics in 'deconstruct 2xco2 switches' section
!
!  log_exp = 20  2x   CO2
!  log_exp = 21  4x   CO2
!  log_exp = 22  10x  CO2
!  log_exp = 23  0.5x CO2
!  log_exp = 24  0x   CO2
!
!  log_exp = 25  CO2-wave 30yrs-period
!  log_exp = 26  2xCO2 30yrs followed by 70yrs CO2-ctrl
!  log_exp = 27  solar constant +27W/m2 (~2xCO2 warming)
!  log_exp = 28  11yrs solar cycle
!
!  log_exp = 30  paleo solar 231Kyr before present & CO2=200ppm
!  log_exp = 31  paleo solar 231Kyr before present
!  log_exp = 32  paleo CO2=200ppm 231Kyr before present
!
!  log_exp = 35  solar radiation obliquity changes
!  log_exp = 36  solar radiation eccentricity changes
!  log_exp = 37  solar radiation radius changes
!
!  log_exp = 40  partial 2xCO2 Northern hemisphere
!  log_exp = 41  partial 2xCO2 Southern hemisphere
!  log_exp = 42  partial 2xCO2 Tropics
!  log_exp = 43  partial 2xCO2 Extratropics
!  log_exp = 44  partial 2xCO2 Ocean
!  log_exp = 45  partial 2xCO2 Land
!  log_exp = 46  partial 2xCO2 Boreal Winter
!  log_exp = 47  partial 2xCO2 Boreal Summer
!
!  log_exp = 95  IPCC A1B scenario
!  log_exp = 96  IPCC RCP26 scenario
!  log_exp = 97  IPCC RCP45 scenario
!  log_exp = 98  IPCC RCP60 scenario
!  log_exp = 99  IPCC RCP85 scenario
!
!  log_exp = 230 run a climate change experiment with forced boundary conditions
!            (surface temperature, hodrizontal winds and omega) of the CMIP5
!            rcp85 ensemble mean response
!
!  log_exp = 240 & 241 run a El Nino (La Nina) experiment with forced boundary conditions
!            (surface temperature, hodrizontal winds and omega) of the ERA-Interim
!            composite mean response
!
!  log_exp = 100 run model with your own CO2 input file
!
!+++++++++++++++++++++++++++++++++++++++
module mo_numerics
!+++++++++++++++++++++++++++++++++++++++

! numerical parameter
  integer, parameter :: xdim = 96, ydim = 48          ! field dimensions
  integer, parameter :: ndays_yr  = 365               ! number of days per year
  integer, parameter :: dt        = 12*3600           ! time step [s]
  integer, parameter :: dt_crcl   = 0.5*3600          ! time step circulation [s]
  integer, parameter :: ndt_days  = 24*3600/dt        ! number of timesteps per day
  integer, parameter :: nstep_yr  = ndays_yr*ndt_days ! number of timesteps per year
  integer            :: time_flux = 0                 ! length of integration for flux correction [yrs]
  integer            :: time_ctrl = 0                 ! length of integration for control run  [yrs]
  integer            :: time_scnr = 0                 ! length of integration for scenario run [yrs]
  integer            :: ipx       = 1                 ! points for diagonstic print outs
  integer            :: ipy       = 1                 ! points for diagonstic print outs
  integer, parameter, dimension(12) :: jday_mon = (/31,28,31,30,31,30,31,31,30,31,30,31/) ! days per
  real, parameter    :: dlon      = 360./xdim         ! linear increment in lon
  real, parameter    :: dlat      = 180./ydim         ! linear increment in lat

  integer            :: ireal     = 4         	      ! record length for IO (machine dependent)
! 						        ireal = 4 for Mac Book Pro and Ubuntu Linux

  namelist / numerics / time_flux, time_ctrl, time_scnr

end module mo_numerics


!+++++++++++++++++++++++++++++++++++++++
module mo_physics
!+++++++++++++++++++++++++++++++++++++++

  use mo_numerics
  integer  :: log_exp   = 0              ! process control logics for expiments (see header)
! deconstruct mean state (dmc) switches
  integer  :: log_cloud_dmc   = 1              ! process control clouds
  integer  :: log_ocean_dmc   = 1              ! process control ocean
  integer  :: log_atmos_dmc   = 1              ! process control Atmosphere
  integer  :: log_co2_dmc     = 1              ! process control CO2
  integer  :: log_hydro_dmc   = 1              ! process control hydro
  integer  :: log_qflux_dmc   = 1              ! process control qflux corrections
! deconstruct 2xco2 (drsp) switches
  integer  :: log_topo_drsp   = 1              ! process control for topo
  integer  :: log_cloud_drsp  = 1              ! process control for clouds
  integer  :: log_humid_drsp  = 1              ! process control for humidity clim
  integer  :: log_ocean_drsp  = 1              ! process control for ocean
  integer  :: log_hydro_drsp  = 1              ! process control for hydro
! switches that are the same for both deconstructions
  integer  :: log_ice         = 1              ! process control ice-albedo
  integer  :: log_hdif        = 1              ! process control Diffusion of heat
  integer  :: log_hadv        = 1              ! process control Advection of heat
  integer  :: log_vdif        = 1              ! process control Diffusion of vapor
  integer  :: log_vadv        = 1              ! process control Advection of vapor
! switches for the hydrological cycle
  integer  :: log_rain        = 0              ! process control precipitation parameterisation
  integer  :: log_eva         = 0              ! process control evaporation parameterisation
  integer  :: log_conv        = 0              ! process control advection parameterisation
  integer  :: log_clim        = 0              ! process control for reference climatology
! switches for external forcing files
  integer  :: log_tsurf_ext   = 0              ! process control evaporation parameterisation
  integer  :: log_hwind_ext   = 0              ! process control advection parameterisation
  integer  :: log_omega_ext   = 0              ! process control for reference climatology

! parameters for scenarios
  real     :: dradius   = 0.		 ! deviations from actual earth radius in %

! physical parameter (natural constants)
  parameter( pi        = 3.1416 )
  parameter( sig       = 5.6704e-8 )     ! stefan-boltzmann constant [W/m^2/K^4]
  parameter( rho_ocean = 999.1 )         ! density of water at T=15C [kg/m^2]
  parameter( rho_land  = 2600. )         ! density of solid rock [kg/m^2]
  parameter( rho_air   = 1.2 )           ! density of air at 20C at NN
  parameter( cp_ocean  = 4186. )         ! specific heat capacity of water at T=15C [J/kg/K]
  parameter( cp_land   = cp_ocean/4.5 )  ! specific heat capacity of dry land [J/kg/K]
  parameter( cp_air    = 1005. )         ! specific heat capacity of air      [J/kg/K]
  parameter( eps       = 1. )            ! emissivity for IR
  real :: S0_var       = 100.            ! variation of solar constant   [%]

! physical parameter (model values)
  parameter( d_ocean   = 50. )                     ! depth of ocean column [m]
  parameter( d_land    = 2. )                      ! depth of land column  [m]
  parameter( d_air     = 5000. )                   ! depth of air column   [m]
  parameter( cap_ocean = cp_ocean*rho_ocean )      ! heat capacity 1m ocean  [J/K/m^2]
  parameter( cap_land  = cp_land*rho_land*d_land ) ! heat capacity land   [J/K/m^2]
  parameter( cap_air   = cp_air*rho_air*d_air )    ! heat capacity air    [J/K/m^2]
  real :: ct_sens   = 22.5                         ! coupling for sensible heat
  real :: da_ice    = 0.25                         ! albedo diff for ice covered points
  real :: a_no_ice  = 0.1                          ! albedo for non-ice covered points
  real :: a_cloud   = 0.35                         ! albedo for clouds
  real :: Tl_ice1   = 273.15-10.                   ! temperature range of land snow-albedo feedback
  real :: Tl_ice2   = 273.15                       ! temperature range of land snow-albedo feedback
  real :: To_ice1   = 273.15-7.                    ! temperature range of ocean ice-albedo feedback
  real :: To_ice2   = 273.15-1.7                   ! temperature range of ocean ice-albedo feedback
  real :: co_turb   = 5.0                          ! turbolent mixing to deep ocean [W/K/m^2]
  real :: kappa     = 8e5                          ! atmos. diffusion coefficient [m^2/s]
  parameter( ce        = 2e-3  )                   ! laten heat transfer coefficient for ocean
  parameter( cq_latent = 2.257e6 )                 ! latent heat of condensation/evapoartion f water [J/kg]
  parameter( cq_rain   = -0.1/24./3600. )          ! decrease in air water vapor due to rain [1/s]
  parameter( z_air     = 8400. )                   ! scaling height atmos. heat, CO2
  parameter( z_vapor   = 5000. )                   ! scaling height atmos. water vapor diffusion
  parameter( grav      = 9.81  )                   ! gravitational acceleration [m/s^2]
  real :: r_qviwv   = 2.6736e3                     ! regres. factor between viwv and q_air  [kg/m^3]

  ! physical paramter (rainfall)
  real :: c_q, c_rq, c_omega, c_omegastd

! parameter emisivity
  real, dimension(10) :: p_emi = (/9.0721, 106.7252, 61.5562, 0.0179, 0.0028,     &
&                                             0.0570, 0.3462, 2.3406, 0.7032, 1.0662/)

! declare climate fields
  real, dimension(xdim,ydim)          ::  z_topo, glacier,z_ocean
  real, dimension(xdim,ydim,nstep_yr) ::  Tclim, uclim, vclim, omegaclim, omegastdclim, wsclim
  real, dimension(xdim,ydim,nstep_yr) ::  qclim, mldclim, Toclim, cldclim
  real, dimension(xdim,ydim,nstep_yr) ::  TF_correct, qF_correct, ToF_correct, swetclim, dTrad
  real, dimension(ydim,nstep_yr)      ::  sw_solar, sw_solar_ctrl, sw_solar_scnr
  real, dimension(xdim,ydim)          ::  co2_part      = 1.0
  real, dimension(xdim,ydim)          ::  co2_part_scn  = 1.0

! declare anomaly fields for enso and climate change
  real, dimension(xdim,ydim,nstep_yr) ::   Tclim_anom_enso     = 0.
  real, dimension(xdim,ydim,nstep_yr) ::   uclim_anom_enso     = 0.
  real, dimension(xdim,ydim,nstep_yr) ::   vclim_anom_enso     = 0.
  real, dimension(xdim,ydim,nstep_yr) ::   omegaclim_anom_enso = 0.
  real, dimension(xdim,ydim,nstep_yr) ::   wsclim_anom_enso    = 0.
  real, dimension(xdim,ydim,nstep_yr) ::   Tclim_anom_cc       = 0.
  real, dimension(xdim,ydim,nstep_yr) ::   uclim_anom_cc       = 0.
  real, dimension(xdim,ydim,nstep_yr) ::   vclim_anom_cc       = 0.
  real, dimension(xdim,ydim,nstep_yr) ::   omegaclim_anom_cc   = 0.
  real, dimension(xdim,ydim,nstep_yr) ::   wsclim_anom_cc      = 0.

! declare constant fields
  real, dimension(xdim,ydim)          ::  cap_surf
  integer jday, ityr

! Mike: declare some program constants
  real, dimension(xdim, ydim)         :: wz_air, wz_vapor
  real, dimension(xdim,ydim,nstep_yr) :: uclim_m, uclim_p
  real, dimension(xdim,ydim,nstep_yr) :: vclim_m, vclim_p

  real :: t0, t1, t2

  namelist / physics / log_exp, ct_sens, da_ice, a_no_ice, a_cloud, co_turb, kappa, 	&
&                      p_emi, Tl_ice1, Tl_ice2, To_ice1, To_ice2, r_qviwv,          	&
&		                   log_cloud_dmc, log_ocean_dmc, log_atmos_dmc, log_co2_dmc,        &
&                      log_hydro_dmc, log_qflux_dmc, 					                          &
&                      log_topo_drsp, log_cloud_drsp, log_humid_drsp, log_hydro_drsp,   &
&                      log_ocean_drsp, log_ice, log_hdif, log_hadv, log_vdif, log_vadv, &
& 		                 S0_var, dradius, log_rain, log_eva, log_conv, log_clim,          &
&                      log_tsurf_ext, log_hwind_ext, log_omega_ext

end module mo_physics

!+++++++++++++++++++++++++++++++++++++++
module mo_diagnostics
!+++++++++++++++++++++++++++++++++++++++

  USE mo_numerics,    ONLY: xdim, ydim

 ! declare diagnostic fields
  real, dimension(xdim,ydim)          :: Tsmn, Tamn, qmn, swmn, lwmn, qlatmn, qsensmn, &
&                                        ftmn, fqmn, icmn, Tomn

! declare output fields
  real, dimension(xdim,ydim)          :: Tmm, Tamm, Tomm, qmm, icmm, prmm, evamm, qcrclmm
  real, dimension(xdim,ydim,12)       :: Tmn_ctrl, Tamn_ctrl, Tomn_ctrl
  real, dimension(xdim,ydim,12)       :: qmn_ctrl, icmn_ctrl, prmn_ctrl, evamn_ctrl, qcrclmn_ctrl

end module mo_diagnostics

!+++++++++++++++++++++++++++++++++++++++
subroutine greb_model
!+++++++++++++++++++++++++++++++++++++++
!   climate model main loop

  use mo_numerics
  use mo_physics
  use mo_diagnostics

! declare temporary fields
  real, dimension(xdim,ydim) :: Ts0, Ts1, Ta0, Ta1, To0, To1, q0, q1,       &
&                               ts_ini, ta_ini, q_ini, to_ini

! open output files
  open(101,file='control.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(102,file='scenario.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(103,file='scenario.gmean.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal)

  dTrad = -0.16*Tclim -5. ! offset Tatmos-rad

! set ocean depth
  z_ocean=0
  do i=1,nstep_yr
     where(mldclim(:,:,i).gt.z_ocean) z_ocean = mldclim(:,:,i)
  end do
  z_ocean = 3.0*z_ocean

! decon mean state switch
  if (log_cloud_dmc ==  0) cldclim = 0.0
  if( log_hydro_dmc ==  0) qclim   = 0.0

! decon2xco2 switch
  if (log_topo_drsp  ==  0) where(z_topo > 1.) z_topo = 1.0     ! sens. exp. constant topo
  if (log_cloud_drsp ==  0) cldclim = 0.7                       ! sens. exp. constant cloud cover
  if (log_humid_drsp ==  0) qclim   = 0.0052                    ! sens. exp. constant water vapor
  if (log_ocean_drsp ==  0) mldclim = d_ocean                   ! sens. exp. no deep ocean

! heat capacity global [J/K/m^2]
  where (z_topo  > 0.) cap_surf = cap_land
  where (z_topo <= 0.) cap_surf = cap_ocean*mldclim(:,:,1)

! decon mean state switch
  if (log_ocean_dmc == 0) cap_surf = cap_land

! initialize fields
  Ts_ini   = Tclim(:,:,nstep_yr)                          ! initial value temp. surf
  Ta_ini   = Ts_ini                                       ! initial value atm. temp.
  To_ini   = Toclim(:,:,nstep_yr)                         ! initial value temp. surf
  q_ini    = qclim(:,:,nstep_yr)                          ! initial value atmos water vapor

  CO2_ctrl = 340.0
! decon mean state switch
  if (log_co2_dmc == 0) CO2_ctrl = 0.

  if (log_exp .ge. 95 .and. log_exp .le. 100 )  CO2_ctrl = 280.  ! IPCC scenarios

  sw_solar = sw_solar_ctrl

  ! define some program constants
  wz_air   = exp(-z_topo/z_air)
  wz_vapor = exp(-z_topo/z_vapor)
  where (uclim(:,:,:) >= 0.0)
     uclim_m = uclim
     uclim_p = 0.0
  elsewhere
     uclim_m = 0.0
     uclim_p = uclim
  end where
  where (vclim(:,:,:) >= 0.0)
     vclim_m = vclim
     vclim_p = 0.0
  elsewhere
     vclim_m = 0.0
     vclim_p = vclim
  end where

  ! initialize the rainfall parameterisation
  select case( log_rain )
  case(-1) ! Original GREB
      c_q=1.; c_rq= 0.; c_omega=0.; c_omegastd=0.
    case(1) ! Adding relative humidity (rq)
      c_q=-1.391649; c_rq=3.018774; c_omega= 0.; c_omegastd=0.
    case(2) ! Adding omega
      c_q=0.862162; c_rq=0.; c_omega=-29.02096; c_omegastd=0.
    case(3) ! Adding rq and omega
      c_q=-0.2685845; c_rq=1.4591853; c_omega=-26.9858807; c_omegastd=0.
    case(0) ! Best GREB
        c_q=-1.88; c_rq=2.25; c_omega=-17.69; c_omegastd=59.07 ! Rainfall parameters (ERA-Interim)
      if (log_clim == 1) then
        c_q=-1.27; c_rq=1.99; c_omega=-16.54; c_omegastd=21.15 ! Rainfall parameters (NCEP)
      end if
  end select

! compute Q-flux corrections
  if ( log_exp .ne. 1 ) then
     print*,'% flux correction ', CO2_ctrl
     1001 format (A4,     T8, A10,   T20, A10,    T32, A15,         T50, A12,      T65, A12,     T80, A15) !TB
     print 1001, "YEAR", "CO2[ppm]", "SW[W/m^2]", "global mean[C]", "Trop Pac[C]", "Hamburg[C]", "North Pole[C]" !TB
     call qflux_correction(CO2_ctrl, Ts_ini, Ta_ini, q_ini, To_ini)

  else if (log_exp .eq. 1 .and. log_qflux_dmc .eq. 1) then
! ONLY FOR WRITING QFLUX FILES
!       print*,'% flux correction ', CO2_ctrl
!       call qflux_correction(CO2_ctrl, Ts_ini, Ta_ini, q_ini, To_ini) ! only for writing files
!       ! write Q-flux corrections
!       open(31,file='Tsurf_flux_correction.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
!       open(32,file='vapour_flux_correction.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
!       open(33,file='Tocean_flux_correction.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
!       do irec=1, nstep_yr
!          write(31,rec=irec)  TF_correct(:,:,irec)
!          write(32,rec=irec)  qF_correct(:,:,irec)
!          write(33,rec=irec)  ToF_correct(:,:,irec)
!       end do
!       stop
       ! read Q-flux corrections
       open(104,file='../input/Tsurf_flux_correction.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
       open(105,file='../input/vapour_flux_correction.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
       open(106,file='../input/Tocean_flux_correction.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
       do irec=1, nstep_yr
          read(104,rec=irec)  TF_correct(:,:,irec)
          read(105,rec=irec)  qF_correct(:,:,irec)
          read(106,rec=irec)  ToF_correct(:,:,irec)
       end do
  else ! => if (log_exp .eq. 1 .and. log_qflux_dmc .eq. 0) then
    TF_correct=0.
    qF_correct=0.
    ToF_correct=0.
  end if

! control run
  print*,'% CONTROL RUN CO2=',CO2_ctrl,'  time=', time_ctrl,'yr'
  print 1001, "YEAR", "CO2[ppm]", "SW[W/m^2]", "global mean[C]", "Trop Pac[C]", "Hamburg[C]", "North Pole[C]" !TB
  Ts1 = Ts_ini; Ta1 = Ta_ini; To1 = To_ini; q1 = q_ini;                   ! initialize fields
  year=1970; mon=1; irec=0; Tmm=0.; Tamm=0.; qmm=0.; apmm=0.;
  do it=1, time_ctrl*nstep_yr                                             ! main time loop
    call time_loop(it, isrec, year, CO2_ctrl, irec, mon, 101, Ts1, Ta1, q1, To1, Ts0,Ta0, q0, To0 )
    Ts1=Ts0; Ta1=Ta0; q1=q0; To1=To0
    if (log_exp .eq. 1 .and. mod(it,nstep_yr) .eq. 0) year=year+1
  end do



! scenario run
if ( log_exp .ne. 1 .or. time_scnr .ne. 0 ) then
  if( log_exp .eq. 30  ) sw_solar = sw_solar_scnr ! paleo 231 kyr bp
  if( log_exp .eq. 31  ) sw_solar = sw_solar_scnr ! paleo 231 kyr bp
  if( log_exp .eq. 35  ) sw_solar = sw_solar_scnr ! change obliquity
  if( log_exp .eq. 36  ) sw_solar = sw_solar_scnr ! change eccentricity
  if( log_exp .eq. 37  ) then ! change solar constant as function of radius
     radius = 1+0.01*(dradius)
     print*,'Solar radius [AU] = ', radius
     rS0 = (1/radius)**2
     sw_solar = rS0*sw_solar
  end if
  if ( log_exp .eq. 230 ) then ! change boundary conditions for Climate Change forcing
     Tclim      = Tclim + Tclim_anom_cc
     uclim      = uclim + uclim_anom_cc
     vclim      = vclim + vclim_anom_cc
     omegaclim  = omegaclim + omegaclim_anom_cc
     wsclim     = wsclim + wsclim_anom_cc
  end if
  if ( log_exp .eq. 240 .or. log_exp .eq. 241 ) then ! change boundary conditions for ENSO forcing
     Tclim      = Tclim + Tclim_anom_enso
     uclim      = uclim + uclim_anom_enso
     vclim      = vclim + vclim_anom_enso
     omegaclim  = omegaclim + omegaclim_anom_enso
     wsclim     = wsclim + wsclim_anom_enso
  end if

  print*,'% SCENARIO EXP: ',log_exp,'  time=', time_scnr,'yr'
  print 1001, "YEAR", "CO2[ppm]", "SW[W/m^2]", "global mean[C]", "Trop Pac[C]", "Hamburg[C]", "North Pole[C]" !TB
  Ts1 = Ts_ini; Ta1 = Ta_ini; q1 = q_ini; To1 = To_ini                     ! initialize fields
  year=1950.; CO2=340.0; mon=1; irec=0; Tmm=0.; Tamm=0.; qmm=0.; apmm=0.;
  if (log_exp .ge. 35 .and. log_exp .le. 37) year=1.

  do it=1, time_scnr*nstep_yr                                              ! main time loop
     call forcing(it, year, CO2, Ts1)
     call time_loop(it,isrec, year, CO2, irec, mon, 102, Ts1, Ta1, q1, To1, Ts0,Ta0, q0, To0 )

     Ts1=Ts0; Ta1=Ta0; q1=q0; To1=To0
     if (mod(it,nstep_yr) == 0) year=year+1
  end do

end if !( log_exp .ne. 1 )

end subroutine

!+++++++++++++++++++++++++++++++++++++++
subroutine time_loop(it, isrec, year, CO2, irec, mon, ionum, Ts1, Ta1, q1, To1, Ts0,Ta0, q0, To0)
!+++++++++++++++++++++++++++++++++++++++
! main time loop

  use mo_numerics
  use mo_physics

  real, dimension(xdim,ydim):: Ts1, Ta1, q1, To1, Ts0,Ta0, q0, To0, sw,       &
&                              ice_cover, Q_sens, Q_lat, Q_lat_air, dq_eva,   &
&                              dq_rain, dTa_crcl, dq_crcl, dq, dT_ocean, dTo, &
&                              LW_surf, LWair_down, LWair_up, em

  jday = mod((it-1)/ndt_days,ndays_yr)+1  ! current calendar day in year
  ityr = mod((it-1),nstep_yr)+1           ! time step in year

  call tendencies(CO2, Ts1, Ta1, To1, q1, ice_cover, SW, LW_surf, Q_lat,      &
&                    Q_sens, Q_lat_air, dq_eva, dq_rain, dq_crcl,             &
&                    dTa_crcl, dT_ocean, dTo, LWair_down, LWair_up, em)

  Tmin_limit = 40 ! no very low Tsurf/Tatmoss;  numerical stability

  ! surface temperature
  Ts0  = Ts1  +dT_ocean +dt*( SW +LW_surf -LWair_down +Q_lat +Q_sens +TF_correct(:,:,ityr)) / cap_surf
  where(Ts0 .le. Tmin_limit )     Ts0 = Tmin_limit ! no very low Tsurf;  numerical stability
  ! air temperature
  Ta0  = Ta1 +dTa_crcl +dt*( LWair_up +LWair_down -em*LW_surf +Q_lat_air -Q_sens )/cap_air
  where(Ta0 .le. Tmin_limit )     Ta0 = Tmin_limit ! no very low Tatmos;  numerical stability
  ! deep ocean temperature
  To0  = To1 +dTo +ToF_correct(:,:,ityr)
  ! air water vapor
  dq = dt*(dq_eva+dq_rain) +dq_crcl + qF_correct(:,:,ityr)
  where(dq .le. -q1 ) dq = -0.9*q1 ! no negative q;  numerical stability
  where(dq .gt.   0.020 ) dq = 0.020   ! no hugh q;  numerical stability

! decon mean state switch
  if( log_hydro_dmc ==  0)    dq   = 0.0

  q0 = q1 + dq
  ! sea ice heat capacity
  call seaice(Ts0)
  ! write output
  call output(it, ionum, irec, mon, ts0, ta0, to0, q0, ice_cover, dq_rain, dq_eva, dq_crcl)
  ! diagnostics: annual means plots
  call diagonstics(it, year, CO2, ts0, ta0, to0, q0, ice_cover, sw, lw_surf, q_lat, q_sens)

end subroutine time_loop

!+++++++++++++++++++++++++++++++++++++++
subroutine tendencies(CO2, Ts1, Ta1, To1, q1, ice_cover, SW, LW_surf, Q_lat, Q_sens, Q_lat_air,  &
&                     dq_eva, dq_rain, dq_crcl, dTa_crcl, dT_ocean, dTo, LWair_down, LWair_up, em)
!+++++++++++++++++++++++++++++++++++++++

  use mo_numerics
  use mo_physics

! declare temporary fields
  real, dimension(xdim,ydim) :: Ts1, Ta1, To1, q1, ice_cover, sw, LWair_up,    &
&                               LWair_down, em, Q_sens, Q_lat, Q_lat_air,      &
&                               dq_eva, dq_rain, dTa_crcl, dq_crcl, LW_surf,   &
&                               dT_ocean, dTo

!$omp parallel sections
!$omp section
    ! SW radiation model
    call SWradiation(Ts1, sw, ice_cover)
!$omp section
    ! LW radiation model
    call LWradiation(Ts1, Ta1, q1, CO2, LW_surf, LWair_up, LWair_down, em)
    ! sensible heat flux
    Q_sens = ct_sens*(Ta1-Ts1)
!$omp section
! decon mean state switch
    if (log_atmos_dmc == 0) Q_sens = 0.

    ! hydro. model
    call hydro(Ts1, q1, Q_lat, Q_lat_air, dq_eva, dq_rain)
    ! atmos. circulation
!$omp section
    call circulation(Ta1, dTa_crcl, z_air, wz_air)       ! air temp
!$omp section
    call circulation( q1,  dq_crcl, z_vapor, wz_vapor)   ! atmos water vapor
!$omp section
    ! deep ocean interaction
    call deep_ocean(Ts1, To1, dT_ocean, dTo)
!$omp end parallel sections

end subroutine tendencies

!+++++++++++++++++++++++++++++++++++++++
subroutine  qflux_correction(CO2_ctrl, Ts1, Ta1, q1, To1)
!+++++++++++++++++++++++++++++++++++++++
!              compute heat flux correction values

  USE mo_numerics
  USE mo_physics

! declare temporary fields
  real, dimension(xdim,ydim) :: Ts0, Ts1, Ta0, Ta1, To0, To1, q0, q1, sw, ice_cover, &
&                               Q_sens, Q_lat, Q_lat_air, dq_eva, dq_rain, LW_surf,  &
&                               LWair_down, LWair_up, em, dTa_crcl, dq_crcl, dTs,    &
&                               dTa, dq, T_error, dT_ocean, dTo

! time loop
  do it=1, time_flux*ndt_days*ndays_yr
     jday = mod((it-1)/ndt_days,ndays_yr)+1  ! current calendar day in year
     ityr = mod((it-1),nstep_yr)+1           ! time step in year
     call tendencies(CO2_ctrl, Ts1, Ta1, To1, q1, ice_cover, SW, LW_surf, Q_lat,     &
&                    Q_sens, Q_lat_air, dq_eva, dq_rain, dq_crcl, dTa_crcl,          &
&                    dT_ocean, dTo, LWair_down, LWair_up, em)

    ! surface temperature without heat flux correction
    dTs = dt*( sw +LW_surf -LWair_down +Q_lat +Q_sens) / cap_surf
    Ts0  = Ts1 +dTs +dT_ocean
    ! air temperature
    dTa = dt*( LWair_up +LWair_down -em*LW_surf +Q_lat_air -Q_sens)/cap_air
    Ta0  = Ta1 + dTa +dTa_crcl
    ! deep ocean temperature without heat flux correction
    To0  = To1 +dTo
   ! air water vapor without flux correction
    dq = dt*(dq_eva+dq_rain)
    q0 = q1 +dq +dq_crcl
   ! heat flux correction Tsurf
    T_error              = Tclim(:,:,ityr) -Ts0 ! error relative to Tclim
    TF_correct(:,:,ityr) = T_error*cap_surf/dt  ! heat flux in [W/m^2]
    ! surface temperature with heat flux correction
    Ts0  = Ts1 +dTs +dT_ocean +TF_correct(:,:,ityr)*dt/ cap_surf
   ! heat flux correction deep ocean
    ToF_correct(:,:,ityr) = Toclim(:,:,ityr) -To0  ! heat flux in [K/dt]
    ! deep ocean temperature with heat flux correction
    To0  = To1 +dTo +ToF_correct(:,:,ityr)
    ! water vapor flux correction
    qF_correct(:,:,ityr) = qclim(:,:,ityr) -q0
    ! air water vapor with flux correction
    q0 = q1 + dq +dq_crcl + qF_correct(:,:,ityr)
    ! sea ice heat capacity
    call seaice(Ts0)
    ! diagnostics: annual means plots
    call diagonstics(it, 0.0, CO2_ctrl, ts0, ta0, to0, q0, ice_cover, sw, lw_surf, q_lat, q_sens)
    ! memory
    Ts1=Ts0; Ta1=Ta0; q1=q0;  To1=To0;

  end do
  1003 format ("On global average a heat flux correction of ", F8.2," W/m^2") !TB
  1004 format ("and a water vapour correction of ", F8.4, " g/kg is applied each time step") !TB
  print 1003, gmean(sum(TF_correct,3)/nstep_yr) !TB
  print 1004, gmean(sum(qF_correct,3)/nstep_yr)*100 !TB
end subroutine qflux_correction

!+++++++++++++++++++++++++++++++++++++++
subroutine SWradiation(Tsurf, sw, ice_cover)
!+++++++++++++++++++++++++++++++++++++++
!    SW radiation model

  USE mo_numerics,    ONLY: xdim, ydim
  USE mo_physics,     ONLY: ityr, sw_solar,da_ice, a_no_ice, a_cloud, z_topo  &
&                         , Tl_ice1, Tl_ice2, To_ice1, To_ice2, glacier       &
&                         , cldclim, log_exp, log_atmos_dmc, log_ice, S0_var

! declare temporary fields
  real, dimension(xdim,ydim)  :: Tsurf, sw, albedo, ice_cover, a_surf, a_atmos

! atmos albedo
  a_atmos=cldclim(:,:,ityr)*a_cloud

! isnow/ice diagnostic only
  ice_cover=0.0
  where(z_topo >= 0. .and. Tsurf <= Tl_ice1) ice_cover = 1.0
  where(z_topo >= 0. .and. Tsurf >  Tl_ice1 .and. Tsurf < Tl_ice2 ) &
&       ice_cover = (1-(Tsurf-Tl_ice1)/(Tl_ice2-Tl_ice1))
  where(z_topo < 0.  .and. Tsurf <= To_ice1) ice_cover = 1.0
  where(z_topo < 0.  .and. Tsurf >  To_ice1 .and. Tsurf < To_ice2 ) &
&       ice_cover = (1-(Tsurf-To_ice1)/(To_ice2-To_ice1))

! surface albedo
! Land:  ice -> albedo linear function of T_surf
   where(z_topo >= 0. .and. Tsurf <= Tl_ice1) a_surf = a_no_ice+da_ice   ! ice
   where(z_topo >= 0. .and. Tsurf >= Tl_ice2) a_surf = a_no_ice          ! no ice
   where(z_topo >= 0. .and. Tsurf > Tl_ice1 .and. Tsurf < Tl_ice2 ) &
&       a_surf = a_no_ice +da_ice*(1-(Tsurf-Tl_ice1)/(Tl_ice2-Tl_ice1))
! Ocean: ice -> albedo/heat capacity linear function of T_surf
  where(z_topo < 0. .and. Tsurf <= To_ice1) a_surf = a_no_ice+da_ice      ! ice
  where(z_topo < 0. .and. Tsurf >= To_ice2) a_surf = a_no_ice             ! no ice
  where(z_topo < 0. .and. Tsurf > To_ice1 .and. Tsurf < To_ice2 ) &
&       a_surf = a_no_ice+da_ice*(1-(Tsurf-To_ice1)/(To_ice2-To_ice1))

! glacier -> no albedo changes
  where(glacier > 0.5) a_surf = a_no_ice+da_ice

! dmc & decon2xco2 switch
  if (log_ice  ==    0) a_surf = a_no_ice

! SW flux
  albedo=a_surf+a_atmos-a_surf*a_atmos
  forall (i=1:xdim)
     sw(i,:)=0.01*S0_var*SW_solar(:,ityr)*(1-albedo(i,:))
  end forall

end subroutine SWradiation


!+++++++++++++++++++++++++++++++++++++++
subroutine LWradiation(Tsurf, Tair, q, CO2, LWsurf, LWair_up, LWair_down, em)
!+++++++++++++++++++++++++++++++++++++++
! new approach with LW atmos

  USE mo_numerics,    ONLY: xdim, ydim
  USE mo_physics,     ONLY: sig, eps, qclim, cldclim, z_topo, jday, ityr,         &
&                           r_qviwv, z_air, z_vapor, dTrad, p_emi, log_exp,       &
&                           log_atmos_dmc, co2_part

! declare temporary fields
  real, dimension(xdim,ydim)  :: Tsurf, Tair, q, LWsurf, LWair, e_co2, e_cloud,   &
&                                LWair_up, LWair_down, e_vapor, em


  e_co2   = exp(-z_topo/z_air)*co2_part*CO2   ! CO2
  e_vapor = exp(-z_topo/z_air)*r_qviwv*q      ! water vapor
  e_cloud = cldclim(:,:,ityr)                 ! clouds

! total
  em      = p_emi(4)*log( p_emi(1)*e_co2 +p_emi(2)*e_vapor +p_emi(3) ) +p_emi(7)   &
&          +p_emi(5)*log( p_emi(1)*e_co2   +p_emi(3) )                             &
&          +p_emi(6)*log( p_emi(2)*e_vapor +p_emi(3) )
  em      = (p_emi(8)-e_cloud)/p_emi(9)*(em-p_emi(10))+p_emi(10)

  LWsurf      = -sig*Tsurf**4
  LWair_down  = -em*sig*(Tair+dTrad(:,:,ityr))**4
  LWair_up    = LWair_down

! decon mean state switch
  if( log_atmos_dmc == 0) LWair_down = 0

end subroutine LWradiation


!+++++++++++++++++++++++++++++++++++++++
subroutine hydro(Tsurf, q, Qlat, Qlat_air, dq_eva, dq_rain)
!+++++++++++++++++++++++++++++++++++++++
!    hydrological model for latent heat and water vapor

  USE mo_numerics,    ONLY: xdim, ydim
  USE mo_physics,     ONLY: rho_air, uclim, vclim, z_topo, swetclim, ityr,   &
&                           ce, cq_latent, cq_rain, z_air, r_qviwv, log_exp, &
&                           log_atmos_dmc, log_hydro_dmc, log_hydro_drsp,    &
&                           omegaclim, omegastdclim, wsclim,  wz_vapor,      &
&                           c_q, c_rq, c_omega, c_omegastd, log_eva, log_rain      ! Rainfall parameters

! declare temporary fields
  real, dimension(xdim,ydim)  :: Tsurf, Tskin, q, Qlat, Qlat_air, qs, dq_eva, &
&                                dq_rain, abswind, rq

  Qlat=0.; Qlat_air=0.; dq_eva=0.; dq_rain=0.

! decon mean state switch
  if( log_atmos_dmc == 0) return
  if( log_hydro_dmc == 0) return
! decon2xco2 switch
  if(log_hydro_drsp == 0) return

  abswind = sqrt(uclim(:,:,ityr)**2 +vclim(:,:,ityr)**2)
  where(z_topo > 0. ) abswind = sqrt(abswind**2 +2.0**2) ! land
  where(z_topo < 0. ) abswind = sqrt(abswind**2 +3.0**2) ! ocean

! saturated humiditiy (max. air water vapor)
  qs = 3.75e-3*exp(17.08085*(Tsurf-273.15)/(Tsurf-273.15+234.175));
  qs = qs*exp(-z_topo/z_air) ! scale qs by topography

! relative humidity
  rq = q/qs

! latent heat flux surface
  if ( log_eva == -1 ) then
    abswind = sqrt(uclim(:,:,ityr)**2 +vclim(:,:,ityr)**2)
    where(z_topo > 0. ) abswind = sqrt(abswind**2 + 2.0**2) !< land turbulent wind
    where(z_topo < 0. ) abswind = sqrt(abswind**2 + 3.0**2) !< ocean turbulent wind
    Qlat   = (q-qs)*abswind*cq_latent*rho_air*ce*swetclim(:,:,ityr) ! latend heat flux
  else if ( log_eva == 1 ) then
    where(z_topo > 0. ) abswind = sqrt(abswind**2 + 144.**2) ! land turbulent wind
    where(z_topo < 0. ) abswind = sqrt(abswind**2 + 7.1**2) ! ocean turbulent wind
    where(z_topo > 0. )  Qlat   = (q-qs)*abswind*cq_latent*rho_air*0.04*ce*swetclim(:,:,ityr) ! latend heat flux land
    where(z_topo <= 0. ) Qlat   = (q-qs)*abswind*cq_latent*rho_air*0.73*ce*swetclim(:,:,ityr) ! latend heat flux ocean
  else if ( log_eva == 2 ) then
    abswind = wsclim(:,:,ityr) ! use the wind speed climatology
    where(z_topo > 0. )  abswind = sqrt(abswind**2 + 9.0**2) ! land turbulent wind
    where(z_topo <= 0. ) abswind = sqrt(abswind**2 + 4.0**2) ! ocean turbulent wind
    where(z_topo > 0. )  Qlat   = (q-qs)*abswind*cq_latent*rho_air*0.56*ce*swetclim(:,:,ityr) ! latend heat flux land
    where(z_topo <= 0. ) Qlat   = (q-qs)*abswind*cq_latent*rho_air*0.79*ce*swetclim(:,:,ityr) ! latend heat flux ocean
  else if ( log_eva == 0 ) then
    where(z_topo > 0. )  Tskin = Tsurf + 5. ! skin temperature land
    where(z_topo <= 0. ) Tskin = Tsurf + 1. ! skin temperature ocean
    qs = 3.75e-3*exp(17.08085*(Tskin-273.15)/(Tskin-273.15+234.175)) ! re-calculate saturation pressure
    qs = qs*exp(-z_topo/z_air) ! scale qs by topography
    where(z_topo > 0. ) abswind = sqrt(wsclim(:,:,ityr)**2 + 11.5**2) ! land turbulent wind
    where(z_topo <= 0. ) abswind = sqrt(wsclim(:,:,ityr)**2 + 5.4**2) ! ocean turbulent wind
    where(z_topo > 0. )  Qlat   = (q-qs)*abswind*cq_latent*rho_air*0.25*ce*swetclim(:,:,ityr) ! latend heat flux land
    where(z_topo <= 0. ) Qlat   = (q-qs)*abswind*cq_latent*rho_air*0.58*ce*swetclim(:,:,ityr) ! latend heat flux ocean
  end if
! change in water vapor
  dq_eva  = -Qlat/cq_latent/r_qviwv  ! evaporation

! precipitation -> Eq. 11 in Stassen et al 2019
! Parameters in unused terms are set to zero
  dq_rain = (c_q + c_rq*rq + c_omega*omegaclim(:,:,ityr) + c_omegastd*omegastdclim(:,:,ityr))*cq_rain*q
  if (log_rain == 1) then
    where(dq_rain >= -0.0015 / (wz_vapor * r_qviwv * 86400.)) dq_rain = -0.0015 / (wz_vapor * r_qviwv * 86400.) !Avoid negative rainfall (dq_rain is negative means positive rainfall!)
  end if

! latent heat flux atmos
  Qlat_air = -dq_rain*cq_latent*r_qviwv

end subroutine hydro

!+++++++++++++++++++++++++++++++++++++++
subroutine seaice(Tsurf)
!+++++++++++++++++++++++++++++++++++++++
!              SW radiation model

  USE mo_numerics,    ONLY: xdim, ydim
  USE mo_physics,     ONLY: ityr, z_topo, cap_surf, cap_land, cap_ocean, &
&                           log_exp, To_ice1, To_ice2, glacier, mldclim, &
&                           log_ice, log_ocean_dmc

! declare temporary fields
  real, dimension(xdim,ydim)  :: Tsurf

! decon mean state switch
  if( log_ocean_dmc == 0) return

  where(z_topo < 0. .and. Tsurf <= To_ice1) cap_surf = cap_land                    ! sea ice
  where(z_topo < 0. .and. Tsurf >= To_ice2) cap_surf = cap_ocean*mldclim(:,:,ityr) ! open ocean
  where(z_topo < 0. .and. Tsurf > To_ice1 .and. Tsurf < To_ice2 ) &
&       cap_surf = cap_land + (cap_ocean*mldclim(:,:,ityr)-cap_land)     &
&                            /(To_ice2-To_ice1)*(Tsurf-To_ice1)

! dmc & decon2xco2 switch
  if( log_ice == 0 ) then
     where(z_topo > 0. ) cap_surf = cap_land                     ! sea ice
     where(z_topo < 0. ) cap_surf = cap_ocean*mldclim(:,:,ityr)  ! open ocean
  end if

! glacier -> no sea ice change
  where(glacier > 0.5) cap_surf = cap_land                       ! ice sheet

end subroutine seaice

!+++++++++++++++++++++++++++++++++++++++
subroutine deep_ocean(Ts, To, dT_ocean, dTo)
!+++++++++++++++++++++++++++++++++++++++
!              deep ocean model

  USE mo_numerics,    ONLY: xdim, ydim, nstep_yr, dt
  USE mo_physics,     ONLY: ityr, z_topo, mldclim, log_exp, To_ice2,     &
&                           cap_ocean, co_turb, z_ocean, log_ocean_dmc, log_ocean_drsp

! declare temporary fields
  real, dimension(xdim,ydim)  :: Ts, To, dT_ocean, dTo, dmld, Tx
  dT_ocean = 0.0;  dTo     = 0.0

! decon mean state switch
  if ( log_ocean_dmc == 0 )  return
! decon2xco2 switch
  if ( log_ocean_drsp == 0 ) return

  if (ityr >  1) dmld = mldclim(:,:,ityr)-mldclim(:,:,ityr-1)
  if (ityr == 1) dmld = mldclim(:,:,ityr)-mldclim(:,:,nstep_yr)

! entrainment & detrainment
  where ( z_topo < 0 .and. Ts >= To_ice2 .and. dmld < 0)     &
&       dTo      = -dmld/(z_ocean-mldclim(:,:,ityr))*(Ts-To)
  where ( z_topo < 0 .and. Ts >= To_ice2 .and. dmld > 0)     &
&       dT_ocean =  dmld/mldclim(:,:,ityr)*(To-Ts)

 c_effmix = 0.5
 dTo      = c_effmix*dTo
 dT_ocean = c_effmix*dT_ocean

! turbulent mixing
  Tx = max(To_ice2,Ts)
 where ( z_topo < 0 ) dTo      = dTo      + dt*co_turb*(Tx-To)/(cap_ocean*(z_ocean-mldclim(:,:,ityr)))
 where ( z_topo < 0 ) dT_ocean = dT_ocean + dt*co_turb*(To-Tx)/(cap_ocean*mldclim(:,:,ityr))

end subroutine deep_ocean

!+++++++++++++++++++++++++++++++++++++++
subroutine circulation(X_in, dX_crcl, h_scl, wz)
!+++++++++++++++++++++++++++++++++++++++
! circulation with shorter time step

  USE mo_numerics,  ONLY: xdim, ydim, dt, dt_crcl
  USE mo_physics,   ONLY: z_vapor, z_air, log_exp, log_atmos_dmc, &
&                         log_vdif, log_vadv, log_hdif, log_hadv, log_conv
  implicit none

  real, dimension(xdim,ydim), intent(in)  :: X_in, wz
  real,                       intent(in)  :: h_scl
  real, dimension(xdim,ydim), intent(out) :: dX_crcl

  real, dimension(xdim,ydim) :: X, dx_diffuse, dx_advec, dx_conv
  integer time, tt

  dX_crcl    = 0.0

! decon mean state switch
  if (log_atmos_dmc == 0 ) return

  dx_diffuse = 0.0
  dx_advec   = 0.0
  dx_conv    = 0.0
  time=max(1,nint(float(dt)/dt_crcl))

  X = X_in;
  do tt=1, time   ! time loop circulation
! dmc & decon2xco2 switch
     if (log_vdif == 1 .and. h_scl .eq. z_vapor) call diffusion(X, dx_diffuse, h_scl, wz)
     if (log_vadv == 1 .and. h_scl .eq. z_vapor) call advection(X, dx_advec, h_scl, wz)
     if (log_conv == 0 .and. h_scl .eq. z_vapor) call convergence(X, dx_conv)
     if (log_hdif == 1 .and. h_scl .eq. z_air)   call diffusion(X, dx_diffuse, h_scl, wz)
     if (log_hadv == 1 .and. h_scl .eq. z_air)   call advection(X, dx_advec, h_scl, wz)
     X = X + dx_diffuse + dx_advec + dx_conv
  end do           ! time loop
  dX_crcl = X - X_in

end subroutine circulation

!+++++++++++++++++++++++++++++++++++++++
subroutine diffusion(T1, dX_diffuse,h_scl, wz)
!+++++++++++++++++++++++++++++++++++++++
!    diffusion

  USE mo_numerics,   ONLY: xdim, ydim, dt, dlon, dlat, dt_crcl
  USE mo_physics,    ONLY: pi, z_topo, log_exp, kappa, z_vapor
  implicit none

  real, dimension(xdim,ydim), intent(in)  :: T1, wz
  real                      , intent(in)  :: h_scl
  real, dimension(xdim,ydim), intent(out) :: dX_diffuse

  integer :: i
  integer, dimension(ydim)   :: ilat = (/(i,i=1,ydim)/)
  real, dimension(ydim)      :: lat, dxlat, ccx
  real, dimension(xdim)      :: T1h, dTxh
  real, dimension(xdim,ydim) :: dTx, dTy

  real    :: deg, dd, dx, dy, dyy, ccy, ccx2
  integer :: j, k, km1, kp1, jm1, jm2, jm3, jp1, jp2, jp3
  integer :: time2, dtdff2, tt2

  deg = 2.*pi*6.371e6/360.;   ! length of 1deg latitude [m]
  dx = dlon; dy=dlat; dyy=dy*deg
  lat = dlat*ilat-dlat/2.-90.;  dxlat=dx*deg*cos(2.*pi/360.*lat)
  ccy=kappa*dt_crcl/dyy**2
  ccx=kappa*dt_crcl/dxlat**2

     ! latitudinal
     do k=1, ydim
        km1=k-1;  kp1=k+1
        if ( k>=2 .and. k<=ydim-1)   dTy(:,k)=ccy*(                                      &
&                         wz(:,km1)*(T1(:,km1)-T1(:,k)) +wz(:,kp1)*(T1(:,kp1)-T1(:,k)) )
        if ( k==1 )                  dTy(:,k)=ccy*wz(:,kp1)*(-T1(:,k)+T1(:,kp1))
        if ( k==ydim )               dTy(:,k)=ccy*wz(:,km1)*(T1(:,km1)-T1(:,k))
        ! longitudinal
        if ( dxlat(k) > 2.5e5) then  ! unitl 25degree
           j = 1
           jp1 = j+1; jp2 = j+2; jp3 = j+3; jm1 = xdim; jm2 = xdim-1; jm3 = xdim-2
           dTx(j,k)=ccx(k)*(                                                           &
&            10*( wz(jm1,k)*(T1(jm1,k)-T1(j,k))   +wz(jp1,k)*(T1(jp1,k) -T1(j,k))   )  &
&            +4*( wz(jm2,k)*(T1(jm2,k)-T1(jm1,k)) +wz(jm1,k)*(T1(j,k)   -T1(jm1,k)) )  &
&            +4*( wz(jp1,k)*(T1(j,k)  -T1(jp1,k)) +wz(jp2,k)*(T1(jp2,k) -T1(jp1,k)) )  &
&            +1*( wz(jm3,k)*(T1(jm3,k)-T1(jm2,k)) +wz(jm2,k)*(T1(jm1,k) -T1(jm2,k)) )  &
&            +1*( wz(jp2,k)*(T1(jp1,k)-T1(jp2,k)) +wz(jp3,k)*(T1(jp3,k) -T1(jp2,k)) ) )/20.
           j = 2
           jp1 = j+1; jp2 = j+2; jp3 = j+3; jm1 = j-1; jm2 = xdim; jm3 = xdim-1
           dTx(j,k)=ccx(k)*(                                                           &
&            10*( wz(jm1,k)*(T1(jm1,k)-T1(j,k))   +wz(jp1,k)*(T1(jp1,k) -T1(j,k))   )  &
&            +4*( wz(jm2,k)*(T1(jm2,k)-T1(jm1,k)) +wz(jm1,k)*(T1(j,k)   -T1(jm1,k)) )  &
&            +4*( wz(jp1,k)*(T1(j,k)  -T1(jp1,k)) +wz(jp2,k)*(T1(jp2,k) -T1(jp1,k)) )  &
&            +1*( wz(jm3,k)*(T1(jm3,k)-T1(jm2,k)) +wz(jm2,k)*(T1(jm1,k) -T1(jm2,k)) )  &
&            +1*( wz(jp2,k)*(T1(jp1,k)-T1(jp2,k)) +wz(jp3,k)*(T1(jp3,k) -T1(jp2,k)) ) )/20.
           j = 3
           jp1 = j+1; jp2 = j+2; jp3 = j+3; jm1 = j-1; jm2 = j-2; jm3 = xdim
           dTx(j,k)=ccx(k)*(                                                           &
&            10*( wz(jm1,k)*(T1(jm1,k)-T1(j,k))   +wz(jp1,k)*(T1(jp1,k) -T1(j,k))   )  &
&            +4*( wz(jm2,k)*(T1(jm2,k)-T1(jm1,k)) +wz(jm1,k)*(T1(j,k)   -T1(jm1,k)) )  &
&            +4*( wz(jp1,k)*(T1(j,k)  -T1(jp1,k)) +wz(jp2,k)*(T1(jp2,k) -T1(jp1,k)) )  &
&            +1*( wz(jm3,k)*(T1(jm3,k)-T1(jm2,k)) +wz(jm2,k)*(T1(jm1,k) -T1(jm2,k)) )  &
&            +1*( wz(jp2,k)*(T1(jp1,k)-T1(jp2,k)) +wz(jp3,k)*(T1(jp3,k) -T1(jp2,k)) ) )/20.
           do j=4, xdim-3              ! longitudinal
              jm1=j-1; jp1=j+1; jm2=j-2; jp2=j+2; jm3=j-3; jp3=j+3
              ! 3.order solution: stable unitl 84degree (dx=2.5degree, a=5e5)
              dTx(j,k)=ccx(k)*(                                                           &
&               10*( wz(jm1,k)*(T1(jm1,k)-T1(j,k))   +wz(jp1,k)*(T1(jp1,k) -T1(j,k))   )  &
&               +4*( wz(jm2,k)*(T1(jm2,k)-T1(jm1,k)) +wz(jm1,k)*(T1(j,k)   -T1(jm1,k)) )  &
&               +4*( wz(jp1,k)*(T1(j,k)  -T1(jp1,k)) +wz(jp2,k)*(T1(jp2,k) -T1(jp1,k)) )  &
&               +1*( wz(jm3,k)*(T1(jm3,k)-T1(jm2,k)) +wz(jm2,k)*(T1(jm1,k) -T1(jm2,k)) )  &
&               +1*( wz(jp2,k)*(T1(jp1,k)-T1(jp2,k)) +wz(jp3,k)*(T1(jp3,k) -T1(jp2,k)) ) )/20.
           end do
           j = xdim-2
           jm1 = j-1; jm2 = j-2; jm3 = j-3; jp1 = j+1; jp2 = j+2; jp3 = 1;
           dTx(j,k)=ccx(k)*(                                                           &
&            10*( wz(jm1,k)*(T1(jm1,k)-T1(j,k))   +wz(jp1,k)*(T1(jp1,k) -T1(j,k))   )  &
&            +4*( wz(jm2,k)*(T1(jm2,k)-T1(jm1,k)) +wz(jm1,k)*(T1(j,k)   -T1(jm1,k)) )  &
&            +4*( wz(jp1,k)*(T1(j,k)  -T1(jp1,k)) +wz(jp2,k)*(T1(jp2,k) -T1(jp1,k)) )  &
&            +1*( wz(jm3,k)*(T1(jm3,k)-T1(jm2,k)) +wz(jm2,k)*(T1(jm1,k) -T1(jm2,k)) )  &
&            +1*( wz(jp2,k)*(T1(jp1,k)-T1(jp2,k)) +wz(jp3,k)*(T1(jp3,k) -T1(jp2,k)) ) )/20.
           j = xdim-1
           jm1 = j-1; jm2 = j-2; jm3 = j-3; jp1 = j+1; jp2 = 1; jp3 = 2
           dTx(j,k)=ccx(k)*(                                                           &
&            10*( wz(jm1,k)*(T1(jm1,k)-T1(j,k))   +wz(jp1,k)*(T1(jp1,k) -T1(j,k))   )  &
&            +4*( wz(jm2,k)*(T1(jm2,k)-T1(jm1,k)) +wz(jm1,k)*(T1(j,k)   -T1(jm1,k)) )  &
&            +4*( wz(jp1,k)*(T1(j,k)  -T1(jp1,k)) +wz(jp2,k)*(T1(jp2,k) -T1(jp1,k)) )  &
&            +1*( wz(jm3,k)*(T1(jm3,k)-T1(jm2,k)) +wz(jm2,k)*(T1(jm1,k) -T1(jm2,k)) )  &
&            +1*( wz(jp2,k)*(T1(jp1,k)-T1(jp2,k)) +wz(jp3,k)*(T1(jp3,k) -T1(jp2,k)) ) )/20.
           j = xdim
           jm1 = j-1; jm2 = j-2; jm3 = j-3; jp1 = 1; jp2 = 2; jp3 = 3
           dTx(j,k)=ccx(k)*(                                                           &
&            10*( wz(jm1,k)*(T1(jm1,k)-T1(j,k))   +wz(jp1,k)*(T1(jp1,k) -T1(j,k))   )  &
&            +4*( wz(jm2,k)*(T1(jm2,k)-T1(jm1,k)) +wz(jm1,k)*(T1(j,k)   -T1(jm1,k)) )  &
&            +4*( wz(jp1,k)*(T1(j,k)  -T1(jp1,k)) +wz(jp2,k)*(T1(jp2,k) -T1(jp1,k)) )  &
&            +1*( wz(jm3,k)*(T1(jm3,k)-T1(jm2,k)) +wz(jm2,k)*(T1(jm1,k) -T1(jm2,k)) )  &
&            +1*( wz(jp2,k)*(T1(jp1,k)-T1(jp2,k)) +wz(jp3,k)*(T1(jp3,k) -T1(jp2,k)) ) )/20.
        else  ! high resolution -> smaller time steps
            dd=max(1,nint(dt_crcl/(1.*dxlat(k)**2/kappa))); dtdff2=dt_crcl/dd
            time2=max(1,nint(float(dt_crcl)/float(dtdff2)))
            ccx2=kappa*dtdff2/dxlat(k)**2
            T1h=T1(:,k)
            do tt2=1, time2      ! additional time loop
              j = 1
              jp1 = j+1; jp2 = j+2; jp3 = j+3; jm1 = xdim; jm2 = xdim-1; jm3 = xdim-2
              dTxh(j) = ccx2*(                                                         &
&                10*( wz(jm1,k)*(T1h(jm1)-T1h(j))   +wz(jp1,k)*(T1h(jp1) -T1h(j))   )  &
&                +4*( wz(jm2,k)*(T1h(jm2)-T1h(jm1)) +wz(jm1,k)*(T1h(j)   -T1h(jm1)) )  &
&                +4*( wz(jp1,k)*(T1h(j)  -T1h(jp1)) +wz(jp2,k)*(T1h(jp2) -T1h(jp1)) )  &
&                +1*( wz(jm3,k)*(T1h(jm3)-T1h(jm2)) +wz(jm2,k)*(T1h(jm1) -T1h(jm2)) )  &
&                +1*( wz(jp2,k)*(T1h(jp1)-T1h(jp2)) +wz(jp3,k)*(T1h(jp3) -T1h(jp2)) ) )/20.
              j = 2
              jp1 = j+1; jp2 = j+2; jp3 = j+3; jm1 = j-1; jm2 = xdim; jm3 = xdim-1
              dTxh(j) = ccx2*(                                                         &
&                10*( wz(jm1,k)*(T1h(jm1)-T1h(j))   +wz(jp1,k)*(T1h(jp1) -T1h(j))   )  &
&                +4*( wz(jm2,k)*(T1h(jm2)-T1h(jm1)) +wz(jm1,k)*(T1h(j)   -T1h(jm1)) )  &
&                +4*( wz(jp1,k)*(T1h(j)  -T1h(jp1)) +wz(jp2,k)*(T1h(jp2) -T1h(jp1)) )  &
&                +1*( wz(jm3,k)*(T1h(jm3)-T1h(jm2)) +wz(jm2,k)*(T1h(jm1) -T1h(jm2)) )  &
&                +1*( wz(jp2,k)*(T1h(jp1)-T1h(jp2)) +wz(jp3,k)*(T1h(jp3) -T1h(jp2)) ) )/20.
              j = 3
              jp1 = j+1; jp2 = j+2; jp3 = j+3; jm1 = j-1; jm2 = j-2; jm3 = xdim;
              dTxh(j) = ccx2*(                                                         &
&                10*( wz(jm1,k)*(T1h(jm1)-T1h(j))   +wz(jp1,k)*(T1h(jp1) -T1h(j))   )  &
&                +4*( wz(jm2,k)*(T1h(jm2)-T1h(jm1)) +wz(jm1,k)*(T1h(j)   -T1h(jm1)) )  &
&                +4*( wz(jp1,k)*(T1h(j)  -T1h(jp1)) +wz(jp2,k)*(T1h(jp2) -T1h(jp1)) )  &
&                +1*( wz(jm3,k)*(T1h(jm3)-T1h(jm2)) +wz(jm2,k)*(T1h(jm1) -T1h(jm2)) )  &
&                +1*( wz(jp2,k)*(T1h(jp1)-T1h(jp2)) +wz(jp3,k)*(T1h(jp3) -T1h(jp2)) ) )/20.
                do j=4, xdim-3     ! longitudinal
                    jm1=j-1; jp1=j+1; jm2=j-2; jp2=j+2; jm3=j-3; jp3=j+3
                    dTxh(j)=ccx2*(                                                           &
&                      10*( wz(jm1,k)*(T1h(jm1)-T1h(j))   +wz(jp1,k)*(T1h(jp1) -T1h(j))   )  &
&                      +4*( wz(jm2,k)*(T1h(jm2)-T1h(jm1)) +wz(jm1,k)*(T1h(j)   -T1h(jm1)) )  &
&                      +4*( wz(jp1,k)*(T1h(j)  -T1h(jp1)) +wz(jp2,k)*(T1h(jp2) -T1h(jp1)) )  &
&                      +1*( wz(jm3,k)*(T1h(jm3)-T1h(jm2)) +wz(jm2,k)*(T1h(jm1) -T1h(jm2)) )  &
&                      +1*( wz(jp2,k)*(T1h(jp1)-T1h(jp2)) +wz(jp3,k)*(T1h(jp3) -T1h(jp2)) ) )/20.

                end do           ! longitudinal
              j = xdim-2
              jm1 = j-1; jm2 = j-2; jm3 = j-3; jp1 = j+1; jp2 = j+2; jp3 = 1
              dTxh(j) = ccx2*(                                                         &
&                10*( wz(jm1,k)*(T1h(jm1)-T1h(j))   +wz(jp1,k)*(T1h(jp1) -T1h(j))   )  &
&                +4*( wz(jm2,k)*(T1h(jm2)-T1h(jm1)) +wz(jm1,k)*(T1h(j)   -T1h(jm1)) )  &
&                +4*( wz(jp1,k)*(T1h(j)  -T1h(jp1)) +wz(jp2,k)*(T1h(jp2) -T1h(jp1)) )  &
&                +1*( wz(jm3,k)*(T1h(jm3)-T1h(jm2)) +wz(jm2,k)*(T1h(jm1) -T1h(jm2)) )  &
&                +1*( wz(jp2,k)*(T1h(jp1)-T1h(jp2)) +wz(jp3,k)*(T1h(jp3) -T1h(jp2)) ) )/20.
              j = xdim-1
              jm1 = j-1; jm2 = j-2; jm3 = j-3; jp1 = j+1; jp2 = 1; jp3 = 2
              dTxh(j) = ccx2*(                                                         &
&                10*( wz(jm1,k)*(T1h(jm1)-T1h(j))   +wz(jp1,k)*(T1h(jp1) -T1h(j))   )  &
&                +4*( wz(jm2,k)*(T1h(jm2)-T1h(jm1)) +wz(jm1,k)*(T1h(j)   -T1h(jm1)) )  &
&                +4*( wz(jp1,k)*(T1h(j)  -T1h(jp1)) +wz(jp2,k)*(T1h(jp2) -T1h(jp1)) )  &
&                +1*( wz(jm3,k)*(T1h(jm3)-T1h(jm2)) +wz(jm2,k)*(T1h(jm1) -T1h(jm2)) )  &
&                +1*( wz(jp2,k)*(T1h(jp1)-T1h(jp2)) +wz(jp3,k)*(T1h(jp3) -T1h(jp2)) ) )/20.
              j = xdim
              jm1 = j-1; jm2 = j-2; jm3 = j-3; jp1 = 1; jp2 = 2; jp3 = 3
              dTxh(j) = ccx2*(                                                         &
&                10*( wz(jm1,k)*(T1h(jm1)-T1h(j))   +wz(jp1,k)*(T1h(jp1) -T1h(j))   )  &
&                +4*( wz(jm2,k)*(T1h(jm2)-T1h(jm1)) +wz(jm1,k)*(T1h(j)   -T1h(jm1)) )  &
&                +4*( wz(jp1,k)*(T1h(j)  -T1h(jp1)) +wz(jp2,k)*(T1h(jp2) -T1h(jp1)) )  &
&                +1*( wz(jm3,k)*(T1h(jm3)-T1h(jm2)) +wz(jm2,k)*(T1h(jm1) -T1h(jm2)) )  &
&                +1*( wz(jp2,k)*(T1h(jp1)-T1h(jp2)) +wz(jp3,k)*(T1h(jp3) -T1h(jp2)) ) )/20.
                where(dTxh .le. -T1h ) dTxh = -0.9*T1h ! no negative q;  numerical stability
                T1h=T1h+dTxh
            end do               ! additional time loop
            dTx(:,k)=T1h-T1(:,k)
        end if
    end do          ! y-loop
    dX_diffuse = wz * (dTx + dTy);

end subroutine diffusion

!+++++++++++++++++++++++++++++++++++++++
subroutine advection(T1, dX_advec,h_scl, wz)
!+++++++++++++++++++++++++++++++++++++++
!    advection after DD

  USE mo_numerics, ONLY: xdim, ydim, dt, dlon, dlat, dt_crcl
  USE mo_physics,  ONLY: pi, z_topo, uclim, vclim, ityr, z_vapor, log_exp
  USE mo_physics,  ONLY: uclim_m, uclim_p, vclim_m, vclim_p
  implicit none

  real, dimension(xdim,ydim), intent(in)  :: T1, wz
  real                      , intent(in)  :: h_scl
  real, dimension(xdim,ydim), intent(out) :: dX_advec

  integer :: i
  integer, dimension(ydim):: ilat = (/(i,i=1,ydim)/)
  real, dimension(ydim) :: lat, dxlat, ccx
  real, dimension(xdim) :: T1h, dTxh
  real, dimension(xdim,ydim) :: ddx, T, dTx, dTy
  integer time2, dtdff2, tt2

  real    :: deg, dx, dy, dd, dyy, ccy, ccx2
  integer :: j, k, km1, km2, kp1, kp2, jm1, jm2, jm3, jp1, jp2, jp3

  deg = 2.*pi*6.371e6/360.;   ! length of 1deg latitude [m]
  dx = dlon; dy=dlat; dyy=dy*deg
  lat = dlat*ilat-dlat/2.-90.;  dxlat=dx*deg*cos(2.*pi/360.*lat)
  ccy=dt_crcl/dyy/2.
  ccx=dt_crcl/dxlat/2.

     ! latitudinal
     k=1
     kp1=k+1; kp2=k+2
     do j = 1, xdim
        dTy(j,k) = ccy * (                                                        &
&                     vclim_p(j,k,ityr)*( wz(j,kp1)*(T1(j,k)-T1(j,kp1))           &
&                                        +wz(j,kp2)*(T1(j,k)-T1(j,kp2)) ) )/3.
     end do
     k=2
     km1=k-1; kp1=k+1; kp2=k+2
     do j = 1, xdim
        dTy(j,k) = ccy * (                                                        &
&                    -vclim_m(j,k,ityr)*( wz(j,km1)*(T1(j,k)-T1(j,km1)))          &
&                   + vclim_p(j,k,ityr)*( wz(j,kp1)*(T1(j,k)-T1(j,kp1))           &
&                                        +wz(j,kp2)*(T1(j,k)-T1(j,kp2)) )/3. )
     end do
     do k=3, ydim-2
        km1=k-1; kp1=k+1; km2=k-2; kp2=k+2
        do j = 1, xdim
           dTy(j,k) = ccy * (                                                     &
&                       -vclim_m(j,k,ityr)*( wz(j,km1)*(T1(j,k)-T1(j,km1))        &
&                                           +wz(j,km2)*(T1(j,k)-T1(j,km2)) )      &
&                      + vclim_p(j,k,ityr)*( wz(j,kp1)*(T1(j,k)-T1(j,kp1))        &
&                                           +wz(j,kp2)*(T1(j,k)-T1(j,kp2)) ) )/3.
        end do
     end do
     k=ydim-1
     km1=k-1; kp1=k+1; km2=k-2
     do j = 1, xdim
        dTy(j,k) = ccy * (                                                        &
&                    -vclim_m(j,k,ityr)*( wz(j,km1)*(T1(j,k)-T1(j,km1))           &
&                                        +wz(j,km2)*(T1(j,k)-T1(j,km2)) )/3.      &
&                   + vclim_p(j,k,ityr)*( wz(j,kp1)*(T1(j,k)-T1(j,kp1)) ) )
     end do
     k=ydim
     km1=k-1; km2=k-2
     do j = 1, xdim
        dTy(j,k) = ccy * (                                                        &
&                    -vclim_m(j,k,ityr)*( wz(j,km1)*(T1(j,k)-T1(j,km1))           &
&                                        +wz(j,km2)*(T1(j,k)-T1(j,km2)) ) )/3.
     end do

     ! longitudinal
     do k=1, ydim
        if ( dxlat(k) > 2.5e5) then  ! unitl 25degree
           j = 1
           jm1 = xdim; jm2 = xdim-1; jp1 = j+1; jp2 = j+2
           dTx(j,k)= ccx(k) * (                                                      &
&                      -uclim_m(j,k,ityr)*( wz(jm1,k)*(T1(j,k)-T1(jm1,k))            &
&                                          +wz(jm2,k)*(T1(j,k)-T1(jm2,k)) )          &
&                     + uclim_p(j,k,ityr)*( wz(jp1,k)*(T1(j,k)-T1(jp1,k))            &
&                                          +wz(jp2,k)*(T1(j,k)-T1(jp2,k)) ) )/3.
           j = 2
           jm1 = j-1; jm2 = xdim; jp1 = j+1; jp2 = j+2
           dTx(j,k)= ccx(k) * (                                                      &
&                      -uclim_m(j,k,ityr)*( wz(jm1,k)*(T1(j,k)-T1(jm1,k))            &
&                                          +wz(jm2,k)*(T1(j,k)-T1(jm2,k)) )          &
&                     + uclim_p(j,k,ityr)*( wz(jp1,k)*(T1(j,k)-T1(jp1,k))            &
&                                          +wz(jp2,k)*(T1(j,k)-T1(jp2,k)) ) )/3.
           do j=3, xdim-2              ! longitudinal
                jm1=j-1; jp1=j+1; jm2=j-2; jp2=j+2
                dTx(j,k)= ccx(k) * (                                                  &
&                           -uclim_m(j,k,ityr)*( wz(jm1,k)*(T1(j,k)-T1(jm1,k))        &
&                                               +wz(jm2,k)*(T1(j,k)-T1(jm2,k)) )      &
&                          + uclim_p(j,k,ityr)*( wz(jp1,k)*(T1(j,k)-T1(jp1,k))        &
&                                               +wz(jp2,k)*(T1(j,k)-T1(jp2,k)) ) )/3.
           end do
           j = xdim-1
           jm1 = j-1; jm2 = j-2; jp1 = j+1; jp2 = 1
           dTx(j,k)= ccx(k) * (                                                      &
&                      -uclim_m(j,k,ityr)*( wz(jm1,k)*(T1(j,k)-T1(jm1,k))            &
&                                          +wz(jm2,k)*(T1(j,k)-T1(jm2,k)) )          &
&                     + uclim_p(j,k,ityr)*( wz(jp1,k)*(T1(j,k)-T1(jp1,k))            &
&                                          +wz(jp2,k)*(T1(j,k)-T1(jp2,k)) ) )/3.
           j = xdim
           jm1 = j-1; jm2 = j-2; jp1 = 1; jp2 = 2
           dTx(j,k)= ccx(k) * (                                                      &
&                      -uclim_m(j,k,ityr)*( wz(jm1,k)*(T1(j,k)-T1(jm1,k))            &
&                                          +wz(jm2,k)*(T1(j,k)-T1(jm2,k)) )          &
&                     + uclim_p(j,k,ityr)*( wz(jp1,k)*(T1(j,k)-T1(jp1,k))            &
&                                          +wz(jp2,k)*(T1(j,k)-T1(jp2,k)) ) )/3.

        else  ! high resolution -> smaller time steps
            dd=max(1,nint(dt_crcl/(dxlat(k)/10.0/1.))); dtdff2=dt_crcl/dd
            time2=max(1,nint(float(dt_crcl)/float(dtdff2)))
            ccx2=dtdff2/dxlat(k)/2
            T1h=T1(:,k)
            do tt2=1, time2      ! additional time loop
                j = 1
                jm1=xdim; jm2=xdim-1; jm3=xdim-2; jp1=j+1; jp2=j+2; jp3=j+3
                dTxh(j)= ccx2 * (                                                              &
&                          -uclim_m(j,k,ityr)*( 10*wz(jm1,k)*(T1h(j)   - T1h(jm1) )            &
&                                               +4*wz(jm2,k)*(T1h(jm1) - T1h(jm2) )            &
&                                               +1*wz(jm3,k)*(T1h(jm2) - T1h(jm3) ) )          &
&                         + uclim_p(j,k,ityr)*( 10*wz(jp1,k)*(T1h(j)   - T1h(jp1) )            &
&                                               +4*wz(jp2,k)*(T1h(jp1) - T1h(jp2) )            &
&                                               +1*wz(jp3,k)*(T1h(jp2) - T1h(jp3) ) ) ) /20.
                j = 2
                jm1=j-1; jm2=xdim; jm3=xdim-1; jp1=j+1; jp2=j+2; jp3=j+3
                dTxh(j)= ccx2 * (                                                              &
&                          -uclim_m(j,k,ityr)*( 10*wz(jm1,k)*(T1h(j)   - T1h(jm1) )            &
&                                               +4*wz(jm2,k)*(T1h(jm1) - T1h(jm2) )            &
&                                               +1*wz(jm3,k)*(T1h(jm2) - T1h(jm3) ) )          &
&                         + uclim_p(j,k,ityr)*( 10*wz(jp1,k)*(T1h(j)   - T1h(jp1) )            &
&                                               +4*wz(jp2,k)*(T1h(jp1) - T1h(jp2) )            &
&                                               +1*wz(jp3,k)*(T1h(jp2) - T1h(jp3) ) ) ) /20.
                j = 3
                jm1=j-1; jm2=j-2; jm3=xdim; jp1=j+1; jp2=j+2; jp3=j+3
                dTxh(j)= ccx2 * (                                                              &
&                          -uclim_m(j,k,ityr)*( 10*wz(jm1,k)*(T1h(j)   - T1h(jm1) )            &
&                                               +4*wz(jm2,k)*(T1h(jm1) - T1h(jm2) )            &
&                                               +1*wz(jm3,k)*(T1h(jm2) - T1h(jm3) ) )          &
&                         + uclim_p(j,k,ityr)*( 10*wz(jp1,k)*(T1h(j)   - T1h(jp1) )            &
&                                               +4*wz(jp2,k)*(T1h(jp1) - T1h(jp2) )            &
&                                               +1*wz(jp3,k)*(T1h(jp2) - T1h(jp3) ) ) ) /20.
                do j=4, xdim-3     ! longitudinal
                    jm1=j-1; jp1=j+1; jm2=j-2; jp2=j+2; jm3=j-3; jp3=j+3
                    dTxh(j)= ccx2 * (                                                          &
&                            -uclim_m(j,k,ityr)*( 10*wz(jm1,k)*(T1h(j)   - T1h(jm1) )          &
&                                                 +4*wz(jm2,k)*(T1h(jm1) - T1h(jm2) )          &
&                                                 +1*wz(jm3,k)*(T1h(jm2) - T1h(jm3) ) )        &
&                           + uclim_p(j,k,ityr)*( 10*wz(jp1,k)*(T1h(j)   - T1h(jp1) )          &
&                                                 +4*wz(jp2,k)*(T1h(jp1) - T1h(jp2) )          &
&                                                 +1*wz(jp3,k)*(T1h(jp2) - T1h(jp3) ) ) ) /20.
                end do           ! longitudinal
                j = xdim-2
                jm1=j-1; jm2=j-2; jm3=j-3; jp1=xdim-1; jp2=xdim-1; jp3=1
                dTxh(j)= ccx2 * (                                                              &
&                          -uclim_m(j,k,ityr)*( 10*wz(jm1,k)*(T1h(j)   - T1h(jm1) )            &
&                                               +4*wz(jm2,k)*(T1h(jm1) - T1h(jm2) )            &
&                                               +1*wz(jm3,k)*(T1h(jm2) - T1h(jm3) ) )          &
&                         + uclim_p(j,k,ityr)*( 10*wz(jp1,k)*(T1h(j)   - T1h(jp1) )            &
&                                               +4*wz(jp2,k)*(T1h(jp1) - T1h(jp2) )            &
&                                               +1*wz(jp3,k)*(T1h(jp2) - T1h(jp3) ) ) ) /20.
                j = xdim-1
                jm1=j-1; jm2=j-2; jm3=j-3; jp1=xdim; jp2=1; jp3=2
                dTxh(j)= ccx2 * (                                                              &
&                          -uclim_m(j,k,ityr)*( 10*wz(jm1,k)*(T1h(j)   - T1h(jm1) )            &
&                                               +4*wz(jm2,k)*(T1h(jm1) - T1h(jm2) )            &
&                                               +1*wz(jm3,k)*(T1h(jm2) - T1h(jm3) ) )          &
&                         + uclim_p(j,k,ityr)*( 10*wz(jp1,k)*(T1h(j)   - T1h(jp1) )            &
&                                               +4*wz(jp2,k)*(T1h(jp1) - T1h(jp2) )            &
&                                               +1*wz(jp3,k)*(T1h(jp2) - T1h(jp3) ) ) ) /20.
                j = xdim
                jm1=j-1; jm2=j-2; jm3=j-3; jp1=1; jp2=2; jp3=3
                dTxh(j)= ccx2 * (                                                              &
&                          -uclim_m(j,k,ityr)*( 10*wz(jm1,k)*(T1h(j)   - T1h(jm1) )            &
&                                               +4*wz(jm2,k)*(T1h(jm1) - T1h(jm2) )            &
&                                               +1*wz(jm3,k)*(T1h(jm2) - T1h(jm3) ) )          &
&                         + uclim_p(j,k,ityr)*( 10*wz(jp1,k)*(T1h(j)   - T1h(jp1) )            &
&                                               +4*wz(jp2,k)*(T1h(jp1) - T1h(jp2) )            &
&                                               +1*wz(jp3,k)*(T1h(jp2) - T1h(jp3) ) ) ) /20.
                where(dTxh .le. -T1h ) dTxh = -0.9*T1h ! no negative q;  numerical stability
                T1h = T1h + dTxh
            end do               ! additional time loop
            dTx(:,k) = T1h - T1(:,k)
        end if
    end do          ! y-loop
    dX_advec = dTx + dTy;

end subroutine advection

!+++++++++++++++++++++++++++++++++++++++
subroutine convergence(T1, div)
!+++++++++++++++++++++++++++++++++++++++
! Calculates divergence (convergence) of a given field (i.e. spec. hum.) when omega is known
! Eq. 18 in Stassen et al 2019
  use mo_numerics, only: xdim, ydim, nstep_yr, dlon, dlat, dt_crcl
  use mo_physics,  only: ityr, rho_air, grav, pi, z_vapor, omegaclim
  implicit none

  real, dimension(xdim,ydim), intent(in)  :: T1
  real, dimension(xdim,ydim), intent(out) :: div

  real    :: w
  integer :: i, j

  do j=1,ydim
    do i=1,xdim
      !< Vertical velocity omega (Pa/s) to m/s
      w = -omegaclim(i,j,ityr) / (rho_air*grav)
      !< Convergence
      div(i,j) = T1(i,j) * w * dt_crcl / z_vapor * 2.5
    end do
  end do

end subroutine convergence

!+++++++++++++++++++++++++++++++++++++++
subroutine forcing(it, year, CO2, Tsurf)
!+++++++++++++++++++++++++++++++++++++++

  USE mo_numerics,    ONLY: xdim, ydim, ndays_yr, ndt_days, nstep_yr
  USE mo_physics,     ONLY: log_exp, sw_solar, sw_solar_ctrl, sw_solar_scnr,     &
&                           co2_part, co2_part_scn, z_topo, ityr, Tclim
  USE mo_diagnostics,  ONLY: icmn_ctrl

  ! input fields
  real, dimension(xdim,ydim)  :: Tsurf

  ! declare variables
  real, dimension(xdim,ydim) 	:: icmn_ctrl1

  ! calculate annual mean ice cover (for log_exp 44 & 45)
  icmn_ctrl1 = 0.0
  do i=1,12
  icmn_ctrl1 = icmn_ctrl1 + icmn_ctrl(:,:,i)
  end do
  icmn_ctrl1 = icmn_ctrl1/12.


  if( log_exp .eq. 10 ) CO2 = 2*340.

  if( log_exp .eq. 20 ) CO2 = 2*340.
  if( log_exp .eq. 21 ) CO2 = 4*340.
  if( log_exp .eq. 22 ) CO2 = 10*340.
  if( log_exp .eq. 23 ) CO2 = 0.5*340.
  if( log_exp .eq. 24 ) CO2 = 0*340.

  if( log_exp .eq. 25 ) CO2 = 340.+170.+170.*cos(2*3.14* (year-13.)/30.) ! sinus wave
  if( log_exp .eq. 26 ) CO2 = 2*340. ! step function
  if( log_exp .eq. 26 .and. year >= 1980 ) CO2 = 340. ! step function
  if( log_exp .eq. 27  ) CO2      = 340.
  if( log_exp .eq. 27  ) dS0      = 27.0 ! [W/m2]
  if( log_exp .eq. 27  ) sw_solar = (1365+dS0)/1365*sw_solar_ctrl ! change S0 by dS0
  if( log_exp .eq. 28  ) CO2      = 340.
  if( log_exp .eq. 28  ) dS0      = 1.0 ! amplitude solar cycle [W/m2]
  if( log_exp .eq. 28  ) sw_solar = (1365+dS0*sin(2*3.14* (year)/11.))/1365*sw_solar_ctrl ! 11yrs cycle

  if( log_exp .eq. 30  ) CO2      = 200. ! paleo 231 kyr BP
  if( log_exp .eq. 30  ) sw_solar = sw_solar_scnr !
  if( log_exp .eq. 31  ) CO2      = 340.
  if( log_exp .eq. 31  ) sw_solar = sw_solar_scnr
  if( log_exp .eq. 32  ) CO2      = 200.

! partial scenarios
! 2xCO2 NH
  if( log_exp .eq. 40 )  CO2      = 2*340.
  if( log_exp .eq. 40 )  co2_part(:,1:24)     = 0.5
! 2xCO2 SH
  if( log_exp .eq. 41 )  CO2      = 2*340.
  if( log_exp .eq. 41 )  co2_part(:,25:48)    = 0.5
! 2xCO2 tropics
  if( log_exp .eq. 42 )  CO2      = 2*340.
  if( log_exp .eq. 42 )  co2_part(:,1:15)     = 0.5
  if( log_exp .eq. 42 )  co2_part(:,33:48)    = 0.5
  if( log_exp .eq. 42 )  co2_part(4:4:96,33)  = 1.0
  if( log_exp .eq. 42 )  co2_part(4:4:96,15)  = 1.0
! 2xCO2 extra tropics
  if( log_exp .eq. 43 )  CO2      = 2*340.
  if( log_exp .eq. 43 )  co2_part(:,16:32)    = 0.5
  if( log_exp .eq. 43 )  co2_part(4:4:96,32)  = 1.0
  if( log_exp .eq. 43 )  co2_part(4:4:96,16)  = 1.0
! 2xCO2 ocean
  if( log_exp .eq. 44 )  CO2      = 2*340.
  if( log_exp .eq. 44 )  where (z_topo  > 0.) co2_part = 0.5
  if( log_exp .eq. 44 )  where (icmn_ctrl1  >= 0.5 ) co2_part = 0.5
! 2xCO2 land/ice
  if( log_exp .eq. 45 )  CO2      = 2*340.
  if( log_exp .eq. 45 )  where (z_topo  <= 0.) co2_part = 0.5
  if( log_exp .eq. 45 )  where (icmn_ctrl1  >= 0.5 ) co2_part = 1.0
! 2xCO2 winter
  if( log_exp .eq. 46 )  CO2      = 340.
  if( log_exp .eq. 46 .and. mod(it,nstep_yr) .le. 181)  CO2 = 2*340.
  if( log_exp .eq. 46 .and. mod(it,nstep_yr) .ge. 547)  CO2 = 2*340.
! 2xCO2 summer
  if( log_exp .eq. 47 )  CO2      = 2*340.
  if( log_exp .eq. 47 .and. mod(it,nstep_yr) .le. 181)  CO2 = 340.
  if( log_exp .eq. 47 .and. mod(it,nstep_yr) .ge. 547)  CO2 = 340.

! IPCC A1B scenario
  if( log_exp .eq. 95 ) then
     CO2_1950=310.;  CO2_2000=370.;  CO2_2050=520.
     if (year <= 2000.) CO2=CO2_1950 + 60./50.*(year-1950.)
     if (year  > 2000. .and. year <= 2050.) CO2=CO2_2000 + 150./50.*(year-2000.)
     if (year  > 2050. .and. year <= 2100.) CO2=CO2_2050 + 180./50.*(year-2050.)
  end if

! IPCC RCP scenarios
  if( log_exp .ge. 96 .and. log_exp .le. 99 ) then ! IPCC RCP scenarios
      if(mod(it,nstep_yr) .eq. 1) read (26,*) t, CO2
  end if

! own CO2 scenario
  if( log_exp .eq. 100 ) then
      if(mod(it,nstep_yr) .eq. 1) read (26,*) t, CO2
  end if

! Forced Climate Change run
  if( log_exp .eq. 230 ) Tsurf = Tclim(:,:,ityr) ! Keep temp on external boundary condition

! Forced ENSO run
  if( log_exp .eq. 240 .or. log_exp .eq. 241 ) Tsurf = Tclim(:,:,ityr)  ! Keep temp on external boundary condition


end subroutine forcing

!+++++++++++++++++++++++++++++++++++++++
subroutine diagonstics(it, year, CO2, ts0, ta0, to0, q0, ice_cover, sw, lw_surf, q_lat, q_sens)
!+++++++++++++++++++++++++++++++++++++++
!    diagonstics plots

  USE mo_numerics,    ONLY: ndays_yr, xdim, ydim, ipx ,ipy, ndt_days, nstep_yr
  USE mo_physics,     ONLY: ityr, TF_correct, qF_correct, cap_surf, Tclim
  use mo_diagnostics

! declare temporary fields
  real, dimension(xdim,ydim)  :: Ts0, Ta0, To0, q0, sw, ice_cover, Q_sens, Q_lat,  LW_surf

  ! diagnostics: annual means
  tsmn=tsmn+Ts0; tamn=tamn+ta0; tomn=tomn+to0; qmn=qmn+q0; icmn=icmn+ice_cover
  swmn=swmn+sw;  lwmn=lwmn+LW_surf; qlatmn=qlatmn+q_lat; qsensmn=qsensmn+Q_sens;
  ftmn=ftmn+TF_correct(:,:,ityr); fqmn=fqmn+qF_correct(:,:,ityr);
  if ( ityr == nstep_yr ) then
     tsmn    = tsmn/nstep_yr;      tamn = tamn/nstep_yr;    tomn = tomn/nstep_yr;
     qmn     = qmn/nstep_yr;
     icmn    = icmn/nstep_yr;      swmn = swmn/nstep_yr;    lwmn = lwmn/nstep_yr;
     qlatmn  = qlatmn/nstep_yr; qsensmn = qsensmn/nstep_yr; ftmn = ftmn/nstep_yr;
     fqmn    = fqmn/nstep_yr;
     1000 format (F5.0, T8, F10.1, T20, F10.2, T37, F10.6, T52, F10.6, T67, F10.6, T85, F10.6) !TB
     print 1000, year, CO2, gmean(swmn), gmean(tsmn)-273.15, tsmn(48,24)-273.15,tsmn(4,39)-273.15, tsmn(1,48)-273.15 !TB
     tsmn=0.; tamn=0.; qmn=0.; icmn=0.; swmn=0.;        ! reset annual mean values
     lwmn=0.; qlatmn=0.; qsensmn=0.; ftmn=0.; fqmn=0.;  ! reset annual mean values
  end if

end subroutine diagonstics

!+++++++++++++++++++++++++++++++++++++++
subroutine output(it, iunit, irec, mon, ts0, ta0, to0, q0, ice_cover, dq_rain, dq_eva, dq_crcl)
!+++++++++++++++++++++++++++++++++++++++
!    write output

  USE mo_numerics,     ONLY: xdim, ydim, jday_mon, ndt_days, nstep_yr, time_scnr &
&                          , time_ctrl, ireal, dt
  USE mo_physics,      ONLY: jday, log_exp, r_qviwv, wz_vapor
  use mo_diagnostics,  ONLY: Tmm, Tamm, Tomm, qmm, icmm, prmm, evamm, qcrclmm &
&                          , Tmn_ctrl, Tamn_ctrl, Tomn_ctrl, qmn_ctrl, icmn_ctrl, prmn_ctrl, evamn_ctrl, qcrclmn_ctrl

  ! declare temporary fields
  real, dimension(xdim,ydim)  :: Ts0, Ta0, To0, q0, ice_cover, dq_rain, dq_eva, dq_crcl

  ! diagnostics: monthly means
  Tmm=Tmm+Ts0; Tamm=Tamm+ta0; Tomm=Tomm+to0; qmm=qmm+q0; icmm=icmm+ice_cover;
  prmm=prmm+dq_rain*(r_qviwv*wz_vapor);          ! kg/m2/s
  evamm=evamm+dq_eva*(r_qviwv*wz_vapor);         ! kg/m2/s
  qcrclmm=qcrclmm+dq_crcl*(r_qviwv*wz_vapor)/dt; ! kg/m2/s

! control output
  if (       jday == sum(jday_mon(1:mon))                   &
&      .and. it/float(ndt_days) == nint(it/float(ndt_days)) &
&      .and. iunit == 101 ) then
     ndm=jday_mon(mon)*ndt_days
     if (it/float(ndt_days)  > 365*(time_ctrl-1)) then
     if (log_exp .eq. 1 .or. log_exp .eq. 230 ) then
     irec=irec+1;
     write(iunit,rec=8*irec-7)  Tmm/ndm
     write(iunit,rec=8*irec-6)  Tamm/ndm
     write(iunit,rec=8*irec-5)  Tomm/ndm
     write(iunit,rec=8*irec-4)   qmm/ndm
     write(iunit,rec=8*irec-3)  icmm/ndm
     write(iunit,rec=8*irec-2)  prmm/ndm
     write(iunit,rec=8*irec-1) evamm/ndm
     write(iunit,rec=8*irec) qcrclmm/ndm
     else
     Tmn_ctrl(:,:,mon)     =   Tmm/ndm
     Tamn_ctrl(:,:,mon)    =  Tamm/ndm
     Tomn_ctrl(:,:,mon)    =  Tomm/ndm
     qmn_ctrl(:,:,mon)     =   qmm/ndm
     icmn_ctrl(:,:,mon)    =  icmm/ndm
     prmn_ctrl(:,:,mon)    =  prmm/ndm
     evamn_ctrl(:,:,mon)   = evamm/ndm
     qcrclmn_ctrl(:,:,mon) = qcrclmm/ndm
     end if
     end if
     Tmm=0.; Tamm=0.;Tomm=0.; qmm=0.; icmm=0.; prmm=0.; evamm=0.; qcrclmm=0.;
     mon=mon+1; if (mon==13) mon=1
  end if

! scenario output
  if (       jday == sum(jday_mon(1:mon))                   &
&      .and. it/float(ndt_days) == nint(it/float(ndt_days)) &
&      .and. iunit == 102 ) then
     ndm=jday_mon(mon)*ndt_days
     irec=irec+1;
     if (log_exp .ge. 35 .and. log_exp .le. 37 ) then
     write(iunit,rec=8*irec-7)   Tmm/ndm
     write(iunit,rec=8*irec-6)  Tamm/ndm
     write(iunit,rec=8*irec-5)  Tomm/ndm
     write(iunit,rec=8*irec-4)   qmm/ndm
     write(iunit,rec=8*irec-3)  icmm/ndm
     write(iunit,rec=8*irec-2)  prmm/ndm
     write(iunit,rec=8*irec-1) evamm/ndm
     write(iunit,rec=8*irec) qcrclmm/ndm

     else
     write(iunit,rec=8*irec-7)   Tmm/ndm  -Tmn_ctrl(:,:,mon)
     write(iunit,rec=8*irec-6)  Tamm/ndm -Tamn_ctrl(:,:,mon)
     write(iunit,rec=8*irec-5)  Tomm/ndm -Tomn_ctrl(:,:,mon)
     write(iunit,rec=8*irec-4)   qmm/ndm -qmn_ctrl(:,:,mon)
     write(iunit,rec=8*irec-3)  icmm/ndm -icmn_ctrl(:,:,mon)
     write(iunit,rec=8*irec-2)  prmm/ndm -prmn_ctrl(:,:,mon)
     write(iunit,rec=8*irec-1) evamm/ndm -evamn_ctrl(:,:,mon)
     write(iunit,rec=8*irec) qcrclmm/ndm -qcrclmn_ctrl(:,:,mon)

     iyrec=floor(float((irec-1)/12))
     write(103,rec=               12*iyrec+mon) gmean(Tmm/ndm  -Tmn_ctrl(:,:,mon))
     write(103,rec=1*12*time_scnr+12*iyrec+mon) gmean(Tamm/ndm -Tamn_ctrl(:,:,mon))
     write(103,rec=2*12*time_scnr+12*iyrec+mon) gmean(Tomm/ndm -Tomn_ctrl(:,:,mon))
     write(103,rec=3*12*time_scnr+12*iyrec+mon) gmean(qmm/ndm -qmn_ctrl(:,:,mon))
     write(103,rec=4*12*time_scnr+12*iyrec+mon) gmean(icmm/ndm -icmn_ctrl(:,:,mon))
     write(103,rec=5*12*time_scnr+12*iyrec+mon) gmean(prmm/ndm -prmn_ctrl(:,:,mon))
     write(103,rec=6*12*time_scnr+12*iyrec+mon) gmean(evamm/ndm -evamn_ctrl(:,:,mon))
     write(103,rec=7*12*time_scnr+12*iyrec+mon) gmean(qcrclmm/ndm -qcrclmn_ctrl(:,:,mon))
     end if
     Tmm=0.; Tamm=0.;Tomm=0.; qmm=0.; icmm=0.; prmm=0.; evamm=0.; qcrclmm=0.;
     mon=mon+1; if (mon==13) mon=1
  end if

end subroutine output

!TB
!+++++++++++++++++++++++++++++++++++++++
function gmean(data)
!+++++++++++++++++++++++++++++++++++++++

use mo_numerics,		ONLY: xdim, ydim, dlat

! declare variables
real, dimension(xdim,ydim) 	:: data, w
real, dimension(ydim)		:: lat

do i=1,ydim
lat(i) = -90-(dlat*0.5)+i*dlat
end do
do i=1,xdim
w(i,:) = cos(2.*3.14159*lat/360.)
end do

gmean = sum(data*w)/sum(w)

end function
