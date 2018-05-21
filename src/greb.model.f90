!
!----------------------------------------------------------
!   The Globally Resolved Energy Balance (GREB) Model
!----------------------------------------------------------
!
!   Authors; Dietmar Dommenget, Janine Flöter and Christian Stassen
!            with numerical opitmizations by Micheal Rezny
!
!   Reference: - Conceptual Understanding of Climate Change with a Globally Resolved Energy Balance Model
!                by Dietmar Dommenget and Janine Flöter, submitted to Climate Dynamics 2010.
!              - A Hydrological Cycle Model for the Globally Resolved Energy Balance Model (GREB)
!                by Christian Stassen, Dietmar Dommenget & Nicholas Loveday
!
!
!  input fields: The GREB model needs the following fields to be specified before
!                the main subroutine greb_model is called:
!
!   z_topo(xdim,ydim):            topography (<0 are ocean points) [m]
!  glacier(xdim,ydim):            glacier mask ( >0.5 are glacier points )
!    Tclim(xdim,ydim,nstep_yr):   mean Tsurf                       [K]
!    uclim(xdim,ydim,nstep_yr):   mean zonal wind speed            [m/s]
!    vclim(xdim,ydim,nstep_yr):   mean meridional wind speed       [m/s]
!    omega_clim(xdim,ydim,nstep_yr):   mean vertical velocity (pressure) [hPa/s]
!omega_std_clim(xdim,ydim,nstep_yr):   standard dev. vertical velocity (pressure) [hPa/s]
!  abswind(xdim,ydim,nstep_yr):   wind speed  [m/s]
!    qclim(xdim,ydim,nstep_yr):   mean atmospheric humidity        [kg/kg]
!  mldclim(xdim,ydim,nstep_yr):   mean ocean mixed layer depth     [m]
!   Toclim(xdim,ydim,nstep_yr):   mean deep ocean temperature      [K]
! swetclim(xdim,ydim,nstep_yr):   soil wetnees, fraction of total  [0-1]
! sw_solar(ydim,nstep_yr):        24hrs mean solar radiation       [W/m^2]
!
!
!+++++++++++++++++++++++++++++++++++++++
module mo_numerics
!+++++++++++++++++++++++++++++++++++++++
  implicit none

! numerical parameter
  integer, parameter :: xdim = 96, ydim = 48          ! field dimensions
  integer, parameter :: ndays_yr  = 365               ! number of days per year
  integer, parameter :: dt        = 12*3600           ! time step [s]
  integer            :: dt_crcl   = 0.5*3600          ! time step circulation [s]
  integer, parameter :: ndt_days  = 24*3600/dt        ! number of timesteps per day
  integer, parameter :: nstep_yr  = ndays_yr*ndt_days ! number of timesteps per year
  integer            :: time_flux = 0                 ! length of integration for flux correction [yrs]
  integer            :: time_ctrl = 0                 ! length of integration for control run  [yrs]
  integer            :: time_scnr = 0                 ! length of integration for scenario run [yrs]
  integer            :: ipx       = 1                 ! points for diagonstic print outs
  integer            :: ipy       = 1                 ! points for diagonstic print outs
  integer, parameter, dimension(12) :: jday_mon = (/31,28,31,30,31,30,31,31,30,31,30,31/) ! days per
  real,    parameter :: dlon      = 360./xdim         ! linear increment in lon
  real,    parameter :: dlat      = 180./ydim         ! linear increment in lat
  integer            :: ireal     = 4         ! record length for IO (machine dependent)
! 												ireal = 4 for Mac Book Pro

  namelist / numerics / time_flux, time_ctrl, time_scnr

end module mo_numerics

!+++++++++++++++++++++++++++++++++++++++
module mo_physics
!+++++++++++++++++++++++++++++++++++++++
  use mo_numerics
  implicit none

  !< Main switches
  integer          :: log_exp    = 0               !< process control logics for sens. exp.
  integer          :: log_rain   = 0               !< process control logics for rainfall model
  integer          :: log_eva    = 0               !< process control logics for evaporation model
  integer          :: log_qflux  = 0               !< process control logics for qflux corrections
  integer          :: log_adv_scheme = 0           !< process control logic for different advection schemes in GREB
  integer          :: log_output = 1               !< process control logics for output
  integer          :: log_clim   = 1               !< process control logic for reference climatology

  !< Forcing switches
  logical          :: log_forceTclim = .FALSE.     !< Logical switch to force surface temperature

  ! physical parameter (natural constants)
  real, parameter :: pi        = 3.1416
  real, parameter :: deg2rad   = pi/180.
  real, parameter :: sig       = 5.6704e-8       ! stefan-boltzmann constant [w/m^2/k^4]
  real, parameter :: rho_ocean = 999.1           ! density of water at t=15c [kg/m^2]
  real, parameter :: rho_land  = 2600.           ! density of solid rock [kg/m^2]
  real, parameter :: rho_air   = 1.2             ! density of air at 20C at NN
  real, parameter :: cp_ocean  = 4186.           ! specific heat capacity of water at T=15C [J/kg/K]
  real, parameter :: cp_land   = cp_ocean/4.5    ! specific heat capacity of dry land [J/kg/K]
  real, parameter :: cp_air    = 1005.           ! specific heat capacity of air      [J/kg/K]
  real, parameter :: eps       = 1.              ! emissivity for IR
  real, parameter :: grav      = 9.81            ! gravitational acceleration [m/s^2]
  real, parameter :: R_d       = 287.04          ! specific gas constant dry air [J/K/kg]

  ! physical parameter (model values)
  real, parameter :: d_ocean   = 50.                       ! depth of ocean column [m]
  real, parameter :: d_land    = 2.                        ! depth of land column  [m]
  real, parameter :: d_air     = 5000.                     ! depth of air column   [m]
  real, parameter :: cap_ocean = cp_ocean*rho_ocean        ! heat capacity 1m ocean  [J/K/m^2]
  real, parameter :: cap_land  = cp_land*rho_land*d_land   ! heat capacity land   [J/K/m^2]
  real, parameter :: cap_air   = cp_air*rho_air*d_air      ! heat capacity air    [J/K/m^2]
  real, parameter :: ct_sens   = 22.5                      ! coupling for sensible heat
  real, parameter :: da_ice    = 0.25                      ! albedo diff for ice covered points
  real, parameter :: a_no_ice  = 0.1                       ! albedo for non-ice covered points
  real, parameter :: a_cloud   = 0.35                      ! albedo for clouds
  real, parameter :: Tl_ice1   = 273.15-10.                ! temperature range of land snow-albedo feedback
  real, parameter :: Tl_ice2   = 273.15                    ! temperature range of land snow-albedo feedback
  real, parameter :: To_ice1   = 273.15-7.                 ! temperature range of ocean ice-albedo feedback
  real, parameter :: To_ice2   = 273.15-1.7                ! temperature range of ocean ice-albedo feedback
  real, parameter :: co_turb   = 5.0                       ! turbolent mixing to deep ocean [W/K/m^2]
  real, parameter :: kappa     = 8e5                       ! atmos. diffusion coefficient [m^2/s]
  real, parameter :: ce        = 2e-3                      ! latent heat transfer coefficient for ocean
  real, parameter :: cq_latent = 2.257e6                   ! latent heat of condensation/evapoartion f water [J/kg]
  real, parameter :: cq_rain   = -0.1/24./3600.            ! decrease in air water vapor due to rain [1/s]
  real, parameter :: z_air     = 8400.                     ! scaling height atmos. heat, CO2
  real, parameter :: z_vapor   = 5000.                     ! scaling height atmos. water vapor diffusion
  real            :: r_qviwv   = 2.6736e3                  ! regres. factor between viwv and q_air  [kg/m^3]

  ! Rainfall Parameters
  real :: c_q, c_rq, c_omega, c_omegastd

  ! parameter emisivity
  real, parameter, dimension(10) :: p_emi = (/9.0721, 106.7252, 61.5562, 0.0179, 0.0028,     &
  &                                             0.0570, 0.3462, 2.3406, 0.7032, 1.0662/)

  ! declare climate fields
  real, dimension(xdim,ydim)          ::  z_topo, glacier,z_ocean
  real, dimension(xdim,ydim,nstep_yr) ::  Tclim, uclim, vclim, omega_clim, omega_std_clim, abswind_clim, &
                                      &   qclim, mldclim, Toclim, cldclim
  real, dimension(xdim,ydim,nstep_yr) ::  TF_correct, qF_correct, ToF_correct, swetclim, dTrad
  real, dimension(ydim,nstep_yr)      ::  sw_solar

  ! declare constant fields
  real, dimension(xdim,ydim)          ::  cap_surf
  integer jday, ityr

  ! declare some program constants
  real, dimension(xdim, ydim)         :: wz_air, wz_vapor
  real, dimension(xdim,ydim,nstep_yr) :: uclim_m, uclim_p
  real, dimension(xdim,ydim,nstep_yr) :: vclim_m, vclim_p

  real :: t0, t1, t2

  namelist / physics / log_exp, log_rain, log_eva, log_adv_scheme, log_qflux, log_clim

end module mo_physics

!+++++++++++++++++++++++++++++++++++++++
module mo_diagnostics
!+++++++++++++++++++++++++++++++++++++++

  use mo_numerics,    only: xdim, ydim

  ! declare diagnostic fields
  real, dimension(xdim,ydim)          :: Tsmn, Tamn, qmn, swmn, lwmn, qlatmn, qsensmn, &
  &                                        ftmn, fqmn, amn, Tomn

  ! declare output fields
  real, dimension(xdim,ydim)          :: Tmm, Tamm, Tomm, qmm, apmm, Prmm  !< monthly means

  real, dimension(xdim,ydim)          :: dqmm, dqrainmm, dqevamm, &  !< monthly means of water vapor tendency (dqmm)
  & dqcrclmm, dq_advmm, dq_diffmm

end module mo_diagnostics

!+++++++++++++++++++++++++++++++++++++++
subroutine greb_model
!+++++++++++++++++++++++++++++++++++++++
!   climate model main loop

  ! external modules to be used by the main routine
  use mo_numerics
  use mo_physics
  use mo_diagnostics

  ! declare temporary fields
  real, dimension(xdim,ydim) :: Ts0, Ts1, Ta0, Ta1, To0, To1, q0, q1,       &
  &                               ts_ini, ta_ini, q_ini, to_ini
  integer            :: year, mon

  ! Naming of the output files
  character(len=4)   :: clog_exp, clog_rain, clog_eva, clog_qflux
  character(len=120) :: csetting, fout_ctrl, fout_scen, fout_flux, add_name
  character(len=4)   :: year_strt, year_last, fformat
  logical            :: log_scen                                  !< Scenario run yes/no

  ! Experiment details to name
  write(clog_exp,'(i4)') log_exp
  write(clog_rain,'(i4)') log_rain
  write(clog_eva,'(i4)') log_eva
  write(clog_qflux,'(i4)') log_qflux

  csetting = ".Exp" // TRIM(ADJUSTL(clog_exp)) // ".P" // TRIM(ADJUSTL(clog_rain)) // &
  & '.E' // TRIM(ADJUSTL(clog_eva)) // '.F' // TRIM(ADJUSTL(clog_qflux))

  !< Output name
  add_name = ''
  if (log_output == 0)     add_name = trim(adjustl(add_name)) //'.' // 'direct'
  if (log_adv_scheme == 0) add_name = trim(adjustl(add_name)) //'.' // 'origadv'
  if (log_adv_scheme == 1) add_name = trim(adjustl(add_name)) //'.' // 'convergence'
  if (log_clim == 0) add_name = trim(adjustl(add_name)) // '.' // 'origclim'
  if (log_clim == 1) add_name = trim(adjustl(add_name)) // '.' // 'newclim'

  dTrad = -0.16*Tclim -5. ! offset Tatmos-rad

  ! set ocean depth
  z_ocean=0
  do i=1,nstep_yr
    where(mldclim(:,:,i).gt.z_ocean) z_ocean = mldclim(:,:,i)
  end do
  z_ocean = 3.0*z_ocean

  if (log_exp ==  1) where(z_topo > 1.) z_topo = 1.0      ! sens. exp. constant topo
  if (log_exp <=  2) cldclim = 0.7                        ! sens. exp. constant cloud cover
  if (log_exp <=  3) qclim   = 0.0052                     ! sens. exp. constant water vapor
  if (log_exp <=  9) mldclim = d_ocean                    ! sens. exp. no deep ocean
  if (log_exp == 11) mldclim = d_ocean                    ! sens. exp. no deep ocean

  ! heat capacity global [J/K/m^2]
  where (z_topo  > 0.) cap_surf = cap_land
  where (z_topo <= 0.) cap_surf = cap_ocean*mldclim(:,:,1)

  ! initialize fields
  Ts_ini   = Tclim(:,:,nstep_yr)                          ! initial value temp. surf
  Ta_ini   = Ts_ini                                       ! initial value atm. temp.
  To_ini   = Toclim(:,:,nstep_yr)                         ! initial value temp. surf
  q_ini    = qclim(:,:,nstep_yr)                          ! initial value atmos water vapor

  CO2_ctrl = 280.
  if ( log_exp == 12 .or. log_exp == 13  ) CO2_ctrl = 298.  ! A1B scenario

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

  !< initialize the rainfall parameterisation
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
      c_q=-1.88; c_rq=2.25; c_omega=-17.69; c_omegastd=59.07 !< Rainfall parameters (ERA-Interim)
      if (log_clim == 0) then
        c_q=-1.27; c_rq=1.99; c_omega=-16.54; c_omegastd=21.15 !< Rainfall parameters (NCEP)
      end if
  end select

  ! compute Q-flux corrections
  print*, ''
  print*,'% flux correction ', CO2_ctrl
  1001 format (A4,     T8, A10,   T20, A10,    T32, A15,         T50, A12,      T65, A12,     T80, A15, T100, A20) !TB
  print 1001, "YEAR", "CO2[ppm]", "SW[W/m^2]", "global mean[C]", "Trop Pac[C]", "Hamburg[C]", &
  & "North Pole[C]", 'Water vapour [kg/kg]' !TB

  write(year_strt,'(i4.4)')  0
  write(year_last,'(i4.4)')  0
  fout_flux = '../output/fluxcorr'//trim(adjustl(csetting))// '.' // &
  & trim(adjustl(year_strt))//'-'//trim(adjustl(year_last))//trim(adjustl(add_name))//trim(adjustl(fformat))
  open(20, file=fout_flux, ACCESS='DIRECT',FORM='UNFORMATTED', STATUS='REPLACE', RECL=ireal*xdim*ydim)
  irec=0
  call qflux_correction(fout_flux, 20, irec, CO2_ctrl, Ts_ini, Ta_ini, q_ini, To_ini)

  ! control run
  print*, ''
  print*,'% CONTROL RUN CO2=',CO2_ctrl,'  time=', time_ctrl,'yr'
  print 1001, "YEAR", "CO2[ppm]", "SW[W/m^2]", "global mean[C]", "Trop Pac[C]", "Hamburg[C]", &
  & "North Pole[C]", 'Water vapour [kg/kg]'
  Ts1 = Ts_ini; Ta1 = Ta_ini; To1 = To_ini; q1 = q_ini;                   ! initialize fields
  mon=1; year=1970; irec=0; Tmm=0.; Tamm=0.; qmm=0.; apmm=0.;
  log_scen=.FALSE.
  write(year_strt,'(i4.4)')  year
  write(year_last,'(i4.4)')  year+time_ctrl-1
  fout_ctrl = '../output/control'//trim(adjustl(csetting))// '.' // &
  & trim(adjustl(year_strt))//'-'//trim(adjustl(year_last))//trim(adjustl(add_name))//trim(adjustl(fformat))
  open(21, file=fout_ctrl, ACCESS='DIRECT',FORM='UNFORMATTED', STATUS='REPLACE', RECL=ireal*xdim*ydim)
  do it=1, time_ctrl*nstep_yr                                             ! main time loop
    call time_loop(it, isrec, fout_ctrl, year, CO2_ctrl, irec, mon, 21, Ts1, Ta1, q1, To1, Ts0,Ta0, q0, To0, log_scen )
    Ts1=Ts0; Ta1=Ta0; q1=q0; To1=To0
  end do

  ! scenario run
  print*, ''
  print*,'% SCENARIO EXP: ',log_exp,'  time=', time_scnr,'yr'
  print 1001, "YEAR", "CO2[ppm]", "SW[W/m^2]", "global mean[C]", "Trop Pac[C]", "Hamburg[C]", &
  & "North Pole[C]", 'Water vapour [kg/kg]'
  Ts1 = Ts_ini; Ta1 = Ta_ini; q1 = q_ini; To1 = To_ini                     ! initialize fields
  year=1970; CO2=CO2_ctrl; mon=1; irec=0; Tmm=0.; Tamm=0.; qmm=0.; apmm=0.;
  write(year_strt,'(i4.4)')  year
  write(year_last,'(i4.4)')  year+time_scnr-1
  fout_scen = '../output/scenario'//trim(adjustl(csetting))// '.' // &
  & trim(adjustl(year_strt))//'-'//trim(adjustl(year_last))//trim(adjustl(add_name))//trim(adjustl(fformat))
  open(22, file=fout_scen, ACCESS='DIRECT',FORM='UNFORMATTED', STATUS='REPLACE', RECL=ireal*xdim*ydim)
  log_scen=.TRUE.

  !< Double CO2 experiments
  if ( log_exp >= 1 .AND. log_exp <= 11 ) then
    CO2 = CO2_ctrl * 2.
  end if

  !< ENSO experiments
  if ( log_exp >= 30 .AND. log_exp <= 39 ) then
    call enso(Tclim)
  end if

  !< Climate change experiments
  if ( log_exp >= 40 .AND. log_exp <= 49 ) then
    call climate_change(Tclim)
  end if

  do it=1, time_scnr*nstep_yr                                              ! main time loop
    !< Call CO2 levels
    if( log_exp .eq. 12 .or.  log_exp .eq. 13 ) call co2_level(it, year, CO2)

    ! sens. exp. SST+1
    if(log_exp >= 14 .and. log_exp <= 16) where (z_topo < 0.0) Ts1 = Tclim(:,:,ityr)+1.0

    !< Time loop
    call time_loop(it,isrec, fout_scen, year, CO2, irec, mon, 22, Ts1, Ta1, q1, To1, Ts0, Ta0, q0, To0, log_scen )

    !< Update variables
    Ts1=Ts0; Ta1=Ta0; q1=q0; To1=To0
    if (mod(it,nstep_yr) == 0) year=year+1
  end do

end subroutine

!+++++++++++++++++++++++++++++++++++++++
subroutine time_loop(it, isrec, fout, year, CO2, irec, mon, ionum, Ts1, Ta1, q1, To1, Ts0,Ta0, q0, To0, log_scen)
!+++++++++++++++++++++++++++++++++++++++
! main time loop
  use mo_numerics
  use mo_physics

  implicit none

  character(len=*), INTENT(in) :: fout
  logical, INTENT(in)          :: log_scen

  integer                      :: it, isrec, year, irec, mon, ionum
  real                         :: CO2

  real, dimension(xdim,ydim):: Ts1, Ta1, q1, To1, Ts0,Ta0, q0, To0, sw,       &
&                              albedo, Q_sens, Q_lat, Q_lat_air, dq_eva,      &
&                              dq_rain, dTa_crcl, dq_crcl, dq, dT_ocean, dTo, &
&                              LW_surf, LWair_down, LWair_up, em,             &
&                              dTs_diab, dTs, dTa_diab

  jday = mod((it-1)/ndt_days,ndays_yr)+1  ! current calendar day in year
  ityr = mod((it-1),nstep_yr)+1           ! time step in year

  call tendencies(CO2, Ts1, Ta1, To1, q1, albedo, SW, LW_surf, Q_lat,   &
  &                    Q_sens, Q_lat_air, dq_eva, dq_rain, dq_crcl,       &
  &                    dTa_crcl, dT_ocean, dTo, LWair_down, LWair_up, em)

  ! surface temperature
  dTs_diab = dt*( SW +LW_surf -LWair_down +Q_lat +Q_sens) / cap_surf
  dTs      = dTs_diab + dT_ocean
  Ts0 = Ts1 +dTs +dt*TF_correct(:,:,ityr) / cap_surf
  !< Keep Ts0 prescribed
  if ( log_scen .AND. log_forceTclim )  Ts0 = Tclim(:,:,ityr)
  ! air temperature
  dTa_diab = dt*( LWair_up +LWair_down -em*LW_surf +Q_lat_air -Q_sens )/cap_air
  Ta0 = Ta1 +dTa_crcl +dTa_diab
  ! deep ocean temperature
  To0 = To1 +dTo +ToF_correct(:,:,ityr)
  ! air water vapor
  dq  = dt*(dq_eva+dq_rain) +dq_crcl +qF_correct(:,:,ityr)
  where(dq .le. -1.*q1 ) dq = -0.9*q1 ! no negative q;  numerical stability
  q0  = q1 + dq
  ! sea ice heat capacity
  call seaice(Ts0)
  ! write output
  call output( it, ionum, fout, irec, mon, ts0, ta0, to0, q0, albedo, dq_eva, dq_rain, dq_crcl )
  ! diagnostics: annual means plots
  call diagonstics(it, year, CO2, ts0, ta0, to0, q0, albedo, sw, lw_surf, q_lat, q_sens)

end subroutine time_loop

!+++++++++++++++++++++++++++++++++++++++
subroutine tendencies(CO2, Ts1, Ta1, To1, q1, albedo, SW, LW_surf, Q_lat, Q_sens, Q_lat_air, dq_eva,   &
  &                   dq_rain, dq_crcl, dTa_crcl, dT_ocean, dTo, LWair_down, LWair_up, em)
!+++++++++++++++++++++++++++++++++++++++

  use mo_numerics
  use mo_physics

  ! declare temporary fields
  real, dimension(xdim,ydim) :: Ts1, Ta1, To1, q1, albedo, sw, LWair_up,      &
  &                               LWair_down, em, Q_sens, Q_lat, Q_lat_air,   &
  &                               dq_eva, dq_rain, dq_crcl, dq_adv, dq_diff,  &
  &                               dTa_crcl, dTa_adv, dTa_diff, LW_surf, dT_ocean, dTo

  ! SW radiation model
  call SWradiation(Ts1, sw, albedo)
  ! LW radiation model
  call LWradiation(Ts1, Ta1, q1, CO2, LW_surf, LWair_up, LWair_down, em)
  ! sensible heat flux
  Q_sens = ct_sens*(Ta1-Ts1)
  ! hydro. model
  call hydro(Ts1, q1, Q_lat, Q_lat_air, dq_eva, dq_rain)
  ! atmos. circulation
  call circulation(Ta1, dTa_crcl, dTa_diff, dTa_adv, z_air, wz_air)     ! air temp
  call circulation( q1,  dq_crcl, dq_diff, dq_adv, z_vapor, wz_vapor)   ! atmos water vapor
  ! deep ocean interaction
  call deep_ocean(Ts1, To1, dT_ocean, dTo)

end subroutine tendencies

!+++++++++++++++++++++++++++++++++++++++
subroutine  qflux_correction(fout, iunit, irec, CO2_ctrl, Ts1, Ta1, q1, To1)
!+++++++++++++++++++++++++++++++++++++++
! compute flux correction values

  use mo_numerics
  use mo_physics

  character(len=*)            :: fout
  integer, intent(in)         :: iunit
  integer, intent(inout)      :: irec

  ! declare temporary fields
  real, dimension(xdim,ydim) :: Ts0, Ts1, Ta0, Ta1, To0, To1, q0, q1, sw, albedo,      &
  &                               Q_sens, Q_lat, Q_lat_air, dq_eva, dq_rain, LW_surf,  &
  &                               LWair_down, LWair_up, em, dTa_crcl, dq_crcl, dTs,    &
  &                               dTa_diab, dq, T_error, dT_ocean, dTo, annual_mean
  ! Init
  q0=q1; Ts0=Ts1
  ! time loop
  do it=1, time_flux*ndt_days*ndays_yr

    jday = mod((it-1)/ndt_days,ndays_yr)+1  ! current calendar day in year
    ityr = mod((it-1),nstep_yr)+1           ! time step in year
    call tendencies(CO2_ctrl, Ts1, Ta1, To1, q1, albedo, SW, LW_surf, Q_lat,    &
    &                    Q_sens, Q_lat_air, dq_eva, dq_rain, dq_crcl, dTa_crcl, &
    &                    dT_ocean, dTo, LWair_down, LWair_up, em)

    ! surface temperature without heat flux correction
    dTs = dt*( sw +LW_surf -LWair_down +Q_lat +Q_sens) / cap_surf
    Ts0  = Ts1 +dTs +dT_ocean
    ! air temperature
    dTa_diab = dt*( LWair_up +LWair_down -em*LW_surf +Q_lat_air -Q_sens)/cap_air
    Ta0  = Ta1 + dTa_diab +dTa_crcl
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
    call diagonstics(it, 0, CO2_ctrl, ts0, ta0, to0, q0, albedo, sw, lw_surf, q_lat, q_sens)
    ! memory
    Ts1=Ts0; Ta1=Ta0; q1=q0;  To1=To0;
  end do

  !< Call output
  !< Time varying Qflux corections
  if( log_qflux == 0 ) THEN
    do it=1,nstep_yr
      call output_fluxcorr(it, iunit, fout, irec, TF_correct(:,:,it)/ cap_surf, ToF_correct(:,:,it)/dt, qF_correct(:,:,it)/dt)
    end do
  end if

  !< Annual mean qflux experiments (specific humidity only)
  if( log_qflux == 1 ) THEN
    annual_mean = SUM(qF_correct, DIM=3) / nstep_yr

    do it=1,nstep_yr
      qF_correct(:,:,it) = annual_mean
      call output_fluxcorr(it, iunit, fout, irec, TF_correct(:,:,it)/ cap_surf, ToF_correct(:,:,it)/dt, qF_correct(:,:,it)/dt)
    end do
  end if

  !< Annual mean qflux experiments
  if( log_qflux == 2 ) THEN
    annual_mean = SUM(qF_correct, DIM=3) / nstep_yr

    do it=1,nstep_yr
      qF_correct(:,:,it) = annual_mean
      call output_fluxcorr(it, iunit, fout, irec, TF_correct(:,:,it)/ cap_surf, ToF_correct(:,:,it)/dt, qF_correct(:,:,it)/dt)
    end do

    annual_mean = SUM(TF_correct, DIM=3) / nstep_yr

    do it=1,nstep_yr
      TF_correct(:,:,it) = annual_mean
      call output_fluxcorr(it, iunit, fout, irec, TF_correct(:,:,it)/ cap_surf, ToF_correct(:,:,it)/dt, qF_correct(:,:,it)/dt)
    end do
  end if

  1003 format ("On global average (RMS) a heat flux correction of ", F8.2," W/m^2") !TB
  1004 format ("and a water vapour correction (RMS) of ", F8.4, " g/kg is applied each time step") !TB

  print 1003, grms(sum(TF_correct,3)/nstep_yr) !TB
  print 1004, grms(sum(qF_correct,3)/nstep_yr)*100 !TB

end subroutine qflux_correction

!+++++++++++++++++++++++++++++++++++++++
subroutine hydro( Tsurf, q, Qlat, Qlat_air, dq_eva, dq_rain )
!+++++++++++++++++++++++++++++++++++++++
! hydrological model for latent heat and water vapor
    use mo_numerics
    use mo_physics

    implicit none

    !< declare temporary fields
    real, dimension(xdim,ydim)  :: Tsurf, q, qs, Qlat, Qlat_air, abswind
    real, dimension(xdim,ydim)  :: dq_eva, dq_rain, dq2mmday, Tskin

    real, dimension(xdim,ydim)  :: rq

    !< Init. values to zero
    Qlat=0.; Qlat_air=0.; dq_eva=0.; dq_rain=0.

    !< Maybe return without hydro cycle
    if(log_exp <=  6 .or. log_exp == 13 .or. log_exp == 15) return

    !< saturated humiditiy (max. air water vapor)
    qs = 3.75e-3*exp(17.08085*(Tsurf-273.15)/(Tsurf-273.15+234.175));
    qs = qs*exp(-z_topo/z_air) ! scale qs by topography

    !< relative humidity
    rq = q/qs

    !< Evaporation
    if ( log_eva == -1 ) then
      abswind = sqrt(uclim(:,:,ityr)**2 +vclim(:,:,ityr)**2)

      ! turbulent factor
      where(z_topo > 0. ) abswind = sqrt(abswind**2 + 2.0**2) !< land
      where(z_topo < 0. ) abswind = sqrt(abswind**2 + 3.0**2) !< ocean

      ! Evaporation
      Qlat   = (q-qs)*abswind*cq_latent*rho_air*ce*swetclim(:,:,ityr)

    else if ( log_eva == 1 ) then
      !< Stat. model evaporation
      where(z_topo > 0. ) abswind = sqrt(abswind**2 + 144.**2) ! land
      where(z_topo < 0. ) abswind = sqrt(abswind**2 + 7.1**2) ! ocean

      where(z_topo > 0. )  Qlat   = (q-qs)*abswind*cq_latent*rho_air*0.04*ce*swetclim(:,:,ityr) !< land
      where(z_topo <= 0. ) Qlat   = (q-qs)*abswind*cq_latent*rho_air*0.73*ce*swetclim(:,:,ityr) !< ocean


    else if ( log_eva == 2 ) then
      !< Absolute wind
      abswind = abswind_clim(:,:,ityr)

      !< Stat. model evaporation
      where(z_topo > 0. )  abswind = sqrt(abswind**2 + 9.0**2) ! land
      where(z_topo <= 0. ) abswind = sqrt(abswind**2 + 4.0**2) ! ocean

      where(z_topo > 0. )  Qlat   = (q-qs)*abswind*cq_latent*rho_air*0.56*ce*swetclim(:,:,ityr) !< land
      where(z_topo <= 0. ) Qlat   = (q-qs)*abswind*cq_latent*rho_air*0.79*ce*swetclim(:,:,ityr) !< ocean

    else if ( log_eva == 0 ) then
      !< maybe change skin temperature for evaporation
      where(z_topo > 0. )  Tskin = Tsurf + 5. !< land
      where(z_topo <= 0. ) Tskin = Tsurf + 1. !< ocean

      ! Re-calculate saturation pressure
      qs = 3.75e-3*exp(17.08085*(Tskin-273.15)/(Tskin-273.15+234.175));
      qs = qs*exp(-z_topo/z_air) ! scale qs by topography

      !< Absolute wind (climatology)
      where(z_topo > 0. ) abswind = sqrt(abswind_clim(:,:,ityr)**2 + 11.5**2)
      where(z_topo <= 0. ) abswind = sqrt(abswind_clim(:,:,ityr)**2 + 5.4**2)

      where(z_topo > 0. )  Qlat   = (q-qs)*abswind*cq_latent*rho_air*0.25*ce*swetclim(:,:,ityr) !< land
      where(z_topo <= 0. ) Qlat   = (q-qs)*abswind*cq_latent*rho_air*0.58*ce*swetclim(:,:,ityr) !< ocean

    end if

    ! change in water vapor
    dq_eva  = -Qlat / cq_latent / (r_qviwv * wz_vapor)  !< evaporation

    !< Call rainfall model
    ! Eq. 8 (a-e) in Stassen & Dommenget 2018
    ! Parameters in unused terms are set to zero
    dq_rain = (c_q + c_rq*rq + c_omega*omega_clim(:,:,ityr) + c_omegastd*omega_std_clim(:,:,ityr))*cq_rain*q
    dq2mmday= wz_vapor * r_qviwv * 86400. !< mm/day

    !< Avoid negative rainfall (dq_rain is negative means positive rainfall!)
    !< Set to minimum ymonmean of GPCP
    where(dq_rain >= -0.0015 / dq2mmday) dq_rain = -0.0015 / dq2mmday

    !< latent heat flux atmos
    Qlat_air = -dq_rain*cq_latent*r_qviwv

  end subroutine hydro

!+++++++++++++++++++++++++++++++++++++++
subroutine climate_change(Tsurf)
!+++++++++++++++++++++++++++++++++++++++
! Routine to do climate change experiments (section 4.3 in Stassen et al 2018)
  use mo_physics
  use mo_numerics

  implicit none
  real, dimension(xdim,ydim,nstep_yr) :: Tsurf
  real, dimension(xdim,ydim,nstep_yr) :: Tsurf_anom, uclim_anom, vclim_anom, & ! Dummies for the anomalies
  &                                      omega_clim_anom, abswind_clim_anom

  integer :: i

  !Print
  write(*,*) '                              '
  write(*,*) '       \  \  \   /\           '
  write(*,*) '        \  \  \ /  \          '
  write(*,*) '         \  \  /    \         '
  write(*,*) '          \  \/      \        '
  write(*,*) '           \  |       |       '
  write(*,*) '            \ |  ~ ~  |       '
  write(*,*) '             \|  ~ ~  |       '
  write(*,*) '                              '
  write(*,*) '            Climate Change    '
  write(*,*) '                              '


  ! Initially set anomalies to zero
  Tsurf_anom=0.; uclim_anom=0.; vclim_anom=0.; omega_clim_anom=0.; abswind_clim_anom=0.

  !< If exp 40: CC -> Prescribe Tsurf, uwind, vwind, omega
  if ( log_exp == 40 ) then
    log_forceTclim = .TRUE.

    ! Read in anomalies
    open(24,file='../input/tsurf.response.cmip5.ensmean', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    open(25,file='../input/zonal.wind.response.cmip5.ensmean', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    open(26,file='../input/meridional.wind.response.cmip5.ensmean', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    open(27,file='../input/omega.response.cmip5.ensmean', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    open(28,file='../input/windspeed.response.cmip5.ensmean', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    do i=1,nstep_yr
      read(24,rec=i) Tsurf_anom(:,:,i)
      read(25,rec=i) uclim_anom(:,:,i)
      read(26,rec=i) vclim_anom(:,:,i)
      read(27,rec=i) omega_clim_anom(:,:,i)
      read(28,rec=i) abswind_clim_anom(:,:,i)
    end do

    ! Add anomalies on top
    Tsurf       = Tsurf + Tsurf_anom
    uclim       = uclim + uclim_anom
    vclim       = vclim + vclim_anom
    omega_clim  = omega_clim + omega_clim_anom
    abswind_clim = abswind_clim + abswind_clim_anom

    ! Close files
    close(24, status='keep'); close(25, status='keep'); close(26, status='keep')
    close(27, status='keep'); close(28, status='keep')

  end if !exp 40

  where (uclim(:,:,:) >= 0.0)
    uclim_m = uclim
    uclim_p = 0.0
  else where
    uclim_m = 0.0
    uclim_p = uclim
  end where
  where (vclim(:,:,:) >= 0.0)
    vclim_m = vclim
    vclim_p = 0.0
  else where
    vclim_m = 0.0
    vclim_p = vclim
  end where

end subroutine climate_change

!+++++++++++++++++++++++++++++++++++++++
subroutine enso(Tsurf)
!+++++++++++++++++++++++++++++++++++++++
! Routine to do the El Nino / La Nina experiments (section 4.2 in Stassen et al 2018)
  use mo_physics
  use mo_numerics

  implicit none

  real, dimension(xdim,ydim,nstep_yr), intent(inout) :: Tsurf

  real, dimension(xdim,ydim,nstep_yr) :: Tsurf_anom, uclim_anom, vclim_anom, & ! Dummies for the anomalies
  &                                      omega_clim_anom, abswind_clim_anom

  integer                       :: i

  !Print
  write(*,*) '                                  '
  write(*,*) '                                  '
  write(*,*) '               ~~~~~~~~~~~~~~~~~  '
  write(*,*) '               ~~~~~~~~~~~~~~~~~  '
  write(*,*) '               ~~~~~~~~~~~~~~~~~  '
  write(*,*) '               ~~~~~~~~~~~~~~~~~  '
  write(*,*) '                                  '
  write(*,*) '               El Nino            '
  write(*,*) '                                  '

  ! Initially set anomalies to zero
  Tsurf_anom=0.; uclim_anom=0.; vclim_anom=0.; omega_clim_anom=0.; abswind_clim_anom=0.

  !< If exp 30: El Nino -> Prescribe T_surf, uwnd, vwnd, omega anomalies
  if ( log_exp == 30 ) THEN
    log_forceTclim = .TRUE.
    ! Read in anomalies
    open(24,file='../input/tsurf.response.erainterim.elnino', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    open(25,file='../input/zonal.wind.response.erainterim.elnino', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    open(26,file='../input/meridional.wind.response.erainterim.elnino', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    open(27,file='../input/omega.response.erainterim.elnino', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    open(28,file='../input/windspeed.response.erainterim.elnino', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    do i=1,nstep_yr
      read(24,rec=i) Tsurf_anom(:,:,i)
      read(25,rec=i) uclim_anom(:,:,i)
      read(26,rec=i) vclim_anom(:,:,i)
      read(27,rec=i) omega_clim_anom(:,:,i)
      read(28,rec=i) abswind_clim_anom(:,:,i)
    end do
  end if !exp 30

  !< If exp 31: El Nino -> Prescribe T_surf anomalies only
  if ( log_exp == 31 ) THEN
    log_forceTclim = .TRUE.
    ! Read in anomalies
    open(24,file='../input/tsurf.response.erainterim.elnino', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    do i=1,nstep_yr
      read(24,rec=i) Tsurf_anom(:,:,i)
    end do
  end if !exp 31

  !< If exp 34: La Nina -> Prescribe T_surf, uwnd, vwnd, omega anomalies
  if ( log_exp == 34 ) THEN
    log_forceTclim = .TRUE.
    ! Read in anomalies
    open(24,file='../input/tsurf.response.erainterim.lanina', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    open(25,file='../input/zonal.wind.response.erainterim.lanina', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    open(26,file='../input/meridional.wind.response.erainterim.lanina', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    open(27,file='../input/omega.response.erainterim.lanina', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    open(28,file='../input/windspeed.response.erainterim.lanina', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    do i=1,nstep_yr
      read(24,rec=i) Tsurf_anom(:,:,i)
      read(25,rec=i) uclim_anom(:,:,i)
      read(26,rec=i) vclim_anom(:,:,i)
      read(27,rec=i) omega_clim_anom(:,:,i)
      read(28,rec=i) abswind_clim_anom(:,:,i)
    end do
  end if !exp 34

  !< If exp 35: La Nina -> Prescribe T_surf anomalies only
  if ( log_exp == 35 ) THEN
    log_forceTclim = .TRUE.
    ! Read in anomalies
    open(24,file='../input/tsurf.response.erainterim.elanina', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    do i=1,nstep_yr
      read(24,rec=i) Tsurf_anom(:,:,i)
    end do
  end if !exp 35

  ! Add anomalies on top (some anomalies might be zero)
  Tsurf       = Tsurf + Tsurf_anom
  uclim       = uclim + uclim_anom
  vclim       = vclim + vclim_anom
  omega_clim  = omega_clim + omega_clim_anom
  abswind_clim = abswind_clim + abswind_clim_anom

  ! Close files
  close(24, status='keep'); close(25, status='keep'); close(26, status='keep')
  close(27, status='keep'); close(28, status='keep')

  where (uclim(:,:,:) >= 0.0)
    uclim_m = uclim
    uclim_p = 0.0
  else where
    uclim_m = 0.0
    uclim_p = uclim
  end where
  where (vclim(:,:,:) >= 0.0)
    vclim_m = vclim
    vclim_p = 0.0
  else where
    vclim_m = 0.0
    vclim_p = vclim
  end where

end subroutine enso

!+++++++++++++++++++++++++++++++++++++++
subroutine SWradiation(Tsurf, sw, albedo)
!+++++++++++++++++++++++++++++++++++++++
!    SW radiation model

  use mo_numerics,    only: xdim, ydim
  use mo_physics,     only: ityr, sw_solar,da_ice, a_no_ice, a_cloud, z_topo,  &
  &                         Tl_ice1, Tl_ice2, To_ice1, To_ice2, glacier,       &
  &                         cldclim, log_exp

  ! declare temporary fields
  real, dimension(xdim,ydim)  :: Tsurf, sw, albedo, a_surf, a_atmos

  ! atmos albedo
  do i=1,xdim
    a_atmos(i,:) = cldclim(i,:,ityr)*a_cloud
  end do

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

  if (log_exp <= 5) a_surf = a_no_ice

  ! SW flux
  albedo=a_surf+a_atmos-a_surf*a_atmos
  forall (i=1:xdim)
    sw(i,:)=SW_solar(:,ityr)*(1-albedo(i,:))
  end forall

end subroutine SWradiation

!+++++++++++++++++++++++++++++++++++++++
subroutine LWradiation(Tsurf, Tair, q, CO2, LWsurf, LWair_up, LWair_down, em)
!+++++++++++++++++++++++++++++++++++++++
! new approach with LW atmos

  use mo_numerics,    only: xdim, ydim
  use mo_physics,     only: sig, eps, qclim, cldclim, z_topo, jday, ityr,         &
  &                           r_qviwv, z_air, z_vapor, dTrad, p_emi, log_exp

  ! declare temporary fields
  real, dimension(xdim,ydim)  :: Tsurf, Tair, q, LWsurf, LWair, e_co2, e_cloud,   &
  &                                LWair_up, LWair_down, e_vapor, em


  e_co2   = exp(-z_topo/z_air)*CO2         ! CO2
  e_vapor = exp(-z_topo/z_air)*r_qviwv*q   ! water vapor
  e_cloud = cldclim(:,:,ityr)              ! clouds

  if(log_exp == 11) e_vapor = exp(-z_topo/z_air)*r_qviwv*qclim(:,:,ityr)     ! sens. exp. linear-function

  ! total
  em      = p_emi(4)*log( p_emi(1)*e_co2 +p_emi(2)*e_vapor +p_emi(3) ) +p_emi(7)   &
  &          +p_emi(5)*log( p_emi(1)*e_co2   +p_emi(3) )                             &
  &          +p_emi(6)*log( p_emi(2)*e_vapor +p_emi(3) )
  em      = (p_emi(8)-e_cloud)/p_emi(9)*(em-p_emi(10))+p_emi(10)
  if(log_exp == 11)  em = em +0.022/(0.15*24.)*r_qviwv*(q-qclim(:,:,ityr)) ! sens. exp. linear-function

  LWsurf      = -sig*Tsurf**4
  LWair_down  = -em*sig*(Tair+dTrad(:,:,ityr))**4
  LWair_up    = LWair_down

end subroutine LWradiation

!+++++++++++++++++++++++++++++++++++++++
subroutine seaice(Tsurf)
!+++++++++++++++++++++++++++++++++++++++
!    Sea ice model

  use mo_numerics,    only: xdim, ydim
  use mo_physics,     only: ityr, z_topo, cap_surf, cap_land, cap_ocean, &
  &                           log_exp, To_ice1, To_ice2, glacier, mldclim

  ! declare temporary fields
  real, dimension(xdim,ydim)  :: Tsurf

  where(z_topo < 0. .and. Tsurf <= To_ice1) cap_surf = cap_land                    ! sea ice
  where(z_topo < 0. .and. Tsurf >= To_ice2) cap_surf = cap_ocean*mldclim(:,:,ityr) ! open ocean
  where(z_topo < 0. .and. Tsurf > To_ice1 .and. Tsurf < To_ice2 ) &
    &       cap_surf = cap_land + (cap_ocean*mldclim(:,:,ityr)-cap_land)     &
    &                            /(To_ice2-To_ice1)*(Tsurf-To_ice1)

    if( log_exp <= 5 ) then
      where(z_topo > 0. ) cap_surf = cap_land                     ! sea ice
      where(z_topo < 0. ) cap_surf = cap_ocean*mldclim(:,:,ityr)  ! open ocean
    end if

    ! glacier -> no sea ice change
    where(glacier > 0.5) cap_surf = cap_land                       ! ice sheet

end subroutine seaice

!+++++++++++++++++++++++++++++++++++++++
subroutine deep_ocean(Ts, To, dT_ocean, dTo)
!+++++++++++++++++++++++++++++++++++++++
! deep ocean model

  use mo_numerics,    only: xdim, ydim, nstep_yr, dt
  use mo_physics,     only: ityr, z_topo, mldclim, log_exp, To_ice2,     &
  &                           cap_ocean, co_turb, z_ocean

  ! declare temporary fields
  real, dimension(xdim,ydim)  :: Ts, To, dT_ocean, dTo, dmld, Tx
  dT_ocean = 0.0;  dTo     = 0.0
  if ( log_exp <= 9 .or. log_exp == 11 )   return
  if ( log_exp >= 14 .and. log_exp <= 16 ) return

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
  dTo      = dTo      + dt*co_turb*(Tx-To)/(cap_ocean*(z_ocean-mldclim(:,:,ityr)))
  dT_ocean = dT_ocean + dt*co_turb*(To-Tx)/(cap_ocean*mldclim(:,:,ityr))

end subroutine deep_ocean

!+++++++++++++++++++++++++++++++++++++++
subroutine circulation(X_in, dX_crcl, dX_diff, dX_advec, h_scl, wz)
!+++++++++++++++++++++++++++++++++++++++
! circulation with shorter time step

  use mo_numerics,  only: xdim, ydim, dt, dt_crcl
  use mo_physics,   only: log_exp, z_vapor, z_air, log_adv_scheme

  implicit none

  real, dimension(xdim,ydim), intent(in)  :: X_in, wz
  real,                       intent(in)  :: h_scl
  real, dimension(xdim,ydim), intent(out) :: dX_crcl

  real, dimension(xdim,ydim) :: X, dx_diffcrcl, dx_adveccrcl, dx_convcrcl
  real, dimension(xdim,ydim) :: dX_diff, dX_advec
  integer time, tt

  !< Init.
  X        = 0.
  dX_diff  = 0.
  dX_advec = 0.
  dX_crcl  = 0.
  dx_diffcrcl  = 0.
  dx_adveccrcl = 0.
  dx_convcrcl  = 0.

  if(log_exp  <=  4 ) return
  if(log_exp .eq.  7 .and. h_scl .eq. z_vapor) return
  if(log_exp .eq. 16 .and. h_scl .eq. z_vapor) return
  if(log_exp .eq. 24) return
  if(log_exp .eq. 25 .and. h_scl .eq. z_air) return

  time=max(1,nint(float(dt)/dt_crcl))

  X = X_in;
  if(log_exp .eq. 8 .and. h_scl .eq. z_vapor) then
    do tt=1, time   ! time loop circulation
      call diffusion(X, dx_diffcrcl, h_scl, wz)
      X = X + dx_diffcrcl
    end do           ! time loop
  else if (h_scl .eq. z_vapor .and. log_adv_scheme .eq. 0) then !For water vapour and convergence
    do tt=1, time   ! time loop circulation
      call diffusion(X, dx_diffcrcl, h_scl, wz)
      call advection(X, dx_adveccrcl, h_scl, wz)
      call convergence(X, dx_convcrcl)
      X = X + dx_diffcrcl + dx_adveccrcl + dx_convcrcl
      dX_diff  = dX_diff + dx_diffcrcl
      dX_advec = dX_advec + dx_adveccrcl + dx_convcrcl
    end do           ! time loop
  else
    do tt=1, time   ! time loop circulation
      call diffusion(X, dx_diffcrcl, h_scl, wz)
      call advection(X, dx_adveccrcl, h_scl, wz)
      X = X + dx_diffcrcl + dx_adveccrcl
      dX_diff  = dX_diff + dx_diffcrcl
      dX_advec = dX_advec + dx_adveccrcl
    end do           ! time loop
  end if
  !dX_crcl = X - X_in
  dX_crcl = dX_diff + dX_advec

end subroutine circulation

!+++++++++++++++++++++++++++++++++++++++
subroutine diffusion(T1, dX_diffuse,h_scl, wz)
!+++++++++++++++++++++++++++++++++++++++
!    diffusion

  use mo_numerics,   only: xdim, ydim, dt, dlon, dlat, dt_crcl
  use mo_physics,    only: pi, z_topo, log_exp, kappa, z_vapor
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
subroutine convergence(T1, div)
!+++++++++++++++++++++++++++++++++++++++
! Calculates divergence (convergence) of a given field (i.e. spec. hum.) when omega is known
! Eq. 14 in Stassen & Dommenget 2018
  use mo_numerics, only: xdim, ydim, nstep_yr, dlon, dlat, dt_crcl
  use mo_physics,  only: ityr, rho_air, grav, pi, z_vapor, omega_clim
  implicit none

  real, dimension(xdim,ydim), intent(in)  :: T1
  real, dimension(xdim,ydim), intent(out) :: div

  real    :: w
  integer :: i, j

  do j=1,ydim
    do i=1,xdim
      !< Vertical velocity omega (Pa/s) to m/s
      w = -omega_clim(i,j,ityr) / (rho_air*grav)
      !< Convergence
      div(i,j) = T1(i,j) * w * dt_crcl / z_vapor * 2.5
    end do
  end do

end subroutine convergence

!+++++++++++++++++++++++++++++++++++++++
subroutine advection(T1, dX_advec,h_scl, wz)
!+++++++++++++++++++++++++++++++++++++++
!    advection after DD

  use mo_numerics, only: xdim, ydim, dt, dlon, dlat, dt_crcl
  use mo_physics,  only: pi, z_topo, uclim, vclim, ityr, z_vapor, log_exp
  use mo_physics,  only: uclim_m, uclim_p, vclim_m, vclim_p
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
    if ( abs(lat(k)) <= 25) then  ! unitl 25degree
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
      dd=max(1,nint(dt_crcl/(dxlat(k)/10.0/1.)))
      !dd=int(100*abs(sin(2.*pi/360.*lat(k))))+2
      !print*, lat(k), dd
      dtdff2=dt_crcl/dd
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

  dX_advec = dTx + dTy

end subroutine advection

!+++++++++++++++++++++++++++++++++++++++
subroutine co2_level(it, year, CO2)
!+++++++++++++++++++++++++++++++++++++++

  use mo_numerics,    only: ndays_yr, ndt_days

  integer :: it, year
  real    :: CO2

  CO2_1950=310.;  CO2_2000=370.;  CO2_2050=520.
  if (year <= 2000.) CO2=CO2_1950 + 60./50.*(year-1950)
  if (year  > 2000. .and. year <= 2050.) CO2=CO2_2000 + 150./50.*(year-2000)
  if (year  > 2050. .and. year <= 2100.) CO2=CO2_2050 + 180./50.*(year-2050)

end subroutine co2_level

!+++++++++++++++++++++++++++++++++++++++
subroutine diagonstics(it, year, CO2, ts0, ta0, to0, q0, albedo, sw, lw_surf, q_lat, q_sens)
!+++++++++++++++++++++++++++++++++++++++
!    diagonstics plots

  use mo_numerics
  use mo_physics
  use mo_diagnostics

  ! declare temporary fields
  real, dimension(xdim,ydim)  :: Ts0, Ta0, To0, q0, sw, albedo, Q_sens, Q_lat,  LW_surf
  integer                     :: year
  real                        :: CO2

  ! diagnostics: annual means
  tsmn=tsmn+Ts0; tamn=tamn+ta0; tomn=tomn+to0; qmn=qmn+q0; amn=amn+albedo
  swmn=swmn+sw;  lwmn=lwmn+LW_surf; qlatmn=qlatmn+q_lat; qsensmn=qsensmn+Q_sens;
  ftmn=ftmn+TF_correct(:,:,ityr); fqmn=fqmn+qF_correct(:,:,ityr);
  if ( ityr == nstep_yr ) then
    tsmn    = tsmn/nstep_yr;      tamn = tamn/nstep_yr;    tomn = tomn/nstep_yr;
    qmn     = qmn/nstep_yr;
    amn     = amn/nstep_yr;       swmn = swmn/nstep_yr;    lwmn = lwmn/nstep_yr;
    qlatmn  = qlatmn/nstep_yr; qsensmn = qsensmn/nstep_yr; ftmn = ftmn/nstep_yr;
    fqmn    = fqmn/nstep_yr;
    1000 format (I4, T8, F10.1, T20, F10.2, T37, F10.6, T52, F10.6, T67, F10.6, T85, F10.6, T110, F10.6) !TB
    print 1000, year, CO2, gmean(swmn), gmean(tsmn)-273.15, tsmn(48,24)-273.15,tsmn(4,39)-273.15,  & !TB
    &     tsmn(1,48)-273.15, sum(qmn)/(96*48)
    tsmn=0.; tamn=0.; qmn=0.; amn=0.; swmn=0.;        ! reset annual mean values
    lwmn=0.; qlatmn=0.; qsensmn=0.; ftmn=0.; fqmn=0.; ! reset annual mean values
  end if

end subroutine diagonstics

!+++++++++++++++++++++++++++++++++++++++
subroutine output( it, iunit, fout, irec, mon, ts0, ta0, to0, q0, albedo, dq_eva, dq_rain, dq_crcl)
!+++++++++++++++++++++++++++++++++++++++
!    write output

  use mo_numerics
  use mo_physics
  use mo_diagnostics

  character(len=*)            :: fout

  ! declare temporary fields
  real, dimension(xdim,ydim)  :: Ts0, Ta0, To0, q0, albedo, dq_eva, dq_rain, dq_crcl

  ! diagnostics: monthly means
  Tmm=Tmm+Ts0; Tamm=Tamm+ta0; Tomm=Tomm+to0; qmm=qmm+q0; apmm=apmm+albedo
  Prmm=Prmm+dq_rain*wz_vapor * r_qviwv * 86400. !mm/day
  dqrainmm=dqrainmm+dq_rain; dqevamm=dqevamm+dq_eva; dqcrclmm=dqcrclmm+dq_crcl/dt
  dq_advmm=dq_advmm+dq_adv/dt; dq_diffmm=dq_diffmm+dq_diff/dt

  if (       jday == sum(jday_mon(1:mon))                   &
  &      .and. it/float(ndt_days) == nint(it/float(ndt_days)) ) then
    ndm=jday_mon(mon)*ndt_days

    irec=irec+1; write(iunit,rec=irec) Tmm/ndm
    irec=irec+1; write(iunit,rec=irec) Tamm/ndm
    irec=irec+1; write(iunit,rec=irec) Tomm/ndm
    irec=irec+1; write(iunit,rec=irec) qmm/ndm
    irec=irec+1; write(iunit,rec=irec) apmm/ndm

    !< Hydro variables
    irec=irec+1; write(iunit,rec=irec) -1.*Prmm/ndm
    irec=irec+1; write(iunit,rec=irec) dqrainmm/ndm
    irec=irec+1; write(iunit,rec=irec) dqevamm/ndm
    irec=irec+1; write(iunit,rec=irec) dqcrclmm/ndm

    !< Set the Monthly mean values back to zero
    Tmm=0.;Tamm=0.;Tomm=0.;qmm=0.;apmm=0.;Prmm=0.;dqrainmm=0.;dqevamm=0.;dqcrclmm=0.;dq_advmm=0.;dq_diffmm=0.
    mon=mon+1; if (mon==13) mon=1
  end if

end subroutine output

!+++++++++++++++++++++++++++++++++++++++
subroutine output_fluxcorr( it, iunit, fout, irec, ts_corr, to_corr, q_corr )
!+++++++++++++++++++++++++++++++++++++++
!    write output
  use mo_numerics,        only: xdim, ydim
  use mo_physics,         only: r_qviwv, wz_vapor

  ! declare temporary fields
  real, dimension(xdim,ydim)  :: Ts_corr, To_corr, q_corr
  character(len=*)            :: fout
  integer, intent(in)         :: iunit

  irec=irec+1; write(iunit,rec=irec) ts_corr
  irec=irec+1; write(iunit,rec=irec) to_corr
  irec=irec+1; write(iunit,rec=irec) q_corr

end subroutine output_fluxcorr

!TB
!+++++++++++++++++++++++++++++++++++++++
function gmean(data)
!+++++++++++++++++++++++++++++++++++++++

  use mo_numerics,		only: xdim, ydim, dlat

  ! declare variables
  real, dimension(xdim,ydim) 	:: data, w
  real, dimension(ydim)		:: lat

  do i=1,ydim
    lat(i) = -90+(dlat*0.5)+(i-1)*dlat
  end do
  do i=1,xdim
    w(i,:) = cos(2.*3.14159*lat/360.)
  end do

  gmean = sum(data*w)/sum(w)

end function

!+++++++++++++++++++++++++++++++++++++++
function grms(data)
!+++++++++++++++++++++++++++++++++++++++

  use mo_numerics,		only: xdim, ydim, dlat

  ! declare variables
  real, dimension(xdim,ydim) 	:: data, w
  real, dimension(ydim)		:: lat

  do i=1,ydim
    lat(i) = -90+(dlat*0.5)+(i-1)*dlat
  end do
  do i=1,xdim
    w(i,:) = cos(2.*3.14159*lat/360.)
  end do

  grms = sqrt(sum(data**2.*w)/sum(w))

end function
