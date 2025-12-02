! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Input control for radiation.

! Description:
!   Module containing input switches/settings as used by the radiation code.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: radiation_control

! Method:
!   Switches are initialised to false and read in from the
!   namelist. The module may then be used directly where the switches
!   are needed within the radiation code.

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3

module rad_input_mod

use missing_data_mod, only: imdi, rmdi
use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim
use errormessagelength_mod, only: errormessagelength

use um_types, only: real_umphys

implicit none

! ----------------
! Control options
! ----------------

logical :: l_radiation = .false.    !  F: Turns off radiation code

logical :: l_use_dust = .false.     !  Use mineral dust in rad calculations
logical :: l_use_biogenic = .false. !  Use biogenic aerosol in radiation code

! Use SO4 aerosol from sulphur cycle for direct/indirect effect
! in radiation, the latter for both SW and LW.
logical :: l_use_sulpc_direct = .false.
logical :: l_use_sulpc_indirect_sw = .false.
logical :: l_use_sulpc_indirect_lw = .false.

! Indirect radiative effect of sea-salt
logical :: l_use_seasalt_indirect = .false.

! Direct radiative effect of sea-salt
logical :: l_use_seasalt_direct = .false.

logical :: l_use_soot_direct = .false.   ! direct radiative effects of soot
logical :: l_use_soot_indirect = .false. ! indirect effects of soot

! Use biomass aerosol for direct/indirect effect in radiation.
logical :: l_use_bmass_direct = .false.
logical :: l_use_bmass_indirect = .false.

! Use fossil-fuel organic carbon aerosol for direct/indirect
! effect in radiation
logical :: l_use_ocff_direct = .false.
logical :: l_use_ocff_indirect = .false.

! Use ammonium nitrate aerosol for direct/indirect effect in radiation
logical :: l_use_nitrate_direct = .false.
logical :: l_use_nitrate_indirect = .false.

! Use the same number of cloud droplets for 1st and 2nd indirect effects
logical :: l_consistent_cdnc = .false.

! Use aerosol climatologies in radiation instead of prognostic variables
! Set on a species by species basis
logical :: l_use_arclbiom = .false. ! biomass burning aerosol
logical :: l_use_arclblck = .false. ! black carbon
logical :: l_use_arclsslt = .false. ! sea salt
logical :: l_use_arclsulp = .false. ! sulpahtes
logical :: l_use_arcldust = .false. ! mineral dust
logical :: l_use_arclocff = .false. ! organic carbon (fossil fuel)
logical :: l_use_arcldlta = .false. ! delta aerosol

! Use droplet number from n_drop_pot array:
logical :: l_use_ndrop = .false.

! Use Liu 2008 spectral broadening in calculation of effective radius
logical :: l_use_liu_spec = .false.
! Parameters for tuning the beta function in the Liu (2008) scheme
real(kind=real_umphys) :: aparam = rmdi
real(kind=real_umphys) :: bparam = rmdi

logical :: L_rad_deg = .false.  ! controls the use of spatial degradation
!                                 of radiation calc.

! Logicals for different radiation packages
logical :: l_forcing     = .false. ! Calculate radiative forcings
logical :: l_radiance    = .false. ! Calculate radiances
logical :: l_timestep    = .false. ! Use new timestepping scheme
logical :: l_rad_perturb = .false. ! Use the perturbation version of
                                   ! the radiative time-stepping
! Control integer for different radiation packages used by check_run_radiation()
! routine to set l_forcing, l_radiance, l_timestep and l_rad_perturb
integer :: i_rad_extra_call = imdi

! Use a solar zenith angle correction based on optical depth
logical :: l_rad_szacor = .false.
!                       Needed in glue_rad-rad_ctl3c.
!                       Switch for the solar zenith angle correction to surface
!                       fluxes using the change in optical depth.

! Scale the condensed water content to simulate
! inhomogeneous clouds
logical :: l_inhom_cloud = .false.

! Orography correction to SW radiation
logical :: l_use_orog_corr  = .false.  !  Find gradients from mean orog
logical :: l_use_grad_corr  = .false.  !  Use ancillary X & Y gradients
! Correction for skyview factor in LW and direct SW
logical :: l_use_skyview    = .false.
logical :: l_orog_unfilt    = .false.  !  Use unfiltered ancillary orog
! Control integer used by check_radiaiton to set the orography correction
! and skyview logicals
integer :: i_rad_topography = imdi

! ----------------------------

integer :: h_swbands    = imdi   ! Number of shortwave radiation bands
integer :: h_lwbands    = imdi   ! Number of longwave radiation bands

! Number of advection steps per prognostic/diagnostic SW and LW step.
!'Prognostic' and 'Diagnostic' refer to the frequency of the calls
! to radiation code in the Unified Model.
! In the case of time stepping prognostic and diagnostic refer to the
! slow and fast radiative timestep respectively. In the case of radiative
! forcing they refer to the prognostic and diagnostic calls to radiation.

integer :: a_sw_radstep_diag = imdi
! Number of advection steps per 'fast' SW step (3C)
integer :: a_lw_radstep_diag = imdi
! Number of advection steps per 'fast' LW step (3C)
integer :: a_sw_radstep_prog = imdi
! Number of advection steps per 'slow' LW step (3C)
integer :: a_lw_radstep_prog = imdi
! Number of advection steps per 'slow' LW step (3C)

integer :: i_ozone_int = imdi ! Option for interpolation of ozone


! The following three switches REMOVED from run_radiation NL (ROSE project)
! They are set in check_run_radiation, dependent on settings of cusack_aero and
! cusack_aero_hgt, which have been added to the NL

! True if climatological aerosol is included.
logical :: L_climat_aerosol = .false.

! True to use real boundary layer heights to specify the boundary
! layer aerosol.
logical :: L_clim_aero_hgt = .false.

! Flag to use HadGEM1 setting for climatological aerosols
logical :: L_HadGEM1_Clim_Aero = .false.

! These two switches ADDED to NL (ROSE project)
integer :: cusack_aero     = imdi
integer :: cusack_aero_hgt = imdi

logical :: lrad_ccrad = .false.
!             Allows access to ccrad code and the logicals
!             lrad_ovrlap and lrad_ccw_scav

! Convert zonal mean ozone to field
logical :: lexpand_ozone

! convert zonal mean tpps ozone to field
logical :: lexpand_tpps_ozone

! Tropopause-based Ozone Scheme
logical :: l_use_tpps_ozone = .false.  !  Use TPPS ozone scheme

!-----------------------------------------------------------
! run_radiation namelists
! ----------------------------------------------------------

logical :: l_sec_var  = .false.     ! true if using time varying astronomy

logical :: l_rad_ovrlap   = .false.
!                           Requires l_ccrad=.true.
!                           Allows Convective and LS Cloud to overlap
!                           for radiative impacts.
!                           (Experimental, defaulted to false,
!                           requires a hand-edit to change)
!             (THIS is EXPERIMENTAL and USED FOR DEVELOPMENT only).
!             Current convective/large-scale cloud fractions in the
!             radiation scheme are mutally exclusive. This assumes CCA
!             and the large-scale to overlap with CCA taking dominance.
!             I.E. Large-scale cloud fraction must exceed the convective
!             cloud fraction before having any presence.

logical :: l_rad_ccw_scav = .false.
!                           Requires l_ccrad=.true. .and. l_rad_ovrlap=.true.
!                           Allows Convective Cloud Water (CCW) to
!                           compensate for LS Cloud water in overlapping
!                           LS/CCA fractions.
!                           (Experimental, defaulted to false, requires a
!                           hand-edit to change)
!             (THIS is EXPERIMENTAL and USED FOR DEVELOPMENT only)
!             Allowing the CCA to negate large-scale fractions of lower
!             values means that the large-scale cloud water in the
!             overlapping fraction is lost. This switch will scavenge
!             the large-scale cloud water from the overlapping fraction
!             and combine it with the convective cloud water to
!             conpensate.

logical :: l_rad_use_clim_volc=.false.
!                           If .true. use climatological volcanic
!                           eruption code in climatological aerosol
!                           code


logical :: l_rad_snow_emis = .false.
!          Switch to adjust the emissivity in radiation for snow cover.
!          This should eventually be moved to the surface scheme, but
!          JULES cannot currently cope with distinct emissivities
!          for snow-covered surfaces, so the switch currently acts
!          only in radiation and logically belongs here for the present.
logical :: l_t_land_nosnow = .false.
!          Switch for emissivity of snow used in averaging the
!          surface temperature. Setting this switch to .true. is
!          deprecated and it is included only for historical reasons.
logical :: l_quad_t_coast = .false.
!          Switch for quadratic averaging of the surface temperature at
!          coastal points. .false. is deprecated.
logical :: l_t_rad_solid = .false.
!          Switch to use common soid temperature at coastal points with
!          sea-ice. .true. is deprecated.

logical :: l_t_bdy_surf = .false.
!          Take the temperature of the air just above the surface as
!          the temperature at the surface.

! ------------------------------------------

! number of components of clouds
integer,parameter:: npd_cloud_component=4

integer :: aero_bl_levels = imdi
!                          Common number of layers taken to be
!                          occupied by the boundary-layer
!                          aerosol if the boundary layer
!                          depth is not used to determine the
!                          number separately at each grid-point
!                          In previous versions of the code,
!                          this was taken to be BL_LEVELS

integer :: clim_rad_volc_eruption_year = imdi  ! Climatological volcano
!                                                eruption year

integer :: clim_rad_volc_eruption_month = imdi    ! Climatological volcano
!                                                eruption month

integer :: rad_mcica_sampling = imdi   ! Version of McICA used (was 1)
!             Needed in open_cloud_gen. Selects the version of McICA
!             used to sample the generated cloud:
!                               0 = full sampling
!                               1 = single sampling
!                               2 = optimal sampling

! --------------------------------

real(kind=real_umphys)    :: rad_mcica_sigma = rmdi
                                   ! Normalised cloud condensate standard
!                                    deviation for the cloud generator.
!                                    Needed in open_cloud_gen.

real(kind=real_umphys)    :: two_d_fsd_factor = rmdi
                                   ! Ratio between fsd in 2D (required by UM)
                                   ! and 1D (as parametrized).

logical :: l_fsd_eff_res = .false. ! Use an effective resolution rather than
                                   ! the actual grid-length in the FSD
                                   ! parametrization

integer :: i_cloud_representation = imdi
integer :: i_inhom = imdi
integer :: i_cloud_entrapment = imdi
integer :: i_overlap = imdi
integer :: i_fsd = imdi
integer :: i_cloud_representation_2 = imdi
integer :: i_inhom_2 = imdi
integer :: i_cloud_entrapment_2 = imdi
integer :: i_overlap_2 = imdi

! Mass Mixing Ratios (MMR) of minor Gases
real(kind=real_umphys) :: co2_mmr    = rmdi ! CO2 concentration (if constant)
real(kind=real_umphys) :: n2ommr     = rmdi ! N2O mmr
real(kind=real_umphys) :: ch4mmr     = rmdi ! CH4 mmr
real(kind=real_umphys) :: c11mmr     = rmdi ! CFC11 mmr
real(kind=real_umphys) :: c12mmr     = rmdi ! CFC12 mmr
real(kind=real_umphys) :: o2mmr      = rmdi ! O2 mmr
real(kind=real_umphys) :: so2mmr     = rmdi ! SO2 mmr
real(kind=real_umphys) :: c113mmr    = rmdi ! CFC113 mmr
real(kind=real_umphys) :: c114mmr    = rmdi ! CFC114 mmr
real(kind=real_umphys) :: hcfc22mmr  = rmdi ! HCFC22 mmr
real(kind=real_umphys) :: hfc125mmr  = rmdi ! HFC125 mmr
real(kind=real_umphys) :: hfc134ammr = rmdi ! HFC134A mmr

! Scaling factors to simulate inhomogeneous cloud.
real(kind=real_umphys)    :: inhom_cloud_sw(npd_cloud_component) = rmdi
real(kind=real_umphys)    :: inhom_cloud_lw(npd_cloud_component) = rmdi


! Decorrelation pressure scale for large scale cloud
real(kind=real_umphys)    :: dp_corr_strat = rmdi

! Decorrelation pressure scale for convective cloud
real(kind=real_umphys)    :: dp_corr_conv  = rmdi


real(kind=real_umphys)    :: clim_rad_volc_eruption_weight = rmdi
! Eruption weighting factor for idealised volcanic aerosol.
! 1.0 is an average 20th century tropical explosive eruption.

real(kind=real_umphys)                                                         &
        :: aeroscl_csk_clim(5) = [ rmdi, rmdi, rmdi, rmdi, rmdi ]
! Scalings for aerosols in Cusack's climatology

! Number of radiation prognostic/diagnostic timesteps per day
integer :: i_sw_radstep_perday_prog = imdi
integer :: i_lw_radstep_perday_prog = imdi
integer :: i_sw_radstep_perday_diag = imdi
integer :: i_lw_radstep_perday_diag = imdi

! Ozone tracer as input to radiation scheme
logical :: l_use_cariolle   = .false.
logical :: l_use_ozoneinrad = .false.

! Use abundances from Burrows & Sharp, ApJ, 1999 (for hot Jupiters)
logical :: l_BS1999_abundances = .false.

! Calculate layer masses using the hydrostatic approximation
logical :: l_hydrostatic_mass = .false.

! Calculate layer heat capacities including the moisture
logical :: l_moist_heat_capacity = .false.

! Create an extra top layer for radiation
logical :: l_extra_top = .false.

! Apply the NLTE correction to heating rates
logical :: l_nlte_corr = .false.

namelist/RUN_Radiation/                                                        &
       cusack_aero, cusack_aero_hgt, aeroscl_csk_clim, co2_mmr,                &
       l_sec_var, inhom_cloud_sw, inhom_cloud_lw, dp_corr_strat,               &
       rad_mcica_sampling, rad_mcica_sigma, two_d_fsd_factor,                  &
       l_fsd_eff_res, dp_corr_conv, aero_bl_levels,                            &
       l_rad_use_clim_volc, clim_rad_volc_eruption_year,                       &
       clim_rad_volc_eruption_month, clim_rad_volc_eruption_weight,            &
       n2ommr, ch4mmr, so2mmr, c11mmr, c12mmr,                                 &
       o2mmr, c113mmr, c114mmr, hcfc22mmr, hfc125mmr, hfc134ammr,              &
       i_cloud_representation, i_cloud_representation_2,                       &
       i_inhom, i_inhom_2, i_overlap, i_overlap_2, i_fsd,                      &
       i_cloud_entrapment, i_cloud_entrapment_2,                               &
       l_rad_snow_emis, l_t_land_nosnow, l_quad_t_coast, l_t_rad_solid,        &
       l_t_bdy_surf,                                                           &
       h_swbands, h_lwbands, i_rad_extra_call, i_rad_topography,               &
       l_radiation, l_rad_deg,                                                 &
       l_rad_szacor, i_sw_radstep_perday_prog,                                 &
       i_lw_radstep_perday_prog,  i_sw_radstep_perday_diag,                    &
       i_lw_radstep_perday_diag, l_use_cariolle, l_use_ozoneinrad,             &
       i_ozone_int, l_consistent_cdnc, l_use_sulpc_direct,                     &
       l_use_sulpc_indirect_sw, l_use_sulpc_indirect_lw,                       &
       l_use_seasalt_direct, l_use_seasalt_indirect, l_use_biogenic,           &
       l_use_nitrate_direct, l_use_nitrate_indirect, l_use_soot_direct,        &
       l_use_soot_indirect, l_use_bmass_direct, l_use_bmass_indirect,          &
       l_use_ocff_direct, l_use_ocff_indirect, l_use_dust,                     &
       l_use_arclbiom, l_use_arclblck, l_use_arcldlta, l_use_arcldust,         &
       l_use_arclocff, l_use_arclsslt, l_use_arclsulp,                         &
       l_BS1999_abundances, l_hydrostatic_mass, l_moist_heat_capacity,         &
       l_extra_top, l_use_liu_spec, aparam, bparam, l_nlte_corr

! Logical variables to control whether different cca type progonostics
! are required by the cloud generator
logical :: l_cca_dp_prog = .false.
logical :: l_cca_md_prog = .false.
logical :: l_cca_sh_prog = .false.

! The number of calls to SW/LW radiation
integer :: n_swcall = imdi
integer :: n_lwcall = imdi

! Parameters for the values that i_rad_extra_call can have:
integer, parameter :: ip_single_call           = 0
integer, parameter :: ip_diagnostic_call       = 1
integer, parameter :: ip_increment_call        = 2
integer, parameter :: ip_radiance_call         = 3

! Timestep number of the first step with a radiation call; assumed to be 1
! in the full UM, but may be set otherwise in the SCM:
integer :: it_rad1 = 1

!DrHook-related parameters
integer(kind=jpim), parameter, private :: zhook_in  = 0
integer(kind=jpim), parameter, private :: zhook_out = 1

character(len=*), parameter, private :: ModuleName='RAD_INPUT_MOD'

contains

subroutine check_run_radiation()

! Description:
!   Subroutine to apply logic controls and set control variables based on the
!   options selected in the run_radiation namelist.

use ereport_mod,  only: ereport
use fsd_parameters_mod, only: ip_fsd_constant, ip_fsd_param,                   &
                              ip_fsd_regime, ip_fsd_regime_no_sh,              &
                              ip_fsd_regime_smooth, ip_fsd_regime_smooth_no_sh,&
                              ip_fsd_boutle, ip_fsd_regime_cca,                &
                              ip_fsd_regime_cca_v8
use cv_run_mod, only: i_convection_vn, i_cv_comorph
use max_calls, only: npd_swcall, npd_lwcall
use chk_opts_mod, only: chk_var, def_src
use umprintmgr, only: newline

implicit none

integer                       :: icode         ! used for ereport
character (len=errormessagelength)   :: cmessage      ! used for ereport
character (len=*), parameter  :: RoutineName = 'CHECK_RUN_RADIATION'

real(kind=jprb)               :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! required by chk_var routine
def_src = ModuleName//':'//RoutineName

! Set logicals for different radiation packages based on a control
! integer choice in namelist
select case (i_rad_extra_call)
case (ip_single_call)
  l_forcing     = .false.
  l_timestep    = .false.
  l_rad_perturb = .false.
  l_radiance    = .false.
  n_lwcall      = 1
  n_swcall      = 1
case (ip_diagnostic_call)
  l_forcing     = .true.
  l_timestep    = .false.
  l_rad_perturb = .false.
  l_radiance    = .false.
  n_lwcall      = npd_lwcall
  n_swcall      = npd_swcall
case (ip_increment_call)
  l_forcing     = .false.
  l_timestep    = .true.
  l_rad_perturb = .true.
  l_radiance    = .false.
  n_lwcall      = 2
  n_swcall      = 2
case (ip_radiance_call)
  l_forcing     = .false.
  l_timestep    = .false.
  l_rad_perturb = .false.
  l_radiance    = .true.
  n_lwcall      = npd_lwcall
  n_swcall      = npd_swcall
case DEFAULT
  write (cmessage,'(A,I1,A)')                                                  &
          'i_rad_extra_call value invalid, default to ',                       &
          ip_single_call, ': single call to radiation'
  icode = -100
  call ereport(RoutineName, icode, cmessage)
end select

! Set logicals for different radiation topography based on a control
! integer choice in namelist
if (i_rad_topography == 0) then
  l_use_orog_corr   = .false.
  l_use_grad_corr   = .false.
  l_use_skyview     = .false.
  l_orog_unfilt     = .false.
else if (i_rad_topography == 1) then
  l_use_orog_corr   = .true.
  l_use_grad_corr   = .false.
  l_use_skyview     = .false.
  l_orog_unfilt     = .false.
else if (i_rad_topography == 2) then
  l_use_orog_corr   = .false.
  l_use_grad_corr   = .true.
  l_use_skyview     = .false.
  l_orog_unfilt     = .false.
else if (i_rad_topography == 3) then
  l_use_orog_corr   = .true.
  l_use_grad_corr   = .false.
  l_use_skyview     = .true.
  l_orog_unfilt     = .false.
else if (i_rad_topography == 4) then
  l_use_orog_corr   = .false.
  l_use_grad_corr   = .true.
  l_use_skyview     = .true.
  l_orog_unfilt     = .true.
else
  write (cmessage,'(A58)') 'i_rad_topography value invalid, default to '       &
                            // '0: flat surface'
  icode = -100
  call ereport(RoutineName, icode, cmessage)
end if

! Warn if duplicate effects selected

if (l_use_arclbiom .and. l_use_bmass_direct) then
  write (cmessage,'(A)') 'arclbiom and bmass_direct should not both be true'
  icode = -110
  call ereport(RoutineName, icode, cmessage)
end if

if (l_use_arclblck .and. l_use_soot_direct) then
  write (cmessage,'(A)') 'arclblck and soot_direct should not both be true'
  icode = -110
  call ereport(RoutineName, icode, cmessage)
end if

if (l_use_arclocff .and. l_use_ocff_direct) then
  write (cmessage,'(A)') 'arclocff and ocff_direct should not both be true'
  icode = -110
  call ereport(RoutineName, icode, cmessage)
end if

if (l_use_arclsslt .and. l_use_seasalt_direct) then
  write (cmessage,'(A)') 'arclsslt and seasalt_direct should not both be true'
  icode = -110
  call ereport(RoutineName, icode, cmessage)
end if

if (l_use_arclsulp .and. l_use_sulpc_direct) then
  write (cmessage,'(A)') 'arclsulp and sulpc_direct should not both be true'
  icode = -110
  call ereport(RoutineName, icode, cmessage)
end if

! Check i_fsd has allowed value
call chk_var( i_fsd, 'i_fsd',                                                  &
              [ ip_fsd_constant, ip_fsd_param,                                 &
                ip_fsd_regime, ip_fsd_regime_no_sh,                            &
                ip_fsd_regime_smooth, ip_fsd_regime_smooth_no_sh,              &
                ip_fsd_boutle, ip_fsd_regime_cca, ip_fsd_regime_cca_v8] )

! Don't allow running comorph convection scheme with i_fsd set to
! use separate deep, shallow, mid scheme outputs
if ( i_convection_vn==i_cv_comorph .and.                                       &
     ( i_fsd==ip_fsd_regime .or. i_fsd==ip_fsd_regime_no_sh .or.               &
       i_fsd==ip_fsd_regime_smooth .or. i_fsd==ip_fsd_regime_smooth_no_sh )    &
   ) then
  write(cmessage,'(A)')                                                        &
    'i_fsd is set to use separate convective cloud arrays from'     //newline//&
    'deep, shallow and mid convection schemes, but i_convection_vn' //newline//&
    'has been set to use the comorph convection scheme, which does' //newline//&
    'not have separate deep, shallow and mid-level schemes.'
  icode=100
  call ereport(RoutineName, icode, cmessage)
end if

! Set logicals for whether cca from different convection types are
! required by the cloud generator
if ((i_fsd == ip_fsd_regime) .or. (i_fsd == ip_fsd_regime_smooth)) then
  l_cca_sh_prog = .true.
  l_cca_dp_prog = .true.
  l_cca_md_prog = .true.
else if ((i_fsd == ip_fsd_regime_no_sh)  .or.                                  &
         (i_fsd == ip_fsd_regime_smooth_no_sh)) then
  l_cca_dp_prog = .true.
  l_cca_md_prog = .true.
end if

! Should also check validity of aerosol radiative effect choices,
! but this needs to be done after run_aerosol has been read.

! Set radiation aerosol switches
if (cusack_aero==2 .or.  cusack_aero==3    ) l_climat_aerosol    = .true.
if (cusack_aero==3 .and. cusack_aero_hgt==1) l_clim_aero_hgt     = .true.
if (cusack_aero==2)                          l_HadGEM1_clim_aero = .true.

! Reset
def_src = ''

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine check_run_radiation


subroutine print_nlist_run_radiation()
use umPrintMgr, only: umPrint
implicit none
character(len=50000) :: lineBuffer
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='PRINT_NLIST_RUN_RADIATION'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

call umPrint('Contents of namelist run_radiation',                             &
    src='rad_input_mod')

write(lineBuffer,*)' cusack_aero = ',cusack_aero
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' cusack_aero_hgt = ',cusack_aero_hgt
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' co2_mmr = ',co2_mmr
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' l_sec_var = ',l_sec_var
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' inhom_cloud_sw = ',inhom_cloud_sw
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' inhom_cloud_lw = ',inhom_cloud_lw
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' dp_corr_strat = ',dp_corr_strat
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' rad_mcica_sampling = ',rad_mcica_sampling
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' rad_mcica_sigma = ',rad_mcica_sigma
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' two_d_fsd_factor = ',two_d_fsd_factor
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' l_fsd_eff_res = ',l_fsd_eff_res
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' dp_corr_conv = ',dp_corr_conv
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' aero_bl_levels = ',aero_bl_levels
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' l_rad_use_clim_volc = ',l_rad_use_clim_volc
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' clim_rad_volc_eruption_year = ',                          &
    clim_rad_volc_eruption_year
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' clim_rad_volc_eruption_month = ',                         &
    clim_rad_volc_eruption_month
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' clim_rad_volc_eruption_weight = ',                        &
    clim_rad_volc_eruption_weight
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' n2ommr = ',n2ommr
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' ch4mmr = ',ch4mmr
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' c11mmr = ',c11mmr
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' c12mmr = ',c12mmr
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' o2mmr = ',o2mmr
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' so2mmr = ',so2mmr
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' c113mmr = ',c113mmr
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' c114mmr = ',c114mmr
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' hcfc22mmr = ',hcfc22mmr
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' hfc125mmr = ',hfc125mmr
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' hfc134ammr = ',hfc134ammr
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' l_use_sulpc_direct = ',l_use_sulpc_direct
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' l_use_sulpc_indirect_sw = ',l_use_sulpc_indirect_sw
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' l_use_sulpc_indirect_lw = ',l_use_sulpc_indirect_lw
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' l_use_seasalt_direct = ',l_use_seasalt_direct
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' l_use_seasalt_indirect = ',l_use_seasalt_indirect
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' l_use_biogenic = ',l_use_biogenic
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*) ' l_use_ocff_direct = ',l_use_ocff_direct
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*) ' l_use_ocff_indirect = ',l_use_ocff_indirect
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*) ' l_use_nitrate_direct = ',l_use_nitrate_direct
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*) ' l_use_nitrate_indirect = ',l_use_nitrate_indirect
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*) ' l_consistent_cdnc = ',l_consistent_cdnc
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*) ' l_use_soot_direct = ',l_use_soot_direct
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*) ' l_use_soot_indirect = ',l_use_soot_indirect
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*) ' l_use_bmass_direct = ',l_use_bmass_direct
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*) ' l_use_bmass_indirect = ',l_use_bmass_indirect
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*) ' l_use_dust = ',l_use_dust
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*) ' l_use_arclbiom = ',l_use_arclbiom
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*) ' l_use_arclblck = ',l_use_arclblck
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*) ' l_use_arcldlta = ',l_use_arcldlta
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*) ' l_use_arcldust = ',l_use_arcldust
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*) ' l_use_arclocff = ',l_use_arclocff
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*) ' l_use_arclsslt = ',l_use_arclsslt
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*) ' l_use_arclsulp = ',l_use_arclsulp
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*) ' l_use_cariolle = ',l_use_cariolle
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*) ' l_bs1999_abundances = ',l_bs1999_abundances
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,'(A,L7)') ' l_hydrostatic_mass = ',l_hydrostatic_mass
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,'(A,L7)') ' l_moist_heat_capacity = ',l_moist_heat_capacity
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,'(A,L7)') ' l_extra_top = ',l_extra_top
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,'(A,L7)') ' l_use_liu_spec = ',l_use_liu_spec
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' aparam = ',aparam
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,*)' bparam = ',bparam
call umPrint(lineBuffer,src='rad_input_mod')
write(lineBuffer,'(A,L7)') ' l_nlte_corr = ',l_nlte_corr
call umPrint(lineBuffer,src='rad_input_mod')

call umPrint('- - - - - - end of namelist - - - - - -',                        &
    src='rad_input_mod')

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine print_nlist_run_radiation


end module rad_input_mod

