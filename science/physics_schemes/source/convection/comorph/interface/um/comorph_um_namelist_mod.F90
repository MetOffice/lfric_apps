! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!  Global data module for switches/options concerned with convection.

module comorph_um_namelist_mod

  ! Description:
  !   Module containing runtime logicals/options used by the comorph code.
  !
  ! Method:
  !   All switches/options which are contained in the &Run_Comorph
  !   namelist in the CNTLATM control file
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: convection_comorph
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !

use missing_data_mod,       only: rmdi, imdi
use yomhook,                only: lhook, dr_hook
use parkind1,               only: jprb, jpim
use errormessagelength_mod, only: errormessagelength

use um_types,               only: real_umphys

implicit none
save

!==============================================================================
! Logical switches in CoMorph namelist
!==============================================================================

logical :: l_core_ent_cmr = .false.       ! include core/mean factor in
                                          ! comorph parcel core dilution
logical :: l_resdep_precipramp = .false.  ! Include grid-length dependence
                                          ! in the parcel radius precip ramp
logical :: l_conv_inc_w = .false.         ! Switch to apply vertical velocity
                                          ! increments from comorph.

!==============================================================================
! Integer options in CoMorph namelist
!==============================================================================

integer :: par_radius_init_method = imdi  ! Switch for sundry enhancements to
                                          ! parcel radius with options:

integer, parameter :: no_dependence       = 0 ! Constant scaling factor
integer, parameter :: rain_dependence     = 1 ! add dependence on surface precip
integer, parameter :: qfacrain_dependence = 2 ! scale precip dependence by 1/q
integer, parameter :: w_dependence        = 3 ! add further dependence on max w
integer, parameter :: linear_qfacrain_dep = 4 ! scale precip dependence by q

! Comorph internal switches
! (allowed values are stored in comorph_constants_mod).

integer :: par_radius_evol_method = imdi  ! Switch for how parcel radius evolves
                                          ! with height in the plume

integer :: autoc_opt = imdi               ! Switch for autoconversion option

integer :: n_dndraft_types = imdi         ! Number of independent downdraft
                                          ! types

!==============================================================================
! Real values in CoMorph namelist
!==============================================================================

! Scaling factor for comorph turbulent parcel initial perturbations to T,q,u,v
real(kind=real_umphys) :: par_gen_pert_fac = rmdi

! comorph non-turbulent parcel initial RH perturbation.
real(kind=real_umphys) :: par_gen_rhpert = rmdi

! Tuning knob for comorph turbulence-based parcel initial radius.
real(kind=real_umphys) :: par_radius_knob = rmdi

! Max tuning knob for comorph turbulence-based parcel initial radius.
real(kind=real_umphys) :: par_radius_knob_max = rmdi

! Precip rate for max tuning knob for comorph turbulence-based parcel
! initial radius.
real(kind=real_umphys) :: par_radius_ppn_max = rmdi

! Reference grid-length for resolution-dependence / m
real(kind=real_umphys) :: dx_ref = rmdi

! Factor for dilution of comorph parcel core, relative to the mean
! entrainment rate.
real(kind=real_umphys) :: core_ent_fac = rmdi

! Imposed minimum allowed precipitation fraction after comorph.
real(kind=real_umphys) :: rain_area_min = rmdi

! Scaling factor for convective cloud fraction
real(kind=real_umphys) :: cf_conv_fac = rmdi

! Drag coefficient for adjusting the parcel winds towards the environment
! (reduces CMT)
real(kind=real_umphys) :: drag_coef_par = rmdi

! Dimensionless constant scaling the initiation mass-sources
real(kind=real_umphys) :: par_gen_mass_fac = rmdi

! Prescribed draft vertical velocity excess, used when
! the vertical momentum equation is disabled / m s-1
real(kind=real_umphys) :: wind_w_fac = rmdi

! Tuning constant for buoyancy-dependent convective fraction:
! Assuming w' = fac * sqrt( buoyancy * radius )
real(kind=real_umphys):: wind_w_buoy_fac = rmdi

! Minimum parcel initial radius
! (assymptotic value above the BL-top; reduced near the surface)
real(kind=real_umphys) :: ass_min_radius = rmdi

! Rate at which minimum parcel radius increases with height near the surface
! min_radius = MIN( fac * height, ass_min_radius )
real(kind=real_umphys) :: min_radius_fac = rmdi

! Scaling factor for turbulence length-scale (scales all parcel radii)
real(kind=real_umphys) :: turb_len_fac = rmdi

! Scaling factor for par_gen core perturbations relative to
! the parcel mean properties (used if l_par_core = .TRUE.)
real(kind=real_umphys) :: par_gen_core_fac = rmdi

! Entrainment:
! Mixing entrainment rate (m-1) = ent_coef / parcel radius
real(kind=real_umphys) :: ent_coef = rmdi

! Power for overlap between liquid and ice cloud fractions
! inside the parcel
! overlap_power => 1 (no overlap)
! overlap_power => 0 (total overlap)
real(kind=real_umphys) :: overlap_power = rmdi
! Note: do not set to exactly zero as this causes a singularity!

! Max and min limits on the core-mean-ratio (determines assumed in-plume
! PDF shape used for detraiment)
real(kind=real_umphys) :: min_cmr = rmdi
real(kind=real_umphys) :: max_cmr = rmdi


! Plume  microphysics parameters

! In-parcel cloud-to-rain autoconversion rate coefficient
real(kind=real_umphys) :: coef_auto = rmdi

! In-parcel cloud-to-rain autoconversion threshold liquid-water content.
real(kind=real_umphys) :: q_cl_auto = rmdi

! Density of rimed ice (used for graupel)
real(kind=real_umphys) :: rho_rim = rmdi

! Temperature-dependent liquid-cloud number concentration slope
! Liquid-cloud number n(T) = n0 exp( ( T - Tmelt ) / tdep_n_cl )
! ( but limited above Tmelt and below T_homnuc)
real(kind=real_umphys) :: tdep_n_cl = rmdi

! Temperature-dependent ice number concentration slope
! Ice number n(T) = n0 exp( -( T - Tmelt ) / tdep_n_cf )
! ( but limited above Tmelt and below T_homnuc)
real(kind=real_umphys) :: tdep_n_cf = rmdi

! Heterogeneous nucleation temeprature / K
! Gradual freezing starts below this
real(kind=real_umphys) :: hetnuc_temp = rmdi

! Ice crystal non-spherical area factor.  This scales up the assumed surface
! area of ice crystals relative to what it would be if they were spheres
! (affects vapour deposition, fall-speed, riming...)
real(kind=real_umphys) :: cf_area_coef = rmdi

! Asymptotic drag coefficient for a sphere at high Reynolds
! number limit
real(kind=real_umphys) :: drag_coef_cond = rmdi

! Coefficent scaling a term added onto the vapour and heat diffusion,
! for additional exchange due to fall-speed ventilation
real(kind=real_umphys) :: vent_factor = rmdi

! Coefficient for reduction of collection efficiency by
! deflection flow around hydrometeors
real(kind=real_umphys) :: col_eff_coef = rmdi

! Number concentrations for liquid-cloud, ice-cloud, rain, snow, graupel
real(kind=real_umphys) :: nconc_cl = rmdi
real(kind=real_umphys) :: nconc_cf = rmdi
real(kind=real_umphys) :: nconc_rain = rmdi
real(kind=real_umphys) :: nconc_snow = rmdi
real(kind=real_umphys) :: nconc_graup = rmdi

!------------------------------------------------------------------------------
! Define namelist &Run_Comorph read in from CNTLATM control file.
!------------------------------------------------------------------------------

namelist/Run_Comorph/                                                          &

! Integers
par_radius_init_method, par_radius_evol_method, autoc_opt,                     &
n_dndraft_types,                                                               &

! Plume model
par_radius_knob, par_radius_knob_max, par_radius_ppn_max, dx_ref,              &
core_ent_fac, rain_area_min, cf_conv_fac, drag_coef_par, par_gen_rhpert,       &
par_gen_mass_fac, wind_w_fac, wind_w_buoy_fac, par_gen_pert_fac,               &
ass_min_radius, min_radius_fac, turb_len_fac, par_gen_core_fac, overlap_power, &
ent_coef, min_cmr, max_cmr,                                                    &

! Plume  microphysics parameters
rho_rim, tdep_n_cl, tdep_n_cf, hetnuc_temp, cf_area_coef, drag_coef_cond,      &
vent_factor, col_eff_coef, q_cl_auto, coef_auto,                               &
nconc_cl, nconc_cf, nconc_rain, nconc_snow, nconc_graup,                       &

! Logical switches
l_core_ent_cmr, l_resdep_precipramp

!------------------------------------------------------------------------------

!DrHook-related parameters
integer(kind=jpim), parameter, private :: zhook_in  = 0
integer(kind=jpim), parameter, private :: zhook_out = 1

character(len=*), parameter, private :: ModuleName='COMORPH_UM_NAMELIST_MOD'
!------------------------------------------------------------------------------

contains

subroutine check_run_comorph()
!------------------------------------------------------------------------------
! Check comorph namelist values - no code at present
!------------------------------------------------------------------------------

use chk_opts_mod, only: chk_var, def_src
use comorph_constants_mod, only: par_radius_evol_const,                        &
                                 par_radius_evol_volume,                       &
                                 par_radius_evol_no_decrease,                  &
                                 par_radius_evol_no_detrain,                   &
                                 autoc_linear,autoc_quadratic

implicit none

character(len=*), parameter :: RoutineName='CHECK_RUN_COMORPH'

real(kind=jprb)                    :: zhook_handle

!---------------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!---------------------------------------------------------------------------
! required by chk_var routine
def_src = ModuleName//':'//RoutineName

! Checking integer switches have allowed values

call chk_var( par_radius_init_method,'par_radius_init_method',                 &
              [no_dependence, rain_dependence, qfacrain_dependence,            &
               w_dependence, linear_qfacrain_dep] )

call chk_var( par_radius_evol_method,'par_radius_evol_method',                 &
              [par_radius_evol_const, par_radius_evol_volume,                  &
               par_radius_evol_no_decrease, par_radius_evol_no_detrain] )

call chk_var( autoc_opt,'autoc_opt',[autoc_linear,autoc_quadratic] )

call chk_var(n_dndraft_types,'n_dndraft_types','[0,1]')

! Checking reals within allowed range - ranges as in meta-data used for GUI

if (l_resdep_precipramp) call chk_var(dx_ref,'dx_ref','[100.0:1000000.0]')

call chk_var(par_gen_mass_fac,'par_gen_mass_fac','[0.01:1.0]')

call chk_var(drag_coef_cond,'drag_coef_cond','[0.2:1.0]')

call chk_var(vent_factor,'vent_factor','[0.0:1.0]')

call chk_var(col_eff_coef,'col_eff_coef','[0.0:10.0]')

call chk_var(hetnuc_temp,'hetnuc_temp','[230.0:273.0]')

call chk_var(cf_area_coef,'cf_area_coef','[1.0:1000.0]')

call chk_var(wind_w_fac,'wind_w_fac','[0.1:10.0]')

call chk_var(wind_w_buoy_fac,'wind_w_buoy_fac','[0.5:2.0]')

call chk_var(ass_min_radius,'ass_min_radius','[0.0:10000.0]')

call chk_var(min_radius_fac,'min_radius_fac','[0.0:1.0E6]')

call chk_var(turb_len_fac,'turb_len_fac','[1.0:100.0]')

call chk_var(par_gen_core_fac,'par_gen_core_fac','[2.0:6.0]')

call chk_var(overlap_power,'overlap_power','[1.0E-6:1.0]')

call chk_var(ent_coef,'ent_coef','[0.1:0.4]')

call chk_var(min_cmr,'min_cmr','[1.0:3.0]')
call chk_var(max_cmr,'max_cmr','[3.0:10.0]')

!---------------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!---------------------------------------------------------------------------
return
end subroutine check_run_comorph

end module comorph_um_namelist_mod
