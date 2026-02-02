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
  !   This code is written to UMDP3 programming standards.
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
                                          ! increments from comorph to the UM.

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

integer :: par_radius_evol_method = imdi  ! Switch for how parcel radius evolves
                                          ! with height in the plume
! (allowed values are stored in comorph_constants_mod).

integer :: autoc_opt = imdi               ! Switch for autoconversion option
! (allowed values are stored in comorph_constants_mod).


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


! Plume  microphysics parameters

! In-parcel cloud-to-rain autoconversion rate coefficient
real(kind=real_umphys) :: coef_auto = rmdi

! In-parcel cloud-to-rain autoconversion threshold liquid-water content.
real(kind=real_umphys) :: q_cl_auto = rmdi

! Density of rimed ice (used for graupel)
real(kind=real_umphys) :: rho_rim = rmdi

! Reciprocal of fac_tdep_n
! Temperature-dependent ice number concentration slope
! The number concentration n(T) will be given by:
! n(T) = n0 exp( fac_tdep_n ( T - Tmelt ) )
! ( but limited above Tmelt and below T_homnuc)
real(kind=real_umphys) :: r_fac_tdep_n = rmdi
! set to 8.18 K, consistent with the UM microphysics.

! Heterogeneous nucleation temeprature / K
! Gradual freezing starts below this
real(kind=real_umphys) :: hetnuc_temp = rmdi

! Asymptotic drag coefficient for a sphere at high Reynolds
! number limit
real(kind=real_umphys) :: drag_coef_cond = rmdi

! Coefficent scaling a term added onto the vapour and heat diffusion,
! for additional exchange due to fall-speed ventilation
real(kind=real_umphys) :: vent_factor = rmdi

! Coefficient for reduction of collection efficiency by
! deflection flow around hydrometeors
real(kind=real_umphys) :: col_eff_coef = rmdi

!------------------------------------------------------------------------------
! Define namelist &Run_Comorph read in from CNTLATM control file.
! Changes made to this list will affect both the Full UM and the SCM
!------------------------------------------------------------------------------

namelist/Run_Comorph/                                                          &

! Integers
par_radius_init_method, par_radius_evol_method, autoc_opt,                     &
n_dndraft_types,                                                               &

! Plume model
par_radius_knob, par_radius_knob_max, par_radius_ppn_max, dx_ref,              &
core_ent_fac, rain_area_min, cf_conv_fac, drag_coef_par, par_gen_rhpert,       &
par_gen_mass_fac, wind_w_fac, wind_w_buoy_fac, par_gen_pert_fac,               &
ass_min_radius, par_gen_core_fac, overlap_power, ent_coef,                     &

! Plume  microphysics parameters
rho_rim, r_fac_tdep_n, hetnuc_temp, drag_coef_cond,                            &
vent_factor, col_eff_coef, q_cl_auto, coef_auto,                               &

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

call chk_var(par_gen_mass_fac,'par_gen_mass_fac','[0.01:1.0]')

call chk_var(drag_coef_cond,'drag_coef_cond','[0.2:1.0]')

call chk_var(vent_factor,'vent_factor','[0.0:1.0]')

call chk_var(col_eff_coef,'col_eff_coef','[0.0:10.0]')

call chk_var(hetnuc_temp,'hetnuc_temp','[230.0:273.0]')

call chk_var(wind_w_fac,'wind_w_fac','[0.1:10.0]')

call chk_var(wind_w_buoy_fac,'wind_w_buoy_fac','[0.5:2.0]')

call chk_var(ass_min_radius,'ass_min_radius','[0.0:10000.0]')

call chk_var(par_gen_core_fac,'par_gen_core_fac','[2.0:6.0]')

call chk_var(overlap_power,'overlap_power','[1.0E-6:1.0]')

call chk_var(ent_coef,'ent_coef','[0.1:0.4]')

if (l_resdep_precipramp) call chk_var(dx_ref,'dx_ref','[100.0:1000000.0]')

!---------------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!---------------------------------------------------------------------------
return
end subroutine check_run_comorph


!---------------------------------------------------------------------------
! Prints the CoMorph namelist and checks values
!---------------------------------------------------------------------------
#if !defined(LFRIC)
subroutine print_nlist_run_comorph()

use umPrintMgr, only: umPrint

implicit none
character(len=50000) :: lineBuffer
real(kind=jprb) :: zhook_handle

character(len=*), parameter :: RoutineName='PRINT_NLIST_RUN_COMORPH'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!---------------------------------------------------------------------------
! Section printing the namelist values
!---------------------------------------------------------------------------

call umPrint('Contents of namelist run_comorph', src=ModuleName)

! Integers

write(lineBuffer,"(A,I0)")' par_radius_init_method = ',par_radius_init_method
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,"(A,I0)")' par_radius_evol_method = ',par_radius_evol_method
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,I0)')' autoc_opt = ',autoc_opt
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,'(A,I0)')' n_dndraft_types = ',n_dndraft_types
call umPrint(lineBuffer,src=ModuleName)

! Reals

write(lineBuffer,"(A,ES14.6)")' par_gen_mass_fac = ',par_gen_mass_fac
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,"(A,ES14.6)")' wind_w_fac = ',wind_w_fac
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,"(A,ES14.6)")' wind_w_buoy_fac = ',wind_w_buoy_fac
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,"(A,ES14.6)")' par_gen_pert_fac = ',par_gen_pert_fac
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,"(A,ES14.6)")' ass_min_radius = ',ass_min_radius
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,"(A,ES14.6)")' par_gen_core_fac = ',par_gen_core_fac
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,"(A,ES14.6)")' overlap_power = ',overlap_power
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,"(A,ES14.6)")' ent_coef = ',ent_coef
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,"(A,ES14.6)")' rho_rim = ',rho_rim
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,"(A,ES14.6)")' r_fac_tdep_n = ',r_fac_tdep_n
call umPrint(lineBuffer,src=ModuleName)

write(lineBuffer,"(A,ES14.6)")' par_radius_knob = ',par_radius_knob
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,"(A,ES14.6)")' par_radius_knob_max = ',par_radius_knob_max
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,"(A,ES14.6)")' par_radius_ppn_max = ',par_radius_ppn_max
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,"(A,ES14.6)")' dx_ref = ',dx_ref
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,"(A,ES14.6)")' core_ent_fac = ',core_ent_fac
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,"(A,ES14.6)")' rain_area_min = ',rain_area_min
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,"(A,ES14.6)")' cf_conv_fac = ',cf_conv_fac
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,"(A,ES14.6)")' drag_coef_par = ',drag_coef_par
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,"(A,ES14.6)")' par_gen_rhpert = ',par_gen_rhpert
call umPrint(lineBuffer,src=ModuleName)


write(lineBuffer,"(A,ES14.6)")' hetnuc_temp = ',hetnuc_temp
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,"(A,ES14.6)")' drag_coef_cond = ',drag_coef_cond
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,"(A,ES14.6)")' vent_factor = ',vent_factor
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,"(A,ES14.6)")' col_eff_coef = ',col_eff_coef
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,"(A,ES14.6)")' q_cl_auto = ',q_cl_auto
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,"(A,ES14.6)")' coef_auto = ',coef_auto
call umPrint(lineBuffer,src=ModuleName)


! Logicals

write(lineBuffer,"(A,L1)")' l_core_ent_cmr = ', l_core_ent_cmr
call umPrint(lineBuffer,src=ModuleName)
write(lineBuffer,"(A,L1)")' l_resdep_precipramp = ', l_resdep_precipramp
call umPrint(lineBuffer,src=ModuleName)

call umPrint('- - - - - - end of namelist - - - - - -', src=ModuleName)
!---------------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine print_nlist_run_comorph



!---------------------------------------------------------------------------
! Reads in the CoMorph namelist
!---------------------------------------------------------------------------
subroutine read_nml_run_comorph(unit_in)

use um_parcore, only: mype

use check_iostat_mod, only: check_iostat

use setup_namelist, only: setup_nml_type

implicit none

integer,intent(in) :: unit_in
integer :: my_comm
integer :: mpl_nml_type
integer :: ErrorStatus
integer :: icode
character(len=errormessagelength) :: iomessage
real(kind=jprb) :: zhook_handle

character(len=*), parameter :: RoutineName='READ_NML_RUN_COMORPH'

! set number of each type of variable in my_namelist type
integer, parameter :: no_of_types = 3
integer, parameter :: n_int = 4
integer, parameter :: n_real = 25
integer, parameter :: n_log = 2

type :: my_namelist
  sequence
  integer :: par_radius_init_method
  integer :: par_radius_evol_method
  integer :: autoc_opt
  integer :: n_dndraft_types
  real(kind=real_umphys) :: par_gen_mass_fac
  real(kind=real_umphys) :: wind_w_fac
  real(kind=real_umphys) :: wind_w_buoy_fac
  real(kind=real_umphys) :: par_radius_knob
  real(kind=real_umphys) :: par_radius_knob_max
  real(kind=real_umphys) :: par_radius_ppn_max
  real(kind=real_umphys) :: dx_ref
  real(kind=real_umphys) :: core_ent_fac
  real(kind=real_umphys) :: rain_area_min
  real(kind=real_umphys) :: cf_conv_fac
  real(kind=real_umphys) :: drag_coef_par
  real(kind=real_umphys) :: par_gen_rhpert
  real(kind=real_umphys) :: par_gen_pert_fac
  real(kind=real_umphys) :: ass_min_radius
  real(kind=real_umphys) :: par_gen_core_fac
  real(kind=real_umphys) :: overlap_power
  real(kind=real_umphys) :: ent_coef
  real(kind=real_umphys) :: rho_rim
  real(kind=real_umphys) :: r_fac_tdep_n
  real(kind=real_umphys) :: hetnuc_temp
  real(kind=real_umphys) :: drag_coef_cond
  real(kind=real_umphys) :: vent_factor
  real(kind=real_umphys) :: col_eff_coef
  real(kind=real_umphys) :: q_cl_auto
  real(kind=real_umphys) :: coef_auto
  logical :: l_core_ent_cmr
  logical :: l_resdep_precipramp
end type my_namelist

type (my_namelist) :: my_nml

!-----------------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

call gc_get_communicator(my_comm, icode)

call setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,                 &
                    n_real_in=n_real, n_log_in=n_log)

if (mype == 0) then

  read(unit=unit_in, nml=RUN_Comorph, iostat=ErrorStatus,                      &
       iomsg=iomessage)
  call check_iostat(errorstatus, "namelist RUN_Comorph", iomessage)

  ! Integers
  my_nml % par_radius_init_method  = par_radius_init_method
  my_nml % par_radius_evol_method  = par_radius_evol_method
  my_nml % autoc_opt               = autoc_opt
  my_nml % n_dndraft_types         = n_dndraft_types
  ! end of integers
  ! Reals
  my_nml % par_gen_mass_fac     = par_gen_mass_fac
  my_nml % wind_w_fac           = wind_w_fac
  my_nml % wind_w_buoy_fac      = wind_w_buoy_fac
  my_nml % par_radius_knob      = par_radius_knob
  my_nml % par_radius_knob_max  = par_radius_knob_max
  my_nml % par_radius_ppn_max   = par_radius_ppn_max
  my_nml % dx_ref               = dx_ref
  my_nml % core_ent_fac         = core_ent_fac
  my_nml % rain_area_min        = rain_area_min
  my_nml % cf_conv_fac          = cf_conv_fac
  my_nml % drag_coef_par        = drag_coef_par
  my_nml % par_gen_rhpert       = par_gen_rhpert
  my_nml % par_gen_pert_fac     = par_gen_pert_fac
  my_nml % ass_min_radius       = ass_min_radius
  my_nml % par_gen_core_fac     = par_gen_core_fac
  my_nml % overlap_power        = overlap_power
  my_nml % ent_coef             = ent_coef
  my_nml % rho_rim              = rho_rim
  my_nml % r_fac_tdep_n         = r_fac_tdep_n
  my_nml % hetnuc_temp          = hetnuc_temp
  my_nml % drag_coef_cond       = drag_coef_cond
  my_nml % vent_factor          = vent_factor
  my_nml % col_eff_coef         = col_eff_coef
  my_nml % q_cl_auto            = q_cl_auto
  my_nml % coef_auto            = coef_auto
  ! end of reals
  ! logicals
  my_nml % l_core_ent_cmr       = l_core_ent_cmr
  my_nml % l_resdep_precipramp  = l_resdep_precipramp
  ! end of logicals

end if

call mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

if (mype /= 0) then
  par_radius_init_method  = my_nml % par_radius_init_method
  par_radius_evol_method  = my_nml % par_radius_evol_method
  autoc_opt               = my_nml % autoc_opt
  n_dndraft_types         = my_nml % n_dndraft_types
  ! end of integers
  par_gen_mass_fac     = my_nml % par_gen_mass_fac
  wind_w_fac           = my_nml % wind_w_fac
  wind_w_buoy_fac      = my_nml % wind_w_buoy_fac
  par_radius_knob      = my_nml % par_radius_knob
  par_radius_knob_max  = my_nml % par_radius_knob_max
  par_radius_ppn_max   = my_nml % par_radius_ppn_max
  dx_ref               = my_nml % dx_ref
  core_ent_fac         = my_nml % core_ent_fac
  rain_area_min        = my_nml % rain_area_min
  cf_conv_fac          = my_nml % cf_conv_fac
  drag_coef_par        = my_nml % drag_coef_par
  par_gen_rhpert       = my_nml % par_gen_rhpert
  par_gen_pert_fac     = my_nml % par_gen_pert_fac
  ass_min_radius       = my_nml % ass_min_radius
  par_gen_core_fac     = my_nml % par_gen_core_fac
  overlap_power        = my_nml % overlap_power
  ent_coef             = my_nml % ent_coef
  rho_rim              = my_nml % rho_rim
  r_fac_tdep_n         = my_nml % r_fac_tdep_n
  hetnuc_temp          = my_nml % hetnuc_temp
  drag_coef_cond       = my_nml % drag_coef_cond
  vent_factor          = my_nml % vent_factor
  col_eff_coef         = my_nml % col_eff_coef
  q_cl_auto            = my_nml % q_cl_auto
  coef_auto            = my_nml % coef_auto
  ! end of reals
  l_core_ent_cmr       = my_nml % l_core_ent_cmr
  l_resdep_precipramp  = my_nml % l_resdep_precipramp
end if

call mpl_type_free(mpl_nml_type,icode)

!------------------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine read_nml_run_comorph
#endif

end module comorph_um_namelist_mod
