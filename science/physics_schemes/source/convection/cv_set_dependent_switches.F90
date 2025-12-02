! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Set variables in module cv_dependent_switch_mod
!
module cv_set_dependent_switches_mod

implicit none

character(len=*),                                                              &
          parameter, private :: ModuleName = 'CV_SET_DEPENDENT_SWITCHES_MOD'
contains

subroutine cv_set_dependent_switches

! Modules used

use cv_dependent_switch_mod, only:                                             &
  l_var_entrain, l_new_det,                                                    &
  sh_on, dp_on, cg_on, md_on,                                                  &
  !mdet_sh_on,
  mdet_dp_on, mdet_md_on,                                                      &
  dp_ent_on, md_ent_on,                                                        &
  sh_sdet_on, dp_sdet_on, cg_sdet_on, md_sdet_on,                              &
  sh_new_termc, dp_new_termc, cg_new_termc, md_new_termc,                      &
  cor_method

use cv_run_mod, only:                                                          &
  icvdiag, adapt, termconv,                                                    &
  i_convection_vn, i_convection_vn_6a, i_cv_comorph,                           &
  l_cv_conserve_check, l_mr_conv, l_wvar_for_conv

use bl_option_mod, only: l_calc_tau_at_p

use cv_param_mod, only: method_en_mx_rho

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! ------------------------------------------------------------------------------
! Description:
!   This routine overrides default values held in cv_dependent_switch_mod
!   with values depending on convection namelist input. Called from readlsta
!   or scm_shell after convection namelist read in.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.1.
!
! ------------------------------------------------------------------------------
! Subroutine arguments - none at present

! Local variables

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CV_SET_DEPENDENT_SWITCHES'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-------------------------------------------------------------------------------
! 1.0 - Variables used in layer_cn  - 5A scheme only, no impact on 4A code
!-------------------------------------------------------------------------------

if (icvdiag == 4 .or. icvdiag == 5) then        ! diurnal cycle diagnosis

  l_var_entrain = .true.      ! Use variable entrainment rate
  l_new_det     = .true.      ! Use new detrainment relationship

end if

!-------------------------------------------------------------------------------
! 2.0 - Variables used in glue_conv and below controlling adaptive options
!-------------------------------------------------------------------------------
! Set flags for adaptive scheme depending on values passed in from namelist
! flags all initialized to zero, so no changes needed if adapt  ==  0

if (adapt  ==  1) then      !HadGEM1a

  md_on = 1
  dp_on = 1
  mdet_dp_on = 1
  mdet_md_on = 1

else if (adapt  ==  2) then  !convection test code

  md_on = 1
  dp_on = 1
  mdet_dp_on = 1
  mdet_md_on = 1
  md_ent_on=1
  dp_ent_on=1

else if (adapt  ==  3) then  !operational

  mdet_dp_on = 1
  dp_on = 1

else if (adapt  ==  4) then  ! Possible HadGEM3 as option 1 plus shallow

  md_on = 1
  dp_on = 1
  mdet_dp_on = 1
  mdet_md_on = 1
  sh_on = 1
  ! mdet_sh_on = 1   Mixing detrainment not set on as has no impact
  cg_on = 1         ! Only used in 5A code

else if (adapt  ==  5) then !adapt det + smoothed forced det for deep + mid

  md_on = 1
  dp_on = 1
  mdet_dp_on = 1
  mdet_md_on = 1
  md_sdet_on=1
  dp_sdet_on=1

else if (adapt  ==  6) then ! as 5 + shallow

  sh_on = 1
  md_on = 1
  dp_on = 1
  ! mdet_sh_on = 1  Mixing detrainment not set on as has no impact
  mdet_md_on = 1
  mdet_dp_on = 1
  sh_sdet_on=1
  md_sdet_on=1
  dp_sdet_on=1

  cg_on = 1          ! Only used in 5A code
  cg_sdet_on=1       ! Only used in 5A code

else if (adapt  ==  7) then ! adapt det + smoothed forced det of theta, q, qcl and qcf for mid and deep

  md_on = 1
  dp_on = 1
  mdet_md_on = 1
  mdet_dp_on = 1
  md_sdet_on=2
  dp_sdet_on=2

else if (adapt  ==  8) then ! as 7 + shallow + congestus

  sh_on = 1
  md_on = 1
  dp_on = 1
  ! mdet_sh_on = 1  Mixing detrainment not set on as has no impact
  mdet_md_on = 1
  mdet_dp_on = 1
  sh_sdet_on=2
  md_sdet_on=2
  dp_sdet_on=2

  cg_on = 1          ! Only used in 5A code
  cg_sdet_on=2       ! Only used in 5A code

end if

! convection termination conditions

if (termconv  ==  0) then

  md_new_termc=0
  dp_new_termc=0

else if (termconv  ==  1) then

  md_new_termc=1
  dp_new_termc=1
  if (adapt == 4 .or. adapt == 6 .or. adapt == 8) then
    sh_new_termc=1          ! use new termination condition
    cg_new_termc=1          ! use new termination condition for congestus
  end if

else if (termconv  ==  2) then
  ! Only available for 6A code
  ! In addition to conditions in termconv=1 also terminates
  ! if the ascent is deeper than that of an undilute parcel ascent
  md_new_termc=2
  dp_new_termc=2
  if (adapt == 4 .or. adapt == 6 .or. adapt == 8) then
    sh_new_termc=2
    cg_new_termc=2
  end if

end if

if (i_convection_vn == i_convection_vn_6a .and. l_cv_conserve_check) then
  ! If the energy correction is switched on in the 6a convection then default
  ! to correcting water and energy on rho-levels in a manner consistent with
  ! the global energy correction.
  cor_method = method_en_mx_rho
end if

!-------------------------------------------------------------------------------
! 3.0 - Settings for the CoMorph convection scheme
!-------------------------------------------------------------------------------
if ( i_convection_vn == i_cv_comorph ) then

  ! CoMorph is written in mixing-ratios, so set flag for convection
  ! using mixing-ratios to true:
  l_mr_conv = .true.

  ! For turbulence-based parcel perturbations, set flags for
  ! calculating turbulent vertical velocity variance,
  ! and the momentum fluxes on the p-grid
  l_calc_tau_at_p = .true.
  l_wvar_for_conv = .true.

end if

!-------------------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine cv_set_dependent_switches
end module cv_set_dependent_switches_mod
