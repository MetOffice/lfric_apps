! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module lat_heat_mod

implicit none

! Allowable values for a flag which sets what type of
! phase-change to do:
integer, parameter :: i_phase_change_none = 0 ! No phase-change
integer, parameter :: i_phase_change_con = 1  ! condensation
integer, parameter :: i_phase_change_evp = -1 ! evaporation
integer, parameter :: i_phase_change_dep = 2  ! deposition
integer, parameter :: i_phase_change_sub = -2 ! sublimation
integer, parameter :: i_phase_change_frz = 3  ! freezing
integer, parameter :: i_phase_change_mlt = -3 ! melting


! Routines to calculate the latent enthalpies of condensation,
! sublimation and fusion (optionally as a function of temperature).
! To switch off temperature dependence, set cp_liq = cp_ice = cp_vap

contains


!----------------------------------------------------------------
! Latent heat of condensation
!----------------------------------------------------------------
subroutine set_l_con( n_points, temperature, L_con )

use comorph_constants_mod, only: L_con_ref, ref_temp_l, cp_vap, cp_liq,        &
                     real_cvprec

implicit none

! Number of points
integer, intent(in) :: n_points

! Temperature (before the phase-change occurs)
real(kind=real_cvprec), intent(in) :: temperature(n_points)

! Latent heat of condensation
real(kind=real_cvprec), intent(out) :: L_con(n_points)

! Loop counter
integer :: ic

do ic = 1, n_points
  L_con(ic) = L_con_ref - ( cp_liq - cp_vap )                                  &
                        * ( temperature(ic) - ref_temp_l )
end do

return
end subroutine set_l_con


!----------------------------------------------------------------
! Latent heat of sublimation
!----------------------------------------------------------------
subroutine set_l_sub( n_points, temperature, L_sub )

use comorph_constants_mod, only: L_sub_ref, ref_temp_l, cp_vap, cp_ice,        &
                     real_cvprec

implicit none

! Number of points
integer, intent(in) :: n_points

! Temperature (before the phase-change occurs)
real(kind=real_cvprec), intent(in) :: temperature(n_points)

! Latent heat of sublimation
real(kind=real_cvprec), intent(out) :: L_sub(n_points)

! Loop counter
integer :: ic

do ic = 1, n_points
  L_sub(ic) = L_sub_ref - ( cp_ice - cp_vap )                                  &
                        * ( temperature(ic) - ref_temp_l )
end do

return
end subroutine set_l_sub


!----------------------------------------------------------------
! Latent heat of fusion
!----------------------------------------------------------------
subroutine set_l_fus( n_points, temperature, L_fus )

use comorph_constants_mod, only: L_fus_ref, ref_temp_l, cp_liq, cp_ice,        &
                     real_cvprec

implicit none

! Number of points
integer, intent(in) :: n_points

! Temperature (before the phase-change occurs)
real(kind=real_cvprec), intent(in) :: temperature(n_points)

! Latent heat of sublimation
real(kind=real_cvprec), intent(out) :: L_fus(n_points)

! Loop counter
integer :: ic

do ic = 1, n_points
  L_fus(ic) = L_fus_ref - ( cp_ice - cp_liq )                                  &
                        * ( temperature(ic) - ref_temp_l )
end do

return
end subroutine set_l_fus


!----------------------------------------------------------------
! Routine to calculate the temperature increment under
! a given water phase-change increment
!----------------------------------------------------------------
subroutine lat_heat_incr( n_points, nc, i_phase_change,                        &
                          cp_tot, temperature,                                 &
                          index_ic, dq_cmpr, dq )

use comorph_constants_mod, only: cp_vap, cp_liq, cp_ice,                       &
                     L_con_ref, L_sub_ref, L_fus_ref, ref_temp_l,              &
                     real_cvprec, newline

use raise_error_mod, only: raise_fatal

implicit none

! Number of points in full arrays
integer, intent(in) :: n_points

! Number of points where phase-change is to be applied
! (only used if optional compression list index_ic is input,
!  but nc is not allowed to be optional as it is used in the
!  declaration of index_ic)
integer, intent(in) :: nc

! Indicator specifying which type of phase-change to do
integer, intent(in) :: i_phase_change

! Properties of the parcel in before and out after the phase-change:
!   Total heat capacity per unit dry-mass
real(kind=real_cvprec), intent(in out) :: cp_tot(n_points)
!   Air temperature
real(kind=real_cvprec), intent(in out) :: temperature(n_points)

! Indices of points where phase-change occurs
integer, optional, intent(in) :: index_ic(nc)

! Mixing ratio increment changed from phase a to phase b
! (alternative argument used for compressed increment)
real(kind=real_cvprec), optional, intent(in) :: dq_cmpr(nc)

! Mixing ratio increment changed from phase a to phase b
! (argument used for full-field increment)
real(kind=real_cvprec), optional, intent(in) :: dq(n_points)

! Specific heat capacities of the 2 water phases
real(kind=real_cvprec) :: cp_a, cp_b

! Reference value of latent heat of the phase-change from a to b,
! valid at the melting point temperature
real(kind=real_cvprec) :: L_ref

! Loop counter
integer :: ic, ic2

character(len=*), parameter :: routinename = "LAT_HEAT_INCR"


! Set properties according to type of phase change occuring
select case ( i_phase_change )

case ( i_phase_change_none )

  ! No phase change; do nothing and exit.
  return

case ( i_phase_change_con )  ! condensation

  ! Phase a is vapour, phase b is liquid
  cp_a = cp_vap
  cp_b = cp_liq
  L_ref = L_con_ref

case ( i_phase_change_evp )  ! evaporation

  ! Phase a is liquid, phase b is vapour
  cp_a = cp_liq
  cp_b = cp_vap
  L_ref = -L_con_ref

case ( i_phase_change_dep )  ! deposition

  ! Phase a is vapour, phase b is ice
  cp_a = cp_vap
  cp_b = cp_ice
  L_ref = L_sub_ref

case ( i_phase_change_sub )  ! sublimation

  ! Phase a is ice, phase b is vapour
  cp_a = cp_ice
  cp_b = cp_vap
  L_ref = -L_sub_ref

case ( i_phase_change_frz )  ! freezing

  ! Phase a is liquid, phase b is ice
  cp_a = cp_liq
  cp_b = cp_ice
  L_ref = L_fus_ref

case ( i_phase_change_mlt )  ! melting

  ! Phase a is ice, phase b is liquid
  cp_a = cp_ice
  cp_b = cp_liq
  L_ref = -L_fus_ref

end select

if ( present(dq) .and. (.not. present(index_ic)) ) then
  ! Phase-change to be applied at all points

  ! Increment the total cp of the air due to the change of some of
  ! the contained water from one phase to another with a different
  ! specific heat capacity
  do ic = 1, n_points
    cp_tot(ic) = cp_tot(ic) + dq(ic) * ( cp_b - cp_a )
  end do

  ! Increment the temperature
  ! Note: this uses the latent heat L (including temperature
  ! dependence) at the temperature before the phase-change occurs,
  ! but uses cp_tot based on the mixing ratios after the
  ! phase-change; this is the correct order of calculation to
  ! conserve moist enthalpy at constant pressure.
  do ic = 1, n_points
    temperature(ic) = temperature(ic)                                          &
     + ( ( L_ref - (cp_b - cp_a)*(temperature(ic) - ref_temp_l) )              &
         / cp_tot(ic) )  *  dq(ic)
  end do

else if ( present(dq) .and. present(index_ic) ) then
  ! Compressed version of exactly the same calculation

  do ic2 = 1, nc
    ic = index_ic(ic2)

    cp_tot(ic) = cp_tot(ic) + dq(ic) * ( cp_b - cp_a )

    temperature(ic) = temperature(ic)                                          &
     + ( ( L_ref - (cp_b - cp_a)*(temperature(ic) - ref_temp_l) )              &
         / cp_tot(ic) )  *  dq(ic)
  end do

else if ( present(dq_cmpr) .and. present(index_ic) ) then
  ! Compressed version where the input increment is already
  ! in compressed form

  do ic2 = 1, nc
    ic = index_ic(ic2)

    cp_tot(ic) = cp_tot(ic) + dq_cmpr(ic2) * ( cp_b - cp_a )

    temperature(ic) = temperature(ic)                                          &
     + ( ( L_ref - (cp_b - cp_a)*(temperature(ic) - ref_temp_l) )              &
         / cp_tot(ic) )  *  dq_cmpr(ic2)
  end do

else

  ! Unworkable combination of inputs present / not present
  call raise_fatal( routinename,                                               &
    "Not allowed: throwing a wobbly!  This routine has been " //               &
    "called with an invalid combination of "         //newline//               &
    "optional input arguments." )

end if

return
end subroutine lat_heat_incr


end module lat_heat_mod
