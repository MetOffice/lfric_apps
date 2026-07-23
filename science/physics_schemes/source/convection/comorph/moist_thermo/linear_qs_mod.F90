! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module linear_qs_mod

implicit none


!----------------------------------------------------------------
! Addresses of linearised qsat fields in a super-array
!----------------------------------------------------------------

! Number of fields contained in a linear_qs super-array
integer, parameter :: n_linear_qs_fields = 5

! Indices of each field in the super-array
integer, parameter :: i_ref_temp = 1
integer, parameter :: i_qsat_liq_ref = 2
integer, parameter :: i_qsat_ice_ref = 3
integer, parameter :: i_dqsatdT_liq = 4
integer, parameter :: i_dqsatdT_ice = 5


contains


!----------------------------------------------------------------
! Subroutine to calculate reference fields for linearised qsat
!----------------------------------------------------------------
subroutine linear_qs_set_ref( n_points,                                        &
                              pressure, linear_qs_super )

use comorph_constants_mod, only: real_cvprec, melt_temp
use set_qsat_mod, only: set_qsat_liq, set_qsat_ice
use set_dqsatdt_mod, only: set_dqsatdt_liq, set_dqsatdt_ice

implicit none

! Number of points
integer, intent(in) :: n_points

! Environment pressure
real(kind=real_cvprec), intent(in) :: pressure(n_points)

! Super-array to store reference temperature and qsat fields
! This routine assumes ref_temp is already set, and just uses
! it to set the qsat and qsat/dT fields.
real(kind=real_cvprec), intent(in out) :: linear_qs_super                      &
                                ( n_points, n_linear_qs_fields )

! Counter and indices for points where above melting-point
integer :: nc
integer :: index_ic(n_points)

! Array to store temperature used in qsat_ice calculation
real(kind=real_cvprec) :: temperature(n_points)

! Loop counters
integer :: ic, ic2


! Set saturation water vapour mixing ratio and dqsat/dT at the
! reference temperature:

! w.r.t. liquid water
call set_qsat_liq( n_points, linear_qs_super(:,i_ref_temp), pressure,          &
                   linear_qs_super(:,i_qsat_liq_ref) )
call set_dqsatdt_liq( n_points,                                                &
                      linear_qs_super(:,i_ref_temp),                           &
                      linear_qs_super(:,i_qsat_liq_ref),                       &
                      linear_qs_super(:,i_dqsatdT_liq) )

! NOTE: many versions of the qsat_ice function don't provide
! values when above freezing (or output values for liquid
! instead of ice).  For CoMorph's linearised qsat calculations,
! any resulting kink in qsat may cause problems when
! extrapolating across the melting threshold.
! To avoid this problem,  where above the melting point, we use
! an extrapolation such that the linearisation should yield the
! correct value of qsat_ice at the melting point.

! Set temperature to ref temp, or melt_temp where above freezing
do ic = 1, n_points
  temperature(ic) = min( linear_qs_super(ic,i_ref_temp),                       &
                         melt_temp )
end do

! Calculate qsat and dqsat/dT w.r.t. ice
call set_qsat_ice( n_points, temperature, pressure,                            &
                   linear_qs_super(:,i_qsat_ice_ref) )
call set_dqsatdt_ice( n_points,                                                &
                      temperature,                                             &
                      linear_qs_super(:,i_qsat_ice_ref),                       &
                      linear_qs_super(:,i_dqsatdT_ice) )

! Find points which are above freezing
nc = 0
do ic = 1, n_points
  if ( linear_qs_super(ic,i_ref_temp) > melt_temp ) then
    nc = nc + 1
    index_ic(nc) = ic
  end if
end do
if ( nc > 0 ) then
  ! At these points, qsat and dqsat/dT w.r.t. ice have been
  ! calculated at Tm instead of Tref.  Extrapolate to find the
  ! value of qsat_ice at Tref such that the linearisation gives
  ! the correct answer at Tm.

  ! qsat(Tm) = qsat(Tref) + dqsat/dT ( Tm - Tref )

  ! => qsat(Tref) = qsat(Tm) + dqsat/dT ( Tref - Tm )
  do ic2 = 1, nc
    ic = index_ic(ic2)
    linear_qs_super(ic,i_qsat_ice_ref)                                         &
      = linear_qs_super(ic,i_qsat_ice_ref)                                     &
      + linear_qs_super(ic,i_dqsatdT_ice)                                      &
        * ( linear_qs_super(ic,i_ref_temp) - melt_temp )
  end do
end if


return
end subroutine linear_qs_set_ref


end module linear_qs_mod
