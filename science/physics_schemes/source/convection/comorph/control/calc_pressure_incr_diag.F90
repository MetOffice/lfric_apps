! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_pressure_incr_diag_mod

implicit none

contains

! Subroutine to calculate the mean pressure that the air now on level
! k2 had at its original level(s) before the convective subsidence
! occurred
! (this routine adds the contribution to that mean from the current
!  level k, and is called multiple times to sum the contributions
!  from all k).
! It performs only a trivial calculation, but is worth having a separate
! subroutine for it because:
! a) It saves duplicating the calculation multiple times in
!    mass_rearrange.F90
! b) The pressure diagnostic is stored as a pointer in the comorph_diags
!    structure; that pointer can be passed into here where it is
!    declared as a normal array, allowing the calculation to vectorise
!    more efficiently.
subroutine calc_pressure_incr_diag( cmpr, lb_p, ub_p, lb_d, ub_d,              &
                                    ij_first, ij_last, k, k2, ij,              &
                                    frac_k2_from_k, pressure,                  &
                                    pressure_incr_env )

use comorph_constants_mod, only: real_cvprec, real_hmprec
use cmpr_type_mod, only: cmpr_type

implicit none

! Structure storing i,j indices of point
type(cmpr_type), intent(in) :: cmpr

! Lower and upper bounds of the 3-D arrays, in case they have halos
integer, intent(in) :: lb_p(3), ub_p(3)
integer, intent(in) :: lb_d(3), ub_d(3)

! First and last ij-indices on the current segment
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last

! Current model-level
integer, intent(in) :: k
! Model-level that fields from level k are being scattered into
integer, intent(in) :: k2(ij_first:ij_last)

! Index list for referencing compressed fields on common grid
integer, intent(in) :: ij( cmpr % n_points )

! Fraction of mass on level k2 to come from level k
! (used as the weight for calculating new mean field values
!  on level k2)
real(kind=real_cvprec), intent(in) :: frac_k2_from_k( cmpr % n_points )

! 3-D array of pressure on thermodynamic levels
real(kind=real_hmprec), intent(in) :: pressure                                 &
                          ( lb_p(1):ub_p(1), lb_p(2):ub_p(2), lb_p(3):ub_p(3) )

! Diagnostic pressure field to be updated; will store the mean
! pressure that the air now on level k2 had at its original level(s)
! before the convective subsidence occurred
! (it later gets converted to a Lagrangian pressure increment, hence the name)
real(kind=real_hmprec), intent(in out) :: pressure_incr_env                    &
                          ( lb_d(1):ub_d(1), lb_d(2):ub_d(2), lb_d(3):ub_d(3) )

! Loop counters
integer :: ic, i, j

do ic = 1, cmpr % n_points
  i = cmpr % index_i(ic)
  j = cmpr % index_j(ic)
  pressure_incr_env(i,j,k2(ij(ic)))                                            &
    = pressure_incr_env(i,j,k2(ij(ic)))                                        &
    + pressure(i,j,k) * real(frac_k2_from_k(ic), real_hmprec)
end do

return
end subroutine calc_pressure_incr_diag

end module calc_pressure_incr_diag_mod
