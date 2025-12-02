! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module mass_rearrange_calc_mod

implicit none

contains

! Subroutine to compute the amount of mass from layer k
! to deposit on layer k2, during the column mass-rearrangement
! calculation.
subroutine mass_rearrange_calc( cmpr_add, n_points_incr_k2,                    &
                                index_ic_incr_k2,                              &
                                ij_first, ij_last, ij,                         &
                                layer_mass_k2, layer_mass_k2_added,            &
                                layer_mass_k, frac_k2_from_k )

use cmpr_type_mod, only: cmpr_type
use comorph_constants_mod, only: real_cvprec, zero

implicit none

! Structure storing compression list info for points where
! mass is to be added
type(cmpr_type), intent(in) :: cmpr_add

! Number of points where the
! k2 level counter is to be incremented
integer, intent(out) :: n_points_incr_k2

! Index list for referencing the cmpr_add compression list
! from the cmpr_incr_k2 compression list, used to
! compress onto points where k2 is incremented
integer, intent(out) :: index_ic_incr_k2( cmpr_add%n_points )

! First and last ij-indices on the current segment
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last

! ij-indices of all points on the segment
integer, intent(in) :: ij(cmpr_add%n_points)

! Target layer-masses on model-level k2
real(kind=real_cvprec), intent(in) :: layer_mass_k2(cmpr_add%n_points)

! Amount of dry-mass so far added to the current model-layer k2.
real(kind=real_cvprec), intent(in out) :: layer_mass_k2_added(ij_first:ij_last)

! Amount of dry-mass from layer k not yet added to any layer
real(kind=real_cvprec), intent(in out) :: layer_mass_k(cmpr_add%n_points)

! Fraction of mass on level k2 to come from level k
! (used as the weight for calculating new mean field values
!  on level k2)
real(kind=real_cvprec), intent(out) :: frac_k2_from_k(cmpr_add % n_points)

! Loop counters
integer :: ic, ic_incr


! Save amount of dry-mass remaining unassigned on level k
! (using the output array as work storage)
do ic = 1, cmpr_add % n_points
  frac_k2_from_k(ic) = layer_mass_k(ic)
end do

! Increment counter for amount of dry-mass added at k2
do ic = 1, cmpr_add % n_points
  layer_mass_k2_added(ij(ic)) = layer_mass_k2_added(ij(ic))                    &
                              + layer_mass_k(ic)
end do

! Store amount of dry-mass from k that doesn't fit on
! level k2 and therefore needs to be added to other layers
do ic = 1, cmpr_add % n_points
  layer_mass_k(ic) = layer_mass_k2_added(ij(ic))                               &
                   - layer_mass_k2(ic)
end do

! Find points where k2 needs to be incremented, and define
! the compression list of points where it does
n_points_incr_k2 = 0
do ic = 1, cmpr_add % n_points
  ! If there is dry-mass on k that doesn't fit on k2...
  if ( layer_mass_k(ic) >= zero ) then
    ! Add this point to the list
    n_points_incr_k2 = n_points_incr_k2 + 1
    index_ic_incr_k2(n_points_incr_k2) = ic
  end if
end do

! At points where k2 will be incremented...
! Set the amount of mass being added to the current k2 level
! = the total amount we were going to add, minus the
!   bit that we're carrying over to the next level instead
do ic_incr = 1, n_points_incr_k2
  ic = index_ic_incr_k2(ic_incr)
  frac_k2_from_k(ic) = frac_k2_from_k(ic) - layer_mass_k(ic)
end do
! At other points, the mass to add at k2 will be all the mass
! that we have at k, which is already stored in frac_k2_from_k.

! At all points where mass is to be added...
! Divide dry-mass to be added by total dry-mass on k2,
! to obtain the fraction of mass used as the weight for
! field properties from layer k in the new means at k2
do ic = 1, cmpr_add % n_points
  frac_k2_from_k(ic) = frac_k2_from_k(ic) / layer_mass_k2(ic)
end do


return
end subroutine mass_rearrange_calc


end module mass_rearrange_calc_mod
