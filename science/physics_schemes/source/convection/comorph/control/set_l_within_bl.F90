! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module set_l_within_bl_mod

implicit none

contains

! Subroutine to set flag for whether the current model-level
! (the theta-level, where the primary thermodynamic variables
! are defined) is within the boundary-layer.
! This is defined based on whether the model-level's lower
! and upper interfaces have any overlap with the boundary-layer
! top height passed in from the boundary-layer scheme.

subroutine set_l_within_bl( cmpr, k,                                           &
                            lb_h, ub_h, height_half, lb_t, ub_t, z_bl_top,     &
                            l_within_bl )

use comorph_constants_mod, only: real_hmprec
use cmpr_type_mod, only: cmpr_type

implicit none

! Structure containing number of points and compression list
! indices for referencing the full 2-D array z_bl_top
type(cmpr_type), intent(in) :: cmpr

! Index of current model-level (the theta-level)
integer, intent(in) :: k

! Full 3-D array of heights of the model-level interfaces
! above the surface (and its array bounds, so we index
! it right if it has halos)
integer, intent(in) :: lb_h(3), ub_h(3)
real(kind=real_hmprec), intent(in) :: height_half                              &
                          ( lb_h(1):ub_h(1), lb_h(2):ub_h(2), lb_h(3):ub_h(3) )


! Full 2-D array of boundary-layer top height above surface
integer, intent(in) :: lb_t(2), ub_t(2)
real(kind=real_hmprec), intent(in) :: z_bl_top                                 &
                               ( lb_t(1):ub_t(1), lb_t(2):ub_t(2) )

! Output flag for whether current level is within the BL
logical, intent(out) :: l_within_bl(cmpr%n_points)

! Loop counters
integer :: ic, i, j


! Test whether the current model-level's lower interface
! is below the BL-top height (level has at least some overlap
! with the BL).
do ic = 1, cmpr % n_points
  i = cmpr % index_i(ic)
  j = cmpr % index_j(ic)
  ! Assuming that height_half(:,:,k) is the lower bound
  ! of theta-level (:,:,k).
  l_within_bl(ic) = height_half(i,j,k) <= z_bl_top(i,j)
end do


return
end subroutine set_l_within_bl

end module set_l_within_bl_mod
