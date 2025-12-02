! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module cfl_limit_scaling_mod

implicit none

contains


! Subroutine to scale the closure scaling down where needed to
! avoid violating the CFL condition on mass entrained from a
! single model-level
subroutine cfl_limit_scaling( n_conv_types, n_conv_layers,                     &
                              n_points_any,                                    &
                              cfl_scaling_k, res_source_k,                     &
                              ij_first, ij_last, index_ic_k,                   &
                              draft_cfl_scaling )

use comorph_constants_mod, only: real_cvprec, nx_full
use res_source_type_mod, only: res_source_type

implicit none

! Number of convection types
integer, intent(in) :: n_conv_types
! Number of distinct convecting layers
integer, intent(in) :: n_conv_layers

! Number of points on the current level with any convection
integer, intent(in) :: n_points_any

! CFL scaling limit for only the current model-level
real(kind=real_cvprec), intent(in) :: cfl_scaling_k                            &
                                      ( n_points_any )

! Structure containing resolved-scale source terms (only needed
! here to get the compression list of where each type / layer
! is active)
type(res_source_type), intent(in) :: res_source_k                              &
                                 ( n_conv_types, n_conv_layers )

! Position indices ( nx*(j-1) + i ) of the first and last
! convecting point on the current segment
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last
! cmpr_any compression list addresses scattered into a
! common array for referencing from res_source compression list
integer, intent(in) :: index_ic_k(ij_first:ij_last)

! Minimum CFL scaling limit over all levels where each
! convection type/layer is active
real(kind=real_cvprec), intent(in out) :: draft_cfl_scaling                    &
         ( ij_first:ij_last, n_conv_types, n_conv_layers )

! Loop counters
integer :: i, j, ij, ic, ic2, i_type, i_layr


! Loop over convection types / layers
do i_layr = 1, n_conv_layers
  do i_type = 1, n_conv_types

    ! If any convection of this type
    if ( res_source_k(i_type,i_layr) % cmpr % n_points > 0 ) then
      do ic2 = 1, res_source_k(i_type,i_layr) % cmpr % n_points
        i = res_source_k(i_type,i_layr) % cmpr % index_i(ic2)
        j = res_source_k(i_type,i_layr) % cmpr % index_j(ic2)
        ij = nx_full*(j-1)+i
        ic = index_ic_k(ij)
        ! Find minimum rescaling factor to avoid exceeding
        ! CFL limit at any level
        draft_cfl_scaling(ij,i_type,i_layr)                                    &
          = min( draft_cfl_scaling(ij,i_type,i_layr),                          &
                 cfl_scaling_k(ic) )
      end do
    end if

  end do
end do


return
end subroutine cfl_limit_scaling


end module cfl_limit_scaling_mod
