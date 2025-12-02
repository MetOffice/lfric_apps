! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module add_conv_cloud_mod

implicit none

contains

! Subroutine to add diagnosed convective cloud contributions for
! a given convective draft on to compressed copies
subroutine add_conv_cloud( n_points,                                           &
                           n_conv_types, n_conv_layers,                        &
                           ij_first, ij_last, index_ic_k,                      &
                           res_source_k, conv_cloud )

use comorph_constants_mod, only: real_cvprec, nx_full
use res_source_type_mod, only: res_source_type
use cloudfracs_type_mod, only: n_convcloud

implicit none

! Number of convecting points on current level (dimensions super)
integer, intent(in) :: n_points
! Number of convection types
integer, intent(in) :: n_conv_types
! Number of distinct convecting layers
integer, intent(in) :: n_conv_layers

! First and last ij-indices on the current segment
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last
! cmpr_any compression list addresses scattered into a
! common array for referencing from res_source compression list
integer, intent(in) :: index_ic_k(ij_first:ij_last)

! Structure containing resolved-scale source terms to be added
type(res_source_type), intent(in) :: res_source_k                              &
                         ( n_conv_types, n_conv_layers )

! Super-array containing convective cloud fields to increment
real(kind=real_cvprec), intent(in out) :: conv_cloud                           &
                                       ( n_points, n_convcloud )

! Store for cmpr_any compression list indices
integer :: ic_any( n_points )

! Loop counters
integer :: i_type, i_layr, i_field, i, j, ij, ic


! Loop over convecting layers and convection types
do i_layr = 1, n_conv_layers
  do i_type = 1, n_conv_types

    ! If there are any resolved-scale source terms
    ! for this convection layer / type at the current level
    if ( res_source_k(i_type,i_layr) % cmpr % n_points > 0 ) then

      ! Extract the cmpr_any compression list index from the grid
      do ic = 1, res_source_k(i_type,i_layr) % cmpr % n_points
        i = res_source_k(i_type,i_layr) % cmpr % index_i(ic)
        j = res_source_k(i_type,i_layr) % cmpr % index_j(ic)
        ij = nx_full*(j-1)+i
        ic_any(ic) = index_ic_k(ij)
      end do

      ! Increment the convective cloud fields
      do i_field = 1, n_convcloud
        do ic = 1, res_source_k(i_type,i_layr) % cmpr % n_points
          conv_cloud(ic_any(ic),i_field)                                       &
            = conv_cloud(ic_any(ic),i_field)                                   &
            + res_source_k(i_type,i_layr)                                      &
              % convcloud_super(ic,i_field)
        end do
      end do

    end if

  end do
end do


return
end subroutine add_conv_cloud

end module add_conv_cloud_mod
