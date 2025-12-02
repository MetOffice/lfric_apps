! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module apply_scaling_mod

implicit none

contains

! Subroutine to apply closure scalings to the resolved-scale
! source terms for a given category of convective draft
subroutine apply_scaling( n_conv_types, n_conv_layers,                         &
                          n_fields_tot, n_points_any,                          &
                          ij_first, ij_last, draft_scaling,                    &
                          res_source_k )

use comorph_constants_mod, only: real_cvprec, nx_full
use res_source_type_mod, only: res_source_type, n_res
use cloudfracs_type_mod, only: n_convcloud

implicit none

! Number of convection types
integer, intent(in) :: n_conv_types
! Number of distinct convecting layers
integer, intent(in) :: n_conv_layers

! Total number of model-fields to update
integer, intent(in) :: n_fields_tot

! Number of points where any type / layer of convection is active
integer, intent(in) :: n_points_any

! First and last ij-indices on the current segment
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last

! Closure scaling for each type / layer on ij-grid
real(kind=real_cvprec), intent(in) :: draft_scaling                            &
               ( ij_first:ij_last, n_conv_types, n_conv_layers )

! Structure containing resolved-scale source terms to be rescaled
type(res_source_type), intent(in out) :: res_source_k                          &
                             ( n_conv_types, n_conv_layers )

! Closure scaling compressed onto convecting points on the
! current level / layer / type
real(kind=real_cvprec) :: scaling_cmpr(n_points_any)

! Loop counters
integer :: i_type, i_layr, i, j, ic, i_field

! Loop over convecting layers and convection types
do i_layr = 1, n_conv_layers
  do i_type = 1, n_conv_types

    ! If there are any resolved-scale source terms
    ! for this convection layer / type at the current level
    if ( res_source_k(i_type,i_layr) % cmpr % n_points > 0 ) then

      ! Compress closure scaling onto active points
      do ic = 1, res_source_k(i_type,i_layr) % cmpr % n_points
        i = res_source_k(i_type,i_layr) % cmpr % index_i(ic)
        j = res_source_k(i_type,i_layr) % cmpr % index_j(ic)
        scaling_cmpr(ic) = draft_scaling                                       &
                           ( nx_full*(j-1)+i, i_type, i_layr )
      end do

      ! Rescale entrainment and detrainment
      do i_field = 1, n_res
        do ic = 1, res_source_k(i_type,i_layr) % cmpr % n_points
          res_source_k(i_type,i_layr) % res_super(ic,i_field)                  &
            = res_source_k(i_type,i_layr) % res_super(ic,i_field)              &
            * scaling_cmpr(ic)
        end do
      end do

      ! Rescale source terms for primary fields
      do i_field = 1, n_fields_tot
        do ic = 1, res_source_k(i_type,i_layr) % cmpr % n_points
          res_source_k(i_type,i_layr) % fields_super(ic,i_field)               &
           =res_source_k(i_type,i_layr)%fields_super(ic,i_field)               &
            * scaling_cmpr(ic)
        end do
      end do

      ! Rescale convective cloud fields
      if ( n_convcloud > 0 ) then
        do i_field = 1, n_convcloud
          do ic = 1, res_source_k(i_type,i_layr) % cmpr%n_points
            res_source_k(i_type,i_layr)                                        &
                                   % convcloud_super(ic,i_field)               &
              = res_source_k(i_type,i_layr)                                    &
                                   % convcloud_super(ic,i_field)               &
               * scaling_cmpr(ic)
          end do
        end do
      end if

    end if !( res_source_k(i_type,i_layr) % cmpr % n_points > 0 )

  end do  ! i_type = 1, n_conv_types
end do  ! i_layr = 1, n_conv_layers

return
end subroutine apply_scaling


end module apply_scaling_mod
