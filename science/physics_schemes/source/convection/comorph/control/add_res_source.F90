! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module add_res_source_mod

implicit none

contains

! Subroutine to add resolved-scale source terms for a given
! convective draft onto compressed copies of the full fields.
subroutine add_res_source( n_points_super, n_fields_tot,                       &
                           n_conv_types, n_conv_layers,                        &
                           ij_first, ij_last, index_ic_k,                      &
                           res_source_k,                                       &
                           fields_k_super, layer_mass_k )

use comorph_constants_mod, only: real_cvprec, comorph_timestep, nx_full
use res_source_type_mod, only: res_source_type, i_ent, i_det

implicit none

! Max number of convective points on one level (dimensions super)
integer, intent(in) :: n_points_super
! Total number of model-fields to update
integer, intent(in) :: n_fields_tot
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

! Super-array containing fields to be incremented
real(kind=real_cvprec), intent(in out) :: fields_k_super                       &
                             ( n_points_super, n_fields_tot )

! Layer-masses to be incremented
real(kind=real_cvprec), intent(in out) :: layer_mass_k                         &
                             ( n_points_super )

! Store for cmpr_any compression list indices
integer :: ic_any( n_points_super )

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

      ! Increment the layer-mass with detrainment - entrainment
      do ic = 1, res_source_k(i_type,i_layr) % cmpr % n_points
        layer_mass_k(ic_any(ic)) = layer_mass_k(ic_any(ic))                    &
                                 + ( res_source_k(i_type,i_layr)               &
                                       % res_super(ic,i_det)                   &
                                   - res_source_k(i_type,i_layr)               &
                                       % res_super(ic,i_ent)                   &
                                   ) * comorph_timestep
      end do

      ! Increment the primary fields with the source terms
      do i_field = 1, n_fields_tot
        do ic = 1, res_source_k(i_type,i_layr) % cmpr % n_points
          fields_k_super(ic_any(ic),i_field)                                   &
                  = fields_k_super(ic_any(ic),i_field)                         &
                  + res_source_k(i_type,i_layr)                                &
                    % fields_super(ic,i_field) * comorph_timestep
        end do
      end do

    end if

  end do
end do


return
end subroutine add_res_source

end module add_res_source_mod
