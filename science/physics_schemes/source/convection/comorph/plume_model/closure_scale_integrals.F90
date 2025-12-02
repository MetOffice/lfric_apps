! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module closure_scale_integrals_mod

implicit none

contains

! Subroutine to apply the convective closure scaling to any vertical
! integral calculations (2D fields) that scale with the mass-flux.
subroutine closure_scale_integrals( ij_first, ij_last,                         &
                                    n_conv_types, n_conv_layers,               &
                                    closure_scaling, fields_2d )

use comorph_constants_mod, only: real_cvprec
use fields_2d_mod, only: n_fields_2d, i_int_m_dz, i_m_ref

implicit none

! i,j indices of the first and last point on the current segment
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last

! Number of convection types
integer, intent(in) :: n_conv_types

! Largest number of convective mass-source layers found in any
! column on the current segment
integer, intent(in) :: n_conv_layers

! Closure scalings for each type / layer on ij-grid
real(kind=real_cvprec), intent(in) :: closure_scaling                          &
                                      ( ij_first:ij_last,                      &
                                        n_conv_types, n_conv_layers )

! Super-array storing 2D work arrays
real(kind=real_cvprec), intent(in out) :: fields_2d                            &
                                          ( ij_first:ij_last, n_fields_2d,     &
                                            n_conv_types, n_conv_layers )

! Loop counters
integer :: ij, i_type, i_layr


! Apply closure scaling to the vertical integrals used to compute
! the mass-flux-weighted CAPE...

do i_layr = 1, n_conv_layers
  do i_type = 1, n_conv_types
    do ij = ij_first, ij_last

      ! Scale integrals of massflux by the massflux closure scaling

      fields_2d(ij,i_m_ref,i_type,i_layr)                                      &
        = fields_2d(ij,i_m_ref,i_type,i_layr)                                  &
        * closure_scaling(ij,i_type,i_layr)

      fields_2d(ij,i_int_m_dz,i_type,i_layr)                                   &
        = fields_2d(ij,i_int_m_dz,i_type,i_layr)                               &
        * closure_scaling(ij,i_type,i_layr)

      ! Note that the mass-flux-weighted CAPE for a single type/layer
      ! has already been normalised by the reference mass-flux m_ref, such that
      ! it is invariant to the closure scaling.
      ! However, the closure scalings matter here
      ! if the integrals are subsequently summed over layers / types
      ! to compute a combined mass-flux-weighted CAPE.

    end do
  end do
end do


return
end subroutine closure_scale_integrals

end module closure_scale_integrals_mod
