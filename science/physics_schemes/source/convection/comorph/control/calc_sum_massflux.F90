! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_sum_massflux_mod

implicit none

contains


! Subroutine to sum the mass-fluxes and initiating masses
! over all convection types and layers on the current level.
! The summed mass-fluxes are needed for numerical safety
! checks within conv_level-step (e.g. to avoid the ensemble
! of updrafts entraining more mass than exists on a level).

subroutine calc_sum_massflux( n_conv_types, n_conv_layers,                     &
                              ij_first, ij_last,                               &
                              par_conv, sum_massflux )

use comorph_constants_mod, only: real_cvprec, zero, nx_full
use parcel_type_mod, only: parcel_type, i_massflux_d

implicit none

! Number of convection types
integer, intent(in) :: n_conv_types
! Number of distinct convecting layers in a column
integer, intent(in) :: n_conv_layers

! ij indices of first & last point on current segment
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last

! Array of pointers to structures containing compression lists
! of parcel properties
! at points where convection is active at the current level
type(parcel_type), intent(in) :: par_conv                                      &
                   ( n_conv_types, n_conv_layers )

! Sum of parcel mass-fluxes over types / layers
! at previous model-level interface
real(kind=real_cvprec), intent(out) :: sum_massflux                            &
                                       ( ij_first:ij_last )

! Loop counters
integer :: i, j, ij, ic, i_type, i_layr


! Initialise sums to zero
do ij = ij_first, ij_last
  sum_massflux(ij) = zero
end do


! Loop over layers and types
do i_layr = 1, n_conv_layers
  do i_type = 1, n_conv_types

    ! If any mass-flux coming up from previous level interface
    if ( par_conv(i_type,i_layr) % cmpr % n_points > 0 ) then
      ! Add contributions from points in parcel compression list
      do ic = 1, par_conv(i_type,i_layr) % cmpr % n_points
        i = par_conv(i_type,i_layr) % cmpr % index_i(ic)
        j = par_conv(i_type,i_layr) % cmpr % index_j(ic)
        ij = nx_full*(j-1)+i
        sum_massflux(ij) = sum_massflux(ij)                                    &
                         + par_conv(i_type,i_layr)                             &
                           % par_super(ic,i_massflux_d)
      end do
    end if

  end do
end do


return
end subroutine calc_sum_massflux


end module calc_sum_massflux_mod
