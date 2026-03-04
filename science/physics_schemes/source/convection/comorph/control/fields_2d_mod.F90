
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module fields_2d_mod

use comorph_constants_mod, only: real_cvprec

implicit none

save

! Flag indicating whether the variables contained in this
! module have been set
logical :: l_init_fields_2d_mod = .false.

!----------------------------------------------------------------
! Addresses in a super-array storing 2D fields used in comorph
!----------------------------------------------------------------

! Total number of 2D fields in use
integer :: n_fields_2d = 0

! Adresses of fields in the super-array
! (need to be set based on switches)

! CAPE
integer :: i_cape = 0
! Vertical integral of mass-flux used to compute reference mass-flux
integer :: i_int_m_dz = 0
! Reference mass-flux used for mass-flux-weighted CAPE
integer :: i_m_ref = 0
! Mass-flux-weighted CAPE
integer :: i_mfw_cape = 0

! Convection termination height
integer :: i_term_height = 0
! Cloud-base height
integer :: i_ccb_height = 0
! Cloud-base area fraction
integer :: i_ccb_area = 0
! Cloud-base mass-flux
integer :: i_ccb_massflux_d = 0
! Cloud-base mean buoyancy
integer :: i_ccb_mean_tv_ex = 0
! Cloud-base core buoyancy
integer :: i_ccb_core_tv_ex = 0

contains

!----------------------------------------------------------------
! Subroutine to set the address of each field in the super-arrays
!----------------------------------------------------------------
subroutine fields_2d_set_addresses()

use comorph_constants_mod, only: l_calc_cape, l_calc_mfw_cape, l_calc_ccb_cct

implicit none

! Set flag to indicate that we don't need to call this routine again!
l_init_fields_2d_mod = .true.

if ( l_calc_cape ) then
  ! CAPE
  i_cape = n_fields_2d + 1
  ! Increment counter for 2D fields
  n_fields_2d = n_fields_2d + 1
end if

if ( l_calc_mfw_cape ) then
  ! Mass-flux-weighted CAPE and vertical integrals required to compute it.
  i_int_m_dz  = n_fields_2d + 1
  i_m_ref     = n_fields_2d + 2
  i_mfw_cape  = n_fields_2d + 3
  ! Increment counter for 2D fields
  n_fields_2d = n_fields_2d + 3
end if

if ( l_calc_ccb_cct ) then
  ! Convective cloud base and top heights and other properties
  i_term_height    = n_fields_2d + 1
  i_ccb_height     = n_fields_2d + 2
  i_ccb_area       = n_fields_2d + 3
  i_ccb_massflux_d = n_fields_2d + 4
  i_ccb_mean_tv_ex = n_fields_2d + 5
  i_ccb_core_tv_ex = n_fields_2d + 6
  n_fields_2d = n_fields_2d + 6
end if

return
end subroutine fields_2d_set_addresses


!----------------------------------------------------------------
! Subroutine to apply convective closure scalings to the 2D fields
!----------------------------------------------------------------
subroutine fields_2d_scaling( ij_first, ij_last,                               &
                              n_conv_types, n_conv_layers,                     &
                              closure_scaling, fields_2d )

use comorph_constants_mod, only: real_cvprec

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
integer :: ij, i_type, i_layr, i_field

! Loop over the 2D fields
do i_field = 1, n_fields_2d

  if ( i_field == i_int_m_dz .or.                                              &
       i_field == i_m_ref .or.                                                 &
       i_field == i_ccb_area .or.                                              &
       i_field == i_ccb_massflux_d ) then
    ! Those fields which are mass-fluxes or convective area fraction must be
    ! scaled by the convective closure scaling.  The other fields are
    ! invarient to the closure scaling

    ! Loop over convection layers, types and points
    do i_layr = 1, n_conv_layers
      do i_type = 1, n_conv_types
        do ij = ij_first, ij_last
          ! Apply closure scaling
          fields_2d(ij,i_field,i_type,i_layr)                                  &
            = fields_2d(ij,i_field,i_type,i_layr)                              &
            * closure_scaling(ij,i_type,i_layr)
        end do
      end do
    end do

  end if

end do  ! i_field = 1, n_fields_2d

return
end subroutine fields_2d_scaling


end module fields_2d_mod
