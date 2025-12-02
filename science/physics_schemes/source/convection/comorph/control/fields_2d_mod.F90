
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

contains

!----------------------------------------------------------------
! Subroutine to set the address of each field in the super-arrays
!----------------------------------------------------------------
subroutine fields_2d_set_addresses()

use comorph_constants_mod, only: l_calc_cape, l_calc_mfw_cape

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

return
end subroutine fields_2d_set_addresses


end module fields_2d_mod
