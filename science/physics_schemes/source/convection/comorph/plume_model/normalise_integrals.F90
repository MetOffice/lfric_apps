! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module normalise_integrals_mod

implicit none

contains

! Subroutine to complete any vertical integral calculations (2D fields)
! performed in the parcel ascents and descents.  Some quantities
! need a final normalisation after the parcel ascents and descents have
! finished; such normalisations are done in here...
subroutine normalise_integrals( ij_first, ij_last,                             &
                                n_conv_types, n_conv_layers,                   &
                                fields_2d )

use comorph_constants_mod, only: real_cvprec, min_float
use fields_2d_mod, only: n_fields_2d, i_int_m_dz, i_m_ref, i_mfw_cape

implicit none

! i,j indices of the first and last point on the current segment
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last

! Number of convection types
integer, intent(in) :: n_conv_types

! Largest number of convective mass-source layers found in any
! column on the current segment
integer, intent(in) :: n_conv_layers

! Super-array storing 2D work arrays
real(kind=real_cvprec), intent(in out) :: fields_2d                            &
                                          ( ij_first:ij_last, n_fields_2d,     &
                                            n_conv_types, n_conv_layers )

! Loop counters
integer :: ij, i_type, i_layr


! Normalise the mass-flux-weighted CAPE integral...

do i_layr = 1, n_conv_layers
  do i_type = 1, n_conv_types
    do ij = ij_first, ij_last

      ! Finish computing the reference mass-flux for each layer and type;
      ! this is a vertical mean of the mass-flux, weighted by itself:
      ! M_ref   =   Int M^2 dz   /   Int M dz
      fields_2d(ij,i_m_ref,i_type,i_layr)                                      &
        = fields_2d(ij,i_m_ref,i_type,i_layr)                                  &
        / max( fields_2d(ij,i_int_m_dz,i_type,i_layr), min_float )

      ! Normalise the integral of mass-flux * buoyancy to finish the
      ! calculation of mass-flux-weighted CAPE
      ! MFW_CAPE   =   Int M b dz   /   M_ref
      !            =   Int M b dz   *   Int M dz   /   Int M^2 dz
      ! (i.e. we have vertically integrated mass-flux times buoyancy,
      !  and we now need to divide by a representative vertical-mean mass-flux
      !  to get back in the right units for CAPE).
      fields_2d(ij,i_mfw_cape,i_type,i_layr)                                   &
        = fields_2d(ij,i_mfw_cape,i_type,i_layr)                               &
        / max( fields_2d(ij,i_m_ref,i_type,i_layr), min_float )

    end do
  end do
end do

return
end subroutine normalise_integrals

end module normalise_integrals_mod
