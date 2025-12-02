!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Computes a correction to horizontal fluxes for consistent transport
!> @details The shifted horizontal fluxes used in consistent conservative
!!          transport equation need an adjustment to correctly capture the
!!          dispersion relation. This kernel computes them.
!!          Only implemented for the lowest-order elements.

module consistent_dispersion_kernel_mod
use argument_mod,            only : arg_type, GH_FIELD,        &
                                    GH_READ, GH_WRITE,         &
                                    GH_INTEGER, GH_REAL,       &
                                    ANY_DISCONTINUOUS_SPACE_2, &
                                    ANY_DISCONTINUOUS_SPACE_3, &
                                    CELL_COLUMN
use fs_continuity_mod,       only : W2
use constants_mod,           only : r_tran, r_def, i_def
use kernel_mod,              only : kernel_type
use reference_element_mod,   only : N

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: consistent_dispersion_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                                    &
       ! NB: This is to be used to write to a continuous W2 field, but using
       ! a discontinuous data pattern, so use discontinuous metadata
       arg_type(GH_FIELD, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), &
       arg_type(GH_FIELD, GH_REAL,    GH_READ,  W2),                        &
       arg_type(GH_FIELD, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), &
       arg_type(GH_FIELD, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), &
       arg_type(GH_FIELD, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3)  &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: consistent_dispersion_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public consistent_dispersion_code

contains

!> @brief Computes a correction to horizontal fluxes for consistent transport
!> @param[in]     nlayers_shifted    Number of layers in the shifted mesh
!> @param[in,out] flux_X_correction  The flux correction to be computed in
!!                                   W2 on the shifted mesh
!> @param[in]     dry_flux_prime     The dry flux in W2 on the prime mesh
!> @param[in]     theta_factor       0.25*dtheta_dz*dz in shifted W2
!> @param[in]     face_selector_ew   2D field indicating which W/E faces
!!                                   to loop over in this column
!> @param[in]     face_selector_ns   2D field indicating which N/S faces
!!                                   to loop over in this column
!> @param[in]     ndf_w2_shifted     Num of DoFs per cell for shifted W2
!> @param[in]     undf_w2_shifted    Num of DoFs per partition for shifted W2
!> @param[in]     map_w2_shifted     Base cell DoF-map for shifted W2
!> @param[in]     ndf_w2_prime       Num of DoFs per cell for prime W2
!> @param[in]     undf_w2_prime      Num of DoFs per partition for prime W2
!> @param[in]     map_w2_prime       Base cell DoF-map for prime W2
!> @param[in]     ndf_w3_2d          Num of DoFs for 2D W3 per cell
!> @param[in]     undf_w3_2d         Num of DoFs for this partition for 2D W3
!> @param[in]     map_w3_2d          Map for 2D W3
subroutine consistent_dispersion_code( nlayers_shifted,    &
                                       flux_X_correction,  &
                                       dry_flux_prime,     &
                                       theta_factor,       &
                                       face_selector_ew,   &
                                       face_selector_ns,   &
                                       ndf_w2_shifted,     &
                                       undf_w2_shifted,    &
                                       map_w2_shifted,     &
                                       ndf_w2_prime,       &
                                       undf_w2_prime,      &
                                       map_w2_prime,       &
                                       ndf_w3_2d,          &
                                       undf_w3_2d,         &
                                       map_w3_2d )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers_shifted
  integer(kind=i_def), intent(in)    :: undf_w2_prime, ndf_w2_prime
  integer(kind=i_def), intent(in)    :: undf_w2_shifted, ndf_w2_shifted
  integer(kind=i_def), intent(in)    :: undf_w3_2d, ndf_w3_2d
  integer(kind=i_def), intent(in)    :: map_w2_prime(ndf_w2_prime)
  integer(kind=i_def), intent(in)    :: map_w2_shifted(ndf_w2_shifted)
  integer(kind=i_def), intent(in)    :: map_w3_2d(ndf_w3_2d)
  real(kind=r_tran),   intent(inout) :: flux_X_correction(undf_w2_shifted)
  real(kind=r_tran),   intent(in)    :: dry_flux_prime(undf_w2_prime)
  real(kind=r_tran),   intent(in)    :: theta_factor(undf_w2_shifted)
  integer(kind=i_def), intent(in)    :: face_selector_ew(undf_w3_2d)
  integer(kind=i_def), intent(in)    :: face_selector_ns(undf_w3_2d)

  ! Internal variables
  integer(kind=i_def) :: j, k, face

  ! -------------------------------------------------------------------------- !
  ! Average dz/4 * dtheta/dz * delta_F to each face
  ! -------------------------------------------------------------------------- !

  ! Loop over horizontal W2 DoFs
  do j = 1, face_selector_ew(map_w3_2d(1)) + face_selector_ns(map_w3_2d(1))
    face = j
    if (j == 3 .and. face_selector_ns(map_w3_2d(1)) == 2                       &
               .and. face_selector_ew(map_w3_2d(1)) == 1) face = N

    ! Faces in bottom layer: zero contribution
    k = 0
    flux_X_correction(map_w2_shifted(face)+k) = 0.0_r_tran

    ! No contributions to bottom or top layers, so loop over internal layers
    do k = 1, nlayers_shifted - 2
      flux_X_correction(map_w2_shifted(face)+k) =                              &
        theta_factor(map_w2_shifted(face)+k)                                   &
        * (dry_flux_prime(map_w2_prime(face)+k) - dry_flux_prime(map_w2_prime(face)+k-1))
    end do

    ! Faces in top layer: zero contribution
    k = nlayers_shifted - 1
    flux_X_correction(map_w2_shifted(face)+k) = 0.0_r_tran
  end do


end subroutine consistent_dispersion_code

end module consistent_dispersion_kernel_mod