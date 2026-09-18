!-------------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------

!> @brief Sets a tracer flux based on dry fluxes and the upwind tracer values
!> @details The tracer flux F_X at W2 points is set to be F_d * m_X, where
!!          F_d is the dry flux at W2 points and m_X is the upwind tracer value
!!          at W2 points.
!!          The kernel is designed for the lowest-order finite element spaces.

module upwind_tracer_flux_kernel_mod
use argument_mod,           only: arg_type,                                    &
                                  GH_FIELD, GH_REAL, GH_READ, GH_INC,          &
                                  CELL_COLUMN
use fs_continuity_mod,       only : W3, W2
use constants_mod,           only : r_def, i_def
use kernel_mod,              only : kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: upwind_tracer_flux_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                          &
      arg_type(GH_FIELD, GH_REAL, GH_INC,  W2),                                &
      arg_type(GH_FIELD, GH_REAL, GH_READ, W2),                                &
      arg_type(GH_FIELD, GH_REAL, GH_READ, W3)                                 &
  /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: upwind_tracer_flux_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public upwind_tracer_flux_code

contains

!> @brief Sets a tracer flux based on dry fluxes and the upwind tracer values
!> @param[in]     nlayers        Number of layers
!> @param[in,out] tracer_flux    Tracer flux to be computed
!> @param[in]     dry_flux       Dry flux at W2 locations
!> @param[in]     tracer_field   Tracer field to be transported
!> @param[in]     ndf_w2         Num DoFs per cell for W2
!> @param[in]     undf_w2        Num DoFs in this partition for W2
!> @param[in]     map_w2         Index of W2 DoFs for lowest cell in column
!> @param[in]     ndf_w3         Num DoFs per cell for W3
!> @param[in]     undf_w3        Num DoFs in this partition for W3
!> @param[in]     map_w3         Index of W3 DoFs for lowest cell in column
subroutine upwind_tracer_flux_code(                                            &
    nlayers,                                                                   &
    tracer_flux,                                                               &
    dry_flux,                                                                  &
    tracer_field,                                                              &
    ndf_w2, undf_w2, map_w2,                                                   &
    ndf_w3, undf_w3, map_w3                                                    &
)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_w2, undf_w2, ndf_w3, undf_w3
  integer(kind=i_def), intent(in)    :: map_w2(ndf_w2)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
  real(kind=r_def),    intent(inout) :: tracer_flux(undf_w2)
  real(kind=r_def),    intent(in)    :: dry_flux(undf_w2)
  real(kind=r_def),    intent(in)    :: tracer_field(undf_w3)

  ! Internal variables
  integer(kind=i_def) :: df
  integer(kind=i_def) :: w2_b_idx, w2_t_idx, w3_b_idx, w3_t_idx
  real(kind=r_def)    :: upwind_switch(nlayers)

  ! Gives the sign associated with the wind direction for each W2 DoF
  real(kind=r_def), parameter :: sign_wind(6) = [                              &
      -1.0_r_def, 1.0_r_def, 1.0_r_def, -1.0_r_def, -1.0_r_def, 1.0_r_def      &
  ]

  w3_b_idx = map_w3(1)
  w3_t_idx = w3_b_idx + nlayers - 1

  do df = 1, ndf_w2
    w2_b_idx = map_w2(df)
    w2_t_idx = w2_b_idx + nlayers - 1

    ! Determine switch for upwind value based on the sign of the dry flux
    upwind_switch(:) = 0.5_r_def*(                                             &
      1.0_r_def + SIGN(1.0_r_def, sign_wind(df)*dry_flux(w2_b_idx:w2_t_idx))   &
    )

    ! Set the tracer flux
    tracer_flux(w2_b_idx:w2_t_idx) = (                                         &
        tracer_flux(w2_b_idx:w2_t_idx)                                         &
        + upwind_switch(:)*tracer_field(w3_b_idx:w3_t_idx)                     &
        * dry_flux(w2_b_idx:w2_t_idx)                                          &
    )
  end do

end subroutine upwind_tracer_flux_code

end module upwind_tracer_flux_kernel_mod
