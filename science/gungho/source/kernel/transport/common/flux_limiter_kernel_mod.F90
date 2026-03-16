!-------------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------

!> @brief Limits a mass flux to ensure that negative values are not produced.
!> @details This kernel limits a W2 flux to ensure that a transported W3 field
!!          will not produce negative values.
!!          All outgoing fluxes F from a cell are scaled by the same factor a,
!!          without considering the incoming fluxes, to ensure that a field f
!!          is updated via:
!!          f(n+1) = f(n) - div(a*F_out) >= 0
!!          This only gives sensible results for Courant numbers less than 1.
!!          The kernel is designed for the lowest-order finite element spaces.

module flux_limiter_kernel_mod
use argument_mod,           only: arg_type,                                    &
                                  GH_FIELD, GH_REAL, GH_READ, GH_READINC,      &
                                  CELL_COLUMN
use fs_continuity_mod,      only: W3, W2
use constants_mod,          only: r_def, i_def, EPS
use kernel_mod,             only: kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: flux_limiter_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                          &
      arg_type(GH_FIELD, GH_REAL, GH_READINC, W2),                             &
      arg_type(GH_FIELD, GH_REAL, GH_READ,    W3),                             &
      arg_type(GH_FIELD, GH_REAL, GH_READ,    W3)                              &
  /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: flux_limiter_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public flux_limiter_code

contains

!> @brief Limits a mass flux to ensure that negative values are not produced.
!> @param[in]     nlayers        Number of layers
!> @param[in,out] flux           Flux to be scaled at W2 locations
!> @param[in]     field          Field at W3 locations that is being transported
!> @param[in]     volume         Volume at W3 locations
!> @param[in]     ndf_w2         Num DoFs per cell for W2
!> @param[in]     undf_w2        Num DoFs in this partition for W2
!> @param[in]     map_w2         Index of W2 DoFs for lowest cell in column
!> @param[in]     ndf_w3         Num DoFs per cell for W3
!> @param[in]     undf_w3        Num DoFs in this partition for W3
!> @param[in]     map_w3         Index of W3 DoFs for lowest cell in column
subroutine flux_limiter_code(                                                  &
    nlayers,                                                                   &
    flux,                                                                      &
    field,                                                                     &
    volume,                                                                    &
    ndf_w2, undf_w2, map_w2,                                                   &
    ndf_w3, undf_w3, map_w3                                                    &
)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_w2, undf_w2, ndf_w3, undf_w3
  integer(kind=i_def), intent(in)    :: map_w2(ndf_w2)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
  real(kind=r_def),    intent(inout) :: flux(undf_w2)
  real(kind=r_def),    intent(in)    :: field(undf_w3)
  real(kind=r_def),    intent(in)    :: volume(undf_w3)

  ! Internal variables
  integer(kind=i_def) :: df
  integer(kind=i_def) :: w2_b_idx, w2_t_idx, w3_b_idx, w3_t_idx
  real(kind=r_def)    :: total_outgoing(nlayers)
  real(kind=r_def)    :: a_scaling(nlayers)
  real(kind=r_def)    :: outgoing_switch(nlayers)
  real(kind=r_def)    :: outgoing_scaling(nlayers)

  ! Gives the sign associated with the divergence operator for each W2 DoF
  real(kind=r_def), parameter :: sign_div(6) = [                               &
      -1.0_r_def, 1.0_r_def, 1.0_r_def, -1.0_r_def, -1.0_r_def, 1.0_r_def      &
  ]

  w3_b_idx = map_w3(1)
  w3_t_idx = w3_b_idx + nlayers - 1

  ! Find the total outgoing flux for each cell
  total_outgoing(:) = 0.0
  do df = 1, ndf_w2
    w2_b_idx = map_w2(df)
    w2_t_idx = w2_b_idx + nlayers - 1
    total_outgoing(:) = total_outgoing(:)                                      &
      + MAX(0.0_r_def, sign_div(df)*flux(w2_b_idx:w2_t_idx))
  end do

  ! Determine scaling factor for each cell
  ! 0 = f(n) - a*total_outgoing/volume, so
  ! a = f(n)*volume/total_outgoing
  ! This needs an upper bound of 1, as fluxes should never be amplified
  a_scaling(:) = MIN(1.0_r_def,                                                &
    field(w3_b_idx:w3_t_idx) * volume(w3_b_idx:w3_t_idx)                       &
    / MAX(total_outgoing(:), EPS)                                              &
  )

  ! Scale outgoing fluxes
  do df = 1, ndf_w2
    w2_b_idx = map_w2(df)
    w2_t_idx = w2_b_idx + nlayers - 1
    ! Define a variable that is 1 for outgoing fluxes and 0 for incoming fluxes
    outgoing_switch(:) = 0.5_r_def * (                                         &
        1.0_r_def - SIGN(1.0_r_def, sign_div(df)*flux(w2_b_idx:w2_t_idx))      &
    )
    ! The scaling should be 1 for incoming fluxes and a_scaling for outgoing
    outgoing_scaling(:) = (                                                    &
        1.0_r_def - outgoing_switch(:) + a_scaling(:)*outgoing_switch(:)       &
    )

    flux(w2_b_idx:w2_t_idx) = flux(w2_b_idx:w2_t_idx) * outgoing_scaling(:)
  end do

end subroutine flux_limiter_code

end module flux_limiter_kernel_mod
