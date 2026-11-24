!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!
! Some of the content of this file has been produced with the assistance of
! Met Office Github Copilot Enterprise
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!> @brief Adds wind metric increments from W3 points to W2H points.
!> @details The increments to the two horizontal components at W3 points are
!!          scaled by the reciprocal multiplicity and incremented onto the
!!          W2H points of the cell, so that W2H points shared by two cells
!!          receive the average of the two. The sign of the xi2 component is
!!          flipped, as S and N face DoFs are positive along minus xi2.
!-----------------------------------------------------------------------------
module wind_metric_increment_to_w2h_kernel_mod

  use argument_mod,      only: arg_type,              &
                               GH_FIELD, GH_REAL,     &
                               GH_READ, GH_INC,       &
                               CELL_COLUMN
  use constants_mod,     only: r_def, i_def
  use fs_continuity_mod, only: W3, W2h
  use kernel_mod,        only: kernel_type
  use reference_element_mod, &
                         only: W, S, E, N

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  type, public, extends(kernel_type) :: wind_metric_increment_to_w2h_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                       &
         arg_type(GH_FIELD, GH_REAL, GH_INC,  W2h),                           & ! uv_w2h
         arg_type(GH_FIELD, GH_REAL, GH_READ, W3),                            & ! d_uv1_w3
         arg_type(GH_FIELD, GH_REAL, GH_READ, W3),                            & ! d_uv2_w3
         arg_type(GH_FIELD, GH_REAL, GH_READ, W2h)                            & ! rmultiplicity
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: wind_metric_increment_to_w2h_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: wind_metric_increment_to_w2h_code

contains

!> @brief Increment W2H winds with the metric increments from W3 points.
!! @param[in]     nlayers        Number of layers
!! @param[in,out] uv_w2h         Horizontal wind in W2H
!! @param[in]     d_uv1_w3       Increment to component along xi1 in W3
!! @param[in]     d_uv2_w3       Increment to component along xi2 in W3
!! @param[in]     rmultiplicity  Reciprocal of the multiplicity of W2H DoFs
!! @param[in]     ndf_w2h        Num DoFs per cell for W2H
!! @param[in]     undf_w2h       Num DoFs in this partition for W2H
!! @param[in]     map_w2h        Map of base cell DoFs for W2H
!! @param[in]     ndf_w3         Num DoFs per cell for W3
!! @param[in]     undf_w3        Num DoFs in this partition for W3
!! @param[in]     map_w3         Map of base cell DoFs for W3
subroutine wind_metric_increment_to_w2h_code(                                &
    nlayers,                                                                 &
    uv_w2h, d_uv1_w3, d_uv2_w3, rmultiplicity,                               &
    ndf_w2h, undf_w2h, map_w2h,                                              &
    ndf_w3, undf_w3, map_w3                                                  &
)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w2h, undf_w2h
  integer(kind=i_def), intent(in) :: ndf_w3, undf_w3

  real(kind=r_def), intent(inout) :: uv_w2h(undf_w2h)
  real(kind=r_def), intent(in)    :: d_uv1_w3(undf_w3)
  real(kind=r_def), intent(in)    :: d_uv2_w3(undf_w3)
  real(kind=r_def), intent(in)    :: rmultiplicity(undf_w2h)

  integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)
  integer(kind=i_def), intent(in) :: map_w3(ndf_w3)

  ! Internal variables
  integer(kind=i_def) :: k

  do k = 0, nlayers - 1
    uv_w2h(map_w2h(W)+k) = uv_w2h(map_w2h(W)+k)                              &
                         + rmultiplicity(map_w2h(W)+k)*d_uv1_w3(map_w3(1)+k)
    uv_w2h(map_w2h(E)+k) = uv_w2h(map_w2h(E)+k)                              &
                         + rmultiplicity(map_w2h(E)+k)*d_uv1_w3(map_w3(1)+k)
    uv_w2h(map_w2h(S)+k) = uv_w2h(map_w2h(S)+k)                              &
                         - rmultiplicity(map_w2h(S)+k)*d_uv2_w3(map_w3(1)+k)
    uv_w2h(map_w2h(N)+k) = uv_w2h(map_w2h(N)+k)                              &
                         - rmultiplicity(map_w2h(N)+k)*d_uv2_w3(map_w3(1)+k)
  end do

end subroutine wind_metric_increment_to_w2h_code

end module wind_metric_increment_to_w2h_kernel_mod
