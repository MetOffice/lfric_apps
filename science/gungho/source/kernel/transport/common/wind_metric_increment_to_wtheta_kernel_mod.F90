!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!
! Some of the content of this file has been produced with the assistance of
! Met Office Github Copilot Enterprise
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!> @brief Adds the vertical wind metric increment from W3 to Wtheta points.
!> @details The increment at each interior Wtheta point is the average of the
!!          increments in the cells above and below. The values at the bottom
!!          and top of the column are copied unchanged.
!-----------------------------------------------------------------------------
module wind_metric_increment_to_wtheta_kernel_mod

  use argument_mod,      only: arg_type,              &
                               GH_FIELD, GH_REAL,     &
                               GH_READ, GH_WRITE,     &
                               CELL_COLUMN
  use constants_mod,     only: r_def, i_def
  use fs_continuity_mod, only: W3, Wtheta
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  type, public, extends(kernel_type) :: wind_metric_increment_to_wtheta_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                       &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, Wtheta),                       & ! w_wt_out
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                       & ! w_wt_in
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W3)                            & ! d_w_w3
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: wind_metric_increment_to_wtheta_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: wind_metric_increment_to_wtheta_code

contains

!> @brief Increment Wtheta winds with the metric increments from W3 points.
!! @param[in]     nlayers   Number of layers
!! @param[in,out] w_wt_out  Vertical wind with the increment added in Wtheta
!! @param[in]     w_wt_in   Vertical wind before the increment in Wtheta
!! @param[in]     d_w_w3    Increment to vertical component in W3
!! @param[in]     ndf_wt    Num DoFs per cell for Wtheta
!! @param[in]     undf_wt   Num DoFs in this partition for Wtheta
!! @param[in]     map_wt    Map of base cell DoFs for Wtheta
!! @param[in]     ndf_w3    Num DoFs per cell for W3
!! @param[in]     undf_w3   Num DoFs in this partition for W3
!! @param[in]     map_w3    Map of base cell DoFs for W3
subroutine wind_metric_increment_to_wtheta_code(                             &
    nlayers,                                                                 &
    w_wt_out, w_wt_in, d_w_w3,                                               &
    ndf_wt, undf_wt, map_wt,                                                 &
    ndf_w3, undf_w3, map_w3                                                  &
)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_wt, undf_wt
  integer(kind=i_def), intent(in) :: ndf_w3, undf_w3

  real(kind=r_def), intent(inout) :: w_wt_out(undf_wt)
  real(kind=r_def), intent(in)    :: w_wt_in(undf_wt)
  real(kind=r_def), intent(in)    :: d_w_w3(undf_w3)

  integer(kind=i_def), intent(in) :: map_wt(ndf_wt)
  integer(kind=i_def), intent(in) :: map_w3(ndf_w3)

  ! Internal variables
  integer(kind=i_def) :: k

  w_wt_out(map_wt(1)) = w_wt_in(map_wt(1))
  do k = 1, nlayers - 1
    w_wt_out(map_wt(1)+k) = w_wt_in(map_wt(1)+k)                             &
                          + 0.5_r_def*( d_w_w3(map_w3(1)+k-1)                &
                                      + d_w_w3(map_w3(1)+k) )
  end do
  w_wt_out(map_wt(1)+nlayers) = w_wt_in(map_wt(1)+nlayers)

end subroutine wind_metric_increment_to_wtheta_code

end module wind_metric_increment_to_wtheta_kernel_mod
