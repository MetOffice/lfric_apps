!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!
! Some of the content of this file has been produced with the assistance of
! Met Office Github Copilot Enterprise
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------!
!> @brief Kernel for scaling a wind increment by the area element on W2.
!> @details Multiplies the W2 wind increment by the area element and applies
!!          the correction to the vertical component due to the slope of the
!!          model levels.
!-----------------------------------------------------------------------------!
module apply_w2_area_and_slope_kernel_mod

  use argument_mod,          only: arg_type, func_type,   &
                                   GH_FIELD, GH_REAL,     &
                                   GH_READ, GH_INC,       &
                                   CELL_COLUMN,           &
                                   GH_DIFF_BASIS,         &
                                   GH_EVALUATOR
  use constants_mod,         only: r_def, i_def
  use fs_continuity_mod,     only: W2, Wchi
  use kernel_mod,            only: kernel_type
  use reference_element_mod, only: B, T, N, E, S, W

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  type, public, extends(kernel_type) :: apply_w2_area_and_slope_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                       &
         arg_type(GH_FIELD, GH_REAL, GH_INC,  W2),                            &
         arg_type(GH_FIELD, GH_REAL, GH_READ, W2),                            &
         arg_type(GH_FIELD, GH_REAL, GH_READ, W2),                            &
         arg_type(GH_FIELD, GH_REAL, GH_READ, Wchi)                           &
         /)
    type(func_type) :: meta_funcs(1) = (/                                     &
         func_type(Wchi, GH_DIFF_BASIS)                                       &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: apply_w2_area_and_slope_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: apply_w2_area_and_slope_code

contains

!> @brief Scale a W2 wind by the area element and correct for the slope.
!! @param[in]     nlayers        Number of layers
!! @param[in,out] scaled_uvw     Output wind scaled by the area element in W2
!! @param[in]     uvw            Input wind in W2
!! @param[in]     da             Area element in W2
!! @param[in]     chi3           Third coordinate field in Wchi
!! @param[in]     ndf_w2         Num DoFs per cell for W2
!! @param[in]     undf_w2        Num DoFs in this partition for W2
!! @param[in]     map_w2         Map of base cell DoFs for W2
!! @param[in]     ndf_wx         Num DoFs per cell for Wchi
!! @param[in]     undf_wx        Num DoFs in this partition for Wchi
!! @param[in]     map_wx         Map of base cell DoFs for Wchi
!! @param[in]     diff_basis_wx  Differential basis functions for Wchi
!!                               evaluated at W2 nodes
subroutine apply_w2_area_and_slope_code(                                     &
    nlayers,                                                                 &
    scaled_uvw, uvw, da, chi3,                                               &
    ndf_w2, undf_w2, map_w2,                                                 &
    ndf_wx, undf_wx, map_wx, diff_basis_wx                                   &
)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
  integer(kind=i_def), intent(in) :: ndf_wx, undf_wx

  real(kind=r_def), intent(inout) :: scaled_uvw(undf_w2)
  real(kind=r_def), intent(in)    :: uvw(undf_w2)
  real(kind=r_def), intent(in)    :: da(undf_w2)
  real(kind=r_def), intent(in)    :: chi3(undf_wx)

  integer(kind=i_def), intent(in) :: map_w2(ndf_w2)
  integer(kind=i_def), intent(in) :: map_wx(ndf_wx)

  real(kind=r_def), intent(in) :: diff_basis_wx(3,ndf_wx,ndf_w2)

  ! Internal variables
  integer(kind=i_def) :: k, df, nm1
  real(kind=r_def)    :: u_av, v_av
  real(kind=r_def)    :: dx_z, dy_z, dz_z

  nm1 = nlayers - 1

  do df = 1, ndf_w2 - 1
    scaled_uvw(map_w2(df):map_w2(df)+nm1) = uvw(map_w2(df):map_w2(df)+nm1)   &
                                           *da(map_w2(df):map_w2(df)+nm1)
  end do

  ! Correction to vertical wind based upon model slope
  scaled_uvw(map_w2(B)) = 0.0_r_def
  do k = 1, nlayers - 1
    dx_z = 0.0_r_def ! delta_x(z) on W point
    dy_z = 0.0_r_def ! delta_y(z) on W point
    dz_z = 0.0_r_def ! delta_z(z) on W point
    do df = 1, ndf_wx
      dx_z = dx_z + chi3(map_wx(df)+k)*diff_basis_wx(1, df, B)
      dy_z = dy_z + chi3(map_wx(df)+k)*diff_basis_wx(2, df, B)
      dz_z = dz_z + chi3(map_wx(df)+k)*diff_basis_wx(3, df, B)
    end do

    u_av = 0.25_r_def*( scaled_uvw(map_w2(W)+k)   + scaled_uvw(map_w2(E)+k)  &
                      + scaled_uvw(map_w2(W)+k-1) + scaled_uvw(map_w2(E)+k-1) )
    v_av = 0.25_r_def*( scaled_uvw(map_w2(S)+k)   + scaled_uvw(map_w2(N)+k)  &
                      + scaled_uvw(map_w2(S)+k-1) + scaled_uvw(map_w2(N)+k-1) )
    scaled_uvw(map_w2(B)+k) = scaled_uvw(map_w2(B)+k)                        &
                            - u_av*dx_z/dz_z + v_av*dy_z/dz_z
  end do
  scaled_uvw(map_w2(T)+nlayers-1) = 0.0_r_def

end subroutine apply_w2_area_and_slope_code

end module apply_w2_area_and_slope_kernel_mod
