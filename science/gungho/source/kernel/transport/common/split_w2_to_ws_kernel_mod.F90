!-------------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------!
!> @brief Kernel to split W2 wind components onto native horizontal/vertical spaces.
!> @details Splits W2 winds into W2h and Wtheta components, removes the area factor
!!          (dA) from the dof values and applies slope correction terms for the
!!          vertical component.
!-----------------------------------------------------------------------------!
module split_w2_to_ws_kernel_mod

  use argument_mod,          only: arg_type, func_type,       &
                                   GH_FIELD, GH_REAL,         &
                                   GH_INTEGER,                &
                                   GH_READ, GH_WRITE,         &
                                   CELL_COLUMN,               &
                                   ANY_SPACE_9,               &
                                   ANY_DISCONTINUOUS_SPACE_1, &
                                   ANY_DISCONTINUOUS_SPACE_2, &
                                   GH_DIFF_BASIS, GH_EVALUATOR
  use constants_mod,         only: r_def, i_def
  use fs_continuity_mod,     only: W2, W2h, Wtheta
  use kernel_mod,            only: kernel_type
  use reference_element_mod, only: W, E, N, S

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  type, public, extends(kernel_type) :: split_w2_to_ws_kernel_type
    private
    type(arg_type) :: meta_args(7) = (/                                       &
         arg_type(GH_FIELD, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! uv
         arg_type(GH_FIELD, GH_REAL,    GH_WRITE, Wtheta),                    & ! w
         arg_type(GH_FIELD, GH_REAL,    GH_READ,  W2),                        & ! uvw
         arg_type(GH_FIELD, GH_REAL,    GH_READ,  ANY_SPACE_9),               & ! chi3
         arg_type(GH_FIELD, GH_REAL,    GH_READ,  W2),                        & ! dA
         arg_type(GH_FIELD, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! face selector ew
         arg_type(GH_FIELD, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_2)  & ! face selector ns
         /)
    type(func_type) :: meta_funcs(1) = (/                                     &
         func_type(ANY_SPACE_9, GH_DIFF_BASIS)                                &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
    integer :: gh_evaluator_targets(1) = (/ Wtheta /)
  contains
    procedure, nopass :: split_w2_to_ws_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: split_w2_to_ws_code

contains

!> @brief Split W2 winds into W2h and Wtheta component fields.
!! @param[in]     nlayers             Number of layers.
!! @param[in,out] uv_w2h              Horizontal wind component in W2h.
!! @param[in,out] w_wt                Vertical wind component in Wtheta.
!! @param[in]     uvw                 Vector wind component in W2.
!! @param[in]     chi3                Third coordinate field in Wchi.
!! @param[in]     da                  Area element in W2.
!! @param[in]     face_selector_ew    Number of east-west face DoFs used per cell.
!! @param[in]     face_selector_ns    Number of north-south face DoFs used per cell.
!! @param[in]     ndf_w2h             Number of degrees of freedom per cell for W2h.
!! @param[in]     undf_w2h            Number of unique degrees of freedom for W2h.
!! @param[in]     map_w2h             DoF map for W2h.
!! @param[in]     ndf_wt              Number of degrees of freedom per cell for Wtheta.
!! @param[in]     undf_wt             Number of unique degrees of freedom for Wtheta.
!! @param[in]     map_wt              DoF map for Wtheta.
!! @param[in]     ndf_w2              Number of degrees of freedom per cell for W2.
!! @param[in]     undf_w2             Number of unique degrees of freedom for W2.
!! @param[in]     map_w2              DoF map for W2.
!! @param[in]     ndf_wx              Number of degrees of freedom per cell for Wchi.
!! @param[in]     undf_wx             Number of unique degrees of freedom for Wchi.
!! @param[in]     map_wx              DoF map for Wchi.
!! @param[in]     diff_basis_wx_on_wt Differential basis functions for Wchi
!!                                    evaluated at Wtheta nodes.
!! @param[in]     ndf_w3_2d           Number of degrees of freedom per cell for 2D W3.
!! @param[in]     undf_w3_2d          Number of unique degrees of freedom for 2D W3.
!! @param[in]     map_w3_2d           DoF map for 2D W3.
subroutine split_w2_to_ws_code(nlayers,                    &
                               uv_w2h, w_wt, uvw,          &
                               chi3, da,                   &
                               face_selector_ew,           &
                               face_selector_ns,           &
                               ndf_w2h, undf_w2h, map_w2h, &
                               ndf_wt, undf_wt, map_wt,    &
                               ndf_w2, undf_w2, map_w2,    &
                               ndf_wx, undf_wx, map_wx,    &
                               diff_basis_wx_on_wt,        &
                               ndf_w3_2d, undf_w3_2d,      &
                               map_w3_2d                   &
                             )

  use sci_coordinate_jacobian_mod, only : pointwise_coordinate_jacobian

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w2h, undf_w2h
  integer(kind=i_def), intent(in) :: ndf_wt, undf_wt
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
  integer(kind=i_def), intent(in) :: ndf_wx, undf_wx
  integer(kind=i_def), intent(in) :: ndf_w3_2d, undf_w3_2d

  integer(kind=i_def), dimension(ndf_w2h),   intent(in) :: map_w2h
  integer(kind=i_def), dimension(ndf_wt),    intent(in) :: map_wt
  integer(kind=i_def), dimension(ndf_w2),    intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_wx),    intent(in) :: map_wx
  integer(kind=i_def), dimension(ndf_w3_2d), intent(in) :: map_w3_2d

  real(kind=r_def),    dimension(undf_w2h),   intent(inout) :: uv_w2h
  real(kind=r_def),    dimension(undf_wt),    intent(inout) :: w_wt
  real(kind=r_def),    dimension(undf_w2),    intent(in)    :: uvw
  real(kind=r_def),    dimension(undf_wx),    intent(in)    :: chi3
  real(kind=r_def),    dimension(undf_w2),    intent(in)    :: da
  integer(kind=i_def), dimension(undf_w3_2d), intent(in)    :: face_selector_ew
  integer(kind=i_def), dimension(undf_w3_2d), intent(in)    :: face_selector_ns

  real(kind=r_def),    dimension(3,ndf_wx,ndf_wt), intent(in) :: diff_basis_wx_on_wt

  ! Internal variables
  integer(kind=i_def) :: k, df, id1, id2
  integer(kind=i_def) :: local_dofs_x(2), local_dofs_y(2)

  real(kind=r_def) :: dx_z, dy_z, dz_z

  local_dofs_x = (/ W, E /)
  local_dofs_y = (/ S, N /)

  do df = 1, face_selector_ew(map_w3_2d(1))
    id1 = map_w2h(local_dofs_x(df))
    id2 = map_w2(local_dofs_x(df))
    uv_w2h(id1:id1+nlayers-1) = uvw(id2:id2+nlayers-1) &
                              /da(id2:id2+nlayers-1)
  end do
  do df = 1, face_selector_ns(map_w3_2d(1))
    id1 = map_w2h(local_dofs_y(df))
    id2 = map_w2(local_dofs_y(df))
    uv_w2h(id1:id1+nlayers-1) = uvw(id2:id2+nlayers-1) &
                              /da(id2:id2+nlayers-1)
  end do
  w_wt(map_wt(1)+1:map_wt(1)+nlayers-1) = uvw(map_w2(5)+1:map_w2(5)+nlayers-1) &
                                         /da(map_w2(5)+1:map_w2(5)+nlayers-1)
  w_wt(map_wt(1)) = 0.0_r_def
  w_wt(map_wt(1)+nlayers) = 0.0_r_def

  ! Correction to vertical wind based upon model slope
  do k = 0, nlayers - 1
    dx_z = 0.0_r_def ! delta_x(z) on W point
    dy_z = 0.0_r_def ! delta_y(z) on W point
    dz_z = 0.0_r_def ! delta_z(z) on W point
    do df = 1, ndf_wx
      dx_z = dx_z + chi3(map_wx(df)+k)*diff_basis_wx_on_wt(1, df, 1)
      dy_z = dy_z + chi3(map_wx(df)+k)*diff_basis_wx_on_wt(2, df, 1)
      dz_z = dz_z + chi3(map_wx(df)+k)*diff_basis_wx_on_wt(3, df, 1)
    end do
    w_wt(map_wt(1)+k) = w_wt(map_wt(1)+k) + 0.5_r_def/da(map_w2(5)+k)    &
                       *(dx_z/dz_z*(uvw(map_w2(1)+k) + uvw(map_w2(3)+k)) &
                       - dy_z/dz_z*(uvw(map_w2(2)+k) + uvw(map_w2(4)+k)) &
                        )
  end do

end subroutine split_w2_to_ws_code

end module split_w2_to_ws_kernel_mod
