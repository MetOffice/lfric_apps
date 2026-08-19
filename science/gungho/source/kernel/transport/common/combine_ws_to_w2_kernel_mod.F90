!-------------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------!
!> @brief Kernel to combine native wind components into a single W2 field.
!> @details Merges horizontal W2h and vertical Wtheta winds into a consistent
!!          W2 representation using face-selector metadata.
!-----------------------------------------------------------------------------!
module combine_ws_to_w2_kernel_mod

  use argument_mod,           only: arg_type,                  &
                                    GH_FIELD, GH_REAL,         &
                                    GH_INTEGER, GH_READ,       &
                                    GH_WRITE, CELL_COLUMN,     &
                                    ANY_DISCONTINUOUS_SPACE_1, &
                                    ANY_DISCONTINUOUS_SPACE_2
  use constants_mod,          only: r_def, i_def
  use fs_continuity_mod,      only: W2, W2h, Wtheta
  use kernel_mod,             only: kernel_type
  use reference_element_mod,  only: W, E, N, S

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  type, public, extends(kernel_type) :: combine_ws_to_w2_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                       &
         arg_type(GH_FIELD, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! uvw
         arg_type(GH_FIELD, GH_REAL,    GH_READ,  W2h),                       & ! uv
         arg_type(GH_FIELD, GH_REAL,    GH_READ,  Wtheta),                    & ! w
         arg_type(GH_FIELD, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! face selector ew
         arg_type(GH_FIELD, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_2)  & ! face selector ns
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: combine_ws_to_w2_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: combine_ws_to_w2_code

contains

!> @brief Combine W2h and Wtheta winds into a W2 wind field.
!! @param[in]     nlayers          Number of layers.
!! @param[in,out] uvw              Combined vector wind component in W2.
!! @param[in]     uv_w2h           Horizontal wind component in W2h.
!! @param[in]     w_wt             Vertical wind component in Wtheta.
!! @param[in]     face_selector_ew Number of east-west face DoFs used per cell.
!! @param[in]     face_selector_ns Number of north-south face DoFs used per cell.
!! @param[in]     ndf_w2           Number of degrees of freedom per cell for W2.
!! @param[in]     undf_w2          Number of unique degrees of freedom for W2.
!! @param[in]     map_w2           DoF map for W2.
!! @param[in]     ndf_w2h          Number of degrees of freedom per cell for W2h.
!! @param[in]     undf_w2h         Number of unique degrees of freedom for W2h.
!! @param[in]     map_w2h          DoF map for W2h.
!! @param[in]     ndf_wt           Number of degrees of freedom per cell for Wtheta.
!! @param[in]     undf_wt          Number of unique degrees of freedom for Wtheta.
!! @param[in]     map_wt           DoF map for Wtheta.
!! @param[in]     ndf_w3_2d        Number of degrees of freedom per cell for 2D W3.
!! @param[in]     undf_w3_2d       Number of unique degrees of freedom for 2D W3.
!! @param[in]     map_w3_2d        DoF map for 2D W3.
subroutine combine_ws_to_w2_code(nlayers,                    &
                                 uvw, uv_w2h, w_wt,          &
                                 face_selector_ew,           &
                                 face_selector_ns,           &
                                 ndf_w2, undf_w2, map_w2,    &
                                 ndf_w2h, undf_w2h, map_w2h, &
                                 ndf_wt, undf_wt, map_wt,    &
                                 ndf_w3_2d, undf_w3_2d,      &
                                 map_w3_2d                   &
                                )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w2h, undf_w2h
  integer(kind=i_def), intent(in) :: ndf_wt, undf_wt
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
  integer(kind=i_def), intent(in) :: ndf_w3_2d, undf_w3_2d

  integer(kind=i_def), dimension(ndf_w2h),   intent(in) :: map_w2h
  integer(kind=i_def), dimension(ndf_wt),    intent(in) :: map_wt
  integer(kind=i_def), dimension(ndf_w2),    intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_w3_2d), intent(in) :: map_w3_2d

  real(kind=r_def),    dimension(undf_w2h),   intent(in)    :: uv_w2h
  real(kind=r_def),    dimension(undf_wt),    intent(in)    :: w_wt
  real(kind=r_def),    dimension(undf_w2),    intent(inout) :: uvw
  integer(kind=i_def), dimension(undf_w3_2d), intent(in)    :: face_selector_ew
  integer(kind=i_def), dimension(undf_w3_2d), intent(in)    :: face_selector_ns

  integer(kind=i_def) :: df, id1, id2
  integer(kind=i_def) :: local_dofs_x(2), local_dofs_y(2)

  local_dofs_x = (/ W, E /)
  local_dofs_y = (/ S, N /)

  do df = 1, face_selector_ew(map_w3_2d(1))
    id1 = map_w2(local_dofs_x(df))
    id2 = map_w2h(local_dofs_x(df))
    uvw(id1:id1+nlayers-1) = uv_w2h(id2:id2+nlayers-1)
  end do
  do df = 1, face_selector_ns(map_w3_2d(1))
    id1 = map_w2(local_dofs_y(df))
    id2 = map_w2h(local_dofs_y(df))
    uvw(id1:id1+nlayers-1) = uv_w2h(id2:id2+nlayers-1)
  end do

  uvw(map_w2(5)) = 0.0_r_def
  uvw(map_w2(5)+1:map_w2(5)+nlayers-1) = w_wt(map_wt(1)+1:map_wt(1)+nlayers-1)
  uvw(map_w2(6)+nlayers-1) = 0.0_r_def

end subroutine combine_ws_to_w2_code

end module combine_ws_to_w2_kernel_mod
