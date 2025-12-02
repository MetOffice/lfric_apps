!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Module for vertical advection forcing of horizontal wind components
module vertadv_wind_kernel_mod

  use argument_mod,      only: arg_type,                &
                               GH_FIELD, GH_SCALAR,     &
                               GH_REAL,                 &
                               GH_READWRITE, GH_READ,   &
                               CELL_COLUMN

  use constants_mod,     only: r_def, i_def, l_def
  use fs_continuity_mod, only: Wtheta, W3
  use kernel_mod,        only: kernel_type

  implicit none

  private

!------------------------------------------------------------------------------
! Public types
!------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
  type, public, extends(kernel_type) :: vertadv_wind_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                     &
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3), &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),     &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta), &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),     &
         arg_type(GH_SCALAR, GH_REAL, GH_READ)              &
         /)
    integer :: operates_on = CELL_COLUMN
contains
    procedure, nopass :: vertadv_wind_code
end type vertadv_wind_kernel_type

!------------------------------------------------------------------------------
! Contained functions/subroutines
!------------------------------------------------------------------------------
public :: vertadv_wind_code

contains
  !----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !----------------------------------------------------------------------------
  !> @brief Calculate vertical advection increments for horizontal wind components
  !> @param[in] nlayers  Number of layers
  !> @param[in,out] field_increment Vertical advection increment
  !> @param[in] field    Field to which vertical advection is applied
  !> @param[in] w_thlev  Vertical wind for vertical advection on theta-levs
  !> @param[in] height   Height coordinate of cell interfaces
  !> @param[in] dt       Model timestep length
  !> @param[in] ndf_w3   Number of degrees of freedom per cell
  !> @param[in] undf_w3  Number of unique degrees of freedom
  !> @param[in] map_w3   Dofmap for the cell at the base of the column
  !> @param[in] ndf_wth  Number of degrees of freedom per cell
  !> @param[in] undf_wth Number of unique degrees of freedom
  !> @param[in] map_wth  Dofmap for the cell at the base of the column
  subroutine vertadv_wind_code( nlayers, field_increment, field,   &
                                w_thlev, height, dt,               &
                                ndf_w3, undf_w3, map_w3,           &
                                ndf_wth, undf_wth, map_wth )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)                    :: nlayers
    integer(kind=i_def), intent(in)                    :: ndf_w3, undf_w3
    integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
    integer(kind=i_def), intent(in)                    :: ndf_wth, undf_wth
    integer(kind=i_def), dimension(ndf_wth),intent(in) :: map_wth

    real(kind=r_def), dimension(undf_w3), intent(inout) :: field_increment
    real(kind=r_def), dimension(undf_w3), intent(in)    :: field
    real(kind=r_def), dimension(undf_wth),intent(in)    :: w_thlev
    real(kind=r_def), dimension(undf_w3), intent(in)    :: height
    real(kind=r_def),                     intent(in)    :: dt

    ! Loop counter over output profile
    integer(kind=i_def) :: k

    ! map(1) + k
    integer(kind=i_def) :: k3, kth

    do k = 0, nlayers - 2
      k3 = map_w3(1) + k
      kth = map_wth(1) + k
      ! use first-order upwind advection
      if (w_thlev(kth+1) < 0.0_r_def) then
        ! subsidence
        field_increment(k3) = - dt * w_thlev(kth+1) *                         &
               ( field(k3+1) - field(k3) ) / ( height(k3+1) - height(k3) )
      else
        ! ascent
        field_increment(k3+1) = - dt * w_thlev(kth+1) *                       &
               ( field(k3+1) - field(k3) ) / ( height(k3+1) - height(k3) )
      end if

    end do

    ! Top level: just copy increment from level below,
    ! to maintain profile gradient
    k3 = map_w3(1) + nlayers - 1
    field_increment(k3) = field_increment(k3-1)

  end subroutine vertadv_wind_code

end module vertadv_wind_kernel_mod

