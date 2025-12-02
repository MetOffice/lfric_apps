!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Module for calculating idealised vertical advection increments
module vertadvect_kernel_mod

  use argument_mod,      only: arg_type,            &
                               GH_FIELD, GH_SCALAR, &
                               GH_REAL,             &
                               GH_WRITE, GH_READ,   &
                               CELL_COLUMN

  use constants_mod,     only: r_def, i_def, l_def
  use fs_continuity_mod, only: Wtheta, W3
  use kernel_mod,        only: kernel_type
  use log_mod,           only: log_event, log_scratch_space,                 &
                               LOG_LEVEL_INFO, LOG_LEVEL_ERROR

  implicit none

  private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
  type, public, extends(kernel_type) :: vertadvect_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                     &
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, Wtheta), &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta), &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),     &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta), &
         arg_type(GH_SCALAR, GH_REAL,    GH_READ)           &
         /)
    integer :: operates_on = CELL_COLUMN
contains
    procedure, nopass :: vertadvect_code
end type vertadvect_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: vertadvect_code

contains
  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  !> @brief Interpolate profile data onto a field
  !> @param[in] nlayers Number of layers
  !> @param[in,out] field_increment vertical advection increment
  !> @param[in] field Field to which vertical advection is applied
  !> @param[in] w_rholev forcing vertical velocity (on w3 (rho)-levels)
  !> @param[in] height Height coordinate of cell centres
  !> @param[in] dt     Model timestep length
  !> @param[in] ndf_wth Number of degrees of freedom per cell
  !> @param[in] undf_wth Number of unique degrees of freedom
  !> @param[in] map_wth Dofmap for the cell at the base of the column
  !> @param[in] ndf_w3 Number of degrees of freedom per cell
  !> @param[in] undf_w3 Number of unique degrees of freedom
  !> @param[in] map_w3 Dofmap for the cell at the base of the column
  subroutine vertadvect_code( nlayers, field_increment, field, &
                              w_rholev, height, dt,            &
                              ndf_wth, undf_wth, map_wth,      &
                              ndf_w3, undf_w3, map_w3 )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)                     :: nlayers
    integer(kind=i_def), intent(in)                     :: ndf_wth, undf_wth
    integer(kind=i_def), dimension(ndf_wth), intent(in) :: map_wth
    integer(kind=i_def), intent(in)                     :: ndf_w3, undf_w3
    integer(kind=i_def), dimension(ndf_w3), intent(in)  :: map_w3

    real(kind=r_def), dimension(undf_wth), intent(inout):: field_increment
    real(kind=r_def), dimension(undf_wth), intent(in)   :: field
    real(kind=r_def), dimension(undf_w3),  intent(in)   :: w_rholev
    real(kind=r_def), dimension(undf_wth), intent(in)   :: height
    real(kind=r_def),                      intent(in)   :: dt

    ! Loop counter over output profile
    integer(kind=i_def) :: k

    ! map(1) + k
    integer(kind=i_def) :: kth, k3

    do k = 0, nlayers-1
      kth = map_wth(1) + k
      k3  = map_w3(1)  + k
      field_increment(kth) = 0.0_r_def

      ! use first-order upwind advection
      if (w_rholev(k3) < 0.0_r_def) then
        ! subsidence
        field_increment(kth) = - dt * ( field(kth+1) - field(kth)  ) * w_rholev(k3) /   &
                                      ( height(kth+1) - height(kth) )
      else
        ! ascent
        field_increment(kth+1) = - dt * ( field(kth+1) - field(kth)  ) * w_rholev(k3) / &
                                        ( height(kth+1) - height(kth) )
      end if

    end do

    ! Top level: just copy increment from level below, to maintain profile gradient
    kth = map_wth(1)+nlayers
    field_increment(kth) = field_increment(kth-1)

  end subroutine vertadvect_code

end module vertadvect_kernel_mod
