!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Calculates the Eulerian horizontal departure distance simply.
!> @details This kernel computes the Eulerian horizontal departure distance
!!          using the advecting departure wind u and time step dt as
!!          \f$\mbox{dist} = x_a - x_d = u dt\f$,
!!          where x_d is the departure point and x_a the arrival point.
!!          The horizontal grid spacing is not taken into account.
!!          This kernel only works with lowest order W2h spaces.

module hori_dep_dist_eulerian_kernel_mod

  use argument_mod,      only : arg_type,              &
                                GH_FIELD, GH_REAL,     &
                                CELL_COLUMN, GH_WRITE, &
                                GH_READ, GH_SCALAR
  use constants_mod,     only : r_tran, i_def
  use fs_continuity_mod, only : W2h
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: hori_dep_dist_eulerian_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/               &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, W2h), & ! dep_dist xy
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W2h), & ! wind
         arg_type(GH_SCALAR, GH_REAL, GH_READ)        & ! dt
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: hori_dep_dist_eulerian_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: hori_dep_dist_eulerian_code

contains

  !> @brief Compute the advective increment in x using PPM for the advective fluxes.
  !> @param[in]     nlayers           Number of layers
  !> @param[in,out] dep_dist_xy       Departure distance
  !> @param[in]     wind              Advecting departure wind
  !> @param[in]     dt                Time step
  !> @param[in]     ndf_w2h           Number of degrees of freedom for W2h per cell
  !> @param[in]     undf_w2h          Number of unique degrees of freedom for W2h
  !> @param[in]     map_w2h           Map for W2h
  !> @param[in]     ndf_w2            Number of degrees of freedom for W2 per cell
  !> @param[in]     undf_w2           Number of unique degrees of freedom for W2
  !> @param[in]     map_w2            Map for W2

  subroutine hori_dep_dist_eulerian_code( nlayers,     &
                                          dep_dist_xy, &
                                          wind,        &
                                          dt,          &
                                          ndf_w2h,     &
                                          undf_w2h,    &
                                          map_w2h )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_w2h
    integer(kind=i_def), intent(in) :: ndf_w2h

    ! Arguments: Maps
    integer(kind=i_def), dimension(ndf_w2h), intent(in) :: map_w2h

    ! Arguments: Fields
    real(kind=r_tran), dimension(undf_w2h), intent(in)    :: wind
    real(kind=r_tran), dimension(undf_w2h), intent(inout) :: dep_dist_xy
    real(kind=r_tran), intent(in)                         :: dt

    integer(kind=i_def) :: k

    do k = 0, nlayers-1
      ! Eulerian departure distance is just dist = x_a - x_d = u dt
      dep_dist_xy( map_w2h(1) + k ) = wind( map_w2h(1) + k )*dt
      dep_dist_xy( map_w2h(3) + k ) = wind( map_w2h(3) + k )*dt
      dep_dist_xy( map_w2h(2) + k ) = wind( map_w2h(2) + k )*dt
      dep_dist_xy( map_w2h(4) + k ) = wind( map_w2h(4) + k )*dt
    end do

  end subroutine hori_dep_dist_eulerian_code

end module hori_dep_dist_eulerian_kernel_mod
