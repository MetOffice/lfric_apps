!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Calculates the Lagrangian trapezoidal horizontal departure distance.
!> @details This kernel computes the Lagrangian horizontal departure distance
!!          using the advecting departure wind u and time step dt as
!!          \f$\mbox{dist} = x_a - x_d = f(u_n,u_np1) dt\f$,
!!          where x_d is the departure point, x_a the arrival point, and
!!          f(u_n,u_np1) is the wind interpolated to the correct point/time using
!!          the trapezoidal approach.
!!          This kernel is only suitable for lowest order W2 spaces and planar geometry.

module hori_dep_dist_trapezoidal_kernel_mod

  use argument_mod,      only : arg_type,              &
                                GH_FIELD, GH_REAL,     &
                                CELL_COLUMN, GH_WRITE, &
                                GH_READ, GH_SCALAR,    &
                                GH_INTEGER, STENCIL,   &
                                CROSS
  use constants_mod,     only : r_tran, i_def
  use fs_continuity_mod, only : W2h
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: hori_dep_dist_trapezoidal_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                  &
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W2h),                 & ! dep_dist
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2h, STENCIL(CROSS)), & ! wind_n
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2h),                 & ! wind_np1
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                       & ! n_iterations
         arg_type(GH_SCALAR, GH_REAL,    GH_READ)                        & ! dt
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: hori_dep_dist_trapezoidal_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: hori_dep_dist_trapezoidal_code

contains

  !> @brief Compute horizontal departure points using the trapezoidal method.
  !>
  !> If the departure point happens to be an integer length which falls on a
  !> cell face then the next cell out along the stencil arm will be needed.
  !> This is true at the extremity of the stencil, so although the departure is
  !> contained wholly within the specified stencil, the kernel will still fail
  !> as it attempts to access beyond the bounds of the array.
  !>
  !> @param[in]     nlayers           Number of layers
  !> @param[in,out] dep_dist          Departure distance
  !> @param[in]     wind_n            Advecting departure wind at time n
  !> @param[in]     stencil_size      Local length of W2 stencil
  !> @param[in]     stencil_map       Dofmap for the W2 stencil
  !> @param[in]     wind_np1          Advecting departure wind at time n+1
  !> @param[in]     n_iterations      Number of iterations
  !> @param[in]     dt                Time step
  !> @param[in]     ndf_w2h           Number of degrees of freedom for W2h per cell
  !> @param[in]     undf_w2h          Number of unique degrees of freedom for W2h
  !> @param[in]     map_w2h           Map for W2h

  subroutine hori_dep_dist_trapezoidal_code( nlayers,      &
                                             dep_dist,     &
                                             wind_n,       &
                                             stencil_size, &
                                             stencil_map,  &
                                             wind_np1,     &
                                             n_iterations, &
                                             dt,           &
                                             ndf_w2h,      &
                                             undf_w2h,     &
                                             map_w2h )

    use departure_points_support_mod, only : interpolate_u_to_x

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_w2h
    integer(kind=i_def), intent(in) :: ndf_w2h
    integer(kind=i_def), intent(in) :: stencil_size

    ! Arguments: Maps
    integer(kind=i_def), dimension(ndf_w2h), intent(in) :: map_w2h
    integer(kind=i_def), dimension(ndf_w2h,stencil_size), intent(in) :: stencil_map

    ! Arguments: Fields
    real(kind=r_tran), dimension(undf_w2h), intent(in)    :: wind_n
    real(kind=r_tran), dimension(undf_w2h), intent(in)    :: wind_np1
    real(kind=r_tran), dimension(undf_w2h), intent(inout) :: dep_dist
    real(kind=r_tran),                      intent(in)    :: dt
    integer(kind=i_def),                    intent(in)    :: n_iterations

    integer(kind=i_def) :: k, jj, ind, stencil_half, stencil_centre, stencil_arm
    real(kind=r_tran)   :: dp_local_x_1
    real(kind=r_tran)   :: dp_local_x_3
    real(kind=r_tran)   :: dp_local_y_2
    real(kind=r_tran)   :: dp_local_y_4
    real(kind=r_tran)   :: u_local_x_1(1:(stencil_size+1)/2)
    real(kind=r_tran)   :: u_local_x_3(1:(stencil_size+1)/2)
    real(kind=r_tran)   :: u_local_y_2(1:(stencil_size+1)/2)
    real(kind=r_tran)   :: u_local_y_4(1:(stencil_size+1)/2)
    real(kind=r_tran)   :: u_n_d_x_1
    real(kind=r_tran)   :: u_n_d_x_3
    real(kind=r_tran)   :: u_n_d_y_2
    real(kind=r_tran)   :: u_n_d_y_4
    real(kind=r_tran)   :: u_np1_a_x_1
    real(kind=r_tran)   :: u_np1_a_x_3
    real(kind=r_tran)   :: u_np1_a_y_2
    real(kind=r_tran)   :: u_np1_a_y_4

    ! Cross stencil has the form
    !     | 5 |
    ! | 2 | 1 | 4 |
    !     | 3 |
    ! The local variables become 1D arrays
    ! | 1 | 2 | 3 |
    !

    ! Get half stencil size and local stencil centre value
    stencil_half   = (stencil_size + 1_i_def) / 2_i_def
    stencil_centre = (stencil_half + 1_i_def) / 2_i_def
    stencil_arm    = (stencil_size - 1_i_def) / 4_i_def

    do k = 0,nlayers-1
      ! Set wind into local stencil
      do jj = 1, stencil_centre
        u_local_x_1(jj) = wind_n(stencil_map(1,stencil_centre+1-jj) + k)
        u_local_x_3(jj) = wind_n(stencil_map(3,stencil_centre+1-jj) + k)
      end do
      do jj = stencil_centre+1, stencil_half
        u_local_x_1(jj) = wind_n(stencil_map(1,stencil_arm+jj) + k)
        u_local_x_3(jj) = wind_n(stencil_map(3,stencil_arm+jj) + k)
      end do
      ! NB: minus sign for y-winds, to align these with the stencil direction
      u_local_y_2(stencil_centre) = -wind_n(stencil_map(2,1) + k)
      u_local_y_4(stencil_centre) = -wind_n(stencil_map(4,1) + k)
      do jj = 1, stencil_centre-1
        u_local_y_2(jj) = -wind_n(stencil_map(2,stencil_arm+stencil_centre+1-jj) + k)
        u_local_y_4(jj) = -wind_n(stencil_map(4,stencil_arm+stencil_centre+1-jj) + k)
      end do
      do jj = stencil_centre+1, stencil_half
        u_local_y_2(jj) = -wind_n(stencil_map(2,2*stencil_arm+jj) + k)
        u_local_y_4(jj) = -wind_n(stencil_map(4,2*stencil_arm+jj) + k)
      end do

      ! Set Eulerian departure distance as first estimate: dist = x_a - x_d = u dt
      dp_local_x_1 = u_local_x_1(stencil_centre)*dt
      dp_local_x_3 = u_local_x_3(stencil_centre)*dt
      dp_local_y_2 = u_local_y_2(stencil_centre)*dt
      dp_local_y_4 = u_local_y_4(stencil_centre)*dt

      ! Trapezoidal method has the form:
      ! x_d^{ind+1} = x_a - dt/2( u_np1(x_a) + u_n(x_d^{ind}) )

      ! Get local u_np1 at arrival point
      u_np1_a_x_1 = wind_np1(map_w2h(1) + k)
      u_np1_a_x_3 = wind_np1(map_w2h(3) + k)
      u_np1_a_y_2 = -wind_np1(map_w2h(2) + k)
      u_np1_a_y_4 = -wind_np1(map_w2h(4) + k)

      do ind = 1, n_iterations
        ! Interpolate u_n to the departure point
        call interpolate_u_to_x(u_n_d_x_1, dp_local_x_1, u_local_x_1, stencil_half, stencil_centre)
        call interpolate_u_to_x(u_n_d_x_3, dp_local_x_3, u_local_x_3, stencil_half, stencil_centre)
        call interpolate_u_to_x(u_n_d_y_2, dp_local_y_2, u_local_y_2, stencil_half, stencil_centre)
        call interpolate_u_to_x(u_n_d_y_4, dp_local_y_4, u_local_y_4, stencil_half, stencil_centre)

        ! Update the departure distance with the correct sign
        dp_local_x_1 = 0.5_r_tran * dt * ( u_np1_a_x_1 + u_n_d_x_1 )
        dp_local_x_3 = 0.5_r_tran * dt * ( u_np1_a_x_3 + u_n_d_x_3 )
        dp_local_y_2 = 0.5_r_tran * dt * ( u_np1_a_y_2 + u_n_d_y_2 )
        dp_local_y_4 = 0.5_r_tran * dt * ( u_np1_a_y_4 + u_n_d_y_4 )
      end do

      dep_dist( map_w2h(1) + k ) = dp_local_x_1
      dep_dist( map_w2h(3) + k ) = dp_local_x_3
      ! NB: minus sign for y-winds, to return these to the same convention as
      ! the rest of the model
      dep_dist( map_w2h(2) + k ) = -dp_local_y_2
      dep_dist( map_w2h(4) + k ) = -dp_local_y_4
    end do

  end subroutine hori_dep_dist_trapezoidal_code

end module hori_dep_dist_trapezoidal_kernel_mod
