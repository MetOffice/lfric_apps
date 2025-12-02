!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Adjoint of coord_transform_mod routines
module adj_coord_transform_mod

  use constants_mod, only: r_def, PI

  implicit none

  private
  public :: adj_cart2sphere_scalar

  contains

  !=============================================================================
  !> @brief Applies the adjoint of cart2sphere_scalar
  !> @param[in]    x x component of location in Cartesian coodinates
  !> @param[in]    y y component of location in Cartesian coodinates
  !> @param[in]    z z component of location in Cartesian coodinates
  !> @param[inout] u u component of flux vector in spherical
  !>                 coordinates, out in Cartesian coordinates
  !> @param[inout] v v component of flux vector in spherical
  !>                 coordinates, out in Cartesian coordinates
  !> @param[inout] w w component of flux vector in spherical
  !>                 coordinates, out in Cartesian coordinates
  subroutine adj_cart2sphere_scalar(x, y, z, u, v, w)

    implicit none

    real(kind=r_def), intent(in)    :: x, y, z
    real(kind=r_def), intent(inout) :: u, v, w

    ! Local
    real(kind=r_def) :: t, r, phi, c1, c2
    real(kind=r_def) :: u_initial, v_initial, w_initial

    t = x**2 + y**2
    r = sqrt(t + z**2)
    phi = 0.5_r_def * PI - acos(z / r)

    ! Note: adjoint formed by taking original subroutine and adjointing
    ! equations by hand, rather than writing line-by-line adjoint
    c1 = r * cos(phi) / t
    c2 = z / (r * sqrt(t))

    u_initial = u
    v_initial = v
    w_initial = w

    u = u_initial * (-y * c1) + v_initial * (-x * c2)     + w_initial * (x / r)
    v = u_initial * (x * c1)  + v_initial * (-y * c2)     + w_initial * (y / r)
    w =                         v_initial * (sqrt(t) / r) + w_initial * (z / r)

  end subroutine adj_cart2sphere_scalar

end module adj_coord_transform_mod
