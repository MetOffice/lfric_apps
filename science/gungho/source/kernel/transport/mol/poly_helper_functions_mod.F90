!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Contains a number of helper functions for use in the poly1d & poly2d
!!        reconstruction kernels used by the MoL advection scheme.
module poly_helper_functions_mod

use constants_mod, only: r_def, i_def, l_def, EPS

implicit none

private

public :: buildadvcoeff
public :: local_distance_2d
public :: local_distance_1d

contains

!=============================================================================!
!> @brief Build the coefficients used to construct polynomial fits for the
!!        advection scheme.
!> @details Build the coefficients used to construct polynomial fits for
!!          advection scheme. The polynomial fit should be exact for
!!          the central cell in the stencil, and the mean square
!!          residual over the other cells should be minimized. This leads to
!!          to a linear problem for the coefficients of the polynomial fit
!!          and some Lagrange multipliers. Finding the inverse of the matrix
!!          involved allows us to save the matrix that gives the polynomial
!!          coefficients in terms of the cell integrals of the advected quantity.
!> @param[in,out] int_monomial Matrix containing the integrals of monomials
!!                             over the stencil
!> @param[in]     ns           Number of cells in the stencil
!> @param[in]     nmonomial    Number of monomials in the polynomial fit
subroutine buildadvcoeff( int_monomial, ns, nmonomial )

  use matrix_invert_mod, only: matrix_invert

  implicit none

  integer(kind=i_def),                           intent(in)    :: ns, nmonomial
  real(kind=r_def),    dimension(ns, nmonomial), intent(inout) :: int_monomial

  integer(kind=i_def) :: i, ne, nm

  real(kind=r_def), allocatable, dimension(:,:) :: l, lt, ltl, m, minv, r, q

  ! Determine the stencil size and the number of cells fitted
  ! exactly; hence allocate space for matrices
  ne = 1_i_def
  nm = nmonomial + ne
  allocate(l(ns,nmonomial), lt(nmonomial,ns), ltl(nmonomial,nmonomial))
  allocate(m(nm,nm), minv(nm,nm), r(nm,ns), q(nm,ns))

  ! Extract the matrix of integrals of monomials
  l = int_monomial(1:ns,1:nmonomial)
  lt = transpose(l)
  ltl = matmul(lt,l)

  ! Build matrix for linear system
  m(1:nmonomial,1:nmonomial) = ltl
  m(1:nmonomial,nmonomial+1:nm) = lt(1:nmonomial,1:ne)
  m(nmonomial+1:nm,1:nmonomial) = l(1:ne,1:nmonomial)
  m(nmonomial+1:nm,nmonomial+1:nm) = 0.0_r_def

  ! Invert matrix
  call matrix_invert(m,minv,nm)

  ! Matrix relating cell integrals to RHS of linear system
  r(1:nmonomial,1:ns) = lt(1:nmonomial,1:ns)
  r(nmonomial+1:nm,1:ns) = 0.0_r_def
  do i = 1,ne
    r(nmonomial+i,i) = 1.0_r_def
  end do

  ! Matrix giving polynomial coefficients in terms of cell integrals
  q = matmul(minv,r)

  ! Save the part we need
  lt = q(1:nmonomial,1:ns)
  l = transpose(lt)
  int_monomial(1:ns,1:nmonomial) = l

  deallocate(l, lt, ltl, m, minv, r, q)

end subroutine buildadvcoeff

!=============================================================================!

  !> @brief Compute the distance between two points x0 & x1
  !!        in terms of a local coordinate.
  !> @param[in] x0        First point (x,y,z)
  !> @param[in] x1        Second point (x,y,z)
  !> @param[in] xn        Normal vector used to define local coordinate direction
  !> @param[in] domain_x  Domain size in x-direction.
  !> @param[in] domain_y  Domain size in y-direction.
  !> @param[in] spherical Switch for spherical or cartesian computation
  !> @return s Distance between x0 and x1
  function local_distance_2d( x0, x1, xn,         &
                              domain_x, domain_y, &
                              spherical ) result(s)

    use cross_product_mod,      only: cross_product

    implicit none

    real(kind=r_def), dimension(3), intent(in) :: x0, x1, xn
    real(kind=r_def),               intent(in) :: domain_x, domain_y
    logical(kind=l_def),            intent(in) :: spherical
    real(kind=r_def), dimension(2)             :: s

    real(kind=r_def), dimension(3)             :: y0, y1, xn1, dx
    real(kind=r_def)                           :: mag, costh, sinth

    if ( spherical ) then
      ! Normalise position vectors
      y0 = x0/sqrt(x0(1)**2 + x0(2)**2 + x0(3)**2)
      y1 = x1/sqrt(x1(1)**2 + x1(2)**2 + x1(3)**2)
      ! Normal to plane containing x0 & x1
      xn1 = cross_product(y0,y1)
      mag = sqrt(xn1(1)**2 + xn1(2)**2 + xn1(3)**2)
      ! Catch error if x0 == x1 and so xn1 == 0 > mag = 0
      if ( mag > EPS ) xn1 = xn1/mag
      ! Angle relative to local x-axis
      sinth = dot_product(y0,cross_product(xn,xn1))
      costh = dot_product(xn,xn1)
      ! Finally obtain local coordinate
      s(1) = asin(mag)*costh
      s(2) = asin(mag)*sinth
    else
      ! Checks if x0 and x1 are across a periodic mesh wrap point
      ! this uses a simple (non-robust) check.
      ! if  the distance is bigger than half the domain size then it has
      ! probably wrapped around so we modify the distance measure
      dx = x1 - x0

      if ( abs(dx(1)) > 0.5_r_def*domain_x ) &
        dx(1) = dx(1) - sign(1.0_r_def,dx(1))*domain_x
      if ( abs(dx(2)) > 0.5_r_def*domain_y ) &
        dx(2) = dx(2) - sign(1.0_r_def,dx(2))*domain_y

      ! Origin of local system is x0 and we only need the horizontal distance
      s(1) = -dx(1)*xn(2) + dx(2)*xn(1)
      s(2) = -dx(2)*xn(2) + dx(1)*xn(1)
    end if

  end function local_distance_2d

  !> @brief Compute the distance between two points x0 & x1
  !!        in terms of a local coordinate.
  !> @param[in] x0        First point (x,y,z)
  !> @param[in] x1        Second point (x,y,z)
  !> @param[in] xn        Normal vector used to define local coordinate direction
  !> @param[in] domain_x  Domain size in x-direction.
  !> @param[in] domain_y  Domain size in y-direction.
  !> @param[in] spherical Switch for spherical or cartesian computation
  !> @return s Distance between x0 and x1
  function local_distance_1d( x0, x1, xn,         &
                              domain_x, domain_y, &
                              spherical ) result(s)

    use cross_product_mod,      only: cross_product

    implicit none

    real(kind=r_def), dimension(3), intent(in) :: x0, x1, xn
    real(kind=r_def),               intent(in) :: domain_x, domain_y
    logical(kind=l_def),            intent(in) :: spherical
    real(kind=r_def)                           :: s

    real(kind=r_def), dimension(3)             :: y0, y1, xn1, dx
    real(kind=r_def)                           :: mag, costh


    if ( spherical ) then
      ! Normalise position vectors
      y0 = x0/sqrt(x0(1)**2 + x0(2)**2 + x0(3)**2)
      y1 = x1/sqrt(x1(1)**2 + x1(2)**2 + x1(3)**2)
      ! Normal to plane containing x0 & x1
      xn1 = cross_product(y0,y1)
      mag = sqrt(xn1(1)**2 + xn1(2)**2 + xn1(3)**2)
      ! Catch error if x0 == x1 and so xn1 == 0 > mag = 0
      if ( mag > EPS ) xn1 = xn1/mag
      ! Angle relative to local x-axis
      costh = dot_product(xn,xn1)
      ! Finally obtain local coordinate
      s = asin(mag)*costh
    else
      ! Checks if x0 and x1 are across a periodic mesh wrap point
      ! this uses a simple (non-robust) check.
      ! if  the distance is bigger than half the domain size then it has
      ! probably wrapped around so we modify the distance measure
      dx = x1 - x0

      if ( abs(dx(1)) > 0.5_r_def*domain_x ) &
        dx(1) = dx(1) - sign(1.0_r_def,dx(1))*domain_x
      if ( abs(dx(2)) > 0.5_r_def*domain_y ) &
        dx(2) = dx(2) - sign(1.0_r_def,dx(2))*domain_y

      ! Origin of local system is x0 and we only need the horizontal distance
      s = -dx(1)*xn(2) + dx(2)*xn(1)
    end if

  end function local_distance_1d

end module poly_helper_functions_mod
