!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Compute the coefficients for reconstructing a 1D vertical upwind
!!        polynomial representation of a tracer field on the centre of a cell.
!> @details Compute the coefficients to reconstruct a tracer field using a high order
!!          1D polynomial fit to the values of the tracer over a given stencil.
!!          The stencil used for the polynomial is centred on the upwind cell for each edge.
!!          A polynomial is used containing all monomials up to the
!!          desired order, i.e. order = 2: 1 + x + x^2.
!!          This is exactly fitted over all cells in the stencil.
!!          The methodology is inspired by that of Thuburn et.al GMD 2014 for
!!          2D reconstructions.
!!          This method is only valid for lowest order elements.
module poly1d_vert_adv_coeffs_kernel_mod

use argument_mod,      only : arg_type, func_type,       &
                              GH_FIELD, GH_SCALAR,       &
                              GH_REAL, GH_INTEGER,       &
                              GH_WRITE, GH_READ,         &
                              ANY_SPACE_1, GH_BASIS,     &
                              ANY_DISCONTINUOUS_SPACE_1, &
                              CELL_COLUMN
use constants_mod,     only : r_def, i_def, r_tran
use kernel_mod,        only : kernel_type

implicit none

private
!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: poly1d_vert_adv_coeffs_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/                                        &
       arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_SPACE_1),               &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                             &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                              &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: poly1d_vert_adv_coeffs_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: poly1d_vert_adv_coeffs_code
contains

!> @brief Compute the coefficients needed for a 1D vertical reconstruction
!!        of a tracer field on vertical faces.
!> @param[in]     nlayers  Number of vertical layers
!> @param[in,out] coeff    Array of fields to store the coefficients for the
!!                         polynomial reconstruction
!> @param[in]     height   Height coordinate
!> @param[in]     ndata    Number of data points per dof location
!> @param[in]     order    Desired polynomial order for advective computations
!> @param[in]     ndf_c    Number of degrees of freedom per cell for the coeff space
!> @param[in]     undf_c   Total number of degrees of freedom for the coeff space
!> @param[in]     map_c    Dofmap for the coeff space
!> @param[in]     ndf_wx   Number of degrees of freedom per cell for the coordinate space
!> @param[in]     undf_wx  Total number of degrees of freedom for the coordinate space
!> @param[in]     map_wx   Dofmap of the coordinate space
subroutine poly1d_vert_adv_coeffs_code( nlayers, &
                                        coeff,   &
                                        height,  &
                                        ndata,   &
                                        order,   &
                                        ndf_c,   &
                                        undf_c,  &
                                        map_c,   &
                                        ndf_wx,  &
                                        undf_wx, &
                                        map_wx )

  use matrix_invert_mod,         only: matrix_invert_lu
  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: order
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndata
  integer(kind=i_def), intent(in) :: ndf_c,  undf_c,  &
                                     ndf_wx, undf_wx

  integer(kind=i_def), dimension(ndf_c),  intent(in) :: map_c
  integer(kind=i_def), dimension(ndf_wx), intent(in) :: map_wx

  real(kind=r_def),  dimension(undf_wx), intent(in)    :: height
  real(kind=r_tran), dimension(undf_c),  intent(inout) :: coeff

  ! Local variables
  integer(kind=i_def)                     :: k, p, df, kmin, kmax, i, j, np
  integer(kind=i_def)                     :: use_upwind, direction, vertical_order
  integer(kind=i_def), dimension(order+1) :: stencil

  real(kind=r_def)                              :: idz
  real(kind=r_def), allocatable, dimension(:)   :: z, alpha, delta
  real(kind=r_def), allocatable, dimension(:,:) :: monomial, inv_monomial

  ! Ensure that we reduce the order if there are only a few layers
  vertical_order = min(order, nlayers-1)

  ! Number of points in stencil
  np = vertical_order + 1
  allocate( z(np), alpha(np), delta(np), monomial(np,np), inv_monomial(np,np) )

  ! If order is odd then we are using an upwind stencil -> use_upwind = 1
  ! For even orders it is zero
  use_upwind = mod(vertical_order, 2_i_def)

  do df = 0, ndata-1
    coeff(map_c(1) + df) = 0.0_r_tran
  end do

  layer_loop: do k = 1, nlayers - 1
    ! Initialise polynomial coefficients to zero
    do df = 0, ndata-1
      coeff(map_c(1) + k*ndata + df) = 0.0_r_tran
    end do

    ! For upwind polynomials (odd order) compute the coefficients
    ! for both positive and negative winds (stencil increments by 1)
    ! if wind > 0 -> direction = 1
    ! if wind < 0 -> direction = 0
    ! For centred polynomials we only need one set of coefficients
    direction_loop: do direction = 0, use_upwind
      ! Compute the stencil of points required
      do p = 0, vertical_order
        stencil(p+1) = k - floor(real(vertical_order,r_def)/2.0_r_def) + p
      end do

      ! Adjust the stencil based upon the wind sign for upwind (odd order)
      ! reconstructions only.
      stencil = stencil - direction
      ! Adjust stencil near boundaries to avoid going out of bounds
      kmin = minval(stencil(1:np))
      if ( kmin < 0 ) stencil = stencil - kmin
      kmax = maxval(stencil(1:np)) - nlayers
      if ( kmax > 0 ) stencil = stencil - kmax
      ! Compute local coordinates (assuming z(k) is origin)
      do p = 1, np
        z(p) = (height(map_wx(1) + stencil(p)) - height(map_wx(1)+k))
      end do
      ! Normalise by average of cell widths around reconstruction point
      ! this term cancels with the implied dz in the mass matrix
      ! when applying the reconstruction so it is important to use the same
      ! value (and not just e.g the upwind spacing)
      idz = 2.0_r_def/( height(map_wx(1) + k + 1) &
                      - height(map_wx(1) + k - 1) )
      z = z*idz

      ! Fit polynomial P = a0 + a1*z + a2*z^2 + a3*z^3 to d
      ! by solving M * [a0, a1, a2, a3]^T = [d(1), d(2), d(3), d(4)]^T
      ! Where monomial matrix is M_ij = z(i)^(j-1)
      ! and d is a delta function such that d_i = 1 in cell i and
      ! solve system for all possible i's

      do i = 1, np
        do j = 1, np
          if (j == 1) then
            monomial(i,j) = 1.0_r_def
          else
            monomial(i,j) = z(i)**(j-1)
          end if
        end do
      end do
      call matrix_invert_lu(monomial, inv_monomial, np)

      do p = 1, np
        delta = 0.0_r_def
        delta(p) = 1.0_r_def

        alpha = matmul(inv_monomial, delta)
        ! dPdz = diff(P,z) = sum_i=1^np ( i*a(i)*z^(i-1) )
        ! evaluated at z = 0 -> dPdz = a(2)
        coeff( map_c(1) + k*ndata + direction*(order+1) + p - 1 ) = real( alpha(2), r_tran )
      end do
    end do direction_loop
  end do layer_loop

  k = nlayers
  do df = 0, ndata-1
    coeff(map_c(1) + k*ndata + df) = 0.0_r_tran
  end do

  deallocate( z, alpha, delta, monomial, inv_monomial )
end subroutine poly1d_vert_adv_coeffs_code

end module poly1d_vert_adv_coeffs_kernel_mod
