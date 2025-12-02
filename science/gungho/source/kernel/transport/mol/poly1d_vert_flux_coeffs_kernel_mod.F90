!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Compute the coefficients for reconstructing a 1D vertical upwind
!!        polynomial representation of a tracer field on the faces of a cell.
!> @details Compute the coefficients of the flux of a tracer density field using a high order
!!          1D polynomial fit to the integrated tracer values over a given stencil.
!!          For even order polynomials the stencil is centred on the upwind cell for each edge.
!!          For odd order polynomials the stencil is centred on the edge the
!!          reconstruction is computed at.
!!          A polynomial is used containing all monomials up to the
!!          desired order, i.e. order = 2: 1 + x + x^2.
!!          This is exactly fitted over all cells in the stencil.
!!          The methodology is inspired by that of Thuburn et.al GMD 2014 for
!!          2D reconstructions.
!!          This method is only valid for lowest order elements.
module poly1d_vert_flux_coeffs_kernel_mod

use argument_mod,      only : arg_type, func_type,       &
                              GH_FIELD, GH_SCALAR,       &
                              GH_REAL, GH_INTEGER,       &
                              GH_WRITE, GH_READ,         &
                              GH_BASIS, CELL_COLUMN,     &
                              GH_QUADRATURE_XYoZ,        &
                              ANY_DISCONTINUOUS_SPACE_1
use constants_mod,     only : r_def, i_def, r_tran
use fs_continuity_mod, only : Wtheta
use kernel_mod,        only : kernel_type

implicit none

private
!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: poly1d_vert_flux_coeffs_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/                                        &
       arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                    &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                             &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                              &
       /)
  type(func_type) :: meta_funcs(1) = (/ &
       func_type(Wtheta, GH_BASIS)      &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_QUADRATURE_XYoZ
contains
  procedure, nopass :: poly1d_vert_flux_coeffs_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: poly1d_vert_flux_coeffs_code
contains

!> @brief Compute the coefficients needed for a 1D vertical reconstruction
!>        of a tracer field on vertical faces.
!> @param[in] nlayers Number of vertical layers
!> @param[in,out] coeff Array of fields to store the coefficients for the polynomial
!!                      reconstruction
!> @param[in] height Physical height above the surface
!> @param[in] ndata Number of data points per dof location
!> @param[in] global_order Desired polynomial order for flux computations
!> @param[in] ndf_c Number of degrees of freedom per cell for the coeff space
!> @param[in] undf_c Total number of degrees of freedom for the coeff space
!> @param[in] map_c Dofmap for the coeff space
!> @param[in] ndf_wx Number of degrees of freedom per cell for the coordinate space
!> @param[in] undf_wx Total number of degrees of freedom for the coordinate space
!> @param[in] map_wx Dofmap of the coordinate space stencil
!> @param[in] basis_wx Basis function of the coordinate space evaluated on
!!                     quadrature points
!!                          quadrature points on the vertical faces
!> @param[in] nqp_h Number of horizontal quadrature points
!> @param[in] nqp_v Number of vertical quadrature points
!> @param[in] wqp_h Weights of horizontal quadrature points
!> @param[in] wqp_v Weights of vertical quadrature points
subroutine poly1d_vert_flux_coeffs_code(nlayers,                    &
                                        coeff,                      &
                                        height,                     &
                                        ndata,                      &
                                        global_order,               &
                                        ndf_c,                      &
                                        undf_c,                     &
                                        map_c,                      &
                                        ndf_wx,                     &
                                        undf_wx,                    &
                                        map_wx,                     &
                                        basis_wx,                   &
                                        nqp_h, nqp_v, wqp_h, wqp_v )

  use matrix_invert_mod, only: matrix_invert_lu

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: global_order
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndata
  integer(kind=i_def), intent(in) :: ndf_c,  undf_c,  &
                                     ndf_wx, undf_wx
  integer(kind=i_def), intent(in) :: nqp_v, nqp_h

  integer(kind=i_def), dimension(ndf_c),  intent(in) :: map_c
  integer(kind=i_def), dimension(ndf_wx), intent(in) :: map_wx

  real(kind=r_def), dimension(undf_wx), intent(in)    :: height
  real(kind=r_tran), dimension(undf_c), intent(inout) :: coeff

  real(kind=r_def), dimension(1,ndf_wx,nqp_h,nqp_v), intent(in) :: basis_wx

  real(kind=r_def), dimension(nqp_h), intent(in) ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in) ::  wqp_v

  ! Local variables
  integer(kind=i_def) :: k, p, df, kmin, kmax, j, np, qp, ik, qh0
  integer(kind=i_def) :: use_upwind, direction, vertical_order

  integer(kind=i_def), dimension(global_order+1) :: stencil

  real(kind=r_def)                              :: z, z0
  real(kind=r_tran)                             :: total_coeff
  real(kind=r_def), allocatable, dimension(:)   :: alpha, delta
  real(kind=r_def), allocatable, dimension(:,:) :: monomial, inv_monomial

  ! Ensure that we reduce the order if there are only a few layers
  vertical_order = min(global_order, nlayers-1)

  ! Number of points in stencil
  np = vertical_order + 1
  allocate( alpha(np), delta(np), monomial(np,np), inv_monomial(np,np) )

  ! If order is even then we are using an upwind stencil -> use_upwind = 1
  ! For odd orders it is zero
  use_upwind = mod(vertical_order+1_i_def, 2_i_def)

  ! Index of quadrature point in the centre of the cell
  ! (this is only true if the number of quadrature points is odd)
  qh0 = (nqp_h+1)/2

  ! Step 1: Build integrals of monomials over all cells in advection stencils

  ! Loop over all layers, coefficients are defined on Wtheta points
  ! so loop 0 (bottom boundary) to nlayers (top boundary)
  layer_loop: do k = 0, nlayers
    ! Initialise polynomial coefficients to zero
    do df = 0, ndata-1
      coeff(map_c(1) + k*ndata + df) = 0.0_r_tran
    end do

    ! Origin of local coordinate system
    z0 = height(map_wx(1)+k)

    ! For upwind polynomials (even order) compute the coefficients
    ! for both positive and negative winds (stencil increments by 1)
    ! if wind < 0 -> direction = 0 (use one more cell above than below)
    ! if wind > 0 -> direction = 1 (use one more cell below than above)
    ! For centred polynomials (odd order) we only need one set of coefficients
    direction_loop: do direction = 0, use_upwind

      ! Compute the stencil of points required
      ! For vertical_order = 2 => stencil = (k-1,k,k+1)
      ! For vertical_order = 3 => stencil = (k-2,k-1,k,k+1)
      do p = 0, vertical_order
        stencil(p+1) = k - floor(real(vertical_order + 1_i_def,r_def)/2.0_r_def) + p
      end do

      ! Adjust the stencil in the upwind direction unless at the top boundary
      if ( k < nlayers ) stencil = stencil - direction

      ! Adjust stencil near boundaries to avoid going out of bounds
      kmin = minval(stencil(1:np))
      if ( kmin < 0 ) stencil = stencil - kmin
      kmax = maxval(stencil(1:np)) - (nlayers - 1)
      if ( kmax > 0 ) stencil = stencil - kmax

      monomial = 0.0_r_def

      quadrature_loop: do qp = 1, nqp_v
        ! Compute local coordinates (assuming z(k) is origin)
        ! and monomial matrix m(i,j) = int( z^(j-1) dz_i)
        do p = 1, np
          z = 0.0_r_def
          do df = 1, ndf_wx
            z = z + height(map_wx(df) + stencil(p))*basis_wx(1,df,qh0,qp)
          end do
          z = z - z0
          do j = 1, np
            monomial(p,j) = monomial(p,j) + wqp_v(qp)*z**(j-1)
          end do
        end do
      end do quadrature_loop
      call matrix_invert_lu(monomial, inv_monomial, np)

      ! Fit polynomial P = a0 + a1*z + a2*z^2 + a3*z^3 to d
      ! by solving M * [a0, a1, a2, a3]^T = [d(1), d(2), d(3), d(4)]^T
      ! Where monomial matrix is M_ij = int (z^(j-1) dz_i)
      ! and d is a delta function such that d_i = 1 in cell i and
      ! solve system for all possible i's
      do p = 1, np
        delta = 0.0_r_def
        delta(p) = 1.0_r_def

        alpha = matmul(inv_monomial, delta)
        ! P(z) = sum_i=1^np ( a(i)*z^(i-1) )
        ! evaluated at z = 0 -> P(0) = a(1)
        ik = map_c(1) + k*ndata + direction*(global_order+1) - 1
        coeff( ik + p ) = real( alpha(1), r_tran )
      end do

      ! Normalise coeffs to ensure a constant remains exacly constant
      p =  map_c(1) + k*ndata + direction*(global_order+1)
      total_coeff = 1.0_r_tran/sum(coeff(p:p+vertical_order))
      coeff(p:p+vertical_order) = coeff(p:p+vertical_order)*total_coeff

    end do direction_loop
  end do layer_loop

  deallocate( alpha, delta, monomial, inv_monomial )

end subroutine poly1d_vert_flux_coeffs_code

end module poly1d_vert_flux_coeffs_kernel_mod
