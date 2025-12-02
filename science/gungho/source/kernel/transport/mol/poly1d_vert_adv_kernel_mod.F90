!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes vertical derivative of a field through fitting a high order 1D
!!        upwind reconstruction.
!> @details Computes the derivative for a tracer density field using a high order
!!          polynomial fit to the tracer values. For odd order polynomials
!!          the stencil is centred on the upwind cell.
!!          Near the boundaries the order of reconstruction may be reduced
!!          if there are not enough points to compute desired order but
!!          is always kept > 1
!!          This method is only valid for lowest order elements.
module poly1d_vert_adv_kernel_mod

use argument_mod,         only : arg_type, func_type,   &
                                 GH_FIELD, GH_SCALAR,   &
                                 GH_REAL, GH_INTEGER,   &
                                 GH_LOGICAL,            &
                                 GH_READWRITE, GH_READ, &
                                 GH_BASIS, CELL_COLUMN, &
                                 GH_EVALUATOR,          &
                                 ANY_DISCONTINUOUS_SPACE_1
use constants_mod,        only : r_tran, i_def, l_def, EPS_R_TRAN
use fs_continuity_mod,    only : W2v, Wtheta
use kernel_mod,           only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: poly1d_vert_adv_kernel_type
  private
  type(arg_type) :: meta_args(7) = (/                                            &
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, Wtheta),                    &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W2v),                       &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                                 &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                                 &
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                                  &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: poly1d_vert_adv_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: poly1d_vert_adv_code

contains

!> @brief Computes the vertical derivative for a tracer field.
!> @param[in]     nlayers      Number of layers
!> @param[in,out] advective    Advective update to increment
!> @param[in]     wind         Wind field
!> @param[in]     tracer       Tracer field to advect
!> @param[in]     coeff        Array of polynomial coefficients for interpolation
!> @param[in]     ndata        Number of data points per dof location
!> @param[in]     global_order Desired order of polynomial reconstruction
!> @param[in]     logspace     If true then perform interpolation in log space
!> @param[in]     ndf_wt       Number of degrees of freedom per cell
!> @param[in]     undf_wt      Number of unique degrees of freedom for the tracer field
!> @param[in]     map_wt       Cell dofmaps for the tracer space
!> @param[in]     ndf_w2v      Number of degrees of freedom per cell
!> @param[in]     undf_w2v     Number of unique degrees of freedom for the flux &
!!                             wind fields
!> @param[in]     map_w2v      Dofmap for the cell at the base of the column
!> @param[in]     ndf_c        Number of degrees of freedom per cell for the coeff space
!> @param[in]     undf_c       Total number of degrees of freedom for the coeff space
!> @param[in]     map_c        Dofmap for the coeff space
subroutine poly1d_vert_adv_code( nlayers,              &
                                 advective,            &
                                 wind,                 &
                                 tracer,               &
                                 coeff,                &
                                 ndata,                &
                                 global_order,         &
                                 logspace,             &
                                 ndf_wt,               &
                                 undf_wt,              &
                                 map_wt,               &
                                 ndf_w2v,              &
                                 undf_w2v,             &
                                 map_w2v,              &
                                 ndf_c,                &
                                 undf_c,               &
                                 map_c)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                     :: nlayers
  integer(kind=i_def), intent(in)                     :: ndf_wt
  integer(kind=i_def), intent(in)                     :: undf_wt
  integer(kind=i_def), intent(in)                     :: ndf_w2v
  integer(kind=i_def), intent(in)                     :: undf_w2v
  integer(kind=i_def), dimension(ndf_w2v), intent(in) :: map_w2v
  integer(kind=i_def), dimension(ndf_wt),  intent(in) :: map_wt
  integer(kind=i_def), intent(in)                     :: ndf_c
  integer(kind=i_def), intent(in)                     :: undf_c
  integer(kind=i_def), dimension(ndf_c),   intent(in) :: map_c
  integer(kind=i_def), intent(in)                     :: ndata
  integer(kind=i_def), intent(in)                     :: global_order

  real(kind=r_tran), dimension(undf_wt),  intent(inout) :: advective
  real(kind=r_tran), dimension(undf_w2v), intent(in)    :: wind
  real(kind=r_tran), dimension(undf_wt),  intent(in)    :: tracer
  real(kind=r_tran), dimension(undf_c),   intent(in)    :: coeff

  logical(kind=l_def), intent(in) :: logspace

  ! Internal variables
  integer(kind=i_def) :: k, kmin, kmax, ij, ik, p
  integer(kind=i_def) :: vertical_order, use_upwind, upwind_offset, upwind

  integer(kind=i_def), dimension(global_order+1) :: stencil

  real(kind=r_tran) :: dpdz
  real(kind=r_tran), dimension(0:nlayers) :: log_tracer

  ij = map_wt(1)

  ! Compute log of tracer. This code should only be used for a positive
  ! quantity, but adding in the abs ensures no errors are thrown
  ! if negative numbers are passed through in redundant calculations
  ! in the halos
  ! TODO #3290: if tracer is zero this could cause problems
  if ( logspace ) then
    do k = 0, nlayers
      log_tracer(k) = log(max(EPS_R_TRAN,abs(tracer(ij+k))))
    end do
  end if

  ! Ensure that we reduce the order if there are only a few layers
  vertical_order = min(global_order, nlayers-1)

  ! If order is odd then we are using an upwind stencil -> use_upwind = 1
  ! For even orders it is zero
  use_upwind = mod(vertical_order, 2_i_def)

  ! Compute dtracer/dz using precomputed weights
  do k = 1, nlayers - 1

    ! Compute the stencil of points required
    do p = 0, vertical_order
      stencil(p+1) = k - floor(real(vertical_order,r_tran)/2.0_r_tran) + p
    end do

    ! Adjust the stencil based upon the wind sign for upwind (odd order)
    ! reconstructions only.
    ! if wind > 0 -> upwind_offset = 1
    ! if wind < 0 -> upwind_offset = 0
    upwind = nint(0.5_r_tran*(1.0_r_tran + sign(1.0_r_tran,wind(map_w2v(1)+k))),i_def)
    upwind_offset = use_upwind*upwind
    stencil = stencil - upwind_offset

    ! Adjust stencil near boundaries to avoid going out of bounds
    kmin = minval(stencil(1:vertical_order+1))
    if ( kmin < 0 ) stencil = stencil - kmin
    kmax = maxval(stencil(1:vertical_order+1)) - nlayers
    if ( kmax > 0 ) stencil = stencil - kmax

    ! Compute the derivative and the advective update
    dpdz = 0.0_r_tran
    if ( logspace ) then
      ! dp/dz = p * d(log(p))/dz
      do p = 1, vertical_order + 1
        ik = p + upwind_offset*(global_order+1) + k*ndata + map_c(1) - 1
        dpdz = dpdz + coeff(ik)*log_tracer(stencil(p))
      end do
      dpdz = tracer(ij + k)*dpdz
    else
      do p = 1, vertical_order + 1
        ik = p + upwind_offset*(global_order+1) + k*ndata + map_c(1) - 1
        dpdz = dpdz + coeff(ik)*tracer(ij + stencil(p))
      end do
    end if
    advective(map_wt(1)+k) = advective(map_wt(1)+k) &
                             + wind(map_w2v(1)+k)*dpdz
  end do

end subroutine poly1d_vert_adv_code

end module poly1d_vert_adv_kernel_mod
