!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Tangent linear for vertical reconstruction.
!> @details The nonlinear model is:
!!          \f$ F = \rho u \f$
!!          The tangent linear model is:
!!          \f$ F = \rho ls_u + ls_\rho u \f$
!!          Near the boundaries the polynomial is no longer upwinded
!!          and reduces to an extrapolation at the boundaries.
module tl_poly1d_vert_w3_reconstruction_kernel_mod

use argument_mod,      only : arg_type, func_type,         &
                              reference_element_data_type, &
                              GH_FIELD, GH_SCALAR,         &
                              GH_REAL, GH_INTEGER,         &
                              GH_LOGICAL,                  &
                              GH_WRITE, GH_READ, GH_BASIS, &
                              CELL_COLUMN,                 &
                              ANY_DISCONTINUOUS_SPACE_1,   &
                              ANY_DISCONTINUOUS_SPACE_2
use constants_mod,     only : r_def, i_def, l_def, EPS
use fs_continuity_mod, only : W3
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: tl_poly1d_vert_w3_reconstruction_kernel_type
  private
  type(arg_type) :: meta_args(7) = (/                                        &
       arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),                        &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),                        &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                             &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                             &
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                              &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: tl_poly1d_vert_w3_reconstruction_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: tl_poly1d_vert_w3_reconstruction_code
public :: tl_poly1d_vert_w3_reconstruction_init
public :: tl_poly1d_vert_w3_reconstruction_final

integer(kind=i_def), private, allocatable, dimension(:,:,:) :: stencil

contains

!> @brief Computes the tangent linear for vertical reconstructions.
!! @param[in]  nlayers      Number of layers
!! @param[in,out] reconstruction  ACTIVE Change in mass reconstruction field
!! @param[in]  tracer       ACTIVE Change in tracer
!! @param[in]  ls_tracer    Linearisation state tracer tracer
!! @param[in]  coeff        Array of polynomial coefficients for interpolation
!! @param[in]  ndata        Number of data points per dof location
!! @param[in]  global_order Desired order of polynomial reconstruction
!! @param[in]  logspace     If true perform interpolation in log space
!! @param[in]  ndf_md      Number of degrees of freedom per cell
!! @param[in]  undf_md     Number of unique degrees of freedom for the
!!                          reconstruction & wind fields
!! @param[in]  map_md      Dofmap for the cell at the base of the column
!! @param[in]  ndf_w3       Number of degrees of freedom per cell
!! @param[in]  undf_w3      Number of unique degrees of freedom for tracer
!! @param[in]  map_w3       Cell dofmaps for the tracer space
!! @param[in]  ndf_c        Number of degrees of freedom per cell for the
!!                          coeff space
!! @param[in]  undf_c       Total number of degrees of freedom for the
!!                          coeff space
!! @param[in]  map_c        Dofmap for the coeff space
subroutine tl_poly1d_vert_w3_reconstruction_code(                &
                                                 nlayers,        &
                                                 reconstruction, &
                                                 tracer,         &
                                                 ls_tracer,      &
                                                 coeff,          &
                                                 ndata,          &
                                                 global_order,   &
                                                 logspace,       &
                                                 ndf_md,         &
                                                 undf_md,        &
                                                 map_md,         &
                                                 ndf_w3,         &
                                                 undf_w3,        &
                                                 map_w3,         &
                                                 ndf_c,          &
                                                 undf_c,         &
                                                 map_c )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                    :: nlayers
  integer(kind=i_def), intent(in)                    :: ndata
  integer(kind=i_def), intent(in)                    :: ndf_w3
  integer(kind=i_def), intent(in)                    :: undf_w3
  integer(kind=i_def), intent(in)                    :: ndf_md
  integer(kind=i_def), intent(in)                    :: undf_md
  integer(kind=i_def), intent(in)                    :: ndf_c
  integer(kind=i_def), intent(in)                    :: undf_c
  integer(kind=i_def), dimension(ndf_md), intent(in) :: map_md
  integer(kind=i_def), dimension(ndf_c),  intent(in) :: map_c
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
  integer(kind=i_def), intent(in)                    :: global_order
  real(kind=r_def), dimension(undf_md), intent(out)  :: reconstruction
  real(kind=r_def), dimension(undf_w3), intent(in)   :: tracer
  real(kind=r_def), dimension(undf_w3), intent(in)   :: ls_tracer
  real(kind=r_def), dimension(undf_c),  intent(in)   :: coeff

  logical(kind=l_def), intent(in) :: logspace

  ! Internal variables
  integer(kind=i_def) :: k, ij, ik, p, f
  integer(kind=i_def) :: vertical_order

  real(kind=r_def) :: new_tracer
  real(kind=r_def) :: ls_new_tracer

  ! Ensure that we reduce the order if there are only a few layers
  vertical_order = min(global_order, nlayers-1)

  ij = map_w3(1)

  if ( logspace ) then

    ! Loop over bottom (f=0) and top (f=1) faces
    do f = 0, 1
      ! Field is reconstructed on a layer-first W3 multidata field
      ! so the index on the bottom face (face=0) is: map_md(1) + 0*nlayers + k
      ! so the index on the top    face (face=1) is: map_md(1) + 1*nlayers + k
      ! i.e for face f index is map_md(1) + f*nlayers + k
      do k = 0, nlayers-1
        ! Compute the tracer reconstructed at W2v points
        ! Interpolate log(tracer)
        ! I.e. polynomial = exp(c_1*log(tracer_1) + c_2*log(tracer_2) + ...)
        !                 = tracer_1**c_1*tracer_2**c_2...
        ! Note that we further take the absolute value before raising to the
        ! fractional power. This code should only be used for a positive
        ! quantity, but adding in the abs ensures no errors are thrown
        ! if negative numbers are passed through in redundant calculations
        ! in the halos
        new_tracer = 0.0_r_def
        ls_new_tracer = 1.0_r_def
        do p = 1, vertical_order + 1
          ik = p + f*(global_order+1) + k*ndata + map_c(1) - 1
          new_tracer = new_tracer &
            + coeff(ik)*tracer(ij + stencil(p,k,f)) / &
              sign(max(EPS,ls_tracer(ij + stencil(p,k,f))),ls_tracer(ij + stencil(p,k,f)))

          ls_new_tracer = ls_new_tracer * &
              max(EPS,abs(ls_tracer(ij + stencil(p,k,f))))**coeff(ik)
        end do

        reconstruction(map_md(1) + f*nlayers + k) = new_tracer * ls_new_tracer

      end do
    end do

  else

    ! Loop over bottom (f=0) and top (f=1) faces
    do f = 0, 1
      ! Field is reconstructed on a layer-first W3 multidata field
      ! so the index on the bottom face (face=0) is: map_md(1) + 0*nlayers + k
      ! so the index on the top    face (face=1) is: map_md(1) + 1*nlayers + k
      ! i.e for face f index is map_md(1) + f*nlayers + k
      do k = 0, nlayers-1

        ! Compute the tracer reconstructed at W2v points
        new_tracer = 0.0_r_def
        do p = 1, vertical_order + 1
          ik = p + f*(global_order+1) + k*ndata + map_c(1) - 1
          new_tracer = new_tracer + coeff(ik)*tracer(ij + stencil(p,k,f))
        end do
        reconstruction(map_md(1) + f*nlayers + k) = new_tracer
      end do
    end do

  end if

end subroutine tl_poly1d_vert_w3_reconstruction_code

!> @brief Computes the offset stencil needed for vertical reconstructions.
!> @param[in] global_order Desired order of reconstruction
!> @param[in] nlayers      Number of layers in the mesh
subroutine tl_poly1d_vert_w3_reconstruction_init(global_order, &
                                                 nlayers)

  implicit none

  integer(kind=i_def), intent(in) :: global_order
  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def) :: k, kmin, kmax, f, p
  integer(kind=i_def) :: vertical_order, use_upwind

  integer(kind=i_def), dimension(global_order+1) :: offset


  ! Ensure that we reduce the order if there are only a few layers
  vertical_order = min(global_order, nlayers-1)

  allocate( stencil(vertical_order+1, 0:nlayers-1, 0:1) )

  ! Compute the stencil offset of points required
  ! For vertical_order = 2 => offset = (-1,0,+1)
  ! For vertical_order = 3 => offset = (-2,-1,0,+1)
  do p = 0, vertical_order
    offset(p+1) = - floor(real(vertical_order+1_i_def,r_def)/2.0_r_def) + p
  end do


  ! If order is even then we are using an upwind stencil -> use_upwind = 1
  ! For odd orders it is zero
  use_upwind = mod(vertical_order+1_i_def, 2_i_def)

  ! Loop over bottom (f=0) and top (f=1) faces
  do f = 0, 1

    ! Field is reconstructed on a layer-first W3 multidata field
    ! so the index on the bottom face (face=0) is: map_md(1) + 0*nlayers + k
    ! so the index on the top    face (face=1) is: map_md(1) + 1*nlayers + k
    ! i.e for face f index is map_md(1) + f*nlayers + k
    do k = 0, nlayers-1

      ! Compute the stencil of points required
      ! For vertical_order = 2 => stencil = (k-1,k,k+1)
      ! For vertical_order = 3 => stencil = (k-2,k-1,k,k+1)
      ! For centred scheme shift to (k-1, k, k+1, k+2) for top face
      ! (when f = 1 and use_upwind = 0)
      do p = 1, vertical_order+1
        stencil(p,k,f) = k + offset(p) + f*(1_i_def-use_upwind)
      end do

      ! Adjust stencil near boundaries to avoid going out of bounds
      kmin = minval(stencil(1:vertical_order+1,k,f))
      if ( kmin < 0 ) stencil(:,k,f) = stencil(:,k,f) - kmin
      kmax = maxval(stencil(1:vertical_order+1,k,f)) - (nlayers-1)
      if ( kmax > 0 ) stencil(:,k,f) = stencil(:,k,f) - kmax
    end do
  end do


end subroutine tl_poly1d_vert_w3_reconstruction_init

!> @brief Frees up the memory for the stencil array.
subroutine tl_poly1d_vert_w3_reconstruction_final()

  implicit none

  if ( allocated(stencil) ) deallocate(stencil)

end subroutine tl_poly1d_vert_w3_reconstruction_final

end module tl_poly1d_vert_w3_reconstruction_kernel_mod
