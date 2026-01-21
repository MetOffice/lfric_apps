!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Adjoint kernel for w3 reconstruction in vertical.
module atl_poly1d_vert_w3_reconstruction_kernel_mod
use argument_mod,     only : arg_type, func_type,         &
                             GH_FIELD,                    &
                             GH_READ, GH_READWRITE,       &
                             GH_WRITE, GH_INC,            &
                             GH_SCALAR, GH_REAL,          &
                             GH_INTEGER, GH_LOGICAL,      &
                             CELL_COLUMN,                 &
                             ANY_DISCONTINUOUS_SPACE_1,   &
                             ANY_DISCONTINUOUS_SPACE_2
use constants_mod,     only : eps, i_def, l_def, r_def
use fs_continuity_mod, only : W3
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: atl_poly1d_vert_w3_reconstruction_kernel_type
type(arg_type) :: meta_args(7) = (/                                        &
     arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), &
     arg_type(GH_FIELD, GH_REAL, GH_READWRITE, W3),                        &
     arg_type(GH_FIELD, GH_REAL, GH_READ, W3),                             &
     arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_2),      &
     arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                             &
     arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                             &
     arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                              &
     /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: atl_poly1d_vert_w3_reconstruction_code
end type atl_poly1d_vert_w3_reconstruction_kernel_type

public :: atl_poly1d_vert_w3_reconstruction_init
public :: atl_poly1d_vert_w3_reconstruction_code
public :: atl_poly1d_vert_w3_reconstruction_final


!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!> @brief Computes the adjoint of tangent linear for vertical reconstructions.
!! @param[in]     nlayers         Number of layers
!! @param[in,out] reconstruction  ACTIVE Change in mass reconstruction field
!! @param[in,out] tracer          ACTIVE Change in tracer
!! @param[in]     ls_tracer       Linearisation state tracer tracer
!! @param[in]     coeff           Array of polynomial coefficients for interpolation
!! @param[in]     ndata           Number of data points per dof location
!! @param[in]     global_order    Desired order of polynomial reconstruction
!! @param[in]     logspace        If true perform interpolation in log space
!! @param[in]     ndf_md          Number of degrees of freedom per cell
!! @param[in]     undf_md         Number of unique degrees of freedom for the
!!                                reconstruction & wind fields
!! @param[in]     map_md          Dofmap for the cell at the base of the column
!! @param[in]     ndf_w3          Number of degrees of freedom per cell
!! @param[in]     undf_w3         Number of unique degrees of freedom for tracer
!! @param[in]     map_w3          Cell dofmaps for the tracer space
!! @param[in]     ndf_c           Number of degrees of freedom per cell for the
!!                                coeff space
!! @param[in]     undf_c          Total number of degrees of freedom for the
!!                                coeff space
!! @param[in]     map_c           Dofmap for the coeff space
!! @param[in]     stencil         Stencil array
subroutine atl_poly1d_vert_w3_reconstruction_code( nlayers,        &
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
                                                   map_c,          &
                                                   stencil )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                     :: nlayers
  integer(kind=i_def), intent(in)                     :: ndata
  integer(kind=i_def), intent(in)                     :: ndf_w3
  integer(kind=i_def), intent(in)                     :: undf_w3
  integer(kind=i_def), intent(in)                     :: ndf_md
  integer(kind=i_def), intent(in)                     :: undf_md
  integer(kind=i_def), intent(in)                     :: ndf_c
  integer(kind=i_def), intent(in)                     :: undf_c
  integer(kind=i_def), dimension(ndf_md), intent(in)  :: map_md
  integer(kind=i_def), dimension(ndf_c), intent(in)   :: map_c
  integer(kind=i_def), dimension(ndf_w3), intent(in)  :: map_w3
  integer(kind=i_def), intent(in)                     :: global_order
  real(kind=r_def), dimension(undf_md), intent(inout) :: reconstruction
  real(kind=r_def), dimension(undf_w3), intent(inout) :: tracer
  real(kind=r_def), dimension(undf_w3), intent(in)    :: ls_tracer
  real(kind=r_def), dimension(undf_c), intent(in)     :: coeff
  logical(kind=l_def), intent(in)                     :: logspace
  ! Added to remove necessity of global variable (which exists in tangent linear kernel)
  ! Necessitates using kernel via PSyKAl lite code
  integer(kind=i_def), intent(in), allocatable, dimension(:,:,:) :: stencil

  ! Internal variables
  integer(kind=i_def) :: k
  integer(kind=i_def) :: ij
  integer(kind=i_def) :: ik
  integer(kind=i_def) :: p
  integer(kind=i_def) :: f
  integer(kind=i_def) :: vertical_order
  real(kind=r_def)    :: new_tracer
  real(kind=r_def)    :: ls_new_tracer

  new_tracer = 0.0_r_def
  vertical_order = MIN(global_order, nlayers - 1)
  ij = map_w3(1)
  if (logspace) then
    do f = 1, 0, -1
      do k = nlayers - 1, 0, -1
        ls_new_tracer = 1.0_r_def
        do p = 1, vertical_order + 1
          ik = f * global_order + f + k * ndata + p + map_c(1) - 1
          ls_new_tracer = ls_new_tracer * MAX(eps, ABS(ls_tracer(ij + stencil(p,k,f)))) ** coeff(ik)
        end do
        new_tracer = new_tracer + ls_new_tracer * reconstruction(map_md(1) + f * nlayers + k)
        reconstruction(map_md(1) + f * nlayers + k) = 0.0_r_def
        do p = vertical_order + 1, 1, -1
          ik = f * global_order + f + k * ndata + p + map_c(1) - 1
          tracer(ij + stencil(p,k,f)) = tracer(ij + stencil(p,k,f)) + coeff(ik) * new_tracer / SIGN(MAX(eps, ls_tracer(ij + &
&stencil(p,k,f))), ls_tracer(ij + stencil(p,k,f)))
        enddo
        new_tracer = 0.0_r_def
      enddo
    enddo
  else
    do f = 1, 0, -1
      do k = nlayers - 1, 0, -1
        new_tracer = new_tracer + reconstruction(map_md(1) + f * nlayers + k)
        reconstruction(map_md(1) + f * nlayers + k) = 0.0_r_def
        do p = vertical_order + 1, 1, -1
          ik = f * global_order + f + k * ndata + p + map_c(1) - 1
          tracer(ij + stencil(p,k,f)) = tracer(ij + stencil(p,k,f)) + coeff(ik) * new_tracer
        enddo
        new_tracer = 0.0_r_def
      enddo
    enddo
  end if

end subroutine atl_poly1d_vert_w3_reconstruction_code

!> @brief Computes the offset stencil needed for vertical reconstructions.
!> @param[in]     global_order Desired order of reconstruction
!> @param[in]     nlayers      Number of layers in the mesh
!> @param[in,out] stencil      Stencil array
subroutine atl_poly1d_vert_w3_reconstruction_init(global_order, &
                                                  nlayers,      &
                                                  stencil)

  implicit none

  integer(kind=i_def), intent(in) :: global_order
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(inout), allocatable, dimension(:,:,:) :: stencil

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


end subroutine atl_poly1d_vert_w3_reconstruction_init

!> @brief Frees up the memory for the stencil array.
subroutine atl_poly1d_vert_w3_reconstruction_final(stencil)

  implicit none

  integer(kind=i_def), intent(inout), allocatable, dimension(:,:,:) :: stencil

  if ( allocated(stencil) ) deallocate(stencil)

end subroutine atl_poly1d_vert_w3_reconstruction_final

end module atl_poly1d_vert_w3_reconstruction_kernel_mod
