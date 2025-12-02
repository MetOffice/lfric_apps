!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!
!> @brief Kernel which computes horizontal cell edges value based on the monotone
!!        Koren scheme.
!> @details The kernel computes the edge values at W1-points from a Wtheta-tracer field.
!!  Ref.   Koren, B. (1993), "A robust upwind discretisation method for advection,
!!         diffusion and source terms", in Vreugdenhil, C.B.; Koren, B. (eds.),
!!         Numerical Methods for Advection/Diffusion Problems, Braunschweig:
!!         Vieweg, p. 117, ISBN 3-528-07645-3.
module polyh_wtheta_koren_kernel_mod

use argument_mod,      only : arg_type, func_type,         &
                              GH_FIELD, GH_REAL,           &
                              GH_READWRITE, GH_READ,       &
                              STENCIL, CROSS,              &
                              CELL_COLUMN,                 &
                              ANY_DISCONTINUOUS_SPACE_1
use constants_mod,     only : r_tran, i_def, l_def, SMALL_R_TRAN, r_tran, EPS_R_TRAN
use fs_continuity_mod, only : Wtheta
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: polyh_wtheta_koren_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                                            &
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta, STENCIL(CROSS))     &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: polyh_wtheta_koren_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: polyh_wtheta_koren_code

contains

!> @brief Computes the horizontal polynomial interpolation of a tracer.
!> @param[in]     nlayers        Number of layers
!> @param[in,out] reconstruction Reconstructed tracer field to compute
!> @param[in]     tracer         Pointwise tracer field to reconstruct
!> @param[in]     stencil_size   Size of the stencil (number of cells)
!> @param[in]     stencil_map    Dofmaps for the stencil
!> @param[in]     ndf_md         Number of degrees of freedom per cell
!> @param[in]     undf_md        Number of unique degrees of freedom for the
!!                               reconstructed field
!> @param[in]     map_md         Dofmap for the cell at the base of the column
!> @param[in]     ndf_wt         Number of degrees of freedom per cell
!> @param[in]     undf_wt        Number of unique degrees of freedom for the tracer field
!> @param[in]     map_wt         Dofmap for the cell at the base of the column for the tracer field
subroutine polyh_wtheta_koren_code(  nlayers,              &
                                     reconstruction,       &
                                     tracer,               &
                                     stencil_size,         &
                                     stencil_map,          &
                                     ndf_md,               &
                                     undf_md,              &
                                     map_md,               &
                                     ndf_wt,               &
                                     undf_wt,              &
                                     map_wt )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                     :: nlayers
  integer(kind=i_def), intent(in)                     :: ndf_wt
  integer(kind=i_def), intent(in)                     :: undf_wt
  integer(kind=i_def), dimension(ndf_wt),  intent(in) :: map_wt
  integer(kind=i_def), intent(in)                     :: ndf_md
  integer(kind=i_def), intent(in)                     :: undf_md
  integer(kind=i_def), dimension(ndf_md),  intent(in) :: map_md
  integer(kind=i_def), intent(in)                     :: stencil_size

  real(kind=r_tran), dimension(undf_md),  intent(inout) :: reconstruction
  real(kind=r_tran), dimension(undf_wt),  intent(in)    :: tracer
  integer(kind=i_def), dimension(ndf_wt,stencil_size), intent(in) :: stencil_map

  ! Internal variables
  integer(kind=i_def), parameter      :: nfaces = 4
  integer(kind=i_def), parameter      :: default_stencil_size = 5
  integer(kind=i_def)                 :: k, df
  real(kind=r_tran)                   :: edge_tracer
  integer(kind=i_def), dimension(3,4) :: point=reshape((/4,1,2,5,1,3,2,1,4,3,1,5/),shape(point))
  real(kind=r_tran)                   :: x, y, r, phi, r1, r2

  ! for order = 2 the cross stencil map is
  !      | 5 |
  !  | 2 | 1 | 4 |
  !      | 3 |
  !if stencil < default_stencil_size use constant reconstrution
  if ( stencil_size < default_stencil_size ) then
    do df = 1, nfaces
      do k = 0, nlayers
        reconstruction(map_md(1) + (df-1)*(nlayers+1) + k ) = tracer(stencil_map(1,1)+k)
      end do
    end do
  else
    do df = 1, nfaces
      do k = 0, nlayers
        x = tracer(stencil_map(1,point(2,df))+k) - tracer(stencil_map(1,point(1,df))+k)
        y = tracer(stencil_map(1,point(3,df))+k) - tracer(stencil_map(1,point(2,df))+k)
        r = (y + SMALL_R_TRAN)/(x + SMALL_R_TRAN)
        r1 = 2.0_r_tran*r
        r2 = ( 1.0_r_tran + r1 )/ 3.0_r_tran
        phi = max (0.0_r_tran, min(r1,r2,2.0_r_tran))
        edge_tracer = tracer(stencil_map(1,point(2,df))+k) + 0.5_r_tran*phi*x
        reconstruction(map_md(1) + (df-1)*(nlayers+1) + k ) = edge_tracer
      end do
    end do
  end if

end subroutine polyh_wtheta_koren_code

end module polyh_wtheta_koren_kernel_mod
