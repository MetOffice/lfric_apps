!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which applies zeroing of reconstruction field
!!        in adj_poly[1,2]d_reconstruction_kernel. Needs to be done
!!        separately to the rest of the adjoint to enable parallel
!!        computation.
module adj_polynd_recon_zeroing_kernel_mod

use argument_mod,      only : arg_type, func_type,       &
                              GH_FIELD, GH_REAL,         &
                              GH_READWRITE, GH_READ,     &
                              CELL_COLUMN,               &
                              ANY_DISCONTINUOUS_SPACE_1, &
                              ANY_DISCONTINUOUS_SPACE_2
use constants_mod,     only : r_tran, i_def
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: adj_polynd_recon_zeroing_kernel_type
  type(arg_type) :: meta_args(2) = (/ &
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),                  & ! reconstruction
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2)                   & ! dummy_ads2
       /)
  integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: adj_polynd_recon_zeroing_code
end type adj_polynd_recon_zeroing_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: adj_polynd_recon_zeroing_code

contains

!> @brief Computes the adjoint of horizontal polynomial interpolation of a tracer.
!> @param[in]     nlayers          Number of layers
!> @param[in,out] reconstruction   Reconstructed tracer field to compute
!> @param[in]     dummy_adcsp2     Dummy tracer field for stencil
!> @param[in]     ndf_md           Number of degrees of freedom per cell
!> @param[in]     undf_md          Number of unique degrees of freedom for the
!!                                 reconstructed field
!> @param[in]     map_md           Dofmap for the cell at the base of the column
!> @param[in]     ndf_ws           Number of degrees of freedom per cell
!> @param[in]     undf_ws          Number of unique degrees of freedom for the tracer field
!> @param[in]     map_ws           Dofmap for the cell at the base of the column for the tracer field
subroutine adj_polynd_recon_zeroing_code( nlayers,          &
                                          reconstruction,   &
                                          dummy_ads2,       &
                                          ndf_md,           &
                                          undf_md,          &
                                          map_md,           &
                                          ndf_ws,           &
                                          undf_ws,          &
                                          map_ws )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                                 :: nlayers
  integer(kind=i_def), intent(in)                                 :: ndf_ws
  integer(kind=i_def), intent(in)                                 :: undf_ws
  integer(kind=i_def), dimension(ndf_ws), intent(in)              :: map_ws
  integer(kind=i_def), intent(in)                                 :: ndf_md
  integer(kind=i_def), intent(in)                                 :: undf_md
  integer(kind=i_def), dimension(ndf_md), intent(in)              :: map_md
  real(kind=r_tran), dimension(undf_md), intent(inout)            :: reconstruction
  real(kind=r_tran), dimension(undf_ws), intent(in)               :: dummy_ads2

  ! Internal variables
  integer(kind=i_def), parameter :: nfaces = 4
  integer(kind=i_def)            :: k
  integer(kind=i_def)            :: f
  integer(kind=i_def)            :: id_r
  integer(kind=i_def)            :: nl

  nl = ndf_ws + nlayers - 2
  do f = nfaces, 1, -1
    id_r = f * nl + f - nl + map_md(1) - 1
    do k = nl, 0, -1
      reconstruction(id_r + k) = 0.0_r_tran
    end do
  end do

end subroutine adj_polynd_recon_zeroing_code

end module adj_polynd_recon_zeroing_kernel_mod
