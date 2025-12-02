!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes adjoint of horizontal tracer values
!>        through fitting a high order 2D upwind reconstruction.
!>        NOTE: Lookup table solution.
module adj_poly2d_recon_lookup_kernel_mod

use argument_mod,      only : arg_type, func_type,       &
                              GH_FIELD, GH_SCALAR,       &
                              GH_REAL, GH_INTEGER,       &
                              GH_READWRITE, GH_READ,     &
                              CELL_COLUMN,               &
                              ANY_DISCONTINUOUS_SPACE_1, &
                              ANY_DISCONTINUOUS_SPACE_2, &
                              ANY_DISCONTINUOUS_SPACE_3, &
                              ANY_DISCONTINUOUS_SPACE_4, &
                              ANY_DISCONTINUOUS_SPACE_5
use constants_mod,     only : r_tran, i_def, l_def
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
!! NOTE: ALWAYS call using the PSy-KAl lite invoke method.
!!       Problem with halo exchange - please see comments in the PSy-lite code.
type, public, extends(kernel_type) :: adj_poly2d_recon_lookup_kernel_type
  type(arg_type) :: meta_args(7) = (/ &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),                  & ! reconstruction
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2),                  & ! tracer
       arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_3),                  & ! lookup_poly2d
       arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_4),                  & ! set_count_poly2d
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_5),                  & ! coeff
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                                                  & ! nsets_max
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                                                   & ! nindices
       /)
  integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: adj_poly2d_recon_lookup_code
end type adj_poly2d_recon_lookup_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: adj_poly2d_recon_lookup_code

contains

!> @brief Computes the adjoint of horizontal polynomial interpolation of a tracer.
!> @param[in]     nlayers          Number of layers
!> @param[in]     reconstruction   Reconstructed tracer field to compute
!> @param[in,out] tracer           Pointwise tracer field to reconstruct
!> @param[in]     lookup_poly2d    Lookup table between base layer read and write indices
!> @param[in]     num_sets_poly2d  Number of sets of lookup indices per cell
!> @param[in]     coeff            Array of polynomial coefficients for interpolation
!> @param[in]     nsets_max        Maximum number of sets of write indices per read index
!> @param[in]     nindices         Number of write indices per set
!> @param[in]     ndf_md           Number of degrees of freedom per cell
!> @param[in]     undf_md          Number of unique degrees of freedom for the
!!                                 reconstructed field
!> @param[in]     map_md           Dofmap for the cell at the base of the column
!> @param[in]     ndf_ws           Number of degrees of freedom per cell
!> @param[in]     undf_ws          Number of unique degrees of freedom for the tracer field
!> @param[in]     map_ws           Dofmap for the cell at the base of the column for the tracer field
!> @param[in]     ndf_lu           Number of degrees of freedom per cell for the lookup table
!> @param[in]     undf_lu          Number of unique degrees of freedom for the lookup table
!> @param[in]     map_lu           Dofmap for the cell at the base of the column for the lookup table
!> @param[in]     ndf_ns           Number of degrees of freedom per cell for the number of index sets
!> @param[in]     undf_ns          Number of unique degrees of freedom for the number of index sets
!> @param[in]     map_ns           Dofmap for the cell at the base of the column for the number of index sets
!> @param[in]     ndf_c            Number of degrees of freedom per cell for the coeff space
!> @param[in]     undf_c           Total number of degrees of freedom for the coeff space
!> @param[in]     map_c            Dofmap for the coeff space
subroutine adj_poly2d_recon_lookup_code( nlayers,          &
                                         reconstruction,   &
                                         tracer,           &
                                         lookup_poly2d,    &
                                         set_count_poly2d,    &
                                         coeff,            &
                                         nsets_max,        &
                                         nindices,         &
                                         ndf_md,           &
                                         undf_md,          &
                                         map_md,           &
                                         ndf_ws,           &
                                         undf_ws,          &
                                         map_ws,           &
                                         ndf_lu,           &
                                         undf_lu,          &
                                         map_lu,           &
                                         ndf_sc,           &
                                         undf_sc,          &
                                         map_sc,           &
                                         ndf_c,            &
                                         undf_c,           &
                                         map_c )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                                 :: nlayers
  integer(kind=i_def), intent(in)                                 :: ndf_ws
  integer(kind=i_def), intent(in)                                 :: undf_ws
  integer(kind=i_def), dimension(ndf_ws), intent(in)              :: map_ws
  integer(kind=i_def), intent(in)                                 :: ndf_md
  integer(kind=i_def), intent(in)                                 :: undf_md
  integer(kind=i_def), dimension(ndf_md), intent(in)              :: map_md
  integer(kind=i_def), intent(in)                                 :: ndf_c
  integer(kind=i_def), intent(in)                                 :: undf_c
  integer(kind=i_def), dimension(ndf_c), intent(in)               :: map_c
  integer(kind=i_def), intent(in)                                 :: ndf_lu
  integer(kind=i_def), intent(in)                                 :: undf_lu
  integer(kind=i_def), dimension(ndf_lu), intent(in)              :: map_lu
  integer(kind=i_def), intent(in)                                 :: ndf_sc
  integer(kind=i_def), intent(in)                                 :: undf_sc
  integer(kind=i_def), dimension(ndf_sc), intent(in)              :: map_sc
  integer(kind=i_def), intent(in)                                 :: nsets_max
  integer(kind=i_def), intent(in)                                 :: nindices
  real(kind=r_tran), dimension(undf_md), intent(in)               :: reconstruction
  real(kind=r_tran), dimension(undf_ws), intent(inout)            :: tracer
  integer(kind=i_def), dimension(undf_lu), intent(in)             :: lookup_poly2d
  integer(kind=i_def), dimension(undf_sc), intent(in)             :: set_count_poly2d
  real(kind=r_tran), dimension(undf_c), intent(in)                :: coeff

  ! Internal variables
  integer(kind=i_def) :: k
  integer(kind=i_def) :: id_r
  integer(kind=i_def) :: id_t, id_t_key
  integer(kind=i_def) :: id_c
  integer(kind=i_def) :: nl
  integer(kind=i_def) :: set
  integer(kind=i_def) :: idx_at_set
  integer(kind=i_def) :: nsets

  nl = ndf_ws + nlayers - 2
  id_t = map_ws(1)
  id_t_key = map_lu(1)

  nsets = set_count_poly2d(map_sc(1))
  do set = 1, nsets
    ! Unpacking the indices from lookup table
    idx_at_set = id_t_key + (set-1)*nindices
    id_c = lookup_poly2d(idx_at_set)
    id_r = lookup_poly2d(idx_at_set + 1)

    do k = nl, 0, -1
      tracer(id_t + k) = tracer(id_t + k) + coeff(id_c) * reconstruction(id_r + k)
    end do
  end do

end subroutine adj_poly2d_recon_lookup_code

end module adj_poly2d_recon_lookup_kernel_mod
