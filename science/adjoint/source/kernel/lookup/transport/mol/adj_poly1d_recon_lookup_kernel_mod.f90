!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes adjoint of horizontal tracer values
!>        through fitting a high order 1D upwind reconstruction.
!>        NOTE: Lookup table solution.
module adj_poly1d_recon_lookup_kernel_mod

use argument_mod,      only : arg_type, func_type,         &
                              reference_element_data_type, &
                              GH_FIELD, GH_SCALAR,         &
                              GH_REAL, GH_INTEGER,         &
                              GH_READWRITE, GH_READ,       &
                              CELL_COLUMN,                 &
                              ANY_DISCONTINUOUS_SPACE_1,   &
                              ANY_DISCONTINUOUS_SPACE_2,   &
                              ANY_DISCONTINUOUS_SPACE_3,   &
                              ANY_DISCONTINUOUS_SPACE_4,   &
                              ANY_DISCONTINUOUS_SPACE_5
use constants_mod,     only : r_tran, i_def
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
!! NOTE: ALWAYS call using the PSy-KAl lite invoke method.
!!       Problem with halo exchange - please see comments in the PSy-lite code.
type, public, extends(kernel_type) :: adj_poly1d_recon_lookup_kernel_type
type(arg_type) :: meta_args(7) = (/                                                            &
     arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),                 & ! reconstruction
     arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2),                 & ! tracer
     arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_3),                 & ! lookup_poly1d
     arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_4),                 & ! num_sets_poly1d
     arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_5),                 & ! coeff
     arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                                                 & ! nsets_max
     arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                                                  & ! nindices
     /)
integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: adj_poly1d_recon_lookup_code
end type adj_poly1d_recon_lookup_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: adj_poly1d_recon_lookup_code

contains
!> @brief Computes the horizontal polynomial interpolation of a tracer.
!> @param[in]     nlayers          Number of layers
!> @param[in]     reconstruction   Reconstructed tracer field to compute
!> @param[in,out] tracer           Pointwise tracer field to reconstruct
!> @param[in]     lookup_poly1d    Lookup table between horizontal read and write indices
!> @param[in]     num_sets_poly1d  Number of sets of lookup indices per cell
!> @param[in]     coeff            Array of polynomial coefficients for interpolation
!> @param[in]     nsets_max        Maximum of sets of write indices per read index for the lookup table
!> @param[in]     nindices         Number of write indices per set for the lookup table
!> @param[in]     ndf_md           Number of degrees of freedom per cell for reconstructed field
!> @param[in]     undf_md          Number of unique degrees of freedom for the
!!                                 reconstructed field
!> @param[in]     map_md           Dofmap for the cell at the base of the column for reconstructed field
!> @param[in]     ndf_ws           Number of degrees of freedom per cell for the tracer
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

subroutine adj_poly1d_recon_lookup_code( nlayers, &
                                         reconstruction, &
                                         tracer, &
                                         lookup_poly1d, &
                                         num_sets_poly1d, &
                                         coeff, &
                                         nsets_max, &
                                         nindices, &
                                         ndf_md, &
                                         undf_md, &
                                         map_md, &
                                         ndf_ws, &
                                         undf_ws, &
                                         map_ws, &
                                         ndf_lu, &
                                         undf_lu, &
                                         map_lu, &
                                         ndf_ns, &
                                         undf_ns, &
                                         map_ns, &
                                         ndf_c, &
                                         undf_c, &
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
  integer(kind=i_def), intent(in)                                 :: ndf_ns
  integer(kind=i_def), intent(in)                                 :: undf_ns
  integer(kind=i_def), dimension(ndf_ns), intent(in)              :: map_ns
  integer(kind=i_def), intent(in)                                 :: nsets_max
  integer(kind=i_def), intent(in)                                 :: nindices
  real(kind=r_tran), dimension(undf_md), intent(in)               :: reconstruction
  real(kind=r_tran), dimension(undf_ws), intent(inout)            :: tracer
  integer(kind=i_def), dimension(undf_lu), intent(in)             :: lookup_poly1d
  integer(kind=i_def), dimension(undf_ns), intent(in)             :: num_sets_poly1d
  real(kind=r_tran), dimension(undf_c), intent(in)                :: coeff

  ! Internal variables
  integer(kind=i_def) :: k
  integer(kind=i_def) :: ijp
  integer(kind=i_def) :: df
  integer(kind=i_def) :: df_t, df_t_key
  integer(kind=i_def) :: nl
  integer(kind=i_def) :: set
  integer(kind=i_def) :: idx_at_set
  integer(kind=i_def) :: nsets

  nl = ndf_ws + nlayers - 2
  df_t = map_ws(1)
  df_t_key = map_lu(1)
  ! Retrieve the number of sets within this cell.
  ! Stored in the last index of the local lookup table array.
  nsets = num_sets_poly1d(map_ns(1))

  do set = 1, nsets
    ! Unpacking the indices from lookup table
    idx_at_set = df_t_key + (set - 1)*nindices
    ijp = lookup_poly1d(idx_at_set)
    df = lookup_poly1d(idx_at_set + 1)

    do k = nl, 0, -1
      tracer(k + df_t) = tracer(k + df_t) + coeff(ijp) * reconstruction(df + k)
    end do
  end do

end subroutine adj_poly1d_recon_lookup_code

end module adj_poly1d_recon_lookup_kernel_mod
