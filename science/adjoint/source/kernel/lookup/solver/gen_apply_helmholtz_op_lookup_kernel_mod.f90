!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes lookup table for apply_helmholtz_op_lookup_kernel

module gen_apply_helmholtz_op_lookup_kernel_mod

use argument_mod,          only : arg_type, func_type,         &
                                  GH_FIELD, GH_SCALAR,         &
                                  GH_REAL, GH_INTEGER,         &
                                  GH_READWRITE, GH_READ,       &
                                  STENCIL, CROSS2D,            &
                                  OWNED_AND_HALO_CELL_COLUMN,  &
                                  ANY_DISCONTINUOUS_SPACE_1,   &
                                  ANY_DISCONTINUOUS_SPACE_2
use constants_mod,         only : r_solver, i_def, l_def
use kernel_mod,            only : kernel_type
use fs_continuity_mod,     only : W3

implicit none

private

!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
! NOTE: Needs to be run serially, due to writing to stencil branches.
!       This is handled by the PSyclone optimisation script within the application.
!       The script looks for `gen_*_lookup_code` - if this naming scheme changes then
!       the script must be updated.
type, public, extends(kernel_type) :: gen_apply_helmholtz_op_lookup_kernel_type
type(arg_type) :: meta_args(6) = (/                                                            &
     arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3,                        STENCIL(CROSS2D)),  & ! dummy_w3
     arg_type(GH_FIELD,  GH_INTEGER, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),                    & ! lookup
     arg_type(GH_FIELD,  GH_INTEGER, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2),                    & ! set_counts
     arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_1, STENCIL(CROSS2D)),  & ! lookup_dummy
     arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_2, STENCIL(CROSS2D)),  & ! set_counts_dummy
     arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                                                     & ! nindices
     /)
integer :: operates_on = OWNED_AND_HALO_CELL_COLUMN
contains
  procedure, nopass :: gen_apply_helmholtz_op_lookup_code
end type gen_apply_helmholtz_op_lookup_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: gen_apply_helmholtz_op_lookup_code
contains

!> @brief Generate the adjoint lookup table for apply_helmholtz_operator_kernel.
!> @param[in]     nlayers      Number of vertical levels to solve over
!> @param[in]     dummy_w3     Dummy W3 field to retrieve dofmaps.
!> @param[in]     smap_w3_size Stencil sizes
!> @param[in]     smap_w3_max  Maximum stencil branch length
!> @param[in]     smap_w3      Stencil dofmap
!> @param[in,out] lookup       Adjoint lookup table.
!> @param[in,out] set_counts   Adjoint lookup table index set counters.
!> @param[in]     lookup_dummy Dummy lookup table to retrieve stencil dofmaps.
!> @param[in]     smap_lu_size Stencil sizes (lookup table)
!> @param[in]     smap_lu_max  Maximum stencil branch length (lookup table)
!> @param[in]     smap_lu      Stencil dofmap (lookup table)
!> @param[in]     set_counts_dummy Dummy set counts field to retrieve stencil dofmaps.
!> @param[in]     smap_sc_size Stencil sizes (set counts)
!> @param[in]     smap_sc_max  Maximum stencil branch length (set counts)
!> @param[in]     smap_sc      Stencil dofmap (set counts)
!> @param[in]     nindices     Number of indices per index set in lookup table.
!> @param[in]     ndf_w3       Number of dofs per cell W3
!> @param[in]     undf_w3      Size of field arrays W3
!> @param[in]     map_w3       W3 dofmap.
!> @param[in]     ndf_lu       Number of dofs per cell (lookup table)
!> @param[in]     undf_lu      Size of field arrays (lookup table)
!> @param[in]     map_lu       Lookup table dofmap.
!> @param[in]     ndf_sc       Number of dofs per cell (set counts)
!> @param[in]     undf_sc      Size of field arrays (set counts)
!> @param[in]     map_sc       Set counts dofmap.
subroutine gen_apply_helmholtz_op_lookup_code( nlayers,          &
                                               dummy_w3,         &
                                               smap_w3_size,     &
                                               smap_w3_max,      &
                                               smap_w3,          &
                                               lookup,           &
                                               set_counts,       &
                                               lookup_dummy,     &
                                               smap_lu_size,     &
                                               smap_lu_max,      &
                                               smap_lu,          &
                                               set_counts_dummy, &
                                               smap_sc_size,     &
                                               smap_sc_max,      &
                                               smap_sc,          &
                                               nindices,         &
                                               ndf_w3,           &
                                               undf_w3,          &
                                               map_w3,           &
                                               ndf_lu,           &
                                               undf_lu,          &
                                               map_lu,           &
                                               ndf_sc,           &
                                               undf_sc,          &
                                               map_sc )

  implicit none

  ! Arguments
  integer(kind=i_def),                                          intent(in) :: nlayers

  integer(kind=i_def),                                          intent(in) :: ndf_w3
  integer(kind=i_def),                                          intent(in) :: undf_w3
  integer(kind=i_def),                      dimension(ndf_w3),  intent(in) :: map_w3
  integer(kind=i_def),                                          intent(in) :: ndf_lu
  integer(kind=i_def),                                          intent(in) :: undf_lu
  integer(kind=i_def),                      dimension(ndf_lu),  intent(in) :: map_lu
  integer(kind=i_def),                                          intent(in) :: ndf_sc
  integer(kind=i_def),                                          intent(in) :: undf_sc
  integer(kind=i_def),                      dimension(ndf_sc),  intent(in) :: map_sc

  real(kind=r_solver),                     dimension(undf_w3),  intent(in) :: dummy_w3

  integer(kind=i_def),                   dimension(undf_lu), intent(inout) :: lookup
  integer(kind=i_def),                   dimension(undf_lu),    intent(in) :: lookup_dummy
  integer(kind=i_def),                   dimension(undf_sc), intent(inout) :: set_counts
  integer(kind=i_def),                   dimension(undf_sc),    intent(in) :: set_counts_dummy

  integer(kind=i_def), dimension(4),                            intent(in) :: smap_lu_size
  integer(kind=i_def),                                          intent(in) :: smap_lu_max
  integer(kind=i_def), dimension(ndf_lu,smap_lu_max,4),         intent(in) :: smap_lu
  integer(kind=i_def), dimension(4),                            intent(in) :: smap_sc_size
  integer(kind=i_def),                                          intent(in) :: smap_sc_max
  integer(kind=i_def), dimension(ndf_sc,smap_sc_max,4),         intent(in) :: smap_sc
  integer(kind=i_def),                                          intent(in) :: smap_w3_max
  integer(kind=i_def), dimension(4),                            intent(in) :: smap_w3_size
  integer(kind=i_def), dimension(ndf_w3,smap_w3_max,4),         intent(in) :: smap_w3

  integer(kind=i_def),                                          intent(in) :: nindices

  integer(kind=i_def), parameter :: nfaces = 4
  ! Internal variables
  integer(kind=i_def) :: f
  integer(kind=i_def) :: stencil
  integer(kind=i_def) :: stencil_sc
  integer(kind=i_def) :: sets_written
  integer(kind=i_def) :: idx_at_set

  do f = 1, nfaces
    stencil = smap_lu(1,2,f)
    stencil_sc = smap_sc(1,2,f)

    sets_written = set_counts(stencil_sc)
    idx_at_set = stencil + sets_written*nindices

    lookup(idx_at_set) = map_w3(1)
    lookup(idx_at_set + 1) = f

    set_counts(stencil_sc) = set_counts(stencil_sc) + 1_i_def
  end do

end subroutine gen_apply_helmholtz_op_lookup_code

end module gen_apply_helmholtz_op_lookup_kernel_mod
