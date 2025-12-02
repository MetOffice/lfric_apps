!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief "Kernel" which computes the adjoint lookup table for adj_poly1d_reconstruction_kernel.
module gen_poly1d_lookup_kernel_mod

use argument_mod,      only : arg_type, func_type,         &
                              GH_FIELD, GH_SCALAR,         &
                              GH_REAL, GH_INTEGER,         &
                              GH_WRITE, &
                              GH_READWRITE, GH_READ,       &
                              GH_REAL, &
                              STENCIL, CROSS,              &
                              OWNED_AND_HALO_CELL_COLUMN,  &
                              ANY_DISCONTINUOUS_SPACE_1,   &
                              ANY_DISCONTINUOUS_SPACE_2,   &
                              ANY_DISCONTINUOUS_SPACE_3,   &
                              ANY_DISCONTINUOUS_SPACE_4,   &
                              ANY_DISCONTINUOUS_SPACE_5
use constants_mod,     only : r_tran, i_def
use kernel_mod,        only : kernel_type

implicit none

private

!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
! NOTE: Needs to be run serially, due to writing to stencil branches.
!       This is handled by the PSyclone optimisation script within the application.
!       The script looks for `gen_*_lookup_code` - if this naming scheme changes then
!       the script must be updated.
type, public, extends(kernel_type) :: gen_poly1d_lookup_kernel_type
type(arg_type) :: meta_args(9) = (/ &
     arg_type(GH_FIELD,  GH_REAL,    GH_READ, ANY_DISCONTINUOUS_SPACE_1), &                 ! dummy_ws
     arg_type(GH_FIELD,  GH_REAL,    GH_READ, ANY_DISCONTINUOUS_SPACE_2), &                 ! dummy_coeff
     arg_type(GH_FIELD,  GH_REAL,    GH_READ, ANY_DISCONTINUOUS_SPACE_3), &                 ! dummy_md
     arg_type(GH_FIELD,  GH_INTEGER, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_4), &            ! lookup_poly1d
     arg_type(GH_FIELD,  GH_INTEGER, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_5), &            ! set_count_poly1d
     arg_type(GH_FIELD,  GH_INTEGER, GH_READ, ANY_DISCONTINUOUS_SPACE_4, STENCIL(CROSS)), & ! lookup_poly1d_dummy
     arg_type(GH_FIELD,  GH_INTEGER, GH_READ, ANY_DISCONTINUOUS_SPACE_5, STENCIL(CROSS)), & ! set_count_poly1d_dummy
     arg_type(GH_SCALAR, GH_INTEGER, GH_READ), &                                            ! order
     arg_type(GH_SCALAR, GH_INTEGER, GH_READ)         &                                     ! nindices
     /)
integer :: operates_on = OWNED_AND_HALO_CELL_COLUMN
contains
  procedure, nopass :: gen_poly1d_lookup_code
end type gen_poly1d_lookup_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: gen_poly1d_lookup_code

contains
!> @brief Computes the lookup table for the horizontal polynomial interpolation of a tracer.
!> @param[in]     nlayers                 Number of layers
!> @param[in]     dummy_ws                Dummy W-space field for retrieving dofmaps.
!> @param[in]     dummy_coeff             Dummy coeffient field, for retrieving dofmaps.
!> @param[in]     dummy_md                Dummy multidata field, for retrieving dofmaps.
!> @param[inout]  lookup_poly1d           Lookup table relating the 0th layer tracer indices
!!                                        to all sets of reconstruction and coefficient indices
!> @param[inout]  set_count_poly1d        Counters relating to the number of sets of indices per cell
!> @param[in]     lookup_poly1d_dummy     Dummy field to retrieve stencil.
!> @param[in]     lu_stencil_map          Dofmaps for the stencil
!> @param[in]     lu_stencil_size         Size of the stencil (number of cells)
!> @param[in]     set_count_poly1d_dummy  Dummy field to retrieve stencil.
!> @param[in]     sc_stencil_map          Dofmaps for the stencil
!> @param[in]     sc_stencil_size         Size of the stencil (number of cells)
!> @param[in]     order                   Desired order of polynomial reconstruction
!> @param[in]     nindices                Number of indices per set
!> @param[in]     ndf_ws                  Number of degrees of freedom per cell for the W-space field.
!> @param[in]     undf_ws                 Number of unique degrees of freedom for the W-space field.
!> @param[in]     map_ws                  Dofmap for the cell at the base of the W-space field.
!> @param[in]     ndf_lu                  Number of degrees of freedom per cell for the lookup table
!> @param[in]     undf_lu                 Number of unique degrees of freedom for the lookup table
!> @param[in]     map_lu                  Dofmap for the cell at the base of the column for the lookup table
!> @param[in]     ndf_sc                  Number of degrees of freedom per cell for the set counter
!> @param[in]     undf_sc                 Number of unique degrees of freedom for the set counter
!> @param[in]     map_sc                  Dofmap for the cell at the base of the column for the set counter
!> @param[in]     ndf_c                   Number of degrees of freedom per cell for the coeff space
!> @param[in]     undf_c                  Total number of degrees of freedom for the coeff space
!> @param[in]     map_c                   Dofmap for the coeff space
!> @param[in]     ndf_md                  Number of degrees of freedom per cell for the multidata space
!> @param[in]     undf_md                 Total number of degrees of freedom for the multidata space
!> @param[in]     map_md                  Dofmap for the multidata space
subroutine gen_poly1d_lookup_code( nlayers,                &
                                   dummy_ws,               &
                                   dummy_coeff,            &
                                   dummy_md,               &
                                   lookup_poly1d,          &
                                   set_count_poly1d,       &
                                   lookup_poly1d_dummy,    &
                                   lu_stencil_size,        &
                                   lu_stencil_map,         &
                                   set_count_poly1d_dummy, &
                                   sc_stencil_size,        &
                                   sc_stencil_map,         &
                                   order,                  &
                                   nindices,               &
                                   ndf_ws,                 &
                                   undf_ws,                &
                                   map_ws,                 &
                                   ndf_c,                  &
                                   undf_c,                 &
                                   map_c,                  &
                                   ndf_md,                 &
                                   undf_md,                &
                                   map_md,                 &
                                   ndf_lu,                 &
                                   undf_lu,                &
                                   map_lu,                 &
                                   ndf_sc,                 &
                                   undf_sc,                &
                                   map_sc )
  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                                     :: nlayers
  integer(kind=i_def), intent(in)                                     :: ndf_lu
  integer(kind=i_def), intent(in)                                     :: undf_lu
  integer(kind=i_def), dimension(ndf_lu), intent(in)                  :: map_lu
  integer(kind=i_def), intent(in)                                     :: ndf_sc
  integer(kind=i_def), intent(in)                                     :: undf_sc
  integer(kind=i_def), dimension(ndf_sc), intent(in)                  :: map_sc
  integer(kind=i_def), intent(in)                                     :: ndf_c
  integer(kind=i_def), intent(in)                                     :: undf_c
  integer(kind=i_def), dimension(ndf_c), intent(in)                   :: map_c
  integer(kind=i_def), intent(in)                                     :: ndf_ws
  integer(kind=i_def), intent(in)                                     :: undf_ws
  integer(kind=i_def), dimension(ndf_ws), intent(in)                  :: map_ws
  integer(kind=i_def), intent(in)                                     :: ndf_md
  integer(kind=i_def), intent(in)                                     :: undf_md
  integer(kind=i_def), dimension(ndf_md), intent(in)                  :: map_md

  integer(kind=i_def), intent(in)                                     :: lu_stencil_size
  integer(kind=i_def), dimension(ndf_lu, lu_stencil_size), intent(in) :: lu_stencil_map
  integer(kind=i_def), intent(in)                                     :: sc_stencil_size
  integer(kind=i_def), dimension(ndf_sc, sc_stencil_size), intent(in) :: sc_stencil_map
  integer(kind=i_def), intent(in)                                     :: order
  integer(kind=i_def), intent(in)                                     :: nindices

  real(kind=r_tran), dimension(undf_c), intent(in)                    :: dummy_coeff
  real(kind=r_tran), dimension(undf_ws), intent(in)                   :: dummy_ws
  real(kind=r_tran), dimension(undf_md), intent(in)                   :: dummy_md
  integer(kind=i_def), dimension(undf_lu), intent(inout)              :: lookup_poly1d
  integer(kind=i_def), dimension(undf_sc), intent(inout)              :: set_count_poly1d
  integer(kind=i_def), dimension(undf_lu), intent(in)                 :: lookup_poly1d_dummy
  integer(kind=i_def), dimension(undf_sc), intent(in)                 :: set_count_poly1d_dummy

  ! Internal variables
  integer(kind=i_def), parameter                      :: nfaces = 4
  integer(kind=i_def)                                 :: f
  integer(kind=i_def)                                 :: p
  integer(kind=i_def)                                 :: face
  integer(kind=i_def)                                 :: stencil, stencil_sc
  integer(kind=i_def)                                 :: stencil_depth
  integer(kind=i_def)                                 :: depth
  integer(kind=i_def)                                 :: face_mod
  integer(kind=i_def)                                 :: ijp
  integer(kind=i_def)                                 :: df
  integer(kind=i_def)                                 :: nl
  integer(kind=i_def)                                 :: idx_at_set
  integer(kind=i_def)                                 :: sets_written
  integer(kind = i_def), dimension(order + 1, nfaces) :: map1d

  nl = ndf_ws + nlayers - 2
  stencil_depth = order / 2
  map1d(1,:) = 1
  do face = 1, nfaces, 1
    depth = 1
    face_mod = stencil_depth * MOD(face + 1, 2)
    do stencil = 2, stencil_depth + 1, 1
      map1d(stencil + depth - 1,face) = face_mod + stencil
      map1d(stencil + depth,face) = face_mod + order + stencil
      depth = depth + 1
    end do
  end do

  do f = 1, nfaces
    df = f * nl + f - nl + map_md(1) - 1
    do p = 1, order + 1
      ijp = f * order + f - order + p + map_c(1) - 2

      stencil = lu_stencil_map(1, map1d(p,f))
      stencil_sc = sc_stencil_map(1, map1d(p,f))

      ! Fetch the number of index sets currently written to this cell.
      sets_written = set_count_poly1d(stencil_sc)
      idx_at_set = stencil + sets_written*nindices

      lookup_poly1d(idx_at_set) = ijp
      lookup_poly1d(idx_at_set + 1) = df

      ! Increment number of index sets written by 1.
      set_count_poly1d(stencil_sc) = set_count_poly1d(stencil_sc) + 1
    end do
  end do

end subroutine gen_poly1d_lookup_code

end module gen_poly1d_lookup_kernel_mod
