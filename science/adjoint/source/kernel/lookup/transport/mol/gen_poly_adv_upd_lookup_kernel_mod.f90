!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief "Kernel" which computes lookup table for horizontal advective update
!!        through fitting a high order upwind reconstruction.


module gen_poly_adv_upd_lookup_kernel_mod

use argument_mod,      only : arg_type, func_type,         &
                              GH_FIELD, GH_SCALAR,         &
                              GH_REAL, GH_INTEGER,         &
                              GH_READWRITE, GH_READ,       &
                              STENCIL, CROSS2D,              &
                              OWNED_AND_HALO_CELL_COLUMN,                 &
                              ANY_DISCONTINUOUS_SPACE_1,   &
                              ANY_DISCONTINUOUS_SPACE_2,   &
                              ANY_DISCONTINUOUS_SPACE_3
use fs_continuity_mod, only : W2, Wtheta
use constants_mod,     only : r_tran, i_def, l_def
use kernel_mod,        only : kernel_type

implicit none

private

!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
! NOTE: Needs to be run serially, due to writing to stencil branches.
!       This is handled by the PSyclone optimisation script within the application.
!       The script looks for `gen_*_lookup_code` - if this naming scheme changes then
!       the script must be updated.
type, public, extends(kernel_type) :: gen_poly_adv_upd_lookup_kernel_type
type(arg_type) :: meta_args(9) = (/                                                            &
     arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta), & ! advective_dummy
     arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1, STENCIL(CROSS2D)), & ! reconstruction_dummy
     arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W2, STENCIL(CROSS2D)), & ! wind_dummy
     arg_type(GH_FIELD,  GH_INTEGER, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & ! lookup_poly_adv_upd
     arg_type(GH_FIELD,  GH_INTEGER, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_3), & ! set_count_poly_adv_upd
     arg_type(GH_FIELD,  GH_INTEGER, GH_READ, ANY_DISCONTINUOUS_SPACE_2, STENCIL(CROSS2D)), & ! lookup_poly_adv_upd_dummy
     arg_type(GH_FIELD,  GH_INTEGER, GH_READ, ANY_DISCONTINUOUS_SPACE_3, STENCIL(CROSS2D)), & ! set_count_poly_adv_upd_dummy
     arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                                              & ! nsets_max
     arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                                               & ! nindices
     /)
integer :: operates_on = OWNED_AND_HALO_CELL_COLUMN
contains
  procedure, nopass :: gen_poly_adv_upd_lookup_code
end type gen_poly_adv_upd_lookup_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: gen_poly_adv_upd_lookup_code
contains

!> @brief Computes the lookup table for poly_adv_update
!> @param[in]     nlayers                       Number of layers
!> @param[in]     advective_dummy               Dummy field to retrieve dofmap.
!> @param[in]     reconstruction_dummy          Dummy field to retrieve dofmap.
!> @param[in]     smap_md_size                  Multidata stencil size.
!> @param[in]     smap_md_size                  Multidata stencil max branch size.
!> @param[in]     smap_md                       Multidata stencil.
!> @param[in]     wind_dummy                    Dummy field to retrieve dofmap.
!> @param[in]     smap_w2_size                  W2 stencil size.
!> @param[in]     smap_w2_size                  W2 stencil max branch size.
!> @param[in]     smap_w2                       W2 stencil.
!> @param[inout]  lookup_poly_adv_upd           Lookup table relating the 0th layer tracer indices
!!                                              to all sets of advective increment, wind and m3_inv indices
!> @param[inout]  set_count_poly_adv_upd        Counters relating to how many index sets there are
!!                                              per cell in the lookup table.
!> @param[in]     lookup_poly_adv_upd_dummy     Dummy field for stencil dofmap.
!> @param[in]     smap_lu_size                  Lookup stencil size.
!> @param[in]     smap_lu_size                  Lookup stencil max branch size.
!> @param[in]     smap_lu                       Lookup stencil.
!> @param[in]     set_count_poly_adv_upd_dummy  Dummy field for stencil dofmap.
!> @param[in]     smap_sc_size                  Set count stencil size.
!> @param[in]     smap_sc_size                  Set count stencil max branch size.
!> @param[in]     smap_sc                       Set count stencil.
!> @param[in]     nsets_max                     Max number of sets of indices
!> @param[in]     nindices                      Number of indices per set
!> @param[in]     ndf_wt                        Number of degrees of freedom per cell (Wtheta)
!> @param[in]     undf_wt                       Number of unique degrees of freedom (Wtheta)
!> @param[in]     map_wt                        Dofmap for Wtheta.
!> @param[in]     ndf_md                        Number of degrees of freedom per cell (multidata)
!> @param[in]     undf_md                       Number of unique degrees of freedom (multidata)
!> @param[in]     map_md                        Dofmap for multidata.
!> @param[in]     ndf_w2                        Number of degrees of freedom per cell W2
!> @param[in]     undf_w2                       Number of unique degrees of freedom (W2)
!> @param[in]     map_w2                        Dofmap for W2.
!> @param[in]     ndf_lu                        Number of degrees of freedom per cell (lookup table)
!> @param[in]     undf_lu                       Number of unique degrees of freedom (lookup table)
!> @param[in]     map_lu                        Dofmap for lookup space.
!> @param[in]     ndf_sc                        Number of degrees of freedom per cell (set counter)
!> @param[in]     undf_sc                       Number of unique degrees of freedom (set counter)
!> @param[in]     map_sc                        Dofmap for set_count space.
subroutine gen_poly_adv_upd_lookup_code( nlayers,                &
                                         advective_dummy,        &
                                         reconstruction_dummy,   &
                                         smap_md_size, &
                                         smap_md_max, &
                                         smap_md, &
                                         wind_dummy, &
                                         smap_w2_size, &
                                         smap_w2_max, &
                                         smap_w2, &
                                         lookup_poly_adv_upd,    &
                                         set_count_poly_adv_upd, &
                                         lookup_poly_adv_upd_dummy, &
                                         smap_lu_size, &
                                         smap_lu_max,            &
                                         smap_lu,                &
                                         set_count_poly_adv_upd_dummy, &
                                         smap_sc_size,           &
                                         smap_sc_max,            &
                                         smap_sc,                &
                                         nsets_max, &
                                         nindices, &
                                         ndf_wt, &
                                         undf_wt, &
                                         map_wt, &
                                         ndf_md, &
                                         undf_md, &
                                         map_md, &
                                         ndf_w2, &
                                         undf_w2, &
                                         map_w2, &
                                         ndf_lu, &
                                         undf_lu, &
                                         map_lu, &
                                         ndf_sc, &
                                         undf_sc, &
                                         map_sc )


  implicit none

  ! Arguments
  integer(kind=i_def),                                          intent(in) :: nlayers

  integer(kind=i_def),                                          intent(in) :: ndf_wt
  integer(kind=i_def),                                          intent(in) :: undf_wt
  integer(kind=i_def),                      dimension(ndf_wt),  intent(in) :: map_wt
  integer(kind=i_def),                                          intent(in) :: ndf_md
  integer(kind=i_def),                                          intent(in) :: undf_md
  integer(kind=i_def),                      dimension(ndf_md),  intent(in) :: map_md
  integer(kind=i_def),                                          intent(in) :: ndf_w2
  integer(kind=i_def),                                          intent(in) :: undf_w2
  integer(kind=i_def),                      dimension(ndf_w2),  intent(in) :: map_w2
  integer(kind=i_def),                                          intent(in) :: ndf_lu
  integer(kind=i_def),                                          intent(in) :: undf_lu
  integer(kind=i_def),                      dimension(ndf_lu),  intent(in) :: map_lu
  integer(kind=i_def),                                          intent(in) :: ndf_sc
  integer(kind=i_def),                                          intent(in) :: undf_sc
  integer(kind=i_def),                      dimension(ndf_sc),  intent(in) :: map_sc

  real(kind=r_tran),                     dimension(undf_wt),    intent(in) :: advective_dummy
  real(kind=r_tran),                     dimension(undf_md),    intent(in) :: reconstruction_dummy
  real(kind=r_tran),                     dimension(undf_w2),    intent(in) :: wind_dummy

  integer(kind=i_def),                   dimension(undf_lu), intent(inout) :: lookup_poly_adv_upd
  integer(kind=i_def),                   dimension(undf_sc), intent(inout) :: set_count_poly_adv_upd
  integer(kind=i_def),                   dimension(undf_lu),    intent(in) :: lookup_poly_adv_upd_dummy
  integer(kind=i_def),                   dimension(undf_sc),    intent(in) :: set_count_poly_adv_upd_dummy

  integer(kind=i_def),                                          intent(in) :: smap_md_max
  integer(kind=i_def), dimension(4),                            intent(in) :: smap_md_size
  integer(kind=i_def), dimension(ndf_md,smap_md_max,4),         intent(in) :: smap_md
  integer(kind=i_def),                                          intent(in) :: smap_w2_max
  integer(kind=i_def), dimension(4),                            intent(in) :: smap_w2_size
  integer(kind=i_def), dimension(ndf_w2,smap_w2_max,4),         intent(in) :: smap_w2

  integer(kind=i_def),                                          intent(in) :: smap_lu_max
  integer(kind=i_def), dimension(4),                            intent(in) :: smap_lu_size
  integer(kind=i_def), dimension(ndf_lu,smap_lu_max,4),         intent(in) :: smap_lu
  integer(kind=i_def),                                          intent(in) :: smap_sc_max
  integer(kind=i_def), dimension(4),                            intent(in) :: smap_sc_size
  integer(kind=i_def), dimension(ndf_sc,smap_sc_max,4),         intent(in) :: smap_sc

  integer(kind=i_def),                                          intent(in) :: nsets_max
  integer(kind=i_def),                                          intent(in) :: nindices

  ! Internal variables
  integer(kind=i_def) :: face, df, df1, df2, idx_at_set
  integer(kind=i_def) :: sets_written, sc_idx, root_idx

  integer(kind=i_def), parameter         :: nfaces = 4
  integer(kind=i_def), dimension(nfaces) :: opposite
  logical(kind=l_def), dimension(nfaces) :: missing_neighbour

  ! For each face of cell, find the index in the neighbouring cell that
  ! corresponds to it.
  ! i.e for no orientation changes opposite = ( 3, 4, 1, 2 )
  ! We use the W2 map to determine these
  ! If there is no neighbour then we ensure the opposite points to
  ! the value on this edge
  opposite = -1
  missing_neighbour = .false.
  do df = 1,nfaces
    df1 = map_w2(df)
    if ( smap_w2_size(df) > 1 ) then
      ! There is a neighbour in direction df so find the
      ! neighboring edge corresponding to edge df
      do df2 = 1, nfaces
        if ( smap_w2(df2,2,df) == df1 ) opposite(df) = df2
      end do
    else
      ! There is no neighbour in direction df so point to itself
      opposite(df) = df
      missing_neighbour(df) = .true.
    end if
  end do

  do face = 1, nfaces
    ! Fetch the number of index sets currently written.
    sc_idx = smap_sc(1,2,face) + opposite(face) - 1
    sets_written = set_count_poly_adv_upd(sc_idx)

    root_idx = smap_lu(1,2,face) + (opposite(face)-1)*nsets_max*nindices
    idx_at_set = root_idx + sets_written*nindices

    lookup_poly_adv_upd(idx_at_set : idx_at_set + nfaces - 1) = map_w2(:)
    lookup_poly_adv_upd(idx_at_set + nfaces) = face

    lookup_poly_adv_upd(idx_at_set + nfaces + 1) = 0.0_i_def
    if (missing_neighbour(face)) lookup_poly_adv_upd(idx_at_set + nfaces + 1) = 1.0_i_def

    lookup_poly_adv_upd(idx_at_set + nfaces + 2) = map_wt(1)

    ! Increase access counter by 1.
    set_count_poly_adv_upd(sc_idx) = set_count_poly_adv_upd(sc_idx) + 1
  end do

end subroutine gen_poly_adv_upd_lookup_code

end module gen_poly_adv_upd_lookup_kernel_mod
