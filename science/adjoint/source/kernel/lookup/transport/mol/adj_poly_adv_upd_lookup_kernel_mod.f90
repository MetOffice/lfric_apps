!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief "Kernel" which computes lookup table for adj_poly_advective_update_kernel.
module adj_poly_adv_upd_lookup_kernel_mod

use argument_mod,      only : arg_type, func_type,   &
                              GH_FIELD, GH_SCALAR,   &
                              GH_READ, GH_READWRITE, &
                              GH_INTEGER, &
                              GH_REAL, GH_WRITE,     &
                              STENCIL, CROSS2D,      &
                              CELL_COLUMN,           &
                              ANY_DISCONTINUOUS_SPACE_1, &
                              ANY_DISCONTINUOUS_SPACE_2, &
                              ANY_DISCONTINUOUS_SPACE_3
use constants_mod,     only : i_def, l_def, r_tran
use fs_continuity_mod, only : W2, Wtheta
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
!! NOTE: ALWAYS call using the PSy-KAl lite invoke method.
!!       Problem with halo exchange - please see comments in the PSy-lite code.
type, public, extends(kernel_type) :: adj_poly_adv_upd_lookup_kernel_type
  type(arg_type) :: meta_args(8) = (/ &
    arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta),                                 & ! advective
    arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),              & ! reconstruction
    arg_type(GH_FIELD, GH_INTEGER, GH_READ,   ANY_DISCONTINUOUS_SPACE_2),              & ! lookup_poly_adv_upd
    arg_type(GH_FIELD, GH_INTEGER, GH_READ,   ANY_DISCONTINUOUS_SPACE_3),              & ! set_count_poly_adv_upd
    arg_type(GH_FIELD, GH_REAL, GH_READ, W2, STENCIL(CROSS2D)),                        & ! wind
    arg_type(GH_FIELD, GH_REAL, GH_READ, W2),                                          & ! wind_dir
    arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                                          & ! nsets_max
    arg_type(GH_SCALAR, GH_INTEGER, GH_READ)/)                                           ! nindices
  integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: adj_poly_adv_upd_lookup_code
end type adj_poly_adv_upd_lookup_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: adj_poly_adv_upd_lookup_code

contains

!> @brief Computes the adjoint of horizontal advective update of a field.
!!        NOTE: Lookup table solution.
!> @param[in]     nlayers                Number of layers
!> @param[in,out] advective              Field containing the horizontal advective
!!                                       update: (u,v).grad_h(reconstruction)
!> @param[in,out] reconstruction         Multidata field containing the edge
!!                                       reconstruction for each cell if it is
!!                                       the upwind cell for that edge
!> @param[in]     lookup_poly_adv_upd    Lookup table between horizontal read and write indices
!> @param[in]     set_count_poly_adv_upd Number of sets of lookup set_count_w3h_adv_up per cell
!> @param[in]     dummy_md               Read only multidata field to get stencil
!> @param[in]     wind                   Passive linearisation state wind field
!> @param[in]     smap_w2_size           Size of the w2stencil map in each direction
!> @param[in]     smap_w2_max            Maximum size of the w2 stencil map
!> @param[in]     smap_w2                Stencil map for the w2 fields
!> @param[in]     wind_dir               Wind field used to determine direction,
!!                                       equal to wind when used in gungho
!!                                       but ls_wind when used in the linear model
!> @param[in]     ndf_wt                 Number of degrees of freedom per cell
!> @param[in]     undf_wt                Number of unique degrees of freedom for the advective field
!> @param[in]     map_wt                 Dofmap for the cell at the base of the column
!> @param[in]     ndf_md                 Number of degrees of freedom per cell
!> @param[in]     undf_md                Number of unique degrees of freedom for the
!!                                       reconstructed field
!> @param[in]     map_md                 Dofmap for the cell at the base of the column
!> @param[in]     ndf_lu                 Number of degrees of freedom per cell for the lookup table
!> @param[in]     undf_lu                Number of unique degrees of freedom for the lookup table
!> @param[in]     map_lu                 Dofmap for the cell at the base of the column for the lookup table
!> @param[in]     ndf_ns                 Number of degrees of freedom per cell for the number of index sets
!> @param[in]     undf_ns                Number of unique degrees of freedom for the number of index sets
!> @param[in]     map_ns                 Dofmap for the cell at the base of the column for the number of index sets
!> @param[in]     ndf_w2                 Number of degrees of freedom per cell
!> @param[in]     undf_w2                Number of unique degrees of freedom for the wind field
!> @param[in]     map_w2                 Dofmap for the cell at the base of the column
subroutine adj_poly_adv_upd_lookup_code( nlayers,                &
                                         advective,              &
                                         reconstruction,         &
                                         lookup_poly_adv_upd,    &
                                         set_count_poly_adv_upd, &
                                         wind,                   &
                                         smap_w2_size,           &
                                         smap_w2_max,            &
                                         smap_w2,                &
                                         wind_dir,               &
                                         nsets_max,              &
                                         nindices,               &
                                         ndf_wt,                 &
                                         undf_wt,                &
                                         map_wt,                 &
                                         ndf_md,                 &
                                         undf_md,                &
                                         map_md,                 &
                                         ndf_lu,                 &
                                         undf_lu,                &
                                         map_lu,                 &
                                         ndf_sc,                 &
                                         undf_sc,                &
                                         map_sc,                 &
                                         ndf_w2,                 &
                                         undf_w2,                &
                                         map_w2 )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                    :: nlayers
  integer(kind=i_def), intent(in)                    :: ndf_wt
  integer(kind=i_def), intent(in)                    :: undf_wt
  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt
  integer(kind=i_def), intent(in)                    :: ndf_md
  integer(kind=i_def), intent(in)                    :: undf_md
  integer(kind=i_def), dimension(ndf_md), intent(in) :: map_md
  integer(kind=i_def), intent(in)                    :: ndf_lu
  integer(kind=i_def), intent(in)                    :: undf_lu
  integer(kind=i_def), dimension(ndf_lu), intent(in) :: map_lu
  integer(kind=i_def), intent(in)                    :: ndf_sc
  integer(kind=i_def), intent(in)                    :: undf_sc
  integer(kind=i_def), dimension(ndf_sc), intent(in) :: map_sc
  integer(kind=i_def), intent(in)                    :: ndf_w2
  integer(kind=i_def), intent(in)                    :: undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2

  integer(kind=i_def), intent(in)                                    :: smap_w2_max
  integer(kind=i_def), dimension(4), intent(in)                      :: smap_w2_size
  integer(kind=i_def), dimension(ndf_w2, smap_w2_max, 4), intent(in) :: smap_w2

  real(kind=r_tran), dimension(undf_wt), intent(inout) :: advective
  real(kind=r_tran), dimension(undf_md), intent(inout) :: reconstruction
  integer(kind=i_def), dimension(undf_lu), intent(in)  :: lookup_poly_adv_upd
  integer(kind=i_def), dimension(undf_sc), intent(in)  :: set_count_poly_adv_upd
  real(kind=r_tran), dimension(undf_w2), intent(in)    :: wind
  real(kind=r_tran), dimension(undf_w2), intent(in)    :: wind_dir
  integer(kind=i_def), intent(in)                      :: nsets_max
  integer(kind=i_def), intent(in)                      :: nindices

  ! Internal variables
  integer(kind=i_def), parameter :: nfaces = 4

  integer(kind=i_def) :: k
  integer(kind=i_def) :: df
  integer(kind=i_def) :: ijp
  integer(kind=i_def) :: df1
  integer(kind=i_def) :: df2

  integer(kind=i_def), dimension(4)      :: direction_dofs
  real(kind=r_tran)                      :: direction
  real(kind=r_tran), dimension(nfaces)   :: v_dot_n
  integer(kind=i_def), dimension(nfaces) :: opposite
  logical(kind=l_def), dimension(nfaces) :: missing_neighbour

  real(kind=r_tran), dimension(4,0:nlayers) :: tracer
  real(kind=r_tran), dimension(2,0:nlayers) :: uv_dir
  real(kind=r_tran)                         :: dtdx
  real(kind=r_tran)                         :: dtdy
  real(kind=r_tran)                         :: tracer_loc

  integer(kind=i_def) :: face, set, idx_at_set, nsets, i_set_counter
  integer(kind=i_def) :: df_key, map_wt_rmt, face_idx_rmt
  integer(kind=i_def), dimension(nfaces) :: map_w2_rmt
  logical(kind=l_def) :: missing_neighbour_rmt


  dtdx = 0.0_r_tran
  dtdy = 0.0_r_tran
  tracer = 0.0_r_tran
  v_dot_n(:) = 1.0_r_tran
  v_dot_n(1) = -1.0_r_tran
  v_dot_n(nfaces) = -1.0_r_tran

  direction_dofs(:) = 1
  direction_dofs(2) = 2
  direction_dofs(4) = 2

  opposite(:) = -1
  missing_neighbour(:) = .false.

  do df = 1, nfaces, 1
    df1 = map_w2(df)
    if (smap_w2_size(df) > 1) then
      do df2 = 1, nfaces, 1
        if (smap_w2(df2,2,df) == df1) then
          opposite(df) = df2
        end if
      end do
    else
      opposite(df) = df
      missing_neighbour(df) = .true.
    end if
  end do


  ! <STENCIL>
  do face = 1, nfaces
    i_set_counter = map_sc(1) + (face-1)
    nsets = set_count_poly_adv_upd(i_set_counter)

    df = map_md(1) + (face-1)*nlayers

    df_key = map_lu(1) + (face-1)*nsets_max*nindices

    do set = 1, nsets
      ! Unpack the indices.
      idx_at_set = df_key + (set - 1)*nindices

      map_w2_rmt(:) = lookup_poly_adv_upd(idx_at_set : idx_at_set + nfaces-1)
      face_idx_rmt = lookup_poly_adv_upd(idx_at_set + nfaces)
      missing_neighbour_rmt = (lookup_poly_adv_upd(idx_at_set + nfaces + 1) == 1.0_i_def)
      map_wt_rmt = lookup_poly_adv_upd(idx_at_set + nfaces + 2)

      ! Compute `direction` for each stencil element.
      k = 0
      uv_dir(1,k) = 0.25_r_tran * wind_dir(map_w2_rmt(1)) + 0.25_r_tran * wind_dir(map_w2_rmt(3))
      uv_dir(2,k) = 0.25_r_tran * wind_dir(map_w2_rmt(2)) + 0.25_r_tran * wind_dir(map_w2_rmt(4))

      do k = 1, nlayers - 1, 1
        uv_dir(1,k) = 0.25_r_tran * wind_dir(k + map_w2_rmt(1)) + 0.25_r_tran * wind_dir(k + map_w2_rmt(3)) + &
                      0.25_r_tran * wind_dir(k + map_w2_rmt(1) - 1) + 0.25_r_tran * wind_dir(k + map_w2_rmt(3) - 1)
        uv_dir(2,k) = 0.25_r_tran * wind_dir(k + map_w2_rmt(2)) + 0.25_r_tran * wind_dir(k + map_w2_rmt(4)) + &
                      0.25_r_tran * wind_dir(k + map_w2_rmt(2) - 1) + 0.25_r_tran * wind_dir(k + map_w2_rmt(4) - 1)
      end do

      k = nlayers
      uv_dir(1,k) = 0.25_r_tran * wind_dir(k + map_w2_rmt(1) - 1) + 0.25_r_tran * wind_dir(k + map_w2_rmt(3) - 1)
      uv_dir(2,k) = 0.25_r_tran * wind_dir(k + map_w2_rmt(2) - 1) + 0.25_r_tran * wind_dir(k + map_w2_rmt(4) - 1)

      do k = nlayers, 0, -1
        direction = uv_dir(direction_dofs(face_idx_rmt),k) * v_dot_n(face_idx_rmt)
        tracer_loc = advective(map_wt_rmt + k) * direction

        if (direction <= 0.0_r_tran .and. .not. missing_neighbour_rmt) then
          ijp = map_md(1) + (face-1)*(nlayers+1)
          reconstruction(ijp + k) = reconstruction(ijp + k) + tracer_loc
        end if
      end do

    end do
  end do

  ! <LOCAL>
  k = 0
  uv_dir(1,k) = 0.25_r_tran * wind_dir(map_w2(1)) + 0.25_r_tran * wind_dir(map_w2(3))
  uv_dir(2,k) = 0.25_r_tran * wind_dir(map_w2(2)) + 0.25_r_tran * wind_dir(map_w2(4))

  do k = 1, nlayers - 1, 1
    uv_dir(1,k) = 0.25_r_tran * wind_dir(k + map_w2(1)) + 0.25_r_tran * wind_dir(k + map_w2(3)) + &
                  0.25_r_tran * wind_dir(k + map_w2(1) - 1) + 0.25_r_tran * wind_dir(k + map_w2(3) - 1)
    uv_dir(2,k) = 0.25_r_tran * wind_dir(k + map_w2(2)) + 0.25_r_tran * wind_dir(k + map_w2(4)) + &
                  0.25_r_tran * wind_dir(k + map_w2(2) - 1) + 0.25_r_tran * wind_dir(k + map_w2(4) - 1)
  end do

  k = nlayers
  uv_dir(1,k) = 0.25_r_tran * wind_dir(k + map_w2(1) - 1) + 0.25_r_tran * wind_dir(k + map_w2(3) - 1)
  uv_dir(2,k) = 0.25_r_tran * wind_dir(k + map_w2(2) - 1) + 0.25_r_tran * wind_dir(k + map_w2(4) - 1)

  do df = nfaces, 1, -1
    do k = nlayers, 0, -1
      direction = uv_dir(direction_dofs(df),k) * v_dot_n(df)
      tracer_loc = advective(map_wt(1) + k) * direction

      if (direction > 0.0_r_tran .OR. missing_neighbour(df)) then
        ijp = map_md(1) + (df-1)*(nlayers+1)
        reconstruction(ijp + k) = reconstruction(ijp + k) + tracer_loc
      end if
    end do
  end do

  ! NOTE: Zeroing of `advective` required for a complete adjoint.
  !       This needs to be done separately to avoid a race condition.

end subroutine adj_poly_adv_upd_lookup_code

end module adj_poly_adv_upd_lookup_kernel_mod
