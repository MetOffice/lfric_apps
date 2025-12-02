!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes adjoint of horizontal advective update
!!        through fitting a high order upwind reconstruction.
module adj_w3h_adv_upd_lookup_kernel_mod

use argument_mod,          only : arg_type,                  &
                                  GH_FIELD, GH_SCALAR,       &
                                  GH_REAL, GH_INTEGER,       &
                                  GH_OPERATOR,               &
                                  GH_WRITE, GH_READ,         &
                                  STENCIL, CROSS2D,          &
                                  ANY_DISCONTINUOUS_SPACE_1, &
                                  ANY_DISCONTINUOUS_SPACE_2, &
                                  ANY_DISCONTINUOUS_SPACE_3, &
                                  ANY_W2, &
                                  CELL_COLUMN
use constants_mod,         only : r_def, i_def, l_def, r_tran
use fs_continuity_mod,     only : W3
use reference_element_mod, only : W, S, E, N
use kernel_mod,            only : kernel_type

implicit none

private
!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
!! NOTE: It is not currently possible to use this as a kernel due to the mix of an
!!       operator and an integer field:
!!
!!          "Parse Error: In the LFRic API a kernel that has an LMA operator argument must
!!           only have field arguments with 'gh_real' data type but kernel
!!           'adj_w3h_adv_upd_lookup_kernel_type' has a field argument with 'gh_integer' data type."
!!
!!       This requires PSyKAl-lite code to invoke.
type, public, extends(kernel_type) :: adj_w3h_adv_upd_lookup_kernel_type
  type(arg_type) :: meta_args(8) = (/                                              &
       arg_type(GH_FIELD,    GH_REAL,    GH_READ,      W3),                        & ! advective_increment
       arg_type(GH_FIELD,    GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! tracer
       arg_type(GH_FIELD,    GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! lookup_w3h_adv_up
       arg_type(GH_FIELD,    GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_3), & ! set_count_w3h_adv_up
       arg_type(GH_FIELD,    GH_REAL,    GH_READ,      ANY_W2, STENCIL(CROSS2D)),  & ! wind
       arg_type(GH_OPERATOR, GH_REAL,    GH_READ,      W3, W3),                    & ! m3_inv
       arg_type(GH_SCALAR,   GH_INTEGER, GH_READ),                                 & ! nsets_max
       arg_type(GH_SCALAR,   GH_INTEGER, GH_READ)                                  & ! nindices
       /)
  integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: adj_w3h_adv_upd_lookup_code
end type adj_w3h_adv_upd_lookup_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: adj_w3h_adv_upd_lookup_code

contains
!> @brief Computes the adjoint horizontal advective update for a tracer in W3.
!> @param[in]     cell                  Horizontal cell index
!> @param[in]     nlayers               Number of layers
!> @param[in]     advective_increment   Advective update field to compute
!> @param[in,out] tracer                Pointwise tracer field to advect stored on cell faces
!> @param[in]     lookup_w3h_adv_upd    Lookup table relating the 0th layer tracer indices
!!                                      to all sets of advective increment, wind and m3_inv indices
!> @param[in]     set_count_w3h_adv_up  Number of sets of lookup indices per cell
!> @param[in]     wind                  Wind field
!> @param[in]     smap_w2_size          Sizes of the stencil map in each direction
!> @param[in]     smap_w2_max           Maximum size of the stencil map
!> @param[in]     smap_w2               Stencil map for the wind space
!> @param[in]     ncell_3d              Total number of cells
!> @param[in]     m3_inv                Inverse mass matrix for W3 space
!> @param[in]     nsets_max             Maximum number of sets of indices for the lookup table
!> @param[in]     nindices              Number of indices per set for the lookup table
!> @param[in]     ndf_w3                Number of degrees of freedom per cell
!> @param[in]     undf_w3               Number of unique degrees of freedom for the
!!                                      advective_update field
!> @param[in]     map_w3                Dofmap for the cell at the base of the column
!> @param[in]     ndf_md                Number of degrees of freedom per cell
!> @param[in]     undf_md               Number of unique degrees of freedom for the
!!                                      tracer field
!> @param[in]     map_md                Dofmap for the cell at the base of the column
!> @param[in]     ndf_lu                Number of degrees of freedom per cell
!> @param[in]     undf_lu               Number of unique degrees of freedom for the
!!                                      lookup table
!> @param[in]     map_lu                Dofmap for the cell at the base of the column
!> @param[in]     ndf_ns                Number of degrees of freedom per cell for the number of index sets
!> @param[in]     undf_ns               Number of unique degrees of freedom for the number of index sets
!> @param[in]     map_ns                Dofmap for the cell at the base of the column for the number of index sets
!> @param[in]     ndf_w2                Number of degrees of freedom per cell for the wind fields
!> @param[in]     undf_w2               Number of unique degrees of freedom for the wind fields
!> @param[in]     map_w2                Dofmap for the cell at the base of the column for the wind fields
subroutine adj_w3h_adv_upd_lookup_code( cell,                &
                                        nlayers,             &
                                        advective_increment, &
                                        tracer,              &
                                        lookup_w3h_adv_upd,  &
                                        set_count_w3h_adv_upd,  &
                                        wind,                &
                                        smap_w2_size,        &
                                        smap_w2_max,         &
                                        smap_w2,             &
                                        ncell_3d,            &
                                        m3_inv,              &
                                        nsets_max,               &
                                        nindices,            &
                                        ndf_w3,              &
                                        undf_w3,             &
                                        map_w3,              &
                                        ndf_md,              &
                                        undf_md,             &
                                        map_md,              &
                                        ndf_lu,              &
                                        undf_lu,             &
                                        map_lu,              &
                                        ndf_sc,              &
                                        undf_sc,             &
                                        map_sc,              &
                                        ndf_w2,              &
                                        undf_w2,             &
                                        map_w2 )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                                  :: nlayers
  integer(kind=i_def), intent(in)                                  :: cell
  integer(kind=i_def), intent(in)                                  :: ncell_3d
  integer(kind=i_def), intent(in)                                  :: nsets_max
  integer(kind=i_def), intent(in)                                  :: nindices

  integer(kind=i_def), intent(in)                                  :: ndf_w3
  integer(kind=i_def), intent(in)                                  :: undf_w3
  integer(kind=i_def), dimension(ndf_w3), intent(in)               :: map_w3
  integer(kind=i_def), intent(in)                                  :: ndf_md
  integer(kind=i_def), intent(in)                                  :: undf_md
  integer(kind=i_def), dimension(ndf_md), intent(in)               :: map_md
  integer(kind=i_def), intent(in)                                  :: ndf_lu
  integer(kind=i_def), intent(in)                                  :: undf_lu
  integer(kind=i_def), dimension(ndf_lu), intent(in)               :: map_lu
  integer(kind=i_def), intent(in)                                  :: ndf_sc
  integer(kind=i_def), intent(in)                                  :: undf_sc
  integer(kind=i_def), dimension(ndf_sc), intent(in)               :: map_sc
  integer(kind=i_def), intent(in)                                  :: ndf_w2
  integer(kind=i_def), intent(in)                                  :: undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in)               :: map_w2

  integer(kind=i_def), intent(in)                                  :: smap_w2_max
  integer(kind=i_def), dimension(4), intent(in)                    :: smap_w2_size
  integer(kind=i_def), dimension(ndf_w2,smap_w2_max,4), intent(in) :: smap_w2
  real(kind=r_tran), dimension(undf_w3), intent(in)                :: advective_increment
  real(kind=r_tran), dimension(undf_w2), intent(in)                :: wind
  real(kind=r_tran), dimension(undf_md), intent(inout)             :: tracer
  integer(kind=i_def), dimension(undf_lu), intent(in)              :: lookup_w3h_adv_upd
  integer(kind=i_def), dimension(undf_sc), intent(in)              :: set_count_w3h_adv_upd
  real(kind=r_def), dimension(ncell_3d,ndf_w3,ndf_w3), intent(in)  :: m3_inv

  ! Internal variables
  integer(kind=i_def), parameter         :: nfaces = 4
  integer(kind=i_def)                    :: k
  integer(kind=i_def)                    :: set, nsets, i_set_counter
  integer(kind=i_def)                    :: idx_at_set
  integer(kind=i_def)                    :: ik, ik_remote
  integer(kind=i_def)                    :: face, face_remote
  integer(kind=i_def)                    :: df, df_key, df_w3_remote
  integer(kind=i_def)                    :: df1
  integer(kind=i_def)                    :: df2
  real(kind=r_tran)                      :: u, u_remote
  real(kind=r_tran)                      :: v, v_remote
  real(kind=r_tran)                      :: dtdx
  real(kind=r_tran)                      :: dtdy
  real(kind=r_tran)                      :: t_E, t_E_remote
  real(kind=r_tran)                      :: t_W, t_W_remote
  real(kind=r_tran)                      :: t_N, t_N_remote
  real(kind=r_tran)                      :: t_S, t_S_remote
  integer(kind=i_def), dimension(nfaces) :: opposite
  logical(kind=l_def), dimension(nfaces) :: missing_neighbour
  logical(kind=i_def)                    :: missing_neighbour_remote

  ! Zeroing of internal active variables
  dtdx = 0.0_r_tran
  dtdy = 0.0_r_tran
  t_N = 0.0_r_tran
  t_E = 0.0_r_tran
  t_S = 0.0_r_tran
  t_W = 0.0_r_tran

  ! For each face of cell, find the index in the neighbouring cell that
  ! corresponds to it.
  ! i.e for no orientation changes opposite = ( 3, 4, 1, 2 )
  ! We use the W2 map to determine these
  ! If there is no neighbour then we ensure the opposite points to
  ! the value on this edge
  opposite = -1
  missing_neighbour = .false.
  do df = 1, nfaces, 1
    df1 = map_w2(df)
    if (smap_w2_size(df) > 1) then
      ! There is a neighbour in direction df so find the
      ! neighboring edge corresponding to edge df
      do df2 = 1, nfaces, 1
        if (smap_w2(df2,2,df) == df1)  opposite(df) = df2
      end do
    else
      ! There is no neighbour in direction df so point to itself
      opposite(df) = df
      missing_neighbour(df) = .true.
    end if
  end do

  do k = nlayers - 1, 0, -1
    u =  0.5_r_tran*( wind(map_w2(1) + k) + wind(map_w2(3) + k) )
    v = -0.5_r_tran*( wind(map_w2(2) + k) + wind(map_w2(4) + k) )
    ik = 1 + k + (cell-1)*nlayers
    dtdx = advective_increment(map_w3(1) + k) * u * real(m3_inv(ik,1,1), r_tran)
    dtdy = advective_increment(map_w3(1) + k) * v * real(m3_inv(ik,1,1), r_tran)
    t_N = dtdy
    t_S = -dtdy
    t_E = dtdx
    t_W = -dtdx

    do face = nfaces, 1, -1
      df = map_md(1) + (face-1)*nlayers
      df_key = map_lu(1) + (face-1)*nsets_max*nindices

      ! Number of sets / acesses done on this cell.
      i_set_counter = map_sc(1) + (face-1)
      nsets = set_count_w3h_adv_upd(i_set_counter)

      do set = 1, nsets
        idx_at_set = df_key + (set - 1)*nindices
        face_remote = lookup_w3h_adv_upd( idx_at_set )
        df1 = lookup_w3h_adv_upd( idx_at_set + 1 )
        df2 = lookup_w3h_adv_upd( idx_at_set + 2 )
        missing_neighbour_remote = ( lookup_w3h_adv_upd( idx_at_set + 3 ) == 1_i_def )
        df_w3_remote = lookup_w3h_adv_upd( idx_at_set + 4 )
        ik_remote = lookup_w3h_adv_upd( idx_at_set + 5 ) + k

        ! Check if neighbours have written to us based on the ls state.
        select case ( face_remote )
          case ( W )
            u_remote = 0.5_r_tran*( wind(df1 + k) + wind(df2 + k) )
            if ( u_remote > 0.0_r_tran .and. .not. missing_neighbour_remote ) then
              t_W_remote = -1.0_r_tran * &
                            advective_increment(df_w3_remote + k) * &
                            u_remote * real(m3_inv(ik_remote, 1, 1), r_tran)
              tracer(df + k) = tracer(df + k) + t_W_remote
            end if
          case ( S )
            v_remote = -0.5_r_tran*( wind(df1 + k) + wind(df2 + k) )
            if ( v_remote > 0.0_r_tran .and. .not. missing_neighbour_remote ) then
              t_S_remote = -1.0_r_tran * &
                            advective_increment(df_w3_remote + k) * &
                            v_remote * real(m3_inv(ik_remote, 1, 1), r_tran)
              tracer(df + k) = tracer(df + k) + t_S_remote
            end if
          case ( E )
            u_remote = 0.5_r_tran*( wind(df1 + k) + wind(df2 + k) )
            if ( u_remote <= 0.0_r_tran .and. .not. missing_neighbour_remote ) then
              t_E_remote = advective_increment(df_w3_remote + k) * &
                           u_remote * real(m3_inv(ik_remote, 1, 1), r_tran)
              tracer(df + k) = tracer(df + k) + t_E_remote
            end if
          case ( N )
            v_remote = -0.5_r_tran*( wind(df1 + k) + wind(df2 + k) )
            if ( v_remote <= 0.0_r_tran .and. .not. missing_neighbour_remote ) then
              t_N_remote = advective_increment(df_w3_remote + k) * &
                           v_remote * real(m3_inv(ik_remote, 1, 1), r_tran)
              tracer(df + k) = tracer(df + k) + t_N_remote
            end if
        end select
      end do

      ! Check if we wrote to ourself
      if ( face == 1_i_def .and. .not. ( u > 0.0_r_tran .and. .not. missing_neighbour(face) ) ) then
        ! West
        tracer(df + k) = tracer(df + k) + t_W
      end if

      if ( face == 2_i_def .and. .not. ( v > 0.0_r_tran .and. .not. missing_neighbour(face) ) ) then
        ! South
        tracer(df + k) = tracer(df + k) + t_S
      end if

      if ( face == 3_i_def .and. .not. ( u <= 0.0_r_tran .and. .not. missing_neighbour(face) ) ) then
        ! East
        tracer(df + k) = tracer(df + k) + t_E
      end if

      if (face == 4_i_def .and. .not. ( v <= 0.0_r_tran .and. .not. missing_neighbour(face) ) ) then
        ! North
        tracer(df + k) = tracer(df + k) + t_N
      end if

    end do ! face

  end do ! k

end subroutine adj_w3h_adv_upd_lookup_code

end module adj_w3h_adv_upd_lookup_kernel_mod
