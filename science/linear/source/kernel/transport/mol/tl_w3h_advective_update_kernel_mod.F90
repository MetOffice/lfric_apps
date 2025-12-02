!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes tangent linear horizontal advective update through fitting a high order
!!        upwind reconstruction. The upwind direction is dictated by the linear wind
!> @details Compute the advective update (u.grad) for a tracer field using a high order
!!          polynomial fit to the integrated tracer values. The stencil used for the
!!          polynomial is centred on the upwind cell.
!!          This method is only valid for lowest order elements.

module tl_w3h_advective_update_kernel_mod

use argument_mod,      only : arg_type,                  &
                              GH_FIELD, GH_REAL,         &
                              GH_OPERATOR,               &
                              GH_READWRITE, GH_READ,     &
                              STENCIL, CROSS2D,          &
                              ANY_DISCONTINUOUS_SPACE_1, &
                              ANY_W2, &
                              CELL_COLUMN
use constants_mod,     only : r_def, i_def, l_def
use fs_continuity_mod, only : W3
use kernel_mod,        only : kernel_type

implicit none

private
!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: tl_w3h_advective_update_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                                                             &
       arg_type(GH_FIELD,    GH_REAL, GH_READWRITE, W3),                                          &
       arg_type(GH_FIELD,    GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1, STENCIL(CROSS2D)), &
       arg_type(GH_FIELD,    GH_REAL, GH_READ,      ANY_W2, STENCIL(CROSS2D)),                    &
       arg_type(GH_FIELD,    GH_REAL, GH_READ,      ANY_W2),                                      &
       arg_type(GH_OPERATOR, GH_REAL, GH_READ,      W3, W3)                                       &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: tl_w3h_advective_update_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: tl_w3h_advective_update_code
contains

!> @brief Computes the horizontal advective update for a tracer in W3.
!> @param[in]     cell                Horizontal cell index
!> @param[in]     nlayers             Number of layers
!> @param[in,out] advective_increment Advective update field to compute
!> @param[in]     tracer              Pointwise tracer field to advect stored on cell faces
!> @param[in]     smap_md_size        Sizes of the stencil map in each direction
!> @param[in]     smap_md_max         Maximum size of the stencil map
!> @param[in]     smap_md             Stencil map for the tracer space
!> @param[in]     ls_wind             Linear wind field
!> @param[in]     smap_w2_size        Sizes of the stencil map in each direction
!> @param[in]     smap_w2_max         Maximum size of the stencil map
!> @param[in]     smap_w2             Stencil map for the wind space
!> @param[in]     wind                Perturbation wind field
!> @param[in]     ncell_3d            Total number of cells
!> @param[in]     m3_inv              Inverse mass matrix for W3 space
!> @param[in]     ndf_w3              Number of degrees of freedom per cell
!> @param[in]     undf_w3             Number of unique degrees of freedom for the
!!                                    advective_update field
!> @param[in]     map_w3              Dofmap for the cell at the base of the column
!> @param[in]     ndf_wd              Number of degrees of freedom per cell
!> @param[in]     undf_wd             Number of unique degrees of freedom for the
!!                                    tracer field
!> @param[in]     map_wd              Dofmap for the cell at the base of the column
!> @param[in]     ndf_w2              Number of degrees of freedom per cell for the wind fields
!> @param[in]     undf_w2             Number of unique degrees of freedom for the wind fields
!> @param[in]     map_w2              Dofmap for the cell at the base of the column for the wind fields
subroutine tl_w3h_advective_update_code( cell,                 &
                                         nlayers,              &
                                         advective_increment,  &
                                         tracer,               &
                                         smap_md_size,         &
                                         smap_md_max,          &
                                         smap_md,              &
                                         ls_wind,              &
                                         smap_w2_size,         &
                                         smap_w2_max,          &
                                         smap_w2,              &
                                         wind,                 &
                                         ncell_3d,             &
                                         m3_inv,               &
                                         ndf_w3,               &
                                         undf_w3,              &
                                         map_w3,               &
                                         ndf_md,               &
                                         undf_md,              &
                                         map_md,               &
                                         ndf_w2,               &
                                         undf_w2,              &
                                         map_w2                &
                                         )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                     :: nlayers, cell, ncell_3d
  integer(kind=i_def), intent(in)                     :: ndf_w3, ndf_w2, ndf_md
  integer(kind=i_def), intent(in)                     :: undf_w3, undf_w2, undf_md
  integer(kind=i_def), dimension(ndf_w3),  intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_w2),  intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_md),  intent(in) :: map_md

  integer(kind=i_def),                                  intent(in) :: smap_md_max
  integer(kind=i_def), dimension(4),                    intent(in) :: smap_md_size
  integer(kind=i_def), dimension(ndf_md,smap_md_max,4), intent(in) :: smap_md
  integer(kind=i_def),                                  intent(in) :: smap_w2_max
  integer(kind=i_def), dimension(4),                    intent(in) :: smap_w2_size
  integer(kind=i_def), dimension(ndf_w2,smap_w2_max,4), intent(in) :: smap_w2

  real(kind=r_def), dimension(undf_w3), intent(inout) :: advective_increment
  real(kind=r_def), dimension(undf_w2), intent(in)    :: ls_wind
  real(kind=r_def), dimension(undf_w2), intent(in)    :: wind
  real(kind=r_def), dimension(undf_md), intent(in)    :: tracer

  real(kind=r_def), dimension(ncell_3d, ndf_w3, ndf_w3), intent(in) :: m3_inv

  ! Internal variables
  integer(kind=i_def) :: k, ik, face, df, df1, df2
  real(kind=r_def)    :: ls_u, ls_v, u, v, dtdx, dtdy
  real(kind=r_def)    :: t_E, t_W, t_N, t_S

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

  ! Horizontal advective update
  ! Reconstruction is stored on a layer first multidata field
  ! with index for face f = (1,2,3,4) = (W,S,E,N): map_md(1) + (f-1)*nlayers + k
  ! each cell contains the values for when it is the upwind cell for each edge
  ! so if u.n > 0 then we set the field to be the value on this edge from this cell
  ! and if u.n < 0 then we set the field to be the value on this edge from a
  ! neighbouring cell
  do k = 0, nlayers - 1

    ! u * dt/dx
    ls_u =  0.5_r_def*( ls_wind(map_w2(1) + k) + ls_wind(map_w2(3) + k) )
    u =  0.5_r_def*( wind(map_w2(1) + k) + wind(map_w2(3) + k) )

    face = 1
    if ( ls_u > 0.0_r_def .and. .not. missing_neighbour(face) ) then
      ! t_W from neighbouring column (if it exists)
      t_W = tracer( smap_md(1,2,face) + (opposite(face)-1)*nlayers + k )
    else
      ! t_W from this column
      t_W = tracer( map_md(1) + k )
    end if

    face = 3
    if ( ls_u <= 0.0_r_def .and. .not. missing_neighbour(face) ) then
      ! t_E from neighbouring column (if it exists)
      t_E = tracer( smap_md(1,2,face) + (opposite(face)-1)*nlayers + k )
    else
      ! t_E from this column
      t_E = tracer( map_md(1) + 2*nlayers + k )
    end if

    dtdx = t_E - t_W

    ! v*dt/dy
    ls_v = -0.5_r_def*( ls_wind(map_w2(2) + k) + ls_wind(map_w2(4) + k) )
    v = -0.5_r_def*( wind(map_w2(2) + k) + wind(map_w2(4) + k) )

    face = 2
    if ( ls_v > 0.0_r_def .and. .not. missing_neighbour(face) ) then
      ! t_S from neighbouring column (if it exists)
      t_S = tracer( smap_md(1,2,face) + (opposite(face)-1)*nlayers + k )
    else
      ! t_S from this column
      t_S = tracer( map_md(1) + nlayers + k )
    end if

    face = 4
    if ( ls_v <= 0.0_r_def .and. .not. missing_neighbour(face) ) then
      ! t_N from neighbouring column (if it exists)
      face = 4
      t_N = tracer( smap_md(1,2,face) + (opposite(face)-1)*nlayers + k )
    else
      ! t_N from this column
      t_N = tracer( map_md(1) + 3*nlayers + k )
    end if

    dtdy = t_N - t_S

    ik = 1 + k + (cell-1)*nlayers
    advective_increment(map_w3(1)+k) = advective_increment(map_w3(1)+k) &
                                     + m3_inv(ik,1,1)*(u*dtdx + v*dtdy)
  end do

end subroutine tl_w3h_advective_update_code

end module tl_w3h_advective_update_kernel_mod
