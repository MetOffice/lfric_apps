!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Apply horizontal Smagorinsky diffusion visc_m * (d2dx2 + d2dy2) to
!>        the components of the momentum equation for lowest order elements.
!>        As in the UM, the Smagorinsky scheme has not been correctly
!>        implemented for the momentum equations - the stress terms have been
!>        implemented purely as diffusion terms. This results in a missing
!>        term as detailed in UMDP28.
!>
module momentum_smagorinsky_kernel_mod

  use argument_mod,      only : arg_type,                    &
                                GH_FIELD, GH_REAL,           &
                                GH_READ, GH_WRITE,           &
                                STENCIL, CROSS, CELL_COLUMN, &
                                ANY_DISCONTINUOUS_SPACE_9
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W2, W1, Wtheta
  use kernel_mod,        only : kernel_type

  implicit none
  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: momentum_smagorinsky_kernel_type
    private
    type(arg_type) :: meta_args(7) = (/                                      &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE,  W2),                       &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,   W2, STENCIL(CROSS)),       &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,   W2, STENCIL(CROSS)),       &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,   W2),                       &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,   W1),                       &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,   Wtheta, STENCIL(CROSS)),   &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,   ANY_DISCONTINUOUS_SPACE_9, &
                                                            STENCIL(CROSS))  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: momentum_smagorinsky_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: momentum_smagorinsky_code

contains

!> @brief Calculates diffusion increment for wind field using horizontal Smagorinsky diffusion
!! @param[in] nlayers Number of layers in the mesh
!! @param[in,out] u_inc Diffusion increment for wind field
!! @param[in] u_n Input wind field
!! @param[in] map_w2_stencil_size Number of cells in the stencil at the base of the column for w2
!! @param[in] map_w2_stencil Array holding the stencil dofmap for the cell at the base of the column for w2
!! @param[in] dx_at_w2 Grid length at the w2 dofs
!! @param[in] map_dx_stencil_size Number of cells in the stencil at the base of the column for w2
!! @param[in] map_dx_stencil Array holding the stencil dofmap for the cell at the base of the column for w2
!! @param[in] height_w2 Height of w3 space levels above surface at cell faces
!! @param[in] height_w1 Height of wtheta levels above surface at cell faces
!! @param[in] visc_m Diffusion coefficient for momentum on wtheta points
!! @param[in] map_wt_stencil_size Number of cells in the stencil at the base of the column for wt
!! @param[in] map_wt_stencil Array holding the stencil dofmap for the cell at the base of the column for wt
!! @param[in] panel_id  The ID number of the current panel
!! @param[in] map_pid_stencil_size Size of the panel ID stencil
!! @param[in] map_pid_stencil Stencil map for the panel ID
!! @param[in] ndf_w2 Number of degrees of freedom per cell for wind space
!! @param[in] undf_w2 Number of unique degrees of freedom for wind space
!! @param[in] map_w2 Array holding the dofmap for the cell at the base of the column for wind space
!! @param[in] ndf_w1 Number of degrees of freedom per cell for w1
!! @param[in] undf_w1 Number of unique degrees of freedom for w1
!! @param[in] map_w1 Array holding the dofmap for the cell at the base of the column for w1
!! @param[in] ndf_wt Number of degrees of freedom per cell for theta space
!! @param[in] undf_wt Number of unique degrees of freedom for theta space
!! @param[in] map_wt Array holding the dofmap for the cell at the base of the column for theta space
!! @param[in] ndf_pid  Number of degrees of freedom per cell for pid space
!! @param[in] undf_pid  Number of unique degrees of freedom for pid space
!! @param[in] map_pid  Cell dofmap for pid space
subroutine momentum_smagorinsky_code( nlayers,                                 &
                                      u_inc,                                   &
                                      u_n,                                     &
                                      map_w2_stencil_size, map_w2_stencil,     &
                                      dx_at_w2,                                &
                                      map_dx_stencil_size, map_dx_stencil,     &
                                      height_w2,                               &
                                      height_w1,                               &
                                      visc_m,                                  &
                                      map_wt_stencil_size, map_wt_stencil,     &
                                      panel_id,                                &
                                      map_pid_stencil_size, map_pid_stencil,   &
                                      ndf_w2, undf_w2, map_w2,                 &
                                      ndf_w1, undf_w1, map_w1,                 &
                                      ndf_wt, undf_wt, map_wt,                 &
                                      ndf_pid, undf_pid, map_pid               &
                                     )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w2, ndf_w1, ndf_wt, ndf_pid
  integer(kind=i_def), intent(in) :: undf_w2, undf_w1, undf_wt, undf_pid
  integer(kind=i_def), intent(in) :: map_w2_stencil_size, map_dx_stencil_size, map_wt_stencil_size, map_pid_stencil_size
  integer(kind=i_def), dimension(ndf_w2,map_w2_stencil_size), intent(in)  :: map_w2_stencil
  integer(kind=i_def), dimension(ndf_w2,map_dx_stencil_size), intent(in)  :: map_dx_stencil
  integer(kind=i_def), dimension(ndf_wt,map_wt_stencil_size), intent(in)  :: map_wt_stencil
  integer(kind=i_def), dimension(ndf_pid,map_pid_stencil_size), intent(in) :: map_pid_stencil
  integer(kind=i_def), dimension(ndf_w2),  intent(in)  :: map_w2
  integer(kind=i_def), dimension(ndf_w1),  intent(in)  :: map_w1
  integer(kind=i_def), dimension(ndf_wt),  intent(in)  :: map_wt
  integer(kind=i_def), dimension(ndf_pid), intent(in)  :: map_pid

  real(kind=r_def), dimension(undf_w2),  intent(inout) :: u_inc
  real(kind=r_def), dimension(undf_w2),  intent(in)    :: u_n, dx_at_w2
  real(kind=r_def), dimension(undf_w2),  intent(in)    :: height_w2
  real(kind=r_def), dimension(undf_w1),  intent(in)    :: height_w1
  real(kind=r_def), dimension(undf_wt),  intent(in)    :: visc_m
  real(kind=r_def), dimension(undf_pid), intent(in)    :: panel_id

  ! Internal variables
  integer(kind=i_def)                      :: k, kp, df
  real(kind=r_def)                         :: d2dx, d2dy
  real(kind=r_def), dimension(1:nlayers-1) :: idx2, idy2
  real(kind=r_def), dimension(0:nlayers-1,4) :: idx2_w2, idy2_w2
  real(kind=r_def), dimension(0:nlayers,4) :: visc_m_w2
  real(kind=r_def)                         :: weight_pl, weight_min
  real(kind=r_def)                         :: visc_m_w2_w3

  integer(kind=i_def) :: true_stencil_map(ndf_w2,map_w2_stencil_size)
  integer(kind=i_def) :: stencil_cell
  integer(kind=i_def) :: cell_panel, stencil_panel, panel_edge
  integer(kind=i_def) :: vec_dir(ndf_w2,map_w2_stencil_size)

  ! Assumed direction for derivatives in this kernel is:
  !  y
  !  ^
  !  |_> x
  !

  ! Layout of dofs for the stencil map
  ! dimensions of map are (ndf, ncell)
  ! Horizontally:
  !
  !   -- 4 --
  !   |     |
  !   1     3
  !   |     |
  !   -- 2 --
  !
  ! df = 5 is in the centre on the bottom face
  ! df = 6 is in the centre on the top face

  ! The layout of the cells in the stencil is:
  !
  !          -----
  !          |   |
  !          | 5 |
  !     ---------------
  !     |    |   |    |
  !     |  2 | 1 |  4 |
  !     ---------------
  !          |   |
  !          | 3 |
  !          -----

  ! If the full stencil isn't available, we must be at the domain edge.
  ! The increment is already 0, so we just exit the routine.
  if (map_w2_stencil_size < 5_i_def) then
    return
  end if

  ! The W2H DoF values change in orientation when we cross over a panel
  ! Vector directions parallel to the boundary (i.e. the winds on faces
  ! perpendicular to the boundary) also flip sign
  ! We need to take this into account by adjusting the stencil map used for
  ! the wind field. Do this by looking at whether the panel changes
  ! for other cells in the stencil
  cell_panel = int(panel_id(map_pid(1)), i_def)

  do stencil_cell = 1, map_w2_stencil_size
    stencil_panel = int(panel_id(map_pid_stencil(1, stencil_cell)), i_def)
    ! Create panel_edge to check whether a panel is changing
    panel_edge = 10*cell_panel + stencil_panel

    select case (panel_edge)
    case (41, 32, 16, 25, 64, 53)
      ! Clockwise rotation of panel
      true_stencil_map(1, stencil_cell) = map_w2_stencil(2, stencil_cell)
      true_stencil_map(2, stencil_cell) = map_w2_stencil(3, stencil_cell)
      true_stencil_map(3, stencil_cell) = map_w2_stencil(4, stencil_cell)
      true_stencil_map(4, stencil_cell) = map_w2_stencil(1, stencil_cell)
      ! Vertical dofs unchanged
      true_stencil_map(5, stencil_cell) = map_w2_stencil(5, stencil_cell)
      true_stencil_map(6, stencil_cell) = map_w2_stencil(6, stencil_cell)
      ! Flip direction of vectors if necessary
      vec_dir(1,stencil_cell) = -1_i_def
      vec_dir(2,stencil_cell) = 1_i_def
      vec_dir(3,stencil_cell) = -1_i_def
      vec_dir(4,stencil_cell) = 1_i_def
    case (14, 23, 61, 52, 46, 35)
      ! Anti-clockwise rotation of panel
      true_stencil_map(1, stencil_cell) = map_w2_stencil(4, stencil_cell)
      true_stencil_map(2, stencil_cell) = map_w2_stencil(1, stencil_cell)
      true_stencil_map(3, stencil_cell) = map_w2_stencil(2, stencil_cell)
      true_stencil_map(4, stencil_cell) = map_w2_stencil(3, stencil_cell)
      ! Vertical dofs unchanged
      true_stencil_map(5, stencil_cell) = map_w2_stencil(5, stencil_cell)
      true_stencil_map(6, stencil_cell) = map_w2_stencil(6, stencil_cell)
      ! Flip direction of vectors if necessary
      vec_dir(1,stencil_cell) = 1_i_def
      vec_dir(2,stencil_cell) = -1_i_def
      vec_dir(3,stencil_cell) = 1_i_def
      vec_dir(4,stencil_cell) = -1_i_def
    case default
      ! Same panel or crossing panel with no rotation, so stencil map is unchanged
      true_stencil_map(:, stencil_cell) = map_w2_stencil(:, stencil_cell)
      vec_dir(:,stencil_cell) = 1_i_def
    end select
  end do

  ! Loop over horizontal faces of the cell for horizontal diffusion
  do df = 1,4

    ! Only calculate face if it's not already been done
    if (u_inc(map_w2(df)+1) == 0.0_r_def) then

      ! Compute horizontal grid spacing
      if (df == 1 .or. df == 3) then
        ! For Dofs 1 & 3, dx is given by the input field,
        ! whilst dy needs to be computed from the 4 neighbouring dy points
        do k = 0, nlayers - 1
          idx2_w2(k,df) = 1.0_r_def/(dx_at_w2(true_stencil_map(df,1)+k))**2
          idy2_w2(k,df) = (4.0_r_def/(dx_at_w2(true_stencil_map(2,1)+k)+       &
                                      dx_at_w2(true_stencil_map(4,1)+k)+       &
                                      dx_at_w2(true_stencil_map(2,df+1)+k)+    &
                                      dx_at_w2(true_stencil_map(4,df+1)+k)))**2
        end do
      else
        ! For Dofs 2 & 4, dy is given by the input field,
        ! whilst dx needs to be computed from the 4 neighbouring dx points
        do k = 0, nlayers - 1
          idx2_w2(k,df) = (4.0_r_def/(dx_at_w2(true_stencil_map(1,1)+k)+       &
                                      dx_at_w2(true_stencil_map(3,1)+k)+       &
                                      dx_at_w2(true_stencil_map(1,df+1)+k)+    &
                                      dx_at_w2(true_stencil_map(3,df+1)+k)))**2
          idy2_w2(k,df) = 1.0_r_def/(dx_at_w2(true_stencil_map(df,1)+k))**2
        end do
      end if

      ! Horizontal interpolation of visc_m to cell face
      do k = 0, nlayers
        visc_m_w2(k,df) = ( visc_m(map_wt_stencil(1,1) + k) + visc_m(map_wt_stencil(1,1+df) + k) ) / 2.0_r_def
      end do

      ! Horizontal velocity diffusion
      do k = 0, nlayers - 1
        kp = k + 1

        ! Vertical interpolation weights:
        weight_pl = (height_w2(map_w2(df) + k) - height_w1(map_w1(df) + k)) /  &
                    (height_w1(map_w1(df) + kp) - height_w1(map_w1(df) + k))
        weight_min = (height_w1(map_w1(df) + kp) - height_w2(map_w2(df) + k)) /&
                     (height_w1(map_w1(df) + kp) - height_w1(map_w1(df) + k))

        ! Vertical interpolation of visc_m from wt to w3 levels
        visc_m_w2_w3 = ( weight_min * visc_m_w2(k,df) ) + ( weight_pl * visc_m_w2(kp,df) )

        ! horizontal diffusion:
        d2dx = (vec_dir(df,2)*u_n(true_stencil_map(df,2) + k)                  &
             - 2.0_r_def*u_n(true_stencil_map(df,1) + k)                       &
             + vec_dir(df,4)*u_n(true_stencil_map(df,4) + k) ) * idx2_w2(k,df)
        d2dy = (vec_dir(df,3)*u_n(true_stencil_map(df,3) + k)                  &
             - 2.0_r_def*u_n(true_stencil_map(df,1) + k)                       &
             + vec_dir(df,5)*u_n(true_stencil_map(df,5) + k) ) * idy2_w2(k,df)
        u_inc(map_w2(df) + k) = visc_m_w2_w3 * (d2dx + d2dy)

      end do

    end if

  end do

  ! Vertical velocity diffusion for this cell
  do k = 1, nlayers - 1

    ! dx and dy are average of cell face values
    ! N.B. not accounting for difference between w3 and wtheta heights
    idx2(k) = (2.0_r_def/(dx_at_w2(true_stencil_map(1,1)+k)+dx_at_w2(true_stencil_map(3,1)+k)))**2
    idy2(k) = (2.0_r_def/(dx_at_w2(true_stencil_map(2,1)+k)+dx_at_w2(true_stencil_map(4,1)+k)))**2

    ! w diffusion:
    d2dx = (u_n(true_stencil_map(5,2) + k)                                     &
         - 2.0_r_def*u_n(true_stencil_map(5,1) + k)                            &
         + u_n(true_stencil_map(5,4) + k) ) * idx2(k)
    d2dy = (u_n(true_stencil_map(5,3) + k)                                     &
         - 2.0_r_def*u_n(true_stencil_map(5,1) + k)                            &
         + u_n(true_stencil_map(5,5) + k) ) * idy2(k)
    u_inc(map_w2(5) + k) = visc_m(map_wt_stencil(1,1) + k) * (d2dx + d2dy)

  end do

end subroutine momentum_smagorinsky_code

end module momentum_smagorinsky_kernel_mod
