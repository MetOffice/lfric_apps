!-------------------------------------------------------------------------------
!(c) Crown copyright 2023 Met Office. All rights reserved.
!The file LICENCE, distributed with this code, contains details of the terms
!under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Compute the SKEB2 divergence and vorticity
module skeb_vort_div_kernel_mod

  use argument_mod,      only: arg_type, GH_FIELD,          &
                               GH_REAL, GH_WRITE, GH_READ,  &
                               CELL_COLUMN, GH_INTEGER,     &
                               GH_SCALAR, STENCIL, CROSS,   &
                               ANY_DISCONTINUOUS_SPACE_9
  use fs_continuity_mod, only: W3, Wtheta, W1, W2
  use constants_mod,     only: r_def, i_def
  use kernel_mod,        only: kernel_type

  implicit none

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: skeb_vort_div_kernel_type
    private
    type(arg_type) :: meta_args(11) = (/                       &
    arg_type(GH_FIELD, GH_REAL, GH_WRITE, W1),                 & ! vorticity
    arg_type(GH_FIELD, GH_REAL, GH_WRITE, W3),                 & ! divergence
    arg_type(GH_FIELD, GH_REAL, GH_READ,  W2, STENCIL(CROSS)), & ! u_n
    arg_type(GH_FIELD, GH_REAL, GH_READ,  W2),                 & ! dx_in_w2
    arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                 & ! wetrho_w3
    arg_type(GH_FIELD, GH_REAL, GH_READ,  WTheta),             & ! wetrho_wth
    arg_type(GH_FIELD, GH_REAL, GH_READ,  WTheta),             & ! w_in_wth
    arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                 & ! rdz_w3
    arg_type(GH_SCALAR, GH_INTEGER, GH_READ ),                 & ! skeb_level_bottom
    arg_type(GH_SCALAR, GH_INTEGER, GH_READ ),                 & ! skeb_level_top
    arg_type(GH_FIELD,   GH_REAL, GH_READ,   ANY_DISCONTINUOUS_SPACE_9, &
                                                       STENCIL(CROSS))  & !panel_id
    /)
    integer :: operates_on = CELL_COLUMN

  contains
    procedure, nopass :: skeb_vort_div_code
  end type skeb_vort_div_kernel_type

  public skeb_vort_div_code
contains

  !> @brief Calculate divergence and vorticity for skeb
  !> @param[in]    nlayers     The number of layers
  !> @param[inout] vorticity   2D Vorticity
  !> @param[inout] divergence  2D Divergence
  !> @param[in]    u_n         Physical wind vector
  !> @param[in]    map_w2_sten_size Size of the stencil for w2 fields
  !> @param[in]    map_w2_sten      Stencil map for w2 fields
  !> @param[in]    dx_in_w2    delta-x and w2 dofs
  !> @param[in]    wetrho_w3   wet density in w3
  !> @param[in]    wetrho_wth  wet density in wth
  !> @param[in]    w_in_wth    vertical velocity in wth
  !> @param[in]    rdz_w3      Spacing of potential temperature space levels above surface
  !> @param[in]    skeb_level_bottom Bottom SKEB level
  !> @param[in]    skeb_level_top    Top SKEB level
  !> @param[in]    panel_id    ID number of the current panel
  !> @param[in]    map_pid_sten_size Size of the panel ID stencil map
  !> @param[in]    map_pid_sten      Stencil map for the panel ID
  !> @param[in]    ndf_w1      Number of DOFs per cell for w1 space
  !> @param[in]    undf_w1     Number of unique DOFs  for w1 space
  !> @param[in]    map_w1      dofmap for the cell at the base of the column for w1 space
  !> @param[in]    ndf_w3      Number of DOFs per cell for density space
  !> @param[in]    undf_w3     Number of unique DOFs  for density space
  !> @param[in]    map_w3      dofmap for the cell at the base of the column for density space
  !> @param[in]    ndf_w2      Number of DOFs per cell for w2 space
  !> @param[in]    undf_w2     Number of unique DOFs  for w2 space
  !> @param[in]    map_w2      dofmap for the cell at the base of the column for w2 space
  !> @param[in]    ndf_wth     Number of DOFs per cell for potential temperature space
  !> @param[in]    undf_wth    Number of unique DOFs for potential temperature space
  !> @param[in]    map_wth     dofmap for the cell at the base of the column for potential temperature space

  subroutine skeb_vort_div_code(nlayers,           &
                                vorticity,         &
                                divergence,        &
                                u_n,               &
                                map_w2_sten_size,  &
                                map_w2_sten,       &
                                dx_at_w2,          &
                                wetrho_w3,         &
                                wetrho_wth,        &
                                w_in_wth,          &
                                rdz_w3,            &
                                skeb_level_bottom, &
                                skeb_level_top,    &
                                panel_id,          &
                                map_pid_sten_size, &
                                map_pid_sten,      &
                                ndf_w1,            &
                                undf_w1,           &
                                map_w1,            &
                                ndf_w3,            &
                                undf_w3,           &
                                map_w3,            &
                                ndf_w2,            &
                                undf_w2,           &
                                map_w2,            &
                                ndf_wth,           &
                                undf_wth,          &
                                map_wth,           &
                                ndf_pid,           &
                                undf_pid,          &
                                map_pid)

    implicit none

    !Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth, ndf_w3, ndf_w2, ndf_w1, ndf_pid
    integer(kind=i_def), intent(in) :: undf_wth, undf_w3, undf_w2, undf_w1, undf_pid
    integer(kind=i_def), intent(in) :: map_w2_sten_size, map_pid_sten_size
    integer(kind=i_def), intent(in),  dimension(ndf_wth)  :: map_wth
    integer(kind=i_def), intent(in),  dimension(ndf_w3)   :: map_w3
    integer(kind=i_def), intent(in),  dimension(ndf_w2)   :: map_w2
    integer(kind=i_def), intent(in),  dimension(ndf_w2, map_w2_sten_size) :: map_w2_sten
    integer(kind=i_def), intent(in),  dimension(ndf_w1)   :: map_w1
    integer(kind=i_def), intent(in),  dimension(ndf_pid)  :: map_pid
    integer(kind=i_def), intent(in),  dimension(ndf_pid, map_pid_sten_size) :: map_pid_sten

    ! Fields
    real(kind=r_def),    intent(inout), dimension(undf_w1)  :: vorticity
    real(kind=r_def),    intent(inout), dimension(undf_w3)  :: divergence
    real(kind=r_def),    intent(in),    dimension(undf_w2)  :: u_n
    real(kind=r_def),    intent(in),    dimension(undf_w2)  :: dx_at_w2
    real(kind=r_def),    intent(in),    dimension(undf_w3)  :: rdz_w3
    real(kind=r_def),    intent(in),    dimension(undf_w3)  :: wetrho_w3
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: wetrho_wth
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: w_in_wth
    real(kind=r_def),    intent(in),    dimension(undf_pid) :: panel_id

    ! Scalars
    integer(kind=i_def), intent(in) :: skeb_level_bottom, skeb_level_top

    integer(kind=i_def) :: true_w2_map(ndf_w2, map_w2_sten_size)
    integer(kind=i_def) :: stencil_cell
    integer(kind=i_def) :: cell_panel, stencil_panel, panel_edge
    integer(kind=i_def) :: vec_dir(ndf_w2, map_w2_sten_size)

    integer(kind=i_def) :: k, df

    ! The W2H DoF values change in orientation when we cross over a panel
    ! Vector directions parallel to the boundary (i.e. the winds on faces
    ! perpendicular to the boundary) also flip sign
    ! We need to take this into account by adjusting the stencil map used for
    ! the wind field. Do this by looking at whether the panel changes
    ! for other cells in the stencil
    cell_panel = int(panel_id(map_pid(1)), i_def)

    do stencil_cell = 1, map_w2_sten_size
      stencil_panel = int(panel_id(map_pid_sten(1, stencil_cell)), i_def)
      ! Create panel_edge to check whether a panel is changing
      panel_edge = 10*cell_panel + stencil_panel

      select case (panel_edge)
      case (41, 32, 16, 25, 64, 53)
        ! Clockwise rotation of panel
        true_w2_map(1, stencil_cell) = map_w2_sten(2, stencil_cell)
        true_w2_map(2, stencil_cell) = map_w2_sten(3, stencil_cell)
        true_w2_map(3, stencil_cell) = map_w2_sten(4, stencil_cell)
        true_w2_map(4, stencil_cell) = map_w2_sten(1, stencil_cell)
        ! Flip direction of vectors if necessary
        vec_dir(1,stencil_cell) = -1_i_def
        vec_dir(2,stencil_cell) = 1_i_def
        vec_dir(3,stencil_cell) = -1_i_def
        vec_dir(4,stencil_cell) = 1_i_def
      case (14, 23, 61, 52, 46, 35)
        ! Anti-clockwise rotation of panel
        true_w2_map(1, stencil_cell) = map_w2_sten(4, stencil_cell)
        true_w2_map(2, stencil_cell) = map_w2_sten(1, stencil_cell)
        true_w2_map(3, stencil_cell) = map_w2_sten(2, stencil_cell)
        true_w2_map(4, stencil_cell) = map_w2_sten(3, stencil_cell)
        ! Flip direction of vectors if necessary
        vec_dir(1,stencil_cell) = 1_i_def
        vec_dir(2,stencil_cell) = -1_i_def
        vec_dir(3,stencil_cell) = 1_i_def
        vec_dir(4,stencil_cell) = -1_i_def
      case default
        ! Same panel or crossing panel with no rotation, so stencil map is unchanged
        true_w2_map(:, stencil_cell) = map_w2_sten(:, stencil_cell)
        vec_dir(:,stencil_cell) = 1_i_def
      end select
    end do

    ! --------------------------------------------------------------------
    ! Compute Vorticity field.
    ! --------------------------------------------------------------------
    !  vort=dv/dx-du/dy=[v(i+1,j)-v(i,j)]/dx-[u(i,j+1)-u(i,j)]/dy
    !  but because beta = -y we actually do -dv/dx-du/dy
    !
    !
    !                       u(3,5)
    !                          |
    !             v(4) --- vort(7)--- v(4,4)
    !                          |
    !                       u(3)
    !
    !                       u(1)
    !                          |
    !            v(2,2) --- vort(5)--- v(2)
    !                          |
    !                       u(1,3)

    do df = 1,3,2

      ! Only calculate this dof if it hasn't already been done
      if (vorticity(map_w1(df+4)+skeb_level_bottom) == 0.0_r_def) then

        do k = skeb_level_bottom, skeb_level_top

          vorticity(map_w1(df+4)+k-1) = ((-1)**((df+1)/2)) *                   &
                          (vec_dir(df+1,df+1)*u_n(true_w2_map(df+1,1)+k-1)-    &
                                              u_n(true_w2_map(df+1,df+1)+k-1)) &
                                        / dx_at_w2(map_w2(df)+k-1) -           &
                                        ((-1)**((df+1)/2)) *                   &
                          (vec_dir(df,df+2)*u_n(true_w2_map(df,df+2)+k-1) -    &
                                            u_n(true_w2_map(df,1)+k-1))        &
                                        / dx_at_w2(map_w2(df+1)+k-1)

        end do

      end if

    end do

    ! --------------------------------------------------------------------
    ! Compute Divergence field.
    ! --------------------------------------------------------------------
    !  div = du/dx + dv/dy = [u(i,j)-u(i-1,j)]/dx + [v(i,j)-v(i,j-1)]/dy =
    !       = - 1/rho x d(rho x w)/dz
    !

    do k = skeb_level_bottom+1, skeb_level_top-1

      divergence(map_w3(1)+k-1) = ( wetrho_wth(map_wth(1)+k) *   &
                                    w_in_wth(map_wth(1)+k) -     &
                                    wetrho_wth(map_wth(1)+k-1) * &
                                    w_in_wth(map_wth(1)+k-1) )   &
                                  * rdz_w3(map_w3(1)+k-1)
      divergence(map_w3(1)+k-1) = (-1.0_r_def) * divergence(map_w3(1)+k-1) &
                                  / wetrho_w3(map_w3(1)+k-1)

    end do

  end subroutine skeb_vort_div_code

end module skeb_vort_div_kernel_mod
