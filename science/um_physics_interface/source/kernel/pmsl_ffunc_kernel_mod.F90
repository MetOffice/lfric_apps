!-------------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Holds code to calculate forcing function of smoothed PMSL

module pmsl_ffunc_kernel_mod

  use argument_mod,         only: arg_type,                  &
                                  GH_FIELD, GH_SCALAR,       &
                                  GH_READ, GH_WRITE,         &
                                  GH_REAL, CELL_COLUMN,      &
                                  ANY_DISCONTINUOUS_SPACE_1, &
                                  ANY_DISCONTINUOUS_SPACE_9, &
                                  STENCIL, CROSS2D
  use fs_continuity_mod,    only: WTHETA, W2
  use constants_mod,        only: r_def, i_def
  use kernel_mod,           only: kernel_type

  implicit none

  private

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: pmsl_ffunc_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                    &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W2),                        & ! dx_at_w2
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1, STENCIL(CROSS2D)), & ! f_vg
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1, STENCIL(CROSS2D)), & ! f_ug
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_9, STENCIL(CROSS2D)), & ! panel_id
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1)  & ! f_func
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: pmsl_ffunc_code
  end type pmsl_ffunc_kernel_type

  public :: pmsl_ffunc_code

contains

  !> @brief Calculate total forcing function
  !> @details Basic formation is based on the old UM diagnostic which is
  !>          described in UM Documentation Paper 80.
  !>          https://code.metoffice.gov.uk/doc/um/latest/papers/umdp_080.pdf
  !>          This kernel calculates the total forcing function from the
  !>          geostrophic winds, against which the PMSL will be relaxed
  !> @param[in]     nlayers       The number of layers
  !> @param[in]     dx_at_w2       cell sizes at w2 dofs
  !> @param[in]     f_vg          'meridional' geostrophic forcing function
  !> @param[in]     smap_vg_size  Size of the stencil map in each direction
  !> @param[in]     vg_len        Max size of the stencil map in any direction
  !> @param[in]     smap_vg       Stencil map
  !> @param[in]     f_ug          'zonal' geostrophic forcing function
  !> @param[in]     smap_ug_size  Size of the stencil map in each direction
  !> @param[in]     ug_len        Max size of the stencil map in any direction
  !> @param[in]     smap_ug       Stencil map
  !> @param[in]     panel_id      The ID number of the current panel
  !> @param[in]     smap_pid_size Size of the stencil map in each direction
  !> @param[in]     pid_len       Max size of the stencil map in any direction
  !> @param[in]     smap_pid      Stencil map
  !> @param[in,out] f_func        total forcing function to relax against PMSL
  !> @param[in]     ndf_w2        Number of degrees of freedom per cell for w2 fields
  !> @param[in]     undf_w2       Number of total degrees of freedom for w2 fields
  !> @param[in]     map_w2        Dofmap for the cell at the base of the column for w2 fields
  !> @param[in]     ndf_2d        Number of degrees of freedom per cell for 2d fields
  !> @param[in]     undf_2d       Number of total degrees of freedom for 2d fields
  !> @param[in]     map_2d        Dofmap for the cell at the base of the column for 2d fields
  subroutine pmsl_ffunc_code(nlayers,                           &
                             dx_at_w2,                          &
                             f_vg,                              &
                             smap_vg_size, vg_len, smap_vg,     &
                             f_ug,                              &
                             smap_ug_size, ug_len, smap_ug,     &
                             panel_id,                          &
                             smap_pid_size, pid_len, smap_pid,  &
                             f_func,                            &
                             ndf_w2, undf_w2, map_w2,           &
                             ndf_2d, undf_2d, map_2d,           &
                             ndf_pid, undf_pid, map_pid)

    implicit none

    ! Arguments added automatically in call to kernel
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in), dimension(ndf_2d)  :: map_2d

    integer(kind=i_def), intent(in) :: ndf_pid, undf_pid
    integer(kind=i_def), intent(in), dimension(ndf_pid)  :: map_pid

    integer(kind=i_def), intent(in) :: vg_len
    integer(kind=i_def), dimension(4), intent(in) :: smap_vg_size
    integer(kind=i_def), dimension(ndf_2d,vg_len,4), intent(in) :: smap_vg

    integer(kind=i_def), intent(in) :: ug_len
    integer(kind=i_def), dimension(4), intent(in) :: smap_ug_size
    integer(kind=i_def), dimension(ndf_2d,ug_len,4), intent(in) :: smap_ug

    integer(kind=i_def), intent(in) :: pid_len
    integer(kind=i_def), dimension(4), intent(in) :: smap_pid_size
    integer(kind=i_def), dimension(ndf_2d,pid_len,4), intent(in) :: smap_pid

    integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
    integer(kind=i_def), dimension(ndf_w2), intent(in)  :: map_w2

    ! Arguments passed explicitly from algorithm
    real(kind=r_def),    intent(in), dimension(undf_w2) :: dx_at_w2
    real(kind=r_def),    intent(in), dimension(undf_2d) :: f_vg
    real(kind=r_def),    intent(in), dimension(undf_2d) :: f_ug
    real(kind=r_def),    intent(in), dimension(undf_pid) :: panel_id
    real(kind=r_def),    intent(inout), dimension(undf_2d) :: f_func

    ! Internal variables
    real(kind=r_def) :: idx, idy
    integer(kind=i_def) :: xp1, xm1, yp1, ym1

    integer(kind=i_def) :: cell_panel, stencil_cell, panel_edge, stencil_panel
    real(kind=r_def),    dimension(undf_2d) :: loc_vg
    real(kind=r_def),    dimension(undf_2d) :: loc_ug

    ! Orientation of geostrophic winds will change when we cross over a panel
    ! Take this into account by creating local versions which are correctly
    ! oriented
    cell_panel = int(panel_id(map_pid(1)), i_def)

    do stencil_cell = 1, 4
      stencil_panel = int(panel_id(smap_pid(1,smap_pid_size(stencil_cell),stencil_cell)), i_def)
      ! Create panel_edge to check whether a panel is changing
      panel_edge = 10*cell_panel + stencil_panel

      select case (panel_edge)
      case (41, 32, 16, 25, 64, 53)
        ! Clockwise rotation of panel
        ! u is v on adjacent panel
        ! v is -u on adjacent panel
        loc_vg(smap_vg(1,smap_vg_size(stencil_cell),stencil_cell)) = &
         -f_ug(smap_ug(1,smap_ug_size(stencil_cell),stencil_cell))
        loc_ug(smap_ug(1,smap_ug_size(stencil_cell),stencil_cell)) = &
          f_vg(smap_vg(1,smap_vg_size(stencil_cell),stencil_cell))
      case (14, 23, 61, 52, 46, 35)
        ! Anti-clockwise rotation of panel
        ! u is -v on adjacent panel
        ! v is u on adjacent panel
        loc_vg(smap_vg(1,smap_vg_size(stencil_cell),stencil_cell)) = &
          f_ug(smap_ug(1,smap_ug_size(stencil_cell),stencil_cell))
        loc_ug(smap_ug(1,smap_ug_size(stencil_cell),stencil_cell)) = &
         -f_vg(smap_vg(1,smap_vg_size(stencil_cell),stencil_cell))
      case default
        ! Same panel or crossing panel with no rotation, so stencil map is unchanged
        loc_vg(smap_vg(1,smap_vg_size(stencil_cell),stencil_cell)) = &
          f_vg(smap_vg(1,smap_vg_size(stencil_cell),stencil_cell))
        loc_ug(smap_ug(1,smap_ug_size(stencil_cell),stencil_cell)) = &
          f_ug(smap_ug(1,smap_ug_size(stencil_cell),stencil_cell))
      end select
    end do

    ! Calculate which cell in the x branch of the stencil to use
    ! This sets the point to use to be the stencil point (2) if it exists,
    ! or the centre point (1) if it doesn't (i.e. we are at a domain edge)
    xp1 = smap_vg_size(3)
    xm1 = smap_vg_size(1)
    ! Calculate which cell in the y branch of the stencil to use
    ! This sets the point to use to be the stencil point (2) if it exists,
    ! or the centre point (1) if it doesn't (i.e. we are at a domain edge)
    yp1 = smap_ug_size(4)
    ym1 = smap_ug_size(2)

    ! Calculate inverse dx values - distance between centres of adjacent cells
    idx = 1.0_r_def / (dx_at_w2(map_w2(1))+dx_at_w2(map_w2(3)))
    idy = 1.0_r_def / (dx_at_w2(map_w2(2))+dx_at_w2(map_w2(4)))

    ! Calulate F function
    ! Because this is obtained from a curl function, there should be a minus
    ! sign applied to the ug term. However, this cancels with the minus
    ! sign that was omitted in the original calculation of the ug term
    f_func(map_2d(1)) = idx* (loc_vg(smap_vg(1,xp1,3)) - loc_vg(smap_vg(1,xm1,1))) &
                      + idy* (loc_ug(smap_ug(1,yp1,4)) - loc_ug(smap_ug(1,ym1,2)))

  end subroutine pmsl_ffunc_code

end module pmsl_ffunc_kernel_mod
