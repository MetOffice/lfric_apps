!-------------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
! Some of the content of this file has been produced with the assistance of
! Anthropic Claude Opus 5 (Claude Code).
!> @brief Height of the W3 levels above the surface.

module height_above_surface_kernel_mod

  use argument_mod,      only: arg_type,            &
                               GH_FIELD,            &
                               GH_READ, GH_WRITE,   &
                               GH_REAL, CELL_COLUMN
  use fs_continuity_mod, only: WTHETA, W3
  use constants_mod,     only: r_def, i_def
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: height_above_surface_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                 &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, W3),     & ! height_agl
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),     & ! height_w3
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA)  & ! height_wth
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: height_above_surface_code
  end type height_above_surface_kernel_type

  public :: height_above_surface_code

contains

  !> @details Heights are held above mean sea level, and the lowest Wtheta
  !!          level is the surface, so subtracting it from the W3 heights
  !!          gives the height of each W3 level above ground. The lowest W3
  !!          level sits part way up the first layer, so the result is
  !!          strictly positive rather than zero at the bottom of the column.
  !> @param[in]     nlayers     The number of layers
  !> @param[in,out] height_agl  Height of the W3 levels above the surface
  !> @param[in]     height_w3   Height of the W3 levels above mean sea level
  !> @param[in]     height_wth  Height of the Wtheta levels above mean sea level
  !> @param[in]     ndf_w3      Number of degrees of freedom per cell for W3
  !> @param[in]     undf_w3     Number of total degrees of freedom for W3
  !> @param[in]     map_w3      Dofmap for the cell at the base of the W3 column
  !> @param[in]     ndf_wth     Number of degrees of freedom per cell for Wtheta
  !> @param[in]     undf_wth    Number of total degrees of freedom for Wtheta
  !> @param[in]     map_wth     Dofmap for the cell at the base of the Wtheta
  !!                            column
  subroutine height_above_surface_code(nlayers,    &
                                       height_agl, &
                                       height_w3,  &
                                       height_wth, &
                                       ndf_w3,     &
                                       undf_w3,    &
                                       map_w3,     &
                                       ndf_wth,    &
                                       undf_wth,   &
                                       map_wth)

    implicit none

    ! Arguments added automatically in call to kernel
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_w3
    integer(kind=i_def), intent(in) :: undf_w3
    integer(kind=i_def), intent(in), dimension(ndf_w3)  :: map_w3
    integer(kind=i_def), intent(in) :: ndf_wth
    integer(kind=i_def), intent(in) :: undf_wth
    integer(kind=i_def), intent(in), dimension(ndf_wth) :: map_wth

    ! Arguments passed explicitly from algorithm
    real(kind=r_def), intent(inout), dimension(undf_w3)  :: height_agl
    real(kind=r_def), intent(in),    dimension(undf_w3)  :: height_w3
    real(kind=r_def), intent(in),    dimension(undf_wth) :: height_wth

    ! Internal variables
    integer(kind=i_def) :: k
    real(kind=r_def)    :: surface

    ! The lowest Wtheta level is the orography
    surface = height_wth(map_wth(1))

    do k = 0, nlayers - 1
      height_agl(map_w3(1)+k) = height_w3(map_w3(1)+k) - surface
    end do

  end subroutine height_above_surface_code

end module height_above_surface_kernel_mod
