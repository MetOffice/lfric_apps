!-------------------------------------------------------------------------------
!(c) Crown copyright 2021 Met Office. All rights reserved.
!The file LICENCE, distributed with this code, contains details of the terms
!under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Apply vertical smoothing (linear ramp/down) for SPT increments
module spt_levels_cap_kernel_mod

  use argument_mod,      only: arg_type, GH_FIELD, &
                               GH_WRITE, GH_REAL,  &
                               GH_SCALAR, GH_INTEGER, &
                               GH_READ, CELL_COLUMN


  use fs_continuity_mod, only: Wtheta
  use constants_mod,     only: r_def, i_def, r_um
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> Metadata describing the kernel to PSyclone
  !>
  type, public, extends(kernel_type) :: spt_levels_cap_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                 &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & !dX
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),      & ! spt_level_bottom
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),      & ! spt_level_top
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),      & ! spt_level_begin_tapering_bottom
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ)       & ! spt_level_begin_tapering_top
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: spt_levels_cap_code
  end type spt_levels_cap_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public spt_levels_cap_code
contains

  !> @brief Apply vertical ramp up/down on SPT buffer verfical levels
  !> @details Apply the vertical ramp up (down) over SPT increments in the
  !>          SPT vertical levels from spt_level_bottom to spt_level_begin_tapering_bottom
  !>          (from spt_level_begin_tapering_top to spt_level_top)
  !> @param[in]     nlayers        The number of layers
  !> @param[in,out] dX             Field with SPT perturbations
  !> @param[in]     ndf_wth        Number of degrees of freedom per cell for wtheta
  !> @param[in]     undf_wth       Number of total degrees of freedom for wtheta
  !> @param[in]     map_wth        Dofmap for the cell at the base of the column
  !> @param[in]     spt_level_bottom      Bottom level of the stochastic scheme
  !> @param[in]     spt_level_top         Top level of the stochastic scheme
  !> @param[in]     spt_level_begin_tapering_bottom        spt_level_begin_tapering_bottom in stochastic_physics_config_mod
  !> @param[in]     spt_level_begin_tapering_top           spt_level_begin_tapering_top in stochastic_physics_config_mod

  subroutine spt_levels_cap_code(nlayers,                         &
                                 dX,                              &
                                 spt_level_bottom,                &
                                 spt_level_top,                   &
                                 spt_level_begin_tapering_bottom, &
                                 spt_level_begin_tapering_top,    &
                                 ndf_wth,                         &
                                 undf_wth,                        &
                                 map_wth                          &
                                 )

    implicit none

    !Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth
    integer(kind=i_def), intent(in) :: undf_wth
    integer(kind=i_def), intent(in), dimension(ndf_wth)  :: map_wth
    integer(kind=i_def), intent(in) :: spt_level_bottom
    integer(kind=i_def), intent(in) :: spt_level_top
    integer(kind=i_def), intent(in) :: spt_level_begin_tapering_bottom
    integer(kind=i_def), intent(in) :: spt_level_begin_tapering_top

    ! field with perturbation
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dX

    real(kind=r_def) :: lev_amp_factor
    integer(kind=i_def) :: k

    ! Set to zero below spt_level_bottom
    do k = 0,spt_level_bottom-1
      dX(map_wth(1) + k) = 0.0_r_def
    end do

    ! Ramp up from spt_level_bottom to spt_level_begin_tapering_bottom
    do k= spt_level_bottom, spt_level_begin_tapering_bottom
      lev_amp_factor=  REAL(k - spt_level_bottom, r_def)/ &
                       REAL(spt_level_begin_tapering_bottom - spt_level_bottom, r_def)

      dX(map_wth(1) + k) = lev_amp_factor*dX(map_wth(1) + k)
    end do

    ! Ramp up from spt_level_begin_tapering_top to spt_level_top
    do k= spt_level_begin_tapering_top, spt_level_top
      lev_amp_factor=REAL(k - spt_level_begin_tapering_top, r_def)/ &
                     REAL(spt_level_top - spt_level_begin_tapering_top, r_def)

      dX(map_wth(1) + k) = lev_amp_factor*dX(map_wth(1) + k)
    end do

    ! Set to zero above  spt_level_top
    do k = spt_level_top+1, nlayers
      dX(map_wth(1) + k) = 0.0_r_def
    end do

  end subroutine spt_levels_cap_code

end module spt_levels_cap_kernel_mod
