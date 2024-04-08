!-------------------------------------------------------------------------------
!(c) Crown copyright 2020 Met Office. All rights reserved.
!The file LICENCE, distributed with this code, contains details of the terms
!under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Set the 3D rh_crit field from its namelist value
module initial_cloud_kernel_mod

  use argument_mod,      only: arg_type,          &
                               GH_FIELD, GH_REAL, &
                               GH_WRITE, CELL_COLUMN
  use fs_continuity_mod, only: Wtheta
  use constants_mod,     only: r_def, i_def
  use kernel_mod,        only: kernel_type

  use cloud_config_mod,  only: rh_crit

  implicit none

  private

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: initial_cloud_kernel_type
    private
    type(arg_type) :: meta_args(1) = (/                &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA) &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: initial_cloud_code
  end type initial_cloud_kernel_type

  public :: initial_cloud_code

contains

  !> @brief Set the 3D rh_crit field from the namelist value
  !> @param[in]     nlayers     The number of layers
  !> @param[in,out] rh_crit_wth Critical relative humidity
  !> @param[in]     ndf_wth     Number of degrees of freedom per cell for wtheta
  !> @param[in]     undf_wth    Number of total degrees of freedom for wtheta
  !> @param[in]     map_wth     Dofmap for the cell at the base of the column
  subroutine initial_cloud_code(nlayers,       &
                                rh_crit_wth,   &
                                ndf_wth,       &
                                undf_wth,      &
                                map_wth)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth
    integer(kind=i_def), intent(in) :: undf_wth
    integer(kind=i_def), intent(in),    dimension(ndf_wth)  :: map_wth
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: rh_crit_wth

    ! Internal variables
    integer(kind=i_def) :: k

    ! Namelist is only specified from level 1 upwards
    ! Set level 0 equal to level 1 value, although this is unused in
    ! the cloud scheme itself
    rh_crit_wth(map_wth(1)) = rh_crit(1)
    do k = 1, nlayers
      rh_crit_wth(map_wth(1) + k) = rh_crit(k)
    end do

  end subroutine initial_cloud_code

end module initial_cloud_kernel_mod
