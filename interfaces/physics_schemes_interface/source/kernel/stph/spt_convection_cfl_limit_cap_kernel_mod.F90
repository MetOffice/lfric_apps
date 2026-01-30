!-------------------------------------------------------------------------------
!(c) Crown copyright 2021 Met Office. All rights reserved.
!The file LICENCE, distributed with this code, contains details of the terms
!under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Remove those points where the CFL is breached in the convection tendencies
module spt_convection_cfl_limit_cap_kernel_mod

  use argument_mod,      only: arg_type, GH_FIELD, &
                               GH_WRITE, GH_REAL,  &
                               GH_SCALAR, GH_INTEGER, &
                               GH_READ, CELL_COLUMN
  use fs_continuity_mod, only: Wtheta
  use constants_mod,     only: r_def, i_def, l_def, &
                               r_second
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> Metadata describing the kernel to PSyclone
  !>
  type, public, extends(kernel_type) :: spt_convection_cfl_limit_cap_kernel_type
    private
    type(arg_type) :: meta_args(7) = (/                 &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! dX_conv_cfl
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! massflux_up
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! fp_spt
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! pressure
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),      & ! spt_level_bottom
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),      & ! spt_level_top
         arg_type(GH_SCALAR, GH_REAL, GH_READ)          & ! dt
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: spt_convection_cfl_limit_cap_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public spt_convection_cfl_limit_cap_code
contains

  !> @brief Remove those points where the CFL is breached in SPT convective Tendencies.
  !> @details Apply the CFL criteria to SPT convective tendencies and
  !>          remove those points breaching the criteria
  !> @param[in]      nlayers      The number of layers
  !> @param[in,out]  dX_conv_cfl  The convection inc. after applying CFL
  !> @param[in]      massflux_up  Convection mass flux
  !> @param[in]      fp_spt       SPT forcing pattern
  !> @param[in]      pressure     Pressure field on theta levels
  !> @param[in]      ndf_wth      Number of degrees of freedom per cell for wtheta
  !> @param[in]      undf_wth     Number of total degrees of freedom for wtheta
  !> @param[in]      map_wth      Dofmap for the cell at the base of the column
  !> @param[in]      spt_level_bottom      Bottom level of the stochastic scheme
  !> @param[in]      spt_level_top         Top level of the stochastic scheme
  !> @param[in]      dt                    Timestep from timestepping_config_mod
  subroutine spt_convection_cfl_limit_cap_code(nlayers,          &
                                               dX_conv_cfl,      &
                                               massflux_up,      &
                                               fp_spt,           &
                                               pressure,         &
                                               spt_level_bottom, &
                                               spt_level_top,    &
                                               dt,               &
                                               ndf_wth,          &
                                               undf_wth,         &
                                               map_wth)

    implicit none

    !Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth
    integer(kind=i_def), intent(in) :: undf_wth
    integer(kind=i_def), intent(in), dimension(ndf_wth)  :: map_wth
    integer(kind=i_def), intent(in) :: spt_level_bottom
    integer(kind=i_def), intent(in) :: spt_level_top
    real(kind=r_second), intent(in) :: dt

    ! Fields perturbations + tendencies
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dX_conv_cfl
    real(kind=r_def), intent(in),    dimension(undf_wth) :: massflux_up
    real(kind=r_def), intent(in),    dimension(undf_wth) :: fp_spt
    real(kind=r_def), intent(in),    dimension(undf_wth) :: pressure

    !Internal variables
    integer(kind=i_def) :: k
    real(kind=r_def) :: cfl_value
    logical(kind=l_def) :: mask_value

    ! Set mask to False
    mask_value = .false.

    ! Set mask to True if CFL criteria is met
    do k = spt_level_bottom+1, spt_level_top-1
      ! The CFL is breached if the massflux scaled by the forcing pattern is higher
      ! than the vertical pressure gadient
      cfl_value=massflux_up(map_wth(1) + k) * (1.0_r_def + fp_spt(map_wth(1) + k))*dt
      if (cfl_value >= abs(pressure(map_wth(1) + k) - pressure(map_wth(1) + k-1) ) ) then
        mask_value = .true.
        exit
      end if
    end do

    ! Remove tendencies where criteria is breached
    if (mask_value) then
      do k = 1, nlayers
        dX_conv_cfl(map_wth(1) + k) = 0.0_r_def
      end do
    end if

  end subroutine spt_convection_cfl_limit_cap_code

end module spt_convection_cfl_limit_cap_kernel_mod
