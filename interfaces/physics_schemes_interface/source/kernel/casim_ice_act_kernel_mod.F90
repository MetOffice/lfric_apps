!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Activate ice number due to convection
!>
module casim_ice_act_kernel_mod

  use argument_mod,       only : arg_type,              &
                                 GH_FIELD, GH_REAL,     &
                                 GH_READ, GH_READWRITE, &
                                 GH_WRITE, CELL_COLUMN

  use constants_mod,      only : r_def, r_double, i_def, i_um, r_um
  use fs_continuity_mod,  only : Wtheta
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: casim_ice_act_kernel_type
    private
    type(arg_type) :: meta_args(6) = (/                   &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),   & ! dms_conv
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),   & ! ni_mphys
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),   & ! theta
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),   & ! exner
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),   & ! rho
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA)    & ! cf_fro
        /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: casim_ice_act_code
  end type

  public :: casim_ice_act_code

contains

  !> @details Determines the ice number of newly detrained ice mass from
  !>          convection, using the relation of Cooper (1986)
  !>          Ice Initiation in Natural Clouds
  !>          Meteor. Monogr. (1986) 43: 29-32.
  !>          https://doi.org/10.1175/0065-9401-21.43.29
  !> @param[in]     nlayers       Number of layers
  !> @param[in]     dms_conv      Cloud snow increment from convection
  !> @param[in,out] ni_mphys      CASIM cloud-ice number concentration
  !> @param[in]     theta         Potential temperature
  !> @param[in]     exner         Exner pressure
  !> @param[in]     rho           Dry density
  !> @param[in]     cf_fro        Frozen cloud fraction
  !> @param[in]     ndf_wth       Number of degrees of freedom per cell for potential temperature space
  !> @param[in]     undf_wth      Number unique of degrees of freedom for potential temperature space
  !> @param[in]     map_wth       Dofmap for the cell at the base of the column for potential temperature space
  subroutine casim_ice_act_code(nlayers,      &
                                dms_conv,     &
                                ni_mphys,     &
                                theta,        &
                                exner,        &
                                rho,          &
                                cf_fro,       &
                                ndf_wth,      &
                                undf_wth,     &
                                map_wth)

    !---------------------------------------
    ! UM modules
    !---------------------------------------
    use thresholds,           only: ni_tidy
    use water_constants_mod,  only: tm

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)     :: nlayers
    integer(kind=i_def), intent(in)     :: ndf_wth
    integer(kind=i_def), intent(in)     :: undf_wth

    integer(kind=i_def), intent(in),    dimension(ndf_wth)  :: map_wth

    real(kind=r_def),    intent(in), dimension(undf_wth) :: dms_conv
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: ni_mphys
    real(kind=r_def),    intent(in), dimension(undf_wth) :: theta
    real(kind=r_def),    intent(in), dimension(undf_wth) :: exner
    real(kind=r_def),    intent(in), dimension(undf_wth) :: rho
    real(kind=r_def),    intent(in), dimension(undf_wth) :: cf_fro

    ! Local variables for the kernel
    integer(i_um) :: k

    real(r_def), parameter :: max_icenumber = 1.0e6_r_def
    real(r_def), parameter :: cooper_a = 5.0_r_def
    real(r_def), parameter :: cooper_b = -0.304_r_def

    real(r_def) :: temperature, cooper_ice_number

    do k = 1, nlayers

      temperature = theta(map_wth(1)+k) * exner(map_wth(1)+k)

      if (temperature < tm) then

        ! If there is ice but no number
        if (dms_conv(map_wth(1) + k) > 0.0_r_def .and. &
             ni_mphys(map_wth(1) + k) <= ni_tidy) then

          ! Use Cooper relation for ice number
          cooper_ice_number = cooper_a * exp(cooper_b * (temperature-tm)) &
                              / rho(map_wth(1) + k)

          ! Limit number to sensible maximum and scale to grid-box mean
          ni_mphys( map_wth(1) + k) = min(max_icenumber, cooper_ice_number) &
                                    * cf_fro(map_wth(1) + k)
        end if
      end if
    end do

    ! Save value of ni_mphys at level 1 for level 0 increment
    ni_mphys( map_wth(1) + 0) = ni_mphys( map_wth(1) + 1)

  end subroutine casim_ice_act_code

end module casim_ice_act_kernel_mod
