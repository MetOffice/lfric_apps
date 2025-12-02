!-------------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Initialise Jules soil fields on soil levels
!> @details Non-standard Surface fields (pseudo-levels) aren't as yet not
!>  implemented in LFRic. As an interim measure Higher-order W3 fields have
!>  been used to mimic psuedo-level field behaviour. This code is written
!>  based on this interim measure and will need to be updated when
!>  suitable infrastructure is available (Ticket #2081)
module initial_soil_kernel_mod

  use argument_mod,  only: arg_type,              &
                           GH_FIELD, GH_REAL,     &
                           GH_WRITE, CELL_COLUMN, &
                           ANY_DISCONTINUOUS_SPACE_1
  use constants_mod, only: r_def, i_def
  use kernel_mod,    only: kernel_type

  use nlsizes_namelist_mod, only: sm_levels
  use ideal_surface_config_mod, only: soil_temperature_in => soil_temperature, &
                                      soil_moisture_in => soil_moisture

  implicit none

  private

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: initial_soil_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                                    &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1)  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: initial_soil_code
  end type initial_soil_kernel_type

  public :: initial_soil_code

contains

  !> @param[in]     nlayers                The number of layers
  !> @param[in,out] soil_temperature       Soil temperature (K)
  !> @param[in,out] soil_moisture          Soil moisture content (kg m-2)
  !> @param[in]     ndf_soil               Number of DOFs per cell for soil levels
  !> @param[in]     undf_soil              Number of total DOFs for soil levels
  !> @param[in]     map_soil               Dofmap for cell for soil levels
  subroutine initial_soil_code(nlayers,                       &
                               soil_temperature,              &
                               soil_moisture,                 &
                               ndf_soil, undf_soil, map_soil)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_soil, undf_soil
    integer(kind=i_def), intent(in) :: map_soil(ndf_soil)

    real(kind=r_def), intent(inout) :: soil_temperature(undf_soil)
    real(kind=r_def), intent(inout) :: soil_moisture(undf_soil)

    integer(i_def) :: k

    ! Copy idealised namelist values into initial prognostic fields

    do k = 1, sm_levels
      soil_temperature(map_soil(1)+k-1) = soil_temperature_in(k)
      soil_moisture(map_soil(1)+k-1) = soil_moisture_in(k)
    end do

  end subroutine initial_soil_code

end module initial_soil_kernel_mod
