!-------------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!
! Section 20 aviation diagnostics kernel.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file currently belongs in section: physics_schemes_interface
! whilst discussions are ongoing about its final location.
!
module aviation_icao_heights_kernel_mod

  use argument_mod,         only: arg_type,            &
                                  gh_field, gh_scalar, gh_logical, &
                                  gh_read, gh_write, gh_integer, &
                                  gh_real, cell_column, &
                                  any_discontinuous_space_1, &
                                  any_discontinuous_space_2
  use constants_mod,        only: r_def, i_def, l_def
  use kernel_mod,           only: kernel_type

  implicit none

  ! the aviation diagnostics kernel type.
  type, extends(kernel_type) :: aviation_icao_heights_kernel_type
    type(arg_type), dimension(2) :: meta_args = (/ &

      ! output icao height.
      arg_type(gh_field, gh_real, gh_write, any_discontinuous_space_1), &

      ! input pressure
      arg_type(gh_field, gh_real, gh_read, any_discontinuous_space_1) &

      /)

    integer :: operates_on = cell_column

    contains
      procedure, nopass :: code => aviation_icao_heights_kernel_code
    end type aviation_icao_heights_kernel_type

contains

  subroutine aviation_icao_heights_kernel_code(nlayers, &
        icao_height, &  ! output icao height.
        pressure_field, &     ! pressure.
        ndf, undf, map)

    ! calculate iaco height from the pressure field.
    ! assumes lowest order w3 data, where ndf is always 1.

    use constants_mod,                  only: rmdi
    use planet_config_mod,              only: gravity, rd
    use science_conversions_mod,        only: feet_to_metres
    use science_aviation_constants_mod, only: mtokft, &
        isa_lapse_ratel, isa_lapse_rateu, isa_press_bot, &
        isa_press_mid, isa_press_top, isa_temp_bot, isa_temp_top, &
        gpm1, gpm2

    implicit none

    ! arguments (kernel)

    ! the number of layers in a column.
    integer(kind=i_def), intent(in) :: nlayers

    ! number of degrees of freedom (columns) in the cell we're processing.
    integer(kind=i_def), intent(in) :: ndf

    ! number of unique degrees of freedom in the fields.
    integer(kind=i_def), intent(in) :: undf

    ! degrees of freedom maps. offsets to the bottom of each column.
    integer(kind=i_def), intent(in), dimension(ndf) :: map


    ! arguments (algorithm)

    ! output icao height.
    real(kind=r_def), intent(inout), dimension(undf) :: icao_height

    ! pressure in pa.
    real(kind=r_def), intent(in), dimension(undf) :: pressure_field


    ! local variables
    integer(kind=i_def) :: df

    real(kind=r_def) :: g_over_rd
    real(kind=r_def) :: zp1
    real(kind=r_def) :: zp2

    real(kind=r_def) :: pressure


    g_over_rd = gravity / rd
    zp1 = isa_lapse_ratel / g_over_rd
    zp2 = isa_lapse_rateu / g_over_rd

    pressure = pressure_field(map(1))

    ! setting a safeguard limit to the lowest pressure to prevent
    ! extremely large icao height values near the top at the
    ! atmosphere.
    if ( (pressure >= 0.0_r_def) .and. (pressure <= 1000.0_r_def) ) then
        pressure = 1000.0_r_def
    end if

    ! pressure must not be greater than surface pressure.
    pressure = min(isa_press_bot, pressure)

    ! missing or invalid data?
    if (pressure == rmdi .or. pressure <= 0.0_r_def) then
        icao_height(map(1)) = rmdi

    ! heights up to 11,000 gpm
    else if ( pressure > isa_press_mid ) then
        pressure = pressure / isa_press_bot
        pressure = 1.0_r_def - pressure**zp1
        icao_height(map(1)) = pressure * isa_temp_bot / isa_lapse_ratel

    ! heights between 11,000 gpm and 20,000 gpm
    else if ( pressure > isa_press_top ) then
        pressure = pressure / isa_press_mid
        pressure = -log(pressure)
        icao_height(map(1)) = gpm1 + pressure * isa_temp_top / g_over_rd

    ! heights above 20,000 gpm
    else
        pressure = pressure / isa_press_top
        pressure = 1.0_r_def - pressure**zp2
        icao_height(map(1)) = gpm2 + pressure * isa_temp_top / isa_lapse_rateu
    end if

  end subroutine aviation_icao_heights_kernel_code

end module aviation_icao_heights_kernel_mod
