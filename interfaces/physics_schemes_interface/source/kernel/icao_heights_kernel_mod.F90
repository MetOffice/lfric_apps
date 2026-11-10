! -----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! -----------------------------------------------------------------------------
!> @brief Calculate icao height.
!>
module icao_heights_kernel_mod

  use argument_mod,         only: arg_type,                 &
                                  gh_field, gh_scalar,      &
                                  gh_read, gh_write,        &
                                  gh_real, cell_column,     &
                                  any_discontinuous_space_1
  use constants_mod,        only: r_def, i_def
  use kernel_mod,           only: kernel_type

  implicit none

  ! The aviation diagnostics kernel type.
  type, extends(kernel_type) :: icao_heights_kernel_type
    type(arg_type), dimension(3) :: meta_args = (/ &

      ! Output icao height.
      arg_type(gh_field, gh_real, gh_write, any_discontinuous_space_1), &

      ! Input pressure
      arg_type(gh_field, gh_real, gh_read, any_discontinuous_space_1), &

      ! g_over_r
      arg_type(gh_scalar, gh_real, gh_read) &

      /)

    integer :: operates_on = cell_column

    contains
      procedure, nopass :: code => icao_heights_kernel_code
    end type icao_heights_kernel_type

contains

  !> @brief Calculate icao height from the pressure field.
  !> @details Assumes lowest order W3 data, where ndf is always 1.
  !> @param[in]     nlayers         The number of layers in a column.
  !> @param[out]    icao_height     Output icao height in kft.
  !> @param[in]     pressure_field  Pressure in pa.
  !> @param[in]     g_over_r        Constant must be passed in.
  !> @param[in]     ndf             Number of DOFs in the cell.
  !> @param[in]     undf            Number of DOFs in the field.
  !> @param[in]     map             Dofmap to the bottom cell.
  subroutine icao_heights_kernel_code(nlayers, &
        icao_height, pressure_field, g_over_r,          &
        ndf, undf, map)

    use constants_mod,                  only: rmdi
    use science_aviation_constants_mod, only: mtokft,             &
        isa_lapse_ratel, isa_lapse_rateu, isa_press_bot,          &
        isa_press_mid, isa_press_top, isa_temp_bot, isa_temp_top, &
        gpm1, gpm2

    implicit none

    ! Arguments (kernel).
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf
    integer(kind=i_def), intent(in) :: undf
    integer(kind=i_def), intent(in), dimension(ndf) :: map

    ! Arguments (algorithm)
    real(kind=r_def), intent(out), dimension(undf) :: icao_height
    real(kind=r_def), intent(in), dimension(undf) :: pressure_field
    real(kind=r_def), intent(in) :: g_over_r


    ! Local variables
    real(kind=r_def) :: zp1
    real(kind=r_def) :: zp2
    real(kind=r_def) :: pressure


    zp1 = isa_lapse_ratel / g_over_r
    zp2 = isa_lapse_rateu / g_over_r

    pressure = pressure_field(map(1))

    ! Setting a safeguard limit to the lowest pressure to prevent
    ! extremely large icao height values near the top at the
    ! atmosphere.
    if ( (pressure >= 0.0_r_def) .and. (pressure <= 1000.0_r_def) ) then
        pressure = 1000.0_r_def
    end if

    ! Pressure must not be greater than surface pressure.
    pressure = min(isa_press_bot, pressure)

    ! Missing or invalid data?
    if (pressure == rmdi) then
        icao_height(map(1)) = rmdi

    ! Heights up to 11,000 gpm
    else if ( pressure > isa_press_mid ) then
        pressure = pressure / isa_press_bot
        pressure = 1.0_r_def - pressure**zp1
        icao_height(map(1)) = pressure * isa_temp_bot / isa_lapse_ratel

    ! Heights between 11,000 gpm and 20,000 gpm
    else if ( pressure > isa_press_top ) then
        pressure = pressure / isa_press_mid
        pressure = -log(pressure)
        icao_height(map(1)) = gpm1 + pressure * isa_temp_top / g_over_r

    ! Heights above 20,000 gpm
    else
        pressure = pressure / isa_press_top
        pressure = 1.0_r_def - pressure**zp2
        icao_height(map(1)) = gpm2 + pressure * isa_temp_top / isa_lapse_rateu
    end if

    ! Convert to kft
    if (icao_height(map(1)) /= rmdi) then
        icao_height(map(1)) = icao_height(map(1)) * mtokft
    end if

  end subroutine icao_heights_kernel_code

end module icao_heights_kernel_mod
