!-------------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
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
MODULE diags_icao_heights_kernel_mod

  USE argument_mod,         ONLY: arg_type,            &
                                  GH_FIELD, GH_SCALAR, GH_LOGICAL, &
                                  GH_READ, GH_WRITE, GH_INTEGER, &
                                  GH_REAL, CELL_COLUMN, &
                                  ANY_DISCONTINUOUS_SPACE_1, &
                                  ANY_DISCONTINUOUS_SPACE_2
  USE constants_mod,        ONLY: r_def, i_def, l_def
  USE kernel_mod,           ONLY: kernel_type

  IMPLICIT NONE

  ! The aviation diagnostics kernel type.
  TYPE, EXTENDS(kernel_type) :: diags_icao_heights_kernel_type
    TYPE(arg_type), DIMENSION(2) :: meta_args = (/ &

      ! Output icao height.
      arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &

      ! Input pressure
      arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1) &

      /)

    INTEGER :: operates_on = cell_column

    CONTAINS
      PROCEDURE, NOPASS :: code => diags_icao_heights_kernel_code
    END TYPE diags_icao_heights_kernel_type

CONTAINS

  SUBROUTINE diags_icao_heights_kernel_code(nlayers, &
        icao_height, &  ! Output icao height.
        pressure_field, &     ! Pressure.
        ndf, undf, map)

    ! Calculate iaco height from the pressure field.
    ! Assumes lowest order W3 data, where ndf is always 1.

    USE constants_mod,                  ONLY: RMDI
    USE planet_config_mod,              ONLY: gravity, rd
    USE science_conversions_mod,        ONLY: feet_to_metres
    USE science_aviation_constants_mod, ONLY: mtokft, &
        isa_Lapse_RateL, isa_Lapse_RateU, isa_Press_Bot, &
        isa_Press_Mid, isa_Press_Top, isa_Temp_Bot, isa_Temp_Top, &
        Gpm1, Gpm2

    IMPLICIT NONE

    ! Arguments (kernel)

    ! The number of layers in a column.
    INTEGER(KIND=i_def), INTENT(IN) :: nlayers

    ! Number of degrees of freedom (columns) in the cell we're processing.
    INTEGER(KIND=i_def), INTENT(IN) :: ndf

    ! Number of unique degrees of freedom in the fields.
    INTEGER(KIND=i_def), INTENT(IN) :: undf

    ! Degrees of freedom maps. Offsets to the bottom of each column.
    INTEGER(KIND=i_def), INTENT(IN), DIMENSION(ndf) :: map


    ! Arguments (algorithm)

    ! Output icao height.
    REAL(KIND=r_def), INTENT(INOUT), DIMENSION(undf) :: icao_height

    ! Pressure in Pa.
    REAL(KIND=r_def), INTENT(IN), DIMENSION(undf) :: pressure_field


    ! Local variables
    INTEGER(KIND=i_def) :: df

    REAL(KIND=r_def) :: g_over_rd
    REAL(KIND=r_def) :: zp1
    REAL(KIND=r_def) :: zp2

    REAL(KIND=r_def) :: pressure


    g_over_rd = gravity / rd
    zp1 = isa_Lapse_RateL / g_over_rd
    zp2 = isa_Lapse_RateU / g_over_rd

    pressure = pressure_field(map(1))

    ! Setting a safeguard limit to the lowest pressure to prevent
    ! extremely large icao height values near the top at the
    ! atmosphere.
    if ( (pressure >= 0.0_r_def) .and. (pressure <= 1000.0_r_def) ) then
        pressure = 1000.0_r_def
    end if

    ! Pressure must not be greater than surface pressure.
    pressure = MIN(isa_Press_bot, pressure)

    ! Missing or invalid data?
    if (pressure == RMDI .or. pressure <= 0.0_r_def) then
        icao_height(map(1)) = RMDI

    ! Heights up to 11,000 GPM
    else if ( pressure > isa_Press_Mid ) then
        pressure = pressure / isa_Press_Bot
        pressure = 1.0_r_def - pressure**zp1
        icao_height(map(1)) = pressure * isa_Temp_Bot / isa_Lapse_RateL

    ! Heights between 11,000 GPM and 20,000 GPM
    else if ( pressure > isa_Press_Top ) then
        pressure = pressure / isa_Press_Mid
        pressure = -LOG(pressure)
        icao_height(map(1)) = Gpm1 + pressure * isa_Temp_Top / g_over_rd

    ! Heights above 20,000 GPM
    else
        pressure = pressure / isa_Press_Top
        pressure = 1.0_r_def - pressure**zp2
        icao_height(map(1)) = Gpm2 + pressure * isa_Temp_Top / isa_Lapse_RateU
    end if

  END SUBROUTINE diags_icao_heights_kernel_code

END MODULE diags_icao_heights_kernel_mod
