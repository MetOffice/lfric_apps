!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Module for printing of IAU related fields.

module printing_iau_mod
  use constants_mod,           only : i_def
  use field_collection_mod,    only : field_collection_type
  use field_mod,               only : field_type
  use log_mod,                 only : log_level
  use print_meanrms_field_mod, only : print_meanrms_field
  use field_minmax_alg_mod,    only : log_field_minmax

  implicit none

  private
  public :: print_minmax_prog,          &
            print_meanrms_prog,         &
            print_minmax_cld,           &
            print_meanrms_cld,          &
            print_minmax_surf

  contains

  !> @brief   print the prognostic_fields min-max
  !> @details print the min-max of the different iau fields
  !> @param[in] prognostic_fields The collection of prognostic fields
  !> @param[in] level             Level of logging. If the configured
  !>                              log_level is less than or equal to
  !>                              level, output will be shown.
  subroutine print_minmax_prog( prognostic_fields, level )

    implicit none

    type( field_collection_type ), intent(in) :: prognostic_fields
    integer( kind = i_def ),       intent(in) :: level

    type( field_type ),            pointer    :: prog_field  => null()

    if( log_level() <= level ) then

      call prognostic_fields % get_field( 'theta', prog_field )
      call log_field_minmax( level, 'theta', prog_field )

      call prognostic_fields % get_field( 'exner', prog_field )
      call log_field_minmax( level, 'exner', prog_field )

      call prognostic_fields % get_field( 'u', prog_field )
      call log_field_minmax( level, 'wind', prog_field )

      call prognostic_fields % get_field( 'm_v', prog_field )
      call log_field_minmax( level, 'm_v', prog_field )

      call prognostic_fields % get_field( 'm_cl', prog_field )
      call log_field_minmax( level, 'm_cl', prog_field )

      call prognostic_fields % get_field( 'm_ci', prog_field )
      call log_field_minmax( level, 'm_ci', prog_field )

      call prognostic_fields % get_field( 'm_r', prog_field )
      call log_field_minmax( level, 'm_r', prog_field )

      call prognostic_fields % get_field( 'm_s', prog_field )
      call log_field_minmax( level, 'm_s', prog_field )

      call prognostic_fields % get_field( 'm_g', prog_field )
      call log_field_minmax( level, 'm_g', prog_field )

      call prognostic_fields % get_field( 'rho', prog_field )
      call log_field_minmax( level, 'dry_rho', prog_field )

      nullify(prog_field)

    end if

  end subroutine print_minmax_prog

  !> @brief   print the prognostic_fields mean/rms
  !> @details print the mean/rms of the different iau fields
  !> @param[in] prognostic_fields The collection of prognostic fields
  !> @param[in] level             Level of logging. If the configured
  !>                              log_level is less than or equal to
  !>                              level, output will be shown.
  subroutine print_meanrms_prog( prognostic_fields, level )

    implicit none

    type( field_collection_type ), intent(in) :: prognostic_fields
    integer( kind = i_def ),       intent(in) :: level

    type( field_type ),            pointer :: prog_field  => null()

    if( log_level() <= level ) then

      call prognostic_fields % get_field( 'theta', prog_field )
      call print_meanrms_field( prog_field, level )

      call prognostic_fields % get_field( 'exner', prog_field )
      call print_meanrms_field( prog_field, level )

      call prognostic_fields % get_field( 'u', prog_field )
      call print_meanrms_field( prog_field, level )

      call prognostic_fields % get_field( 'm_v', prog_field )
      call print_meanrms_field( prog_field, level )

      call prognostic_fields % get_field( 'm_cl', prog_field )
      call print_meanrms_field( prog_field, level )

      call prognostic_fields % get_field( 'm_ci', prog_field )
      call print_meanrms_field( prog_field, level )

      call prognostic_fields % get_field( 'm_r', prog_field )
      call print_meanrms_field( prog_field, level )

      call prognostic_fields % get_field( 'm_s', prog_field )
      call print_meanrms_field( prog_field, level )

      call prognostic_fields % get_field( 'm_g', prog_field )
      call print_meanrms_field( prog_field, level )

      call prognostic_fields % get_field( 'rho', prog_field )
      call print_meanrms_field( prog_field, level )

      nullify( prog_field )

    end if

  end subroutine print_meanrms_prog

  !> @brief Print the cloud_fields min-max.
  !> @param[in] cloud_fields  Collection of cloud fields
  !> @param[in] level         Level of logging. If the configured
  !>                          log_level is less than or equal to
  !>                          level, output will be shown.
  subroutine print_minmax_cld( cloud_fields, level )

    implicit none

    type( field_collection_type ), intent(in) :: cloud_fields
    integer( kind = i_def ),       intent(in) :: level

    type( field_type ),            pointer :: cld_field  => null()

    if( log_level() <= level ) then

      call cloud_fields % get_field( 'area_fraction', cld_field )
      call log_field_minmax( level, 'acf', cld_field )

      call cloud_fields % get_field( 'bulk_fraction', cld_field )
      call log_field_minmax( level, 'bcf', cld_field )

      call cloud_fields % get_field( 'frozen_fraction', cld_field )
      call log_field_minmax( level, 'cff', cld_field )

      call cloud_fields % get_field( 'liquid_fraction', cld_field )
      call log_field_minmax( level, 'cfl', cld_field )

      nullify( cld_field )

    end if

  end subroutine print_minmax_cld

  !> @brief Print the cloud_fields mean and rms.
  !> @param[in] cloud_fields  Collection of cloud fields
  !> @param[in] level         Level of logging. If the configured
  !>                          log_level is less than or equal to
  !>                          level, output will be shown.
  subroutine print_meanrms_cld( cloud_fields, level )

    implicit none

    type( field_collection_type ), intent(in) :: cloud_fields
    integer( kind = i_def ),       intent(in) :: level

    type( field_type ),            pointer :: cld_field  => null()

    if( log_level() <= level ) then

      call cloud_fields % get_field( 'area_fraction', cld_field )
      call print_meanrms_field( cld_field, level )

      call cloud_fields % get_field( 'bulk_fraction', cld_field )
      call print_meanrms_field( cld_field, level )

      call cloud_fields % get_field( 'frozen_fraction', cld_field )
      call print_meanrms_field( cld_field, level )

      call cloud_fields % get_field( 'liquid_fraction', cld_field )
      call print_meanrms_field( cld_field, level )

      nullify( cld_field )

    end if

  end subroutine print_meanrms_cld

  !> @brief Print the land surface fields min-max
  !> @param[in] soil_fields       Collection of soil fields
  !> @param[in] surface_fields    Collection of surface fields
  !> @param[in] snow_fields       Collection of snow fields
  !> @param[in] level             Level of logging. If the configured
  !>                              log_level is less than or equal to
  !>                              level, output will be shown.
  subroutine print_minmax_surf( soil_fields, surface_fields, snow_fields, level )

    implicit none

    type( field_collection_type ), intent(in) :: soil_fields
    type( field_collection_type ), intent(in) :: surface_fields
    type( field_collection_type ), intent(in) :: snow_fields
    integer( kind = i_def ),       intent(in) :: level

    type( field_type ),            pointer :: surf_field  => null()

    if( log_level() <= level ) then

      call soil_fields % get_field( 'soil_moisture', surf_field )
      call log_field_minmax( level, 'soil_moisture', surf_field )

      call soil_fields % get_field( 'soil_temperature', surf_field )
      call log_field_minmax( level, 'soil_temperature', surf_field )

      call surface_fields % get_field( 'tile_temperature', surf_field )
      call log_field_minmax( level, 'tile_temperature', surf_field )

      call snow_fields % get_field( 'snow_layer_temp', surf_field )
      call log_field_minmax( level, 'snow_layer_temp', surf_field )

      nullify( surf_field )

    end if

  end subroutine print_minmax_surf

end module printing_iau_mod
