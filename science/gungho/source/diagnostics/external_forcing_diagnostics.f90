!-------------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!
!> @brief Module for writing external forcing increment diagnostics
!> @details Set of subroutines for writing external forcing increment fields
!
module external_forcing_diagnostics_mod

  use field_mod,                  only: field_type
  use log_mod,                    only: log_event, LOG_LEVEL_INFO
  use constants_mod,              only: l_def
  use timing_mod,                 only: start_timing, stop_timing, tik, LPROF
  use initialise_diagnostics_mod, only: init_diag => init_diagnostic_field
  use physics_mappings_alg_mod,   only: map_physics_winds

  implicit none

  private
  public :: write_forcing_diagnostics

contains

  !> @brief Write external forcing increments when requested
  !> @details If any external forcing increments was requested for output then write it
  !> @param[in]  du_forcing       3D wind increment from external forcing
  !> @param[in]  output_wind_inc  Logical flag to output wind increments
  subroutine write_forcing_diagnostics(du_forcing, output_wind_inc)

    implicit none

    type( field_type ),    intent( in ) :: du_forcing
    logical( kind=l_def ), intent( in ) :: output_wind_inc


    logical( kind=l_def ) :: output_du_force, &
                             output_dv_force, &
                             output_dw_force
    type( field_type ) :: du_force, dv_force, dw_force
    integer(tik)       :: id

    if ( LPROF ) call start_timing( id, 'write_forcing_diagnostics' )

    if ( output_wind_inc ) then
      !
      ! Set wind forcing increment output check flags. NOTE: By default the init_diag function
      ! calls below will also initialise the wind increment fields that are requested.
      !
      output_du_force = init_diag(du_force, 'forcing__du_force')
      output_dv_force = init_diag(dv_force, 'forcing__dv_force')
      output_dw_force = init_diag(dw_force, 'forcing__dw_force')

      !
      ! Output requested wind forcing increments.
      !
      if ( output_du_force .or. output_dv_force .or. output_dw_force) then

        call log_event( 'slow_physics: Output wind forcing increments', &
                        LOG_LEVEL_INFO )

        !
        ! Initialise any still uninitialised wind forcing increment fields. NOTE: The
        ! easterly wind increment is logically treated differently here, as it is
        ! initialised first, and is then subsequently used to initialise the other two
        ! wind increments.
        !
        if ( .not. output_du_force ) then
          if ( output_dv_force ) then
            call du_force%initialise( vector_space = dv_force%get_function_space() )
          else
            call du_force%initialise( vector_space = dw_force%get_function_space() )
          end if
        end if
        if ( .not. output_dv_force ) call dv_force%initialise( vector_space = du_force%get_function_space() )
        if ( .not. output_dw_force ) call dw_force%initialise( vector_space = du_force%get_function_space() )

        !
        ! Remap from 3D wind increment to individual wind forcing increments
        !
        call map_physics_winds(du_force, dv_force, dw_force, du_forcing)

        !
        ! Only write requested wind forcing increments
        !
        if ( output_du_force ) call du_force%write_field('forcing__du_force')
        if ( output_dv_force ) call dv_force%write_field('forcing__dv_force')
        if ( output_dw_force ) call dw_force%write_field('forcing__dw_force')

      end if

    end if

    if ( LPROF ) call stop_timing( id, 'write_forcing_diagnostics' )

  end subroutine write_forcing_diagnostics

end module external_forcing_diagnostics_mod
