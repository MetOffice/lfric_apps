!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Outputs the diagnostics from the gravity-wave miniapp

!> @details Calls the routine that generates diagnostic output for all
!>          fields used by the gravity-wave miniapp. This is only a temporary
!>          hard-coded solution in lieu of a proper dianostic system

module gravity_wave_diagnostics_driver_mod

  use constants_mod,        only : i_def
  use field_mod,            only : field_type
  use mesh_mod,             only : mesh_type
  use model_clock_mod,      only : model_clock_type
  use diagnostics_io_mod,   only : write_scalar_diagnostic, &
                                   write_vector_diagnostic
  implicit none

  private
  public gravity_wave_diagnostics_driver

contains

  !> @brief Outputs the diagnostics from the gravity-wave miniapp
  !> @param [in]    mesh        Primary mesh
  !> @param [inout] state       A collection containing the fields that will
  !>                            be written to diagnostic output
  !> @param [in]    model_clock Time within the model.
  !> @param [in]    W3_project  Flag that determines if vector fields should be
  !>                            projected to W3
  !>
  subroutine gravity_wave_diagnostics_driver( mesh,        &
                                              wind,        &
                                              pressure,    &
                                              buoyancy,    &
                                              model_clock, &
                                              W3_project )

    implicit none
    type( field_type ),      intent(inout)       :: wind
    type( field_type ),      intent(inout)       :: buoyancy
    type( field_type ),      intent(inout)       :: pressure
    type( mesh_type  ),      intent(in), pointer :: mesh
    class(model_clock_type), intent(in)          :: model_clock
    logical,                 intent(in)          :: W3_project

    ! Calculation and output of diagnostics
    call write_vector_diagnostic( 'wind', wind, &
                                  model_clock, mesh, W3_project )
    call write_scalar_diagnostic( 'pressure', pressure, &
                                  model_clock, mesh, W3_project )
    call write_scalar_diagnostic( 'buoyancy', buoyancy, &
                                  model_clock, mesh, W3_project )

  end subroutine gravity_wave_diagnostics_driver

end module gravity_wave_diagnostics_driver_mod
