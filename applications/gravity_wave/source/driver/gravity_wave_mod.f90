!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> Gravity wave parameters.
!>
module gravity_wave_mod

  implicit none

  private

  character(*), public, parameter :: &
    gravity_wave_required_namelists(8) = ['base_mesh             ', &
                                          'planet                ', &
                                          'extrusion             ', &
                                          'initial_temperature   ', &
                                          'initial_wind          ', &
                                          'timestepping          ', &
                                          'formulation           ', &
                                          'gravity_wave_constants']

end module gravity_wave_mod
