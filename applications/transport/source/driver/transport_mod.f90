!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Transport knows what namelists it needs.
!>
module transport_mod

  implicit none

  private

  character(*), public, parameter ::                           &
     transport_required_namelists(8) = ['base_mesh          ', &
                                        'planet             ', &
                                        'extrusion          ', &
                                        'initial_temperature', &
                                        'initial_wind       ', &
                                        'initial_density    ', &
                                        'transport          ', &
                                        'timestepping       ']

end module transport_mod
