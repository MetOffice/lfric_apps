!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> NAME Transport knows what namelists it needs.
!>
module name_transport_mod

  implicit none

  private

  character(*), public, parameter ::                                &
     name_transport_required_namelists(9) = ['base_mesh          ', &
                                             'planet             ', &
                                             'extrusion          ', &
                                             'initial_name_field ', &
                                             'initial_wind       ', &
                                             'initial_density    ', &
                                             'name_options       ', &
                                             'transport          ', &
                                             'timestepping       ']

end module name_transport_mod
