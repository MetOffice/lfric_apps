!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Gungho's requirements for configuration.
!>
!> The job of setting up configuration belongs to the host program but Gungho
!> has requirements for that configuration.
!>
module gungho_mod

  implicit none

  private

  !> Lists namelists which must be present.
  !>
  !> @todo A mechanism whereby this list is dynamic and dependent on values
  !>       within the configuration. i.e. If a certain feature is switched on
  !>       then require the list which configures that feature.
  !>
  character(*), public, parameter ::                       &
      gungho_required_namelists(12) = ['finite_element  ', &
                                       'formulation     ', &
                                       'base_mesh       ', &
                                       'initial_wind    ', &
                                       'planet          ', &
                                       'solver          ', &
                                       'mixed_solver    ', &
                                       'helmholtz_solver', &
                                       'timestepping    ', &
                                       'extrusion       ', &
                                       'transport       ', &
                                       'orography       ']

end module gungho_mod
