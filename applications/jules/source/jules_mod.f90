!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> JULES app requirements for configuration.
!>
!> The JULES app runs the Gungho driver without dynamics, so it does not need
!> the solver or transport namelists required by Gungho. It does need the energy
!> correction namelist, as its absence would switch on energy correction.
!>
module jules_mod

  implicit none

  private

  !> Lists namelists which must be present.
  !>
  character(*), public, parameter ::                        &
      jules_required_namelists(9) = ['finite_element   ', &
                                      'formulation      ', &
                                      'base_mesh        ', &
                                      'initial_wind     ', &
                                      'planet           ', &
                                      'timestepping     ', &
                                      'extrusion        ', &
                                      'orography        ', &
                                      'energy_correction']

end module jules_mod
