!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Shallow water knows what namelists it needs.
!>
module shallow_water_mod

  implicit none

  private

  character(*), public, parameter :: &
    shallow_water_required_namelists(8) &
      = ['base_mesh             ', &
         'files                 ', &
         'formulation           ', &
         'planet                ', &
         'shallow_water_settings', &
         'time                  ', &
         'timestepping          ', &
         'transport             ']

end module shallow_water_mod
