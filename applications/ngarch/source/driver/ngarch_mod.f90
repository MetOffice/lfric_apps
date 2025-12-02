!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> ngarch knows what configuration it needs.
!>
module ngarch_mod

  implicit none

  private

  character(*), public, parameter ::                      &
      ngarch_required_namelists(5) =  [ 'base_mesh     ', &
                                        'extrusion     ', &
                                        'finite_element', &
                                        'partitioning  ', &
                                        'planet        ']

end module ngarch_mod
