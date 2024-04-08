!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> jedi_lfric_tests parameters.
!>
module jedi_lfric_tests_mod

  implicit none

  private

  character(*), public, parameter ::                        &
    jedi_lfric_tests_required_namelists(5) =  [ 'base_mesh     ', &
                                                'extrusion     ', &
                                                'finite_element', &
                                                'partitioning  ', &
                                                'planet        ']

end module jedi_lfric_tests_mod
