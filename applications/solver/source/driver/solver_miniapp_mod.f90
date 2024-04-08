!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> Miniapp skeleton program support functions.
!>
!> Originally these were "block" constructs within the program but neither
!> GNU or Intel Fortran where properly able to cope with that.
!>
module solver_miniapp_mod

  implicit none

  private
  public :: solver_required_namelists

  character(*), parameter :: &
    solver_required_namelists(9) = ['finite_element      ', &
                                    'base_mesh           ', &
                                    'formulation         ', &
                                    'planet              ', &
                                    'extrusion           ', &
                                    'solver              ', &
                                    'partitioning        ', &
                                    'extrusion           ', &
                                    'logging             ']

end module solver_miniapp_mod
