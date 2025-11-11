!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Lists the namelists that lfric2lfric needs to configure
!!          the application.
!>
module lfric2lfric_mod

  implicit none

  private

  character(len=*), public, parameter ::                       &
      lfric2lfric_required_namelists(6) =  [ 'extrusion     ', &
                                             'finite_element', &
                                             'partitioning  ', &
                                             'planet        ', &
                                             'lfric2lfric   ', &
                                             'files         ']

end module lfric2lfric_mod
