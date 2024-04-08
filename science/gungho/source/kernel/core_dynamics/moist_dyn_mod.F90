!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!>  @brief   Indexing for moisture-dependent factors in dynamics
!!
!-------------------------------------------------------------------------------
module moist_dyn_mod

implicit none

! Indexing for moist_dyn field bundle
integer, parameter :: num_moist_factors = 3
integer, parameter :: gas_law           = 1 ! 1 + m_v / epsilon
integer, parameter :: total_mass        = 2 ! 1 + sum (m_x)
integer, parameter :: water             = 3 ! For future development

character(10), parameter :: moist_dyn_names(num_moist_factors) = &
  [character(10) :: 'gas_law', 'total_mass', 'water']

end module moist_dyn_mod
