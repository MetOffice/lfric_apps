! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Timestepping values.

module timestep_mod

! Description:
!              Timestep information
!              Updated at start of ATM_STEP every timestep

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: top_level

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

use missing_data_mod, only: rmdi

implicit none

integer, target :: timestep_number

real       :: timestep      = rmdi ! atmosphere model timestep
real       :: recip_timestep= rmdi ! recip model timestep

! These need the save attribute, because they are set from a subroutine
! call outside of atm_step; timestep_mod then goes out of scope before
! timestepping begins.
real, save :: radiation_tstep_diag ! timestep of fast radiation (3C)
real, save :: radiation_tstep_prog ! timestep of slow radiation (3C)

real       :: chemistry_timestep   ! must be  <=  model timestep
real       :: pos_timestep         ! = +timestep.
real       :: neg_timestep         ! = -timestep.

end module timestep_mod
