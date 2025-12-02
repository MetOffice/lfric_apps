! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Control variables for idealised dry planet suite

! Description:
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised

! Method:
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3

module planet_suite_mod

use tforce_mod, only: tf_none
use trelax_mod, only: tr_none

implicit none

! Items read in through IDEALISED namelist
integer :: tforce_number = tf_none    ! Choice of forcing profile
integer :: trelax_number = tr_none    ! Choice of relaxation timescale
integer :: nsteps_consv_print = 0     ! Frequency of printing of AAM and KE

end module planet_suite_mod
