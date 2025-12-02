! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Module containing variables which are common to the tcs warm routines
!
module tcs_common_warm

use tcs_class_scales, only:                                                    &
   scales_conv

use tcs_class_interface, only:                                                 &
   interface_input

implicit none
!
! Description:
!   This routine holds a common set of variables for the tcs
! warm routines
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.1 programming standards.
!

type(scales_conv)     :: scales     !
type(interface_input) :: cb     ! Cloud base
type(interface_input) :: cb_p1  ! Level above cloud base
type(interface_input) :: cb_m1  ! Level below cloud base
type(interface_input) :: inv    ! Inversion base
type(interface_input) :: inv_p1 ! Level above inv
type(interface_input) :: inv_m1 ! Level below inv

end module tcs_common_warm
