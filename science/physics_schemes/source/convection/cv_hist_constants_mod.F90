! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!  Convection history prognostics - constants

module cv_hist_constants_mod

use um_types, only: real_umphys

implicit none
save

! Description:
!   Module containing constants used in the calculation of the
!   convective history prognostics.
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.1.
!
!------------------------------------------------------------------------------

! Decay period for deep convection flag

real(kind=real_umphys), parameter :: decay_period = 10800.0 ! 3 hours (s)


!------------------------------------------------------------------------------

end module cv_hist_constants_mod
