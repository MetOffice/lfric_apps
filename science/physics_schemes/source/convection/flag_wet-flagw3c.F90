! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate a mask for when condensation is liquid
!
module flag_wet_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'FLAG_WET_MOD'
contains

subroutine flag_wet (np_field, npnts, nlev,                                    &
                     th, exner_layer_centres,                                  &
                     bwater )

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

use cv_run_mod, only:                                                          &
    tice

implicit none

!
! Description : Calculates a mask for when condensation is liquid
!
! Method:
!   If 0.5 * (TK + TK+1) > TICE  then any condensation  in layer k+1 is liquid
!
!   If 0.5 * (TK + TK+1) < TICE  then any condensation  in layer k+1 is ice
!
!  Returns bwater - logical true if liquid rather than ice
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP programming standards vn8.2.
!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------

integer, intent(in) ::                                                         &
  np_field             &  ! Full vector length
 ,npnts                &  ! Vector length
 ,nlev                    ! Number of model layers

real(kind=real_umphys), intent(in) ::                                          &
  th(np_field,nlev)                    & ! Potential temperature (K)
 ,exner_layer_centres(np_field,0:nlev)   ! Exner ratio at layer centres
                                         ! (starting with the surface).
logical, intent(out) ::                                                        &
  bwater(npnts,2:nlev)     ! mask for those points at which condensate
                           ! is liquid

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

integer ::                                                                     &
  i,k               ! Loop counters

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='FLAG_WET'

!----------------------------------------------------------------------
!  Calculate mask
!----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do k=1,nlev-1
  do i=1,npnts

    bwater(i,k+1) = 0.5*(th(i,k)*exner_layer_centres(i,k) +                    &
                         th(i,k+1)*exner_layer_centres(i,k+1))  >  tice

  end do  ! npnts
end do    !  nlev -1

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine flag_wet
end module flag_wet_mod
