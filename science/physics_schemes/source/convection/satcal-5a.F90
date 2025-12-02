! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate saturated temperature
!
module satcal_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'SATCAL_MOD'
contains

subroutine satcal (npnts, th, pk, exk, q_k, the_k, qse_k, qs, thdds)

use planet_constants_mod, only: cp, rv

use water_constants_mod, only: lc,  tm

use cv_derived_constants_mod, only: ls

use qsat_mod, only: qsat

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

!
! Description: Calculate saturated temperature
!
! Method: UM documentation paper 27
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards vn8.2.
!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------
integer, intent(in) ::                                                         &
  npnts                   ! Vector length

real(kind=real_umphys), intent(in) ::                                          &
  th(npnts)          & ! Potential temperature (K)

 ,pk(npnts)          & ! Pressure of layer k (Pa)

 ,exk(npnts)         & ! Exner ratio of layer k

 ,q_k(npnts)         & ! Mixing ratio of layer k (kg/kg)

 ,the_k(npnts)       & ! Environmental potential temperature in layer k
 ,qse_k(npnts)         ! qsat for environment of layer k (kg/kg)

real(kind=real_umphys), intent(out) ::                                         &
  qs(npnts)          & ! Saturated specific  humidity (kg/kg)

 ,thdds(npnts)         ! Saturated environmental potential temperature (K)

!-----------------------------------------------------------------------
! Model constants
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

integer ::                                                                     &
  i, ic            ! Loop counters

real(kind=real_umphys) ::                                                      &
  l                ! Latent heat

real(kind=real_umphys) ::                                                      &
  t_fg(npnts)    & ! Temperature first guess (K)

 ,th_fg(npnts)   & ! Potential temperature first guess (K)

 ,dqbydt           ! First guess at mixing ratio increment (kg/kg/s)

real(kind=real_umphys) ::                                                      &
  cpexk(npnts)     ! cp * exner

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='SATCAL'


!-----------------------------------------------------------------------
! Set initial first guess temperature and theta - based upon
! environmental temperature in layer k
!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do i=1,npnts
  th_fg(i) = the_k(i)
  t_fg(i)  = th_fg(i)*exk(i)
  qs(i)    = qse_k(i)        ! first guess qsat from environment value
  cpexk(i) = cp*exk(i)       ! save CPU by calculating before interation
end do

!----------------------------------------------------------------------
! Do two iterations to find saturation point due to evaporation
!----------------------------------------------------------------------

do ic=1,2

  !----------------------------------------------------------------------
  ! Calculate dqsat/dT for first guess temperature
  !----------------------------------------------------------------------

  !  call dqs_dth(dqbydt,th_fg,qs,exk,npnts)
  ! Cheaper to remove call and do inline making use of l

    !----------------------------------------------------------------------
    ! Calculate updated temperature at saturation
    !----------------------------------------------------------------------

  do i=1,npnts

    if (t_fg(i) >  tm) then
      l=lc
    else
      l=ls     ! lc+lf
    end if
    !----------------------------------------------------------------------
    ! Calculate dqsat/dT for first guess temperature
    !----------------------------------------------------------------------
    dqbydt = l*qs(i)/(exk(i)*rv*th_fg(i)*th_fg(i))

    thdds(i) = (th(i)*cpexk(i) - l*(qs(i)-q_k(i)-th_fg(i)*dqbydt)) /           &
                  (cpexk(i)+l*dqbydt)


    !----------------------------------------------------------------------
    ! Calculate temperature at saturation and update first guess
    !----------------------------------------------------------------------

    th_fg(i) = thdds(i)
    t_fg(i) = th_fg(i)*exk(i)

  end do

  !----------------------------------------------------------------------
  ! Calculate revised saturation mixing ratio at saturation
  !---------------------------------------------------------------------
  call qsat(qs,t_fg,pk,npnts)

end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine satcal
end module satcal_mod
