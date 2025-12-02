! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Evaporation routine
!
! Subroutine Interface:
!
module devap_mod

use um_types, only: real_umphys

implicit none

! Description: Evaporation routine
!              Carries out evaporation and updates precipitation
!
! Method: UM documentation paper 27
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards vn8.2.

character(len=*), parameter, private :: ModuleName = 'DEVAP_MOD'

contains

subroutine devap(npnts, bddwt_km1                                              &
                 ,thdd_k, thdds, qdds, flx_dd_km1, exk, exkm1                  &
                 ,qsatdd, delpkm1, cca, pkm1                                   &
                 ,thdd_km1, qdd_km1, rain, snow)

use planet_constants_mod, only: g, cp, r

use cv_param_mod, only: dd_area_fac

use water_constants_mod, only: lc, lf
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

use evp_mod, only: evp
implicit none

!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------
integer, intent(in) ::                                                         &
  npnts                  ! Vector length

logical, intent(in) ::                                                         &
  bddwt_km1(npnts)       ! Mask for those points in downdraught where
                         ! precipitation is liquid

real(kind=real_umphys), intent(in) ::                                          &
  thdd_k(npnts)        & ! Potential temperature of downdraught in layer k (K)

 ,thdds(npnts)         & ! Potential temperature of saturated downdraught
                         ! in layer (K)
 ,qdds(npnts)          & ! Saturated downdraught mixing ratio of layer
                         ! (kg/kg)
 ,flx_dd_km1(npnts)    & ! Downdraught mass flux in layer k-1 (Pa/s)

 ,exk(npnts)           & ! Exner ratio in layer k
 ,exkm1(npnts)         & ! Exner ratio in layer k-1
 ,qsatdd(npnts)        & ! Saturated downdraught mixing ratio (kg/kg)
 ,delpkm1(npnts)       & ! Change in pressure across layer k-1  (Pa)
 ,cca(npnts)           & ! Convective cloud amount (fraction)
 ,pkm1(npnts)            ! Pressure in layer k-1  (Pa)

real(kind=real_umphys), intent(in out) ::                                      &
  thdd_km1(npnts)      & ! In  Potential temperature of DD in layer k-1(kg/kg)
                         ! Out Potential temperature of DD in layer k-1 after
                         ! evaporation of saturation (K)
 ,qdd_km1(npnts)       & ! In  Mixing ratio of downdraught in layer k-1 (kg/kg)
                         ! Out Mixing ratio of DD in layer k-1 after
                         ! evaporation of saturation (kg/kg)
 ,rain(npnts)          & ! In  Amount of rain
                         ! Out Updated amount of rainfall (kg/m**2/s)
 ,snow(npnts)            ! In  Amount of snow
                         ! Out Updated amount of snowfall (kg/m**2/s)

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='DEVAP'

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

integer :: i               ! Loop counter

logical ::                                                                     &
  bevap(npnts)          & ! Mask for those points at which evaporation
                          ! calculation is to be carried out
 ,bsat(npnts)           & ! Mask for those points which are subsaturated
 ,l_rates_adjusted      & ! True if evaporation / sublimation rates have
                          ! been adjusted to give saturation
 ,full_evap_rain(npnts) & ! True if full rain evap in evp
 ,full_evap_snow(npnts)   ! True if full snow evap in evp

real(kind=real_umphys) ::                                                      &
  tevp(npnts)      & ! Temperature used in evaporation calculation (K)

 ,evap_rain(npnts) & ! Amount of evaporation of rain

 ,sub_snow(npnts)  & ! Amount of snow sublimation

 ,delq(npnts)      & ! Difference in mixing ratio (kg/kg)

 ,delth(npnts)     & ! Increment to downdraught potential temperature in layer
                     ! k-1  due to evaporation (K)
 ,delqe            & ! Increment to downdraught mixing ratio in layer k-1
                     !  due to evaporation (kg/kg)
 ,delths           & ! Saturated potential temperature minus potential
                     ! temperature of downdraught
 ,factor           & ! delths/delth

 ,pincr            & ! Increase in precipitation if parcel supersaturates

 ,rho(npnts)         ! Density of air in parcel

!-----------------------------------------------------------------------
! Check if evaporation possible
!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do i=1,npnts
  delq(i) = qsatdd(i)-qdd_km1(i)

  bevap(i) =((rain(i) >  0.0) .or. (snow(i) >  0.0)) .and. (delq(i) >  0.0)
  bsat(i) = delq(i)  <   0.0

  !-----------------------------------------------------------------------
  ! Calculate temperature used in calculation of evaporation constants
  ! based on temeprature of parcel after unsaturated descent
  !-----------------------------------------------------------------------

  if (bevap(i)) then
    tevp(i) = ((thdd_k(i)*exk(i))+(thdd_km1(i)*exkm1(i)))*0.5
    rho(i) = pkm1(i) / (r*tevp(i))
  end if
end do

!-----------------------------------------------------------------------
! Evaporation calculation - calculate rates for rain and snow
!-----------------------------------------------------------------------

call evp(npnts,1,bevap,dd_area_fac,rain,tevp,cca,rho,delq,delpkm1,             &
         pkm1,evap_rain,full_evap_rain)

call evp(npnts,2,bevap,dd_area_fac,snow,tevp,cca,rho,delq,delpkm1,             &
         pkm1,sub_snow,full_evap_snow)

do i=1,npnts
  if (bevap(i)) then

    !-----------------------------------------------------------------------
    ! Adjust evaporation and sublimation rates back to grid box means
    !-----------------------------------------------------------------------

    evap_rain(i) = evap_rain(i) * cca(i) * dd_area_fac
    sub_snow(i) = sub_snow(i) * cca(i) * dd_area_fac

    !-----------------------------------------------------------------------
    ! Check if parcel supersaturated
    !-----------------------------------------------------------------------

    delth(i) = -((lc*evap_rain(i))+((lc+lf)*sub_snow(i)))*g/                   &
                                        (cp*exkm1(i)*flx_dd_km1(i))
    delqe = (evap_rain(i)+sub_snow(i))*g/flx_dd_km1(i)

    delths = thdds(i)-thdd_km1(i)

    l_rates_adjusted = .false.

    if (delth(i) <  delths) then

      !-----------------------------------------------------------------------
      ! Adjust evaporation and sublimation rates to give saturation
      !-----------------------------------------------------------------------

      l_rates_adjusted = .true.
      factor = delths/delth(i)
      delth(i)  = delths
      delqe  = delqe*factor
      evap_rain(i) = evap_rain(i)*factor
      sub_snow(i)  = sub_snow(i)*factor
    end if

    !-----------------------------------------------------------------------
    ! Update T, q and precipitation
    !-----------------------------------------------------------------------

    rain(i) = rain(i)-evap_rain(i)
    if (rain(i) <  0.0) rain(i)=0.0
    snow(i) = snow(i)-sub_snow(i)
    if (snow(i) <  0.0) snow(i)=0.0
    thdd_km1(i) = thdd_km1(i)+ delth(i)
    qdd_km1(i)  = qdd_km1(i) + delqe

    ! If the call to evp has produced an evaporation rate that evaporates
    ! all the precipitation, and there was no adjustment of the rate to
    ! account for saturation, logic dictates that the precipitation rate
    ! should now be exactly zero. However, because of rounding effects it is
    ! possible to end up with slightly non-zero rates. To get around this,
    ! we directly reset to zero:
    if (.not. l_rates_adjusted) then
      if (full_evap_rain(i)) rain(i) = 0.0
      if (full_evap_snow(i)) snow(i) = 0.0
    end if

    !-----------------------------------------------------------------------
    ! Parcel is supersaturated before evaporation occurs
    ! Bring parcel to saturation and precipitate water
    !-----------------------------------------------------------------------

  else if (bsat(i)) then
    pincr       = (qdd_km1(i)-qdds(i))*flx_dd_km1(i)/g
    qdd_km1(i)  = qdds(i)
    thdd_km1(i) = thdds(i)
    if (bddwt_km1(i)) then
      rain(i) = rain(i)+pincr
    else
      snow(i) = snow(i)+pincr
    end if
  end if
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine devap

end module devap_mod
