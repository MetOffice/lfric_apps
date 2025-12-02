! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Evaporate rain below cloud base if no downdraught
!
module pevp_bcb_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'PEVP_BCB_MOD'
contains

subroutine pevp_bcb (npnts,k,iccb,th,pk,q,qse,delp,rain,snow,                  &
                     dthbydt,dqbydt,exk,timestep,cca)


use cv_param_mod, only: precip_area_fac

use planet_constants_mod, only: g, r, cp

use water_constants_mod, only: lc, lf

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

use evp_mod, only: evp
use satcal_mod, only: satcal
implicit none

!-----------------------------------------------------------------------
! Description:  Evaporate rain below cloud base if no downdraught.
!
! Method: UM documentation paper 27
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards vn8.3.
!-----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------

integer, intent(in) ::                                                         &
  npnts                & ! Vector length
 ,k                      ! present model layer

integer, intent(in) ::                                                         &
  iccb(npnts)            ! convective cloud base layer

real(kind=real_umphys), intent(in) ::                                          &
  th(npnts)            & ! Potential temperature (K)
 ,pk(npnts)            & ! Pressure (Pa)
 ,q(npnts)             & ! Mixing ratio (kg/kg)
 ,qse(npnts)           & ! Environmental qsat Mixing ratio (kg/kg)
 ,delp(npnts)          & ! Change in pressure across layer k-1 (Pa)
 ,exk(npnts)           & ! Exner ratio of layer k
 ,timestep             & ! Convection timestep (s)
 ,cca(npnts)             ! Convective cloud amount

real(kind=real_umphys), intent(in out) ::                                      &
  rain(npnts)            & ! in  Amount of falling rain (kg/m**2/s)
                           ! out Updated amount of falling rain
 ,snow(npnts)            & ! in  Amount of falling snow (kg/m**2/s)
                           ! out Updated amount of falling snow
 ,dthbydt(npnts)         & ! in  Increment to model potential temperature
                           ! out Updated Increment to model potential
                           !     temperature (K/s)
 ,dqbydt(npnts)            ! in  Increment to model mixing ratio (kg/kg/s)
                           ! out Updated Increment to model mixing ratio
                           !     after evaporation below cloud base.

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

integer ::                                                                     &
  i                 ! loop counter

logical ::                                                                     &
  bevap(npnts)          & ! Mask for those points where evaporation occurs
 ,l_rates_adjusted      & ! True if evaporation / sublimation rates have
                          ! been adjusted to give saturation
 ,full_evap_rain(npnts) & ! True if full rain evap in evp
 ,full_evap_snow(npnts)   ! True if full snow evap in evp

real(kind=real_umphys) ::                                                      &
  t(npnts)            & ! Model temperature (K)
 ,evap_rain(npnts)    & ! Amount of evaporation of rain
 ,sub_snow(npnts)     & ! Amount of snow sublimation
 ,delq(npnts)         & ! Change in mixing ratio across layer k  (kg/kg)
 ,ths(npnts)          & ! Saturated parcel potential temperature (K)
 ,qs(npnts)           & ! Saturated parcel parcel mixing ratio (kg/kg)
 ,dthbydt_evp         & ! Increment to potential temperature due to
                        ! evaporation (K)
 ,dqbydt_evp          & ! Increment to mixing ratio due to
                        ! evaporation (kg/kg)
 ,dthbydt_sat         & ! Increment to potential temperature due to
                        ! saturation (K)
 ,factor              & ! dthbydt_sat/dthbydt_evp
 ,rho(npnts)            ! Density of air in parcel

real(kind=real_umphys) ::                                                      &
  precip_area         & ! Area fraction of precip below cloud-base
                        !  = cca * precip_area_fac
 ,rtimestep             ! 1/timestep

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='PEVP_BCB'

!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Evaporate rain in layer k if layer k is below cloud base.
! Calculate moisture sub-saturation
!-----------------------------------------------------------------------

do i=1,npnts
  t(i) = th(i)*exk(i)
  bevap(i) = .false.

  if (k  <   iccb(i)) then
    delq(i) = qse(i)-q(i)

    !-----------------------------------------------------------------------
    ! Check if evaporation possible
    !-----------------------------------------------------------------------

    if ((rain(i) >  0.0 .or. snow(i) >  0.0) .and.                             &
        delq(i)  >   0.0) then

      bevap(i) = .true.
      rho(i) = pk(i) / (r*t(i))
    end if
  end if
end do

!-----------------------------------------------------------------------
! Calculate evaporation
!-----------------------------------------------------------------------

call evp (npnts,1,bevap,precip_area_fac,rain,t,cca,rho,delq,delp,              &
            pk,evap_rain,full_evap_rain)

call evp (npnts,2,bevap,precip_area_fac,snow,t,cca,rho,delq,delp,              &
            pk,sub_snow,full_evap_snow)

!-----------------------------------------------------------------------
! Calculate temperature and mixing ratio if layer brought to saturation
! by evaporation and sublimation
!-----------------------------------------------------------------------

call satcal(npnts,th,pk,exk,q,th,qse,qs,ths)

rtimestep = 1.0/timestep

do i=1,npnts
  if (bevap(i)) then
    dthbydt_evp = -((lc*evap_rain(i))+((lc+lf)*sub_snow(i)))*g/                &
                       (cp*exk(i)*delp(i))
    dqbydt_evp  = (evap_rain(i)+sub_snow(i))*g/delp(i)

    dthbydt_sat = (ths(i)-th(i))*rtimestep

    l_rates_adjusted = .false.

    if (dthbydt_evp <  dthbydt_sat) then

      !---------------------------------------------------------------------
      ! Adjust evaporation and sublimation rates to give saturation
      !---------------------------------------------------------------------

      l_rates_adjusted = .true.
      factor = dthbydt_sat/dthbydt_evp
      dthbydt_evp = dthbydt_sat
      dqbydt_evp = dqbydt_evp*factor
      evap_rain(i) = evap_rain(i)*factor
      sub_snow(i) = sub_snow(i)*factor
    end if

    !---------------------------------------------------------------------
    ! Update increments and rainfall and adjust back to gridbox means
    !---------------------------------------------------------------------

    precip_area = cca(i)*precip_area_fac
    dthbydt(i) = dthbydt(i)+dthbydt_evp*precip_area
    dqbydt(i)  = dqbydt(i) +dqbydt_evp*precip_area
    rain(i) = rain(i)-evap_rain(i)*precip_area
    snow(i) = snow(i)-sub_snow(i)*precip_area

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

  end if
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine pevp_bcb
end module pevp_bcb_mod
