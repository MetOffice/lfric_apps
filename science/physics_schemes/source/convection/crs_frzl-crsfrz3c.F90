! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Change of phase routine where precip crosses a melting or freezing level
!
module crs_frzl_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'CRS_FRZL_MOD'
contains

subroutine crs_frzl (npnts, bddwt_km1, exk, exkm1, flx_dd_km1,                 &
                     thdd_k, thdd_km1, rain, snow, wtrac_conv_old)

use cv_run_mod, only: t_melt_snow, pr_melt_frz_opt

use cv_param_mod, only: t_frez_rain, pr_melt_tdep, pr_melt_frz_tdep

use planet_constants_mod, only: g, cp
use water_constants_mod,  only: lf, tm
use wtrac_conv_mod,       only: l_wtrac_conv
use wtrac_conv_store_mod, only: conv_old_wtrac_type

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

!
! Description: Change of phase routine where precipitation crosses a melting
!              or freezing level (The code of this routine is almost identical
!              to chg_phse. At some future date it would be better to have
!              one common routine if possible.)
!
! Method: UM documentation paper 27
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards 8.2.
!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------
integer, intent(in) ::                                                         &
  npnts                  ! Vector length

logical, intent(in) ::                                                         &
  bddwt_km1(npnts)       ! Mask for those points in downdraught where
                         ! precipitation is liquid in layer k-1

real(kind=real_umphys), intent(in) ::                                          &
  exk(npnts)           & ! Exner ratio in layer k
 ,exkm1(npnts)         & ! Exner ratio in layer k-1
 ,flx_dd_km1(npnts)      ! Downdraught mass flux in layer k-1 (Pa/s)

real(kind=real_umphys), intent(in) ::                                          &
  thdd_k(npnts)          ! Potential temperature of DD in layer k (K)

real(kind=real_umphys), intent(in out) ::                                      &
  thdd_km1(npnts)      & ! In  Potential temperature of DD in layer k-1(K)
                         ! Out Potential temperature of DD in layer k-1
                         ! updated due to change of phase
 ,rain(npnts)          & ! In  Amount of rain descending from k-1 to k-2
                         ! Out Updated amount of rainfall (kg/m**2/s)
 ,snow(npnts)            ! In  Amount of snow descending from k-1 to k-2
                         ! Out Updated amount of snowfall (kg/m**2/s)

type(conv_old_wtrac_type), optional, intent(in out) :: wtrac_conv_old
                         ! Store of amount of precip freezing/melting
                         ! (optional as only used when called from 6A routine)

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

integer :: i          ! Loop counter

real(kind=real_umphys) ::                                                      &
  factor           &  ! Used in the calculation of the change of phase of
                      ! falling precipitation.
 ,precip_limit        ! Limit on amount of precip melting/freezing to avoid
                      ! a phase change
real(kind=real_umphys) ::                                                      &
  precip_frez(npnts)  &  ! Amount of precipitation freezing
                         ! (This field is stored for water tracer use)
 ,precip_melt(npnts)     ! Amount of precipitation melting
                         ! (This field is stored for water tracer use)
real(kind=real_umphys) :: melt_tscale
                      ! e-folding temperature thickness for snow melt
real(kind=real_umphys) :: r_melt_tscale2
                      ! reciprocal squared of e-folding temperature thickness
                      ! for snow melt
real(kind=real_umphys) :: frez_tscale
                      ! e-folding temperature thickness for freezing rain
real(kind=real_umphys) :: r_frez_tscale2
                      ! reciprocal squared of e-folding temperature thickness
                      ! for freezing rain

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CRS_FRZL'
!-----------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!-----------------------------------------------------------------------
! Add latent heating where precip crosses a melting or freezing level
!
!   UM Documentation paper 27
!   Section (11), equation (42)
!-----------------------------------------------------------------------

if (pr_melt_frz_opt == pr_melt_tdep .or.                                       &
    pr_melt_frz_opt == pr_melt_frz_tdep) then
  ! Melt snow at at rate proportional to the temperature
  ! above the melting temperature (usually 0C) and (if switched on) freeze
  ! rain below the melting temperature

  !Melting scale
  melt_tscale   = t_melt_snow - tm
  r_melt_tscale2= 1.0/melt_tscale**2

  !Freezing scale
  frez_tscale   = t_frez_rain - tm
  r_frez_tscale2= 1.0/frez_tscale**2

  do i=1,npnts

    precip_melt(i) = 0.0
    precip_frez(i) = 0.0
    factor = (exkm1(i)*cp*flx_dd_km1(i))/(lf*g)

    if ( (thdd_km1(i)*exkm1(i) > tm) .and. (snow(i) > 0.0) ) then
      ! Melt
      if ( (3.0*melt_tscale) > (thdd_km1(i)*exkm1(i)-tm) ) then
        ! If temperature still below emergency melting level
        !then calculate melting rate
        precip_melt(i)   = snow(i)*(thdd_km1(i)*exkm1(i)-tm)                   &
                      *(thdd_km1(i)*exkm1(i)-thdd_k(i)*exk(i))*r_melt_tscale2
      else
        ! If above emergency melting level melt all snow
        precip_melt(i)   = snow(i)
      end if

      ! Check that precip_melt is not negative.
      ! to prevent refreezing of rain
      precip_melt(i)     = max(0.0,precip_melt(i))

      ! Limit the amount of melting so that is cannot
      ! lower the temperature to below tm.
      precip_limit    = (thdd_km1(i) - tm/exkm1(i))*factor
      precip_melt(i)     = min(precip_limit,precip_melt(i))

      ! Limit the melting to the amount of snow available
      precip_melt(i)     = min(snow(i),precip_melt(i))

      ! Calculate temperature change and update rain and snow
      thdd_km1(i)     = thdd_km1(i) - precip_melt(i)/factor
      rain(i)         = rain(i)     + precip_melt(i)
      snow(i)         = snow(i)     - precip_melt(i)

    end if  ! test on melting

    if (pr_melt_frz_opt == pr_melt_frz_tdep .and.                              &
       (thdd_km1(i)*exkm1(i) < tm) .and. (rain(i) > 0.0) ) then
      ! Calculate the freezing rate. Note no emergency freeze.
      precip_frez(i)     = rain(i)*(tm-thdd_km1(i)*exkm1(i))                   &
                      *(thdd_km1(i)*exkm1(i)-thdd_k(i)*exk(i))*r_frez_tscale2

      ! Check that precip_frez is not negative.
      ! to prevent melting of snow
      precip_frez(i)     = max(0.0,precip_frez(i))

      ! Limit the amount of melting so that it cannot
      ! raise the temperature to above tm.
      precip_limit    = (tm/exkm1(i) - thdd_km1(i))*factor
      precip_frez(i)     = min(precip_limit,precip_frez(i))

      ! Limit the freezing to the amount of rain
      precip_frez(i)     = min(rain(i),precip_frez(i))

      !Calculate temperature increment and update rain and snow
      thdd_km1(i)     = thdd_km1(i) + precip_frez(i)/factor
      rain(i)         = rain(i)     - precip_frez(i)
      snow(i)         = snow(i)     + precip_frez(i)

    end if  ! test on freezing
  end do

else
  ! original code melts snow at freezing level and immediately
  ! freezes any supercooled rain.
  do i=1,npnts

    precip_melt(i) = 0.0
    precip_frez(i) = 0.0
    factor = (exkm1(i)*cp*flx_dd_km1(i))/(lf*g)

    if (.not. bddwt_km1(i) .and. rain(i) >  0.0 .and. thdd_km1(i)              &
         *exkm1(i) <  tm) then
      ! Freeze
      precip_frez(i) = (tm/exkm1(i)-thdd_km1(i))* factor
      precip_frez(i) = min(rain(i),precip_frez(i))
      thdd_km1(i) = thdd_km1(i)+precip_frez(i)/factor
      rain(i) = rain(i)-precip_frez(i)
      snow(i) = snow(i)+precip_frez(i)

    else if (bddwt_km1(i) .and. snow(i) >  0.0) then
      ! Melt
      precip_melt(i) = (thdd_km1(i)-tm/exkm1(i))*factor
      precip_melt(i) = min(snow(i),precip_melt(i))
      thdd_km1(i) = thdd_km1(i)-precip_melt(i)/factor
      rain(i) = rain(i)+precip_melt(i)
      snow(i) = snow(i)-precip_melt(i)

    end if  ! test on melting/freezing

  end do

end if

! If water tracers, save the amount of melt/freeze.  (Note, if
! pr_melt_frz_opt == pr_melt_tdep .or. pr_melt_frz_opt == pr_melt_frz_tdep)
! then there can be melt followed by freeze at the same point.)
if (l_wtrac_conv) then
  do i = 1, npnts
    wtrac_conv_old%precip_frez(i) = precip_frez(i) - precip_melt(i)
  end do
end if


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine crs_frzl

end module crs_frzl_mod
