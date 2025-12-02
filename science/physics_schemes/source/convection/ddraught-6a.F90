! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Downdraught routine
!
module ddraught_6a_mod

use um_types, only: real_umphys

implicit none

! Description: Downdraught routine
!              Convective downdraught based on parcel theory.
!              Carry out dry descent.
!              Calculate subsaturation.
!              Calculate effect on the environment.
!
! Method: UM documentation paper 27
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards vn8.2.

character(len=*), parameter, private :: ModuleName = 'DDRAUGHT_6A_MOD'

contains

subroutine ddraught_6a (npnts, np_full, k, kct, ntra, n_wtrac,                 &
                     l_tracer,                                                 &
                     thdd_k,qdd_k,the_k,                                       &
                     the_km1,qe_k,qe_km1,qse_km1,dthbydt_k,dthbydt_km1,        &
                     dqbydt_k,dqbydt_km1,flx_dd_k,p_km1,delpk,                 &
                     delpkm1,exk,exkm1,deltd,delqd,amdetk,ekm14,               &
                     ekm34,rain,snow,                                          &
                     bdd_start,bddwt_k,bddwt_km1,                              &
                     bdd_on,b_dd_end,                                          &
                     deltrad, cca, ppn_mix_dd,                                 &
                     tradd_k,trae_k,trae_km1,dtrabydt_k,                       &
                     dtrabydt_km1,wtrac_dd2)

use water_constants_mod, only: tm
use cv_run_mod, only: pr_melt_frz_opt
use cv_param_mod, only: pr_melt_tdep, pr_melt_frz_tdep

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

use crs_frzl_mod, only: crs_frzl
use dd_env_6a_mod, only: dd_env_6a
use devap_mod, only: devap
use satcal_mod, only: satcal
use termdd_mod, only: termdd

use qsat_mod, only: qsat

use wtrac_conv_mod,            only: l_wtrac_conv, conv_dd2_wtrac_type
use wtrac_conv_store_mod,      only: conv_old_wtrac_type,                      &
                                     wtrac_alloc_conv_store2
use wtrac_precip_chg_phse_mod, only: wtrac_precip_chg_phse
use wtrac_precip_evap_mod,     only: wtrac_precip_evap

implicit none

!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------

integer, intent(in) ::                                                         &
  npnts                &  ! Vector length
 ,np_full              &  ! Full vector length
 ,k                    &  ! Present model layer
 ,kct                  &  ! Convective cloud top
 ,ntra                 &  ! Number of tracers
 ,n_wtrac                 ! Number of water tracers

logical, intent(in) ::                                                         &
  l_tracer                ! Switch for tracers


real(kind=real_umphys), intent(in) ::                                          &
  the_km1(npnts)        & ! Potential temperature of enviroment in layer k-1(K)
 ,qe_km1(npnts)         & ! Mixing ratio of enviroment in layer k-1 (kg/kg)
 ,qse_km1(npnts)        & ! qsat Mixing ratio of enviroment in layer k-1 (kg/kg)
 ,p_km1(npnts)          & ! Pressure in layer k-1  (Pa)
 ,delpk(npnts)          & ! Change in pressure across layer k (Pa)
 ,delpkm1(npnts)        & ! Change in pressure across layer k-1  (Pa)
 ,exk(npnts)            & ! Exner ratio in layer k
 ,exkm1(npnts)          & ! Exner ratio in layer k-1
 ,amdetk(npnts)         & ! Mixing detrainment rate
 ,ekm14(npnts)          & ! Exner ratio at layer k-1/4
 ,ekm34(npnts)          & ! Exner ratio at layer k-3/4
 ,deltd(npnts)          & ! Cooling necessary to achieve saturation (K)
 ,delqd(npnts)          & ! moistening necessary to achieve saturation (kg/kg)
 ,cca(npnts)            & ! Convective cloud amount (fraction)
 ,ppn_mix_dd(npnts)       ! precipitation mixing ratio (kg/kg)

real(kind=real_umphys), intent(in) ::                                          &
  trae_km1(np_full,ntra)  & ! Tracer content of enviroment in layer k-1 (kg/kg)
 ,deltrad(npnts,ntra)       ! Depletion of environment tracer due to
                            ! downdraught formation (kg/kg)

logical, intent(in out) ::                                                     &
  bdd_on(npnts)            & ! In  Mask for those points where DD has continued
                             !     from layer k+1
                             ! Out Mask for those points where DD continues
                             !     to layer k-1
 ,bddwt_k(npnts)           & ! In Mask for those points in downdraught where
                             ! precipitation is liquid in layer k
 ,bddwt_km1(npnts)         & ! Mask for those points in downdraught where
                             ! precipitation is liquid in layer k-1
 ,bdd_start(npnts)           ! Mask for those points where DD may start
                             ! in layer k-1

real(kind=real_umphys), intent(in out) ::                                      &
  thdd_k(npnts)        & ! In  Potential temperature of downdraught in layer k
                         ! Out Potential temperature reset for next layer (K)
 ,qdd_k(npnts)         & ! In  Mixing ratio of downdraught in layer k
                         ! Out Mixing ratio reset for next layer (kg/kg)
 ,tradd_k(np_full,ntra)& ! In  Downdraught tracer content of layer k
                         ! Out tracer content reset for next layer (kg/kg)
 ,the_k(npnts)         & ! In  Potential temperature of environment in layer k
                         ! Out Potential temperature reset for next layer (K)
 ,qe_k(npnts)          & ! In  Mixing ratio of environment in layer k
                         ! Out Mixing ratio reset for next layer (kg/kg)
 ,trae_k(np_full,ntra) & ! In  environment tracer content of layer k
                         ! Out tracer content reset for next layer (kg/kg)
 ,flx_dd_k(npnts)      & ! In  Downdraught mass flux of layer k
                         ! Out Downdraught mass flux reset for next layer (Pa/s)
 ,rain(npnts)          & ! In  Amount of rain
                         ! Out Updated amount of rainfall (kg/m**2/s)
 ,snow(npnts)          & ! In  Amount of snow
                         ! Out Updated amount of snowfall (kg/m**2/s)
 ,dthbydt_k(npnts)     & ! In  Increment to potential temperature of layer k
                         ! Out Updated increment potential temperature layer k
                         !           (K/s)
 ,dthbydt_km1(npnts)   & ! In  Increment to potential temperature of layer k-1
                         ! Out Updated increment potential temperature layer k-1
                         !           (K/s)
 ,dqbydt_k(npnts)      & ! In  Increment to mixing ratio of layer k
                         ! Out Updated increment mixing ratio layer k (kg/kg)
 ,dqbydt_km1(npnts)      ! In  Increment to mixing ratio  of layer k-1
                         ! Out Updated increment mixing ratio layer k-1 (kg/kg)

real(kind=real_umphys), intent(in out) ::                                      &
  dtrabydt_k(np_full,ntra)   & ! In  Increment to tracer content of layer k
                               ! Out Updated increment tracer content layer k
                               ! (kg/kg)
 ,dtrabydt_km1(np_full,ntra)   ! In  Increment to tracer content of layer k-1
                               ! Out Updated increment tracer content layer k-1
                               ! (kg/kg)

type(conv_dd2_wtrac_type), intent(in out) :: wtrac_dd2(n_wtrac)
                               ! Structure containing 2nd compression water
                               ! tracer fields used for downdraught calculation

logical, intent(out) ::                                                        &
  b_dd_end(npnts)          ! Mask for those points where DD is ending in
                           ! layer k-1



!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

integer ::                                                                     &
  i, ktra, i_wt         ! Loop counters

real(kind=real_umphys) ::                                                      &
  thdd_km1(npnts)     & ! Potential temperature of downdraught in layer k-1 (K)

 ,qdd_km1(npnts)      & ! Downdraught mixing ratio of layer k-1 (kg/kg)

 ,qsatdd(npnts)       & ! Saturated downdraught mixing ratio of layer k-1
                        ! (kg/kg)
 ,tdd_km1(npnts)      & ! Temperature of downdraught in layer k-1 (K)

 ,thdds(npnts)        & ! Potential temperature of saturated downdraught
                        ! in layer k-1 (K)
 ,qdds(npnts)         & ! Saturated downdraught mixing ratio of layer k-1
                        ! (kg/kg)
 ,flx_dd_km1(npnts)   & ! Downdraught mass flux in layer K-1 (Pa/s)

 ,rain_tmp(npnts)     & ! Liquid precipitation store

 ,snow_tmp(npnts)       ! Snow Store

real(kind=real_umphys) ::                                                      &
  tradd_km1(npnts,ntra)   ! Tracer content of downdraught in layer k-1 (kg/kg)

real(kind=real_umphys), allocatable :: qdd_km1_wtrac(:,:)
                                         ! Downdraught water tracer mixing
                                         ! ratio of layer k-1 (kg/kg)
real(kind=real_umphys), allocatable :: rain_tmp_wtrac(:,:)
                                         ! Water tracer rain store
real(kind=real_umphys), allocatable :: snow_tmp_wtrac(:,:)
                                         ! Water tracer snow store

type(conv_old_wtrac_type) :: wtrac_conv_old  ! Store of water values prior to
                                             ! phase change

real(kind=real_umphys) ::                                                      &
  ekm14_plus1           & ! 1+ekm14
 ,ekm34_plus1           & ! 1+ekm34
 ,dnom                  & ! (1+ekm14)*(1+ekm34)
 ,rdnom                   ! 1/dnom

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='DDRAUGHT_6A'

!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Calculate mask for those points in downdraught where precipitation
! is liquid
!
! Store precipitation in layer k in temporary variables.
!-----------------------------------------------------------------------

do i=1,npnts
  if (k  ==  kct .or. bdd_start(i)) then
    bddwt_k(i) = thdd_k(i)*exk(i)  >   tm
  else
    bddwt_k(i) = bddwt_km1(i)
  end if

  rain_tmp(i) = rain(i)
  snow_tmp(i) = snow(i)

  !-----------------------------------------------------------------------
  ! Dry descent from layer k to k-1
  !
  ! Entrainment calculation
  !-----------------------------------------------------------------------

  ekm14_plus1 = 1.0+ekm14(i)
  ekm34_plus1 = 1.0+ekm34(i)
  dnom        = ekm14_plus1 *ekm34_plus1
  rdnom        = 1.0/dnom

  thdd_km1(i) = (thdd_k(i)+(ekm14(i)*the_k(i)) +                               &
                                   ekm14_plus1*ekm34(i)*the_km1(i))*rdnom

  qdd_km1(i) = (qdd_k(i)+(ekm14(i)*qe_k(i)) +                                  &
                                   ekm14_plus1*ekm34(i)*qe_km1(i))*rdnom

  !-----------------------------------------------------------------------
  ! Update mass flux and calculate temperature of layer k-1
  !-----------------------------------------------------------------------

  flx_dd_km1(i) = flx_dd_k(i)*dnom*(1.0-amdetk(i))

  tdd_km1(i)    = thdd_km1(i)*exkm1(i)

end do

!----------------------------------------------------------------------
! Dry descent for tracers
!----------------------------------------------------------------------

if (l_tracer) then

  do ktra=1,ntra
    do i=1,npnts
      ekm14_plus1 = 1.0+ekm14(i)
      dnom        = ekm14_plus1 * (1.0+ekm34(i))

      tradd_km1(i,ktra)=(tradd_k(i,ktra)+(ekm14(i)*trae_k(i,ktra))             &
                             + ekm14_plus1*ekm34(i)*trae_km1(i,ktra))/dnom


    end do
  end do

end if

if (l_wtrac_conv) then
  ! Dry descent from layer k to k-1 for water tracers
  ! (Entrainment calculation)

  ! Allocate water tracer local arrays
  allocate(qdd_km1_wtrac(npnts,n_wtrac))
  allocate(rain_tmp_wtrac(npnts,n_wtrac))
  allocate(snow_tmp_wtrac(npnts,n_wtrac))

  do i_wt=1,n_wtrac
    do i=1,npnts
      ekm14_plus1 = 1.0+ekm14(i)
      ekm34_plus1 = 1.0+ekm34(i)
      dnom        = ekm14_plus1 *ekm34_plus1
      rdnom        = 1.0/dnom

      qdd_km1_wtrac(i,i_wt) =                                                  &
                (wtrac_dd2(i_wt)%qdd_k(i)+(ekm14(i)*wtrac_dd2(i_wt)%qe_k(i))   &
                 + ekm14_plus1*ekm34(i)*wtrac_dd2(i_wt)%qe_km1(i))*rdnom

      rain_tmp_wtrac(i,i_wt) = wtrac_dd2(i_wt)%rain(i)
      snow_tmp_wtrac(i,i_wt) = wtrac_dd2(i_wt)%snow(i)

    end do
  end do
else
  ! Note, no need to allocate rain_tmp_wtrac or snow_tmp_wtrac
  allocate(qdd_km1_wtrac(1,1))
end if

!-----------------------------------------------------------------------
! Calculate subsaturation
! Calculate temperature if brought to saturation
!-----------------------------------------------------------------------
if (pr_melt_frz_opt == pr_melt_tdep .or.                                       &
    pr_melt_frz_opt == pr_melt_frz_tdep) then
  ! No longer calculating saturation temperature as seems to contribute to noise
  ! when using the option to melt snow over a range of temperatures.
  do i=1,npnts
    bddwt_km1(i) = thdd_km1(i)*exkm1(i)  >   tm
  end do
else
  call satcal(npnts,thdd_km1,p_km1,exkm1,qdd_km1,the_km1,qse_km1,              &
            qdds,thdds)


  do i=1,npnts
    bddwt_km1(i) = thdds(i)*exkm1(i)  >   tm
  end do
end if

!-----------------------------------------------------------------------
! Calculate change of phase due to downdraught saturation temperature or
!   downdraught temperature depending on option
!-----------------------------------------------------------------------

if (l_wtrac_conv) then
  ! Store values before phase change for use in water tracer calculations
  call wtrac_alloc_conv_store2(npnts, qdd_km1, rain, snow, wtrac_conv_old)
end if

call crs_frzl(npnts, bddwt_km1, exk, exkm1, flx_dd_km1, thdd_k, thdd_km1,      &
              rain, snow, wtrac_conv_old=wtrac_conv_old)

if (l_wtrac_conv) then
  ! Update water tracers for phase change
  call wtrac_precip_chg_phse(npnts, n_wtrac, rain, snow, wtrac_conv_old,       &
                             wtrac_dd2 = wtrac_dd2)
end if

do i=1,npnts
  tdd_km1(i) = thdd_km1(i)*exkm1(i)
end do

!-----------------------------------------------------------------------
! Recalculate subsaturation temperature
!-----------------------------------------------------------------------

call satcal(npnts,thdd_km1,p_km1,exkm1,qdd_km1,the_km1,qse_km1,                &
                 qdds, thdds)

!-----------------------------------------------------------------------
! Calculate moisture subsaturation
!-----------------------------------------------------------------------
call qsat(qsatdd,tdd_km1,p_km1,npnts)

!-----------------------------------------------------------------------
! Evaporation calculation and adujstment of downdraught temperature
! and moisture
!-----------------------------------------------------------------------

call devap(npnts, bddwt_km1                                                    &
           , thdd_k, thdds, qdds, flx_dd_km1, exk, exkm1                       &
           , qsatdd, delpkm1, cca, p_km1                                       &
           , thdd_km1, qdd_km1, rain, snow)

if (l_wtrac_conv) then
  ! Update water tracers for phase change (rain/snow -> vapour)
  call wtrac_precip_evap (npnts, n_wtrac, flx_dd_km1, qdd_km1, rain, snow,     &
                          thdd_km1, exkm1, qsatdd,                             &
                          wtrac_conv_old, qdd_km1_wtrac, l_evap_call=.false.,  &
                          wtrac_dd2 = wtrac_dd2)
  ! Update water tracers for phase change (vapour -> rain/snow)
  call wtrac_precip_evap (npnts, n_wtrac, flx_dd_km1, qdd_km1, rain, snow,     &
                          thdd_km1, exkm1, qsatdd,                             &
                          wtrac_conv_old, qdd_km1_wtrac, l_evap_call=.true.,   &
                          wtrac_dd2 = wtrac_dd2)
end if

!-----------------------------------------------------------------------
! Check if parcel still negatively buoyant such that downdraught can
! continue to k-1
!-----------------------------------------------------------------------

call termdd(npnts, k, bdd_start                                                &
           , thdd_km1, qdd_km1, the_km1, qe_km1, ppn_mix_dd                    &
           , b_dd_end, bdd_on)

!-----------------------------------------------------------------------
! Calculate the effect on the environment in layer k
!-----------------------------------------------------------------------

call dd_env_6a  (npnts, np_full, ntra, n_wtrac                                 &
                ,l_tracer, b_dd_end, bdd_start, bdd_on                         &
                ,thdd_k, thdd_km1, qdd_k, qdd_km1, the_k, the_km1              &
                ,qe_k, qe_km1, flx_dd_k, flx_dd_km1, delpk, delpkm1            &
                ,deltd, delqd, amdetk, ekm14                                   &
                ,tradd_k, tradd_km1, trae_k, trae_km1, deltrad, qdd_km1_wtrac  &
                ,dthbydt_k, dthbydt_km1, dqbydt_k, dqbydt_km1                  &
                ,dtrabydt_k, dtrabydt_km1, wtrac_dd2)


!-----------------------------------------------------------------------
! Reset downdraught bit vectors
!-----------------------------------------------------------------------

! (Reset water tracers first before bdd_on is reset)
if (l_wtrac_conv) then
  do i_wt = 1, n_wtrac
    do i=1,npnts
      if (.not. bdd_on(i)) then
        wtrac_dd2(i_wt)%rain(i) = rain_tmp_wtrac(i,i_wt)
        wtrac_dd2(i_wt)%snow(i) = snow_tmp_wtrac(i,i_wt)
      end if
    end do
  end do
end if

do i=1,npnts
  bdd_start(i) = .false.
  if (.not. bdd_on(i)) then
    rain(i) = rain_tmp(i)
    snow(i) = snow_tmp(i)
  end if
  if (b_dd_end(i)) bdd_on(i) = .false.
end do


!-----------------------------------------------------------------------
! Switch potential temperature, mixing ratio,  mass flux and
! tracer ready for calculation at next model layer
!-----------------------------------------------------------------------

if (k >  2) then
  do i=1,npnts
    if (bdd_on(i)) then
      thdd_k(i)   = thdd_km1(i)
      qdd_k(i)    = qdd_km1(i)
      flx_dd_k(i) = flx_dd_km1(i)
    end if
  end do

  if (l_tracer) then

    do ktra=1,ntra
      do i=1,npnts
        if (bdd_on(i)) then
          tradd_k(i,ktra) = tradd_km1(i,ktra)
        end if
      end do
    end do

  end if

  if (l_wtrac_conv) then
    do i_wt = 1, n_wtrac
      do i=1,npnts
        if (bdd_on(i)) then
          wtrac_dd2(i_wt)%qdd_k(i) = qdd_km1_wtrac(i,i_wt)
        end if
      end do
    end do
  end if

end if  ! K > 2

! Deallocate water tracer local arrays
if (l_wtrac_conv) then
  deallocate(snow_tmp_wtrac)
  deallocate(rain_tmp_wtrac)
end if

deallocate(qdd_km1_wtrac)


!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine ddraught_6a
end module ddraught_6a_mod
