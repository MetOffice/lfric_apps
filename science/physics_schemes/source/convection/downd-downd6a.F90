! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calls downdraught calculation or checks for change of phase
!
module downd_6a_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'DOWND_6A_MOD'
contains

subroutine downd_6a (npnts,np_full,k,kct,ntra,n_wtrac,nddon_a, iccb,           &
                  l_tracer,bwater_k,                                           &
                  timestep,the_k,the_km1,qe_k,qe_km1,qse_km1,                  &
                  p_km1,delpk,                                                 &
                  delpkm1,exk,exkm1,deltd,delqd,amdetk,ekm14,                  &
                  ekm34,flx_ud_k,cca, trae_k,trae_km1,deltrad,                 &
                  bdd_start,bddwt_k,bddwt_km1,bdd_on,                          &
                  thdd_k,qdd_k,dthbydt_k,dthbydt_km1,dqbydt_k,                 &
                  dqbydt_km1,rain,snow,precip_k,rain_env,snow_env,             &
                  rain_dd,snow_dd,flx_dd_k,                                    &
                  tradd_k,dtrabydt_k,dtrabydt_km1,wtrac_dd,wtrac_ev)

use planet_constants_mod, only: g

use cv_param_mod, only:                                                        &
    ddptef

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

use chg_phse_mod, only: chg_phse
use ddraught_6a_mod, only: ddraught_6a
use pevp_bcb_mod, only: pevp_bcb

use wtrac_conv_mod,            only: l_wtrac_conv, conv_dd_wtrac_type,         &
                                  conv_dd2_wtrac_type, conv_ev_wtrac_type,     &
                                  wtrac_alloc_conv_dd2, wtrac_dealloc_conv_dd2
use wtrac_conv_store_mod,      only: conv_old_wtrac_type,                      &
                                     wtrac_alloc_conv_store2
use wtrac_precip_chg_phse_mod, only: wtrac_precip_chg_phse
use wtrac_precip_evap_mod,     only: wtrac_precip_evap

implicit none
! ------------------------------------------------------------------------------
! Description:
! Calls downdraught calculation.
! Change of phase calculation where no downdraught
!
!  See UM Documentation paper No. 27
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!
! ------------------------------------------------------------------------------
! Subroutine arguments

integer, intent(in) ::                                                         &
  npnts                & ! Number of points
 ,np_full              & ! Full vector length
 ,k                    & ! Present model layer
 ,kct                  & ! Convective cloud top layer
 ,ntra                 & ! Number of tracers
 ,n_wtrac              & ! Number of water tracers
 ,nddon_a                ! Number of points at which downdraught does occur

integer, intent(in) ::                                                         &
  iccb(npnts)            ! Cloud base model level

logical, intent(in) ::                                                         &
  l_tracer               ! Switch for tracers

logical, intent(in) ::                                                         &
  bwater_k(npnts)        ! Mask for points at which condensate is water
                         ! in layer k

real(kind=real_umphys), intent(in) ::                                          &
  timestep               ! timestep

real(kind=real_umphys), intent(in) ::                                          &
  the_k(npnts)      & ! potential temperature of environment in layer k (K)
 ,the_km1(npnts)    & ! Potential temperature of environment in layer k-1 (K)
 ,qe_k(npnts)       & ! environment mixing ratio of layer k (kg/kg)
 ,qe_km1(npnts)     & ! environment mixing ratio of layer k-1 (kg/kg)
 ,qse_km1(npnts)    & ! environment qsat mixing ratio of layer k-1 (kg/kg)
 ,p_km1(npnts)      & ! Pressure of layer k-1 (Pa)
 ,delpk(npnts)      & ! Pressure difference across layer K  (Pa)
 ,delpkm1(npnts)    & ! Pressure difference across layer K-1 (Pa)
 ,exk(npnts)        & ! exner ratio for layer K
 ,exkm1(npnts)      & ! exner ratio for layer K-1
 ,deltd(npnts)      & ! Cooling necessary to achieve saturation (K)
 ,delqd(npnts)      & ! moistening necessary to achieve saturation (kg/kg)
 ,amdetk(npnts)     & ! Mixing detrainment at level k multiplied by
                      ! appropriate layer thickness
 ,ekm14(npnts)      & ! exner ratio at layer k-1/4
 ,ekm34(npnts)      & ! exner ratio at layer k-3/4
 ,flx_ud_k(npnts)   & ! updraught mass flux at layer K
 ,cca(npnts)          ! convective cloud amount

real(kind=real_umphys), intent(in) ::                                          &
  trae_k(npnts,ntra)   & ! tracer of environment in layer k (kg/kg)
 ,trae_km1(npnts,ntra) & ! tracer of environment in layer k-1 (kg/kg)
 ,deltrad(npnts,ntra)    ! Depletion of environment tracer due to

                         ! downdraught formation (kg/kg)

logical, intent(in out) ::                                                     &
  bdd_start(npnts)        & !in Mask for those points where downdraught may
                            !   form in layer k
                            !out Mask for those points where downdraught may
                            !   form in layer k-1
 ,bddwt_k(npnts)          & !   Mask for those points in downdraught where
                            !   precipitation is liquid in layer k
 ,bddwt_km1(npnts)        & !   Mask for those points in downdraught where
                            !   precipitation is liquid in layer k-1
 ,bdd_on(npnts)             !in Mask for those points where DD has continued
                            !   from previous layer
                            !out Mask for those points where downdraught
                            !    continues to layer k-1

real(kind=real_umphys), intent(in out) ::                                      &
  thdd_k(npnts)        & ! in Model potential temperature of downdraught
                         !   in layer k (K)
                         !out Updated potential temperature of downdraught
                         !   in layer k (K)
 ,qdd_k(npnts)         & ! in mixing ratio of downdraught in layer k (kg/kg)
                         !out Updated mixing ratio of downdraught in layer k
 ,dthbydt_k(npnts)     & ! in Increment to model potential temperature of
                         !    layer k (K/s)
                         !out Updated increment to model potential temperature
                         !    of layer k (K/s)
 ,dthbydt_km1(npnts)   & ! in Increment to model potential temperature of
                         !    layer k-1 (K/s)
                         !out Updated increment to model potential temperature
                         !    of layer-1 (K/s)
 ,dqbydt_k(npnts)      & ! in Increment to model mixing ratio of layer k
                         !    (kg/kg/s)
                         !out Updated increment to  mixing ratio of layer k
                         !    (kg/kg/s)
 ,dqbydt_km1(npnts)    & ! in Increment to model mixing ratio of layer k-1
                         !    (kg/kg/s)
                         !out Updated increment to  mixing ratio of layer k-1
                         !    (kg/kg/s)
 ,rain_env(npnts)      & ! Amount of rainfall passing through environment
                         ! (kg/m**2/s)
 ,snow_env(npnts)      & ! Amount of snowfall passing through environment
                         ! (kg/m**2/s)
 ,rain_dd(npnts)       & ! Amount of rainfall passing through
                         ! downdraught (KG/M**2/S)
 ,snow_dd(npnts)       & ! Amount of snowfall passing through
                         ! downdraught (KG/M**2/S)
 ,rain(npnts)          & ! in Initialised rainfall (kg/m**2/s)
                         ! out Surface rainfall (kg/m**2/s)
 ,snow(npnts)          & ! in Initialised snowfall (kg/m**2/s)
                         ! out Surface snowfall (kg/m**2/s)
 ,precip_k(npnts)      & ! Precipitation added when descending from layer
                         ! k to k-1 (kg/m**2/s)
 ,flx_dd_k(npnts)        ! Downdraught initial mass flux (Pa/s)

real(kind=real_umphys), intent(in out) ::                                      &
  tradd_k(npnts,ntra)       & ! Model tracer of downdraught in layer k (kg/kg)
 ,dtrabydt_k(np_full,ntra)  & ! in Increment to model tracer of layer k
                              !    (kg/kg/s)
                              ! out Updated increment to model tracer of layer k
                              !    (kg/kg/s)
 ,dtrabydt_km1(np_full,ntra)  ! in Increment to model tracer of layer k-1
                              !    (kg/kg/s)
                              ! out Updated increment to model tracer of
                              !    layer k-1 (kg/kg/s)

type(conv_dd_wtrac_type), intent(in out) :: wtrac_dd(n_wtrac)
                              ! Structure containing water tracer fields
                              ! used in downdraught calculations

type(conv_ev_wtrac_type), intent(in out) :: wtrac_ev(n_wtrac)
                              ! Structure containing water tracer fields
                              ! used for phase change and evap of precip
                              ! in environment calculations

!-------------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

integer ::                                                                     &
  i,ktra,i_wt     & ! loop counter
 ,nddon             ! Number of points at which downdraught does occur

integer ::                                                                     &
  index1(nddon_a)   !  Index for compressing points

logical ::                                                                     &
  bwork(nddon_a,5)  & !  Work space for 'bit' masks
 ,b_dd_end(npnts)     !  mask for points where downdraught has ended


real(kind=real_umphys) ::                                                      &
  work(nddon_a,38)      ! Work space  for compression

real(kind=real_umphys) ::                                                      &
  qse_km1_c(nddon_a)      ! qsat k-1 for compression

real(kind=real_umphys) ::                                                      &
  tradd_k_c(nddon_a,ntra)  & ! Tracer content in downdraught at layer k
                             ! compressed (kg/kg)
 ,trae_k_c(nddon_a,ntra)   & ! Tracer content of environment layer k
                             ! compressed (kg/kg)
 ,trae_km1_c(nddon_a,ntra) & ! Tracer content of environment layer k-1
                             ! compressed (kg/kg)
 ,dtra_k_c(nddon_a,ntra)   & ! Increment to model tracer in layer k
                             ! compressed (kg/kg)
 ,dtra_km1_c(nddon_a,ntra) & ! Increment to model tracer in layer k-1
                             ! compressed (kg/kg)
 ,deltrad_c(nddon_a,ntra)    ! Depletion of environmen  tracer due to
                             ! downdraught formation compressed

real(kind=real_umphys) ::                                                      &
  ppn_mix_dd_c(nddon_a)     ! Precip mixing ratio in DD

real(kind=real_umphys) ::                                                      &
  factor_env          ! Proportion of rainfall going into DD from falling ppn

real(kind=real_umphys), allocatable :: qe_km1_wtrac(:,:)
                            ! Water tracer mixing ratio of layer K-1

type(conv_dd2_wtrac_type) :: wtrac_dd2(n_wtrac)
                            ! Structure containing water tracer fields
                            ! following 2nd compression used in downdraught
                            ! calculations

type(conv_old_wtrac_type) :: wtrac_conv_old
                            ! Store of water values prior to phase change

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='DOWND_6A'

!-----------------------------------------------------------------------
! Start of main loop
!   Update precipitation and calculate mask for where precipitation
!   is liquid.
!-----------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do i=1,npnts
  b_dd_end(i) = .false.
end do

if (k == kct) then
  do i=1,npnts
    rain_dd(i) = 0.0
    rain_env(i) = 0.0
    snow_dd(i) = 0.0
    snow_env(i) = 0.0
  end do
end if

if (l_wtrac_conv) then
  if (k == kct) then
    do i_wt = 1, n_wtrac
      do i=1,npnts
        wtrac_dd(i_wt)%rain_dd(i)  = 0.0
        wtrac_ev(i_wt)%rain_env(i) = 0.0
        wtrac_dd(i_wt)%snow_dd(i)  = 0.0
        wtrac_ev(i_wt)%snow_env(i) = 0.0
      end do
    end do
  end if
end if

!----------------------------------------------------------------------
! Injection of precipitation from UD at level k
!----------------------------------------------------------------------


!----------------------------------------------------------------------
! Corrected DD
! Transfer all the precipitation created in the UD into the environment
!----------------------------------------------------------------------
do i=1,npnts

  if (bwater_k(i)) then
    rain_env(i) = rain_env(i) + precip_k(i)
  else
    snow_env(i) = snow_env(i) + precip_k(i)
  end if
  precip_k(i) = 0.0
end do

if (l_wtrac_conv) then
  do i_wt = 1, n_wtrac
    do i=1,npnts

      if (bwater_k(i)) then
        wtrac_ev(i_wt)%rain_env(i) =                                           &
                     wtrac_ev(i_wt)%rain_env(i) + wtrac_dd(i_wt)%precip_k(i)
      else
        wtrac_ev(i_wt)%snow_env(i) =                                           &
                     wtrac_ev(i_wt)%snow_env(i) + wtrac_dd(i_wt)%precip_k(i)
      end if
      wtrac_dd(i_wt)%precip_k(i) = 0.0

    end do
  end do

end if

!----------------------------------------------------------------------
! Interaction of downdraught with reserve of precipitation outside
! downdraught

! Based upon continuuity of precipitation mixing ratio within
! downdraught - either after injection of rain from UD in level
! k or with PPN mixing ratio in lowest preciptating layer

! if downdraught increases in mass then water injected
! if downdraught decreases in mass then water is removed

!----------------------------------------------------------------------

do i=1,npnts

  if (bdd_on(i)) then

    if (flx_ud_k(i) > 0.0) then
      factor_env = ddptef*flx_dd_k(i)/flx_ud_k(i)*                             &
      delpk(i)/5000.0
    else
      factor_env = 1.0*delpk(i)/5000.0
    end if
    factor_env = amin1(factor_env,1.0)


    if (factor_env >  0.0) then
      rain_dd(i)  = rain_dd(i) + rain_env(i)*factor_env
      rain_env(i) = rain_env(i) * (1.0-factor_env)
      snow_dd(i)  = snow_dd(i) + snow_env(i)*factor_env
      snow_env(i) = snow_env(i) * (1.0-factor_env)
    else
      rain_env(i) = rain_env(i) - rain_dd(i)*factor_env
      rain_dd(i)  = rain_dd(i) * (1.0+factor_env)
      snow_env(i) = snow_env(i) - snow_dd(i)*factor_env
      snow_dd(i)  = snow_dd(i) * (1.0+factor_env)
    end if

    if (l_wtrac_conv) then
      do i_wt = 1, n_wtrac
        if (factor_env >  0.0) then
          wtrac_dd(i_wt)%rain_dd(i)  = wtrac_dd(i_wt)%rain_dd(i)               &
                                     + wtrac_ev(i_wt)%rain_env(i)*factor_env
          wtrac_ev(i_wt)%rain_env(i) = wtrac_ev(i_wt)%rain_env(i)              &
                                       * (1.0-factor_env)
          wtrac_dd(i_wt)%snow_dd(i)  = wtrac_dd(i_wt)%snow_dd(i)               &
                                     + wtrac_ev(i_wt)%snow_env(i)*factor_env
          wtrac_ev(i_wt)%snow_env(i) = wtrac_ev(i_wt)%snow_env(i)              &
                                       * (1.0-factor_env)
        else
          wtrac_ev(i_wt)%rain_env(i) = wtrac_ev(i_wt)%rain_env(i)              &
                                     - wtrac_dd(i_wt)%rain_dd(i)*factor_env
          wtrac_dd(i_wt)%rain_dd(i)  = wtrac_dd(i_wt)%rain_dd(i)               &
                                       * (1.0+factor_env)
          wtrac_ev(i_wt)%snow_env(i) = wtrac_ev(i_wt)%snow_env(i)              &
                                     - wtrac_dd(i_wt)%snow_dd(i)*factor_env
          wtrac_dd(i_wt)%snow_dd(i)  = wtrac_dd(i_wt)%snow_dd(i)               &
                                       * (1.0+factor_env)
        end if
      end do
    end if  ! l_wtrac_conv

  end if

end do   ! i

!-----------------------------------------------------------------------
! Compress out on basis of bit vector BDDON - those points with a
! downdraught
!-----------------------------------------------------------------------

nddon=0

do i=1,npnts
  if (bdd_on(i)) then
    nddon = nddon+1
    index1(nddon) = i
  end if
end do

if (nddon  >  0) then
  do i=1,nddon
    work(i,1) = thdd_k(index1(i))
    work(i,2) = qdd_k(index1(i))
    work(i,3) = the_k(index1(i))
    work(i,4) = the_km1(index1(i))
    work(i,5) = qe_k(index1(i))
    work(i,6) = qe_km1(index1(i))
    work(i,7) = dthbydt_k(index1(i))
    work(i,8) = dthbydt_km1(index1(i))
    work(i,9) = dqbydt_k(index1(i))
    work(i,10) = dqbydt_km1(index1(i))
    work(i,11) = flx_dd_k(index1(i))
    work(i,12) = p_km1(index1(i))
    work(i,13) = delpk(index1(i))
    work(i,14) = delpkm1(index1(i))
    work(i,15) = exk(index1(i))
    work(i,16) = exkm1(index1(i))
    work(i,17) = deltd(index1(i))
    work(i,18) = delqd(index1(i))
    work(i,19) = amdetk(index1(i))
    work(i,20) = ekm14(index1(i))
    work(i,21) = ekm34(index1(i))
    work(i,22) = rain_dd(index1(i))
    work(i,23) = snow_dd(index1(i))
    work(i,24) = cca(index1(i))
    qse_km1_c(i) = qse_km1(index1(i))

    bwork(i,1) = bdd_start(index1(i))
    bwork(i,2) = bddwt_k(index1(i))
    bwork(i,3) = bddwt_km1(index1(i))
    bwork(i,4) = bdd_on(index1(i))
    bwork(i,5) = b_dd_end(index1(i))
  end do

  do i=1,nddon

    ppn_mix_dd_c(i) = g*(rain_dd(index1(i)) +                                  &
                         snow_dd(index1(i)))/flx_dd_k(index1(i))
  end do

  if (l_tracer) then

    do ktra=1,ntra
      do i=1,nddon
        tradd_k_c(i,ktra) = tradd_k(index1(i),ktra)
        trae_k_c(i,ktra) = trae_k(index1(i),ktra)
        trae_km1_c(i,ktra) = trae_km1(index1(i),ktra)
        dtra_k_c(i,ktra)  = dtrabydt_k(index1(i),ktra)
        dtra_km1_c(i,ktra) = dtrabydt_km1(index1(i),ktra)
        deltrad_c(i,ktra) = deltrad(index1(i),ktra)
      end do
    end do

  end if

  if (l_wtrac_conv) then

    ! Allocate water tracer arrays
    ! Note, these are size nddon (not nddon_a) which is the required size
    ! for passing these arrays into ddraught
    call wtrac_alloc_conv_dd2(nddon, n_wtrac, wtrac_dd2)

    do i_wt = 1, n_wtrac
      do i=1,nddon
        wtrac_dd2(i_wt)%qdd_k(i)      = wtrac_dd(i_wt)%qdd_k(index1(i))
        wtrac_dd2(i_wt)%qe_k(i)       = wtrac_dd(i_wt)%qe_k(index1(i))
        wtrac_dd2(i_wt)%qe_km1(i)     = wtrac_dd(i_wt)%qe_km1(index1(i))
        wtrac_dd2(i_wt)%dqbydt_k(i)   = wtrac_dd(i_wt)%dqbydt_k(index1(i))
        wtrac_dd2(i_wt)%dqbydt_km1(i) = wtrac_ev(i_wt)%dqbydt_km1(index1(i))
        wtrac_dd2(i_wt)%delqd(i)      = wtrac_dd(i_wt)%delqd(index1(i))
        wtrac_dd2(i_wt)%rain(i)       = wtrac_dd(i_wt)%rain_dd(index1(i))
        wtrac_dd2(i_wt)%snow(i)       = wtrac_dd(i_wt)%snow_dd(index1(i))
      end do
    end do
  end if


  !-----------------------------------------------------------------------
  ! Start downdraught calculation
  !-----------------------------------------------------------------------


  call ddraught_6a (nddon,nddon_a,k,kct,ntra,n_wtrac,                          &
                  l_tracer,                                                    &
                  work(1,1),work(1,2),                                         &
                  work(1,3),work(1,4),work(1,5),work(1,6),                     &
                  qse_km1_c,                                                   &
                  work(1,7),work(1,8),work(1,9),work(1,10),                    &
                  work(1,11),work(1,12),work(1,13),work(1,14),                 &
                  work(1,15),work(1,16),work(1,17),work(1,18),                 &
                  work(1,19),work(1,20),work(1,21),work(1,22),                 &
                  work(1,23),bwork(1,1),bwork(1,2),bwork(1,3),                 &
                  bwork(1,4),bwork(1,5),                                       &
                  deltrad_c, work(1,24),ppn_mix_dd_c,                          &
                  tradd_k_c,trae_k_c,trae_km1_c,                               &
                  dtra_k_c,dtra_km1_c,wtrac_dd2)

  !-----------------------------------------------------------------------
  ! Expand requried vectors back to full fields
  !-----------------------------------------------------------------------

  do i=1,nddon
    thdd_k(index1(i)) = work(i,1)
    qdd_k(index1(i)) = work(i,2)
    dthbydt_k(index1(i)) = work(i,7)
    dthbydt_km1(index1(i)) = work(i,8)
    dqbydt_k(index1(i)) = work(i,9)
    dqbydt_km1(index1(i)) = work(i,10)
    flx_dd_k(index1(i)) = work(i,11)
    rain_dd(index1(i)) = work(i,22)
    snow_dd(index1(i)) = work(i,23)

    bdd_start(index1(i)) = bwork(i,1)
    bddwt_k(index1(i)) = bwork(i,2)
    bddwt_km1(index1(i)) = bwork(i,3)
    bdd_on(index1(i)) = bwork(i,4)
    b_dd_end(index1(i)) = bwork(i,5)
  end do

  if (l_tracer) then

    do ktra=1,ntra
      do i=1,nddon
        tradd_k(index1(i),ktra) = tradd_k_c(i,ktra)
        dtrabydt_k(index1(i),ktra) = dtra_k_c(i,ktra)
        dtrabydt_km1(index1(i),ktra) = dtra_km1_c(i,ktra)
      end do
    end do

  end if

  if (l_wtrac_conv) then
    do i_wt = 1, n_wtrac
      do i = 1, nddon
        wtrac_dd(i_wt)%qdd_k(index1(i))      = wtrac_dd2(i_wt)%qdd_k(i)
        wtrac_dd(i_wt)%dqbydt_k(index1(i))   = wtrac_dd2(i_wt)%dqbydt_k(i)
        wtrac_ev(i_wt)%dqbydt_km1(index1(i)) = wtrac_dd2(i_wt)%dqbydt_km1(i)
        wtrac_dd(i_wt)%rain_dd(index1(i))    = wtrac_dd2(i_wt)%rain(i)
        wtrac_dd(i_wt)%snow_dd(index1(i))    = wtrac_dd2(i_wt)%snow(i)
      end do
    end do

    ! Deallocate water tracer compressed arrays
    call wtrac_dealloc_conv_dd2(n_wtrac, wtrac_dd2)

  end if

end if   ! nddon > 0

!-----------------------------------------------------------------------
! Reset precipitation falling through environment if downdraught
! did not form
!-----------------------------------------------------------------------

do i=1,npnts
  if (.not. bdd_on(i) .and. .not. b_dd_end(i)) then
    rain_env(i) = rain_env(i)+rain_dd(i)
    snow_env(i) = snow_env(i)+snow_dd(i)
    rain_dd(i) = 0.0
    snow_dd(i) = 0.0
  end if
end do

if (l_wtrac_conv) then
  allocate(qe_km1_wtrac(npnts, n_wtrac))
  do i_wt = 1, n_wtrac
    do i=1,npnts
      if (.not. bdd_on(i) .and. .not. b_dd_end(i)) then
        wtrac_ev(i_wt)%rain_env(i) = wtrac_ev(i_wt)%rain_env(i)                &
                                   + wtrac_dd(i_wt)%rain_dd(i)
        wtrac_ev(i_wt)%snow_env(i) = wtrac_ev(i_wt)%snow_env(i)                &
                                   + wtrac_dd(i_wt)%snow_dd(i)
        wtrac_dd(i_wt)%rain_dd(i) = 0.0
        wtrac_dd(i_wt)%snow_dd(i) = 0.0
      end if
      ! Set environment value for call to precip_evap_wtrac
      qe_km1_wtrac(i,i_wt) = wtrac_dd(i_wt)%qe_km1(i)
    end do
  end do
end if

!-----------------------------------------------------------------------
! Carry out change of phase calculation for precipitation falling
! through environment
!-----------------------------------------------------------------------

if (l_wtrac_conv) then
  ! Store values before phase change for use in water tracer calculations
  call wtrac_alloc_conv_store2(npnts, qe_km1, rain_env, snow_env,              &
                               wtrac_conv_old)
end if

call chg_phse (npnts,rain_env,snow_env,dthbydt_km1,                            &
                  exk,exkm1,delpkm1,the_k,the_km1,timestep,cca,wtrac_conv_old)

if (l_wtrac_conv) then
  ! Update water tracers for phase change
  call wtrac_precip_chg_phse(npnts, n_wtrac, rain_env, snow_env,               &
                             wtrac_conv_old, wtrac_ev = wtrac_ev)
end if

!-----------------------------------------------------------------------
! Evaporate rain falling through environment if layer k below
! cloud base
!-----------------------------------------------------------------------

call pevp_bcb (npnts,k-1,iccb,the_km1,p_km1,qe_km1,qse_km1,delpkm1,            &
               rain_env,snow_env,dthbydt_km1,dqbydt_km1,                       &
               exkm1,timestep,cca)

if (l_wtrac_conv) then
  ! Update water tracers for phase change
  call wtrac_precip_evap (npnts, n_wtrac, delpkm1, qe_km1,                     &
                          rain_env, snow_env, the_km1, exkm1, qse_km1,         &
                          wtrac_conv_old, qe_km1_wtrac,                        &
                          timestep = timestep, wtrac_ev = wtrac_ev)
  deallocate(qe_km1_wtrac)
end if

!-----------------------------------------------------------------------
! Reset precipitation falling through environment if downdraught
! terminates
!-----------------------------------------------------------------------

do i=1,npnts
  if (b_dd_end(i)) then
    rain_env(i) = rain_env(i)+rain_dd(i)
    snow_env(i) = snow_env(i)+snow_dd(i)
    rain_dd(i) = 0.0
    snow_dd(i) = 0.0
  end if
end do

if (l_wtrac_conv) then
  do i_wt = 1, n_wtrac
    do i=1,npnts
      if (b_dd_end(i)) then
        wtrac_ev(i_wt)%rain_env(i) = wtrac_ev(i_wt)%rain_env(i)                &
                                   + wtrac_dd(i_wt)%rain_dd(i)
        wtrac_ev(i_wt)%snow_env(i) = wtrac_ev(i_wt)%snow_env(i)                &
                                   + wtrac_dd(i_wt)%snow_dd(i)
        wtrac_dd(i_wt)%rain_dd(i) = 0.0
        wtrac_dd(i_wt)%snow_dd(i) = 0.0
      end if
    end do
  end do
end if


!-----------------------------------------------------------------------
! Update rain and snow
!-----------------------------------------------------------------------

if (k == 2) then
  do i=1,npnts
    rain(i) = rain(i)+rain_dd(i)+rain_env(i)
    snow(i) = snow(i)+snow_dd(i)+snow_env(i)
  end do
end if

if (k == 2 .and. l_wtrac_conv) then
  do i_wt = 1, n_wtrac
    do i=1,npnts
      wtrac_dd(i_wt)%rain(i) = wtrac_dd(i_wt)%rain(i) +                        &
                   wtrac_dd(i_wt)%rain_dd(i) + wtrac_ev(i_wt)%rain_env(i)
      wtrac_dd(i_wt)%snow(i) = wtrac_dd(i_wt)%snow(i) +                        &
                   wtrac_dd(i_wt)%snow_dd(i) + wtrac_ev(i_wt)%snow_env(i)
    end do
  end do
end if ! l_wtrac_conv and k==2

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine downd_6a
end module downd_6a_mod
