! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate initial Downdraught massflux.
!
module dd_all_call_6a_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'DD_ALL_CALL_6A_MOD'
contains

subroutine dd_all_call_6a (npnts, npossdd, klev, nlev, trlev, ntra, n_wtrac    &
                       ,kterm, iccb, icct, index1                              &
                       ,l_tracer                                               &
                       ,bwater                                                 &
                       ,exner_layer_centres, exner_layer_boundaries            &
                       ,p_layer_centres, p_layer_boundaries, pstar             &
                       ,recip_pstar, timestep , cca                            &
                       ,thp, qp, the, qe, qse, trap, trae, flx                 &
                       ,precip, dthbydt, dqbydt, dtrabydt                      &
                       ,rain, snow, rain_3d, snow_3d                           &
                       ,wtrac_p, wtrac_e                                       &
                       ,dd_flux, entrain_dwn, detrain_dwn, dt_dd, dq_dd)


use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

use cv_stash_flg_mod, only:                                                    &
  flg_dwn_flx, flg_entr_dwn, flg_detr_dwn

use dd_init_6a_mod, only: dd_init_6a
use downd_6a_mod, only: downd_6a
use flx_init_6a_mod, only: flx_init_6a
use layer_dd_6a_mod, only: layer_dd_6a
use wtrac_conv_mod, only: l_wtrac_conv, conv_e_wtrac_type, conv_p_wtrac_type,  &
                          conv_dd_wtrac_type, conv_ev_wtrac_type,              &
                          wtrac_alloc_conv_dd, wtrac_dealloc_conv_dd,          &
                          wtrac_alloc_conv_ev, wtrac_dealloc_conv_ev


implicit none

!
! Description: Calculate initial Downdraught massflux.
!            Reset en/detrainment rates for Downdraught
!            Compress/expand variables
!            Initialise downdrought
!            Call downdraught routine
!
! Method:
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
  npnts                &  ! Vector length of full arrays
 ,npossdd              &  ! Length of array holding index of points where a
                          ! downdraught is possible.
 ,klev                 &  ! Number of levels (may be model levels or a reduced
                          ! set required here equal to kmax_term+1).
 ,nlev                 &  ! Number of model levels
 ,trlev                &  ! Number of tracer levels
 ,ntra                 &  ! Number of tracers
 ,n_wtrac                 ! Number of water tracers


integer, intent(in) ::                                                         &
  kterm(npnts)         &  ! Convective cloud top layer
 ,iccb(npnts)          &  ! Convective cloud base level (m)
 ,icct(npnts)          &  ! Convective cloud top level (m)
 ,index1(npnts)           ! index of points where downdraught possible

logical, intent(in) ::                                                         &
  l_tracer                ! Switch for tracers

logical, intent(in) ::                                                         &
  bwater(npnts,2:klev)  ! Mask for points at which condensate is liquid

real(kind=real_umphys), intent(in) ::                                          &
  exner_layer_centres(npnts,0:klev)       & ! exner pressure
 ,exner_layer_boundaries(npnts,0:klev)    & ! exner at half level above
                                            ! exner_layer_centres
 ,p_layer_centres(npnts,0:klev+1)         & ! Pressure (Pa)
 ,p_layer_boundaries(npnts,0:klev)        & ! Pressure at half level above
                                            ! p_layer_centres (Pa)
 ,pstar(npnts)                            & ! Surface pressure (Pa)
 ,recip_pstar(npnts)                        ! 1/pstar (Pa)

real(kind=real_umphys), intent(in) ::                                          &
  timestep                    ! timestep

real(kind=real_umphys), intent(in) ::                                          &
  cca(npnts)                & ! 2d convective cloud amount
 ,thp(npnts,klev)           & ! Parcel potential temperature (K)
 ,qp(npnts,klev)            & ! Parcel mixing ratio (kg/kg)
 ,the(npnts,klev)           & ! Model enviromental potential temperature (K)
 ,qe(npnts,klev)            & ! Model enviromental mixing ratio (kg/kg)
 ,qse(npnts,klev)           & ! Model enviromental qsat mixing ratio (kg/kg)
 ,trap(npnts,nlev,ntra)     & ! Parcel tracer (kg/kg)
 ,trae(npnts,trlev,ntra)    & ! Environment tracer (kg/kg)
 ,flx(npnts,klev)             ! updraught mass flux (Pa/s)

real(kind=real_umphys), intent(in out) ::                                      &
  precip(npnts,klev)        & ! precipitation added when descending
                              ! from layer k to k-1 (kg/m**2/s)
 ,dthbydt(npnts,klev)       & ! increment to model potential temperature (K/s)
 ,dqbydt(npnts,klev)        & ! increment to model mixing ratio (kg/kg/s)
 ,dtrabydt(npnts,nlev,ntra) & ! increment to model tracers (kg/kg/s)
 ,rain(npnts)               & ! rainfall at surface (kg/m**2/s)
 ,snow(npnts)               & ! snowfall at surface (kg/m**2/s)
 ,rain_3d(npnts,klev)       & ! rainfall flux  (kg/m**2/s)
 ,snow_3d(npnts,klev)       & ! snowfall flux  (kg/m**2/s)
 ,dt_dd(npnts,klev)         & ! dT/dt from Downdraught (K/s)
 ,dq_dd(npnts,klev)           ! dq/dt from Downdraught (kg/kg/s)

type(conv_p_wtrac_type), intent(in out) :: wtrac_p(n_wtrac)
                                          ! Structure containing parcel
                                          ! water tracer fields

type(conv_e_wtrac_type), intent(in out) :: wtrac_e(n_wtrac)
                                          ! Structure containing environment
                                          ! water tracer fields

real(kind=real_umphys), intent(out) ::                                         &
  entrain_dwn(npnts,klev)   & ! fractional entrainment rate for downdraught
                              ! mass flux
 ,detrain_dwn(npnts,klev)   & ! fractional detrainment rate for downdraught
                              ! mass flux
 ,dd_flux(npnts,klev)         ! Downdraught mass flux (Pa/s)

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

integer ::                                                                     &
  i,k ,i2,ktra,i_wt & ! Loop counters
 ,ndd               & ! Compressed vector length for downdraught calculation
 ,nddon_tmp           ! Number of points with active downdraught


logical ::                                                                     &
  bddi(npnts)         & ! Mask for points where downdraught might occur
 ,bdd_start(npnts)    & ! mask for those points where downdraught is able
                        ! to start from level k
 ,bddwt_k(npnts)      & ! mask for points in downdraught where ppt in
                        ! layer k is liquid
 ,bddwt_km1(npnts)    & ! mask for points in downdraught where ppt in
                        ! layer k-1 is liquid
 ,bdd_on(npnts)         ! mask for those points where DD continues from
                        ! layer k+1

real(kind=real_umphys) ::                                                      &
  flx_strt(npnts)       ! Mass flux at level where downdraught starts (Pa/s)

!-----------------------------------------------------------------------
! Compressed arrays
!-----------------------------------------------------------------------

integer ::                                                                     &
  iccb_c(npossdd) & ! Compressed cloud base level
 ,kmin(npossdd)     ! Freezing level where entrainment rates are increased

logical ::                                                                     &
  bwater_k_c(npossdd)  & ! Compressed mask for those points at which condensate
                         ! is water in layer k
 ,bddi_c(npossdd)      & ! Compressed mask for points where downdraught may
                         ! initiate
 ,bdd_start_c(npossdd) & ! Compressed mask for those points where downdraught
                         ! is able to start from level k
 ,bddwt_k_c(npossdd)   & ! Compressed mask for points in DD where ppt in
                         ! layer k is liquid
 ,bddwt_km1_c(npossdd) & ! Compressed mask for points in DD where ppt in
                         ! layer k-1 is liquid
 ,bdd_on_c(npossdd)      ! Compressed mask for points where DD continues from
                         ! layer K+1

real(kind=real_umphys) ::                                                      &
  exner_km12_c(npossdd) & ! Compressed exner function at layer k
 ,exner_kp12_c(npossdd) & ! Compressed exner function at layer k+1
 ,exner_KM32_c(npossdd) & ! Compressed exner function at layer k-1

 ,pk(npossdd)           & ! Pressure of layer k (Pa)
 ,p_km1(npossdd)        & ! Pressure of layer k-1 (Pa)
 ,exk(npossdd)          & ! exner ratio for layer K
 ,exkm1(npossdd)        & ! exner ratio for layer K-1

 ,delpk(npossdd)        & ! Pressure difference across layer K  (Pa)
 ,delpkm1(npossdd)      & ! Pressure difference across layer K-1 (Pa)

 ,amdetk(npossdd)       & ! Mixing detrainment at level k multiplied by
                          ! appropriate layer thickness
 ,ekm14(npossdd)        & ! exner ratio at layer k-1/4
 ,ekm34(npossdd)        & ! exner ratio at layer k-3/4

 ,precip_k_c(npossdd)   & ! Compressed precipitation added when descending
                          ! from layer K to K-1 (kg/m**2/s)
 ,q_k_c(npossdd)        & ! Compressed parcel mixing ratio of layer K (kg/kg)

 ,th_k_c(npossdd)       & ! Compressed parcel potential temperature of
                          ! layer k (K)
 ,tra_k_c(npossdd,ntra) & ! Compressed parcel tracer in layer K (kg/kg)

 ,pstar_c(npossdd)      & ! Compressed surface pressure (Pa)

 ,recip_pstar_c(npossdd)  ! Reciprocal of comp. pstar array

real(kind=real_umphys) ::                                                      &
  P_layer_centres_c(npossdd,0:klev+1)       & ! Pressure (Pa)
 ,P_layer_boundaries_c(npossdd,0:klev)      & ! Pressure (Pa)
 ,exner_layer_centres_c(npossdd,0:klev)       ! exner

real(kind=real_umphys) ::                                                      &
  dthbydt_K_c(npossdd)    & ! Compressed increment to model potential
                            ! temperature of layer k (K/s)
 ,dthbydt_km1_c(npossdd)  & ! Compressed increment to model potential
                            ! temperature of layer k-1 (K/s)
 ,dqbydt_k_c(npossdd)     & ! Compressed increment to model mixing ratio
                            ! of layer k (kg/kg/s)
 ,dqbydt_km1_c(npossdd)   & ! Compressed increment to model mixing ratio
                            ! of layer k-1(kg/kg/s)
 ,dtra_k_c(npossdd,ntra)  & ! Compressed increment to model tracer of
                            ! layer k (kg/kg/s)
 ,dtra_km1_c(npossdd,ntra)  ! Compressed increment to model tracer of
                            ! layer k-1 (kg/kg/s)

real(kind=real_umphys) ::                                                      &
  deltd(npossdd)          & ! Cooling necessary to achieve saturation (K)

 ,delqd(npossdd)          & ! moistening necessary to achieve saturation (K)

 ,deltrad(npossdd,ntra)   & ! Depletion of environment tracer due to
                            ! downdraught formation (kg/kg)
 ,qdd_k(npossdd)          & ! mixing ratio of downdraught in layer K (kg/kg)

 ,thdd_k(npossdd)         & ! Model potential temperature of downdraught
                            ! in layer k (K)
 ,tradd_k(npossdd,ntra)     ! Model tracer of downdraught in layer K (kg/kg)

real(kind=real_umphys) ::                                                      &
  flx_dd_k(npnts)          & ! Downdraught initial mass flux (Pa/s)
 ,flx_dd_k_c(npossdd)      & ! Compressed downdraught initial mass flux (Pa/s)
 ,qe_k_c(npossdd)          & ! Compressed environment mixing ratio of
                             ! layer k  (kg/kg)
 ,qe_km1_c(npossdd)        & ! Compressed environment mixing
                             ! ratio of layer k-1 (kg/kg)
 ,qse_k_c(npossdd)         & ! Compressed environment qsat mixing ratio of
                             ! layer k  (kg/kg)
 ,qse_km1_c(npossdd)       & ! Compressed environment qsat mixing
                             ! ratio of layer k-1 (kg/kg)
 ,the_k_c(npossdd)         & ! Compressed potential temperature
                             ! of environment in layer k (K)
 ,the_km1_c(npossdd)       & ! Compressed potential temperature
                             ! of environment in layer k-1 (K)
 ,trae_k_c(npossdd,ntra)   & ! Compressed tracer of environment in layer k
                             ! (kg/kg)
 ,trae_km1_c(npossdd,ntra) & ! Compressed tracer of environment in layer k-1
                             ! (kg/kg)
 ,rain_c(npossdd)          & ! Compressed surface rainfall (kg/m**2/s)

 ,snow_c(npossdd)          & ! Compressed surface snowfall (kg/m**2/s)

 ,flx_ud_k_c(npossdd)      & ! updraught mass flux at layer K

 ,rain_env(npossdd)        & ! Amount of rainfall passing through environment
                             ! (kg/m**2/s)
 ,snow_env(npossdd)        & ! Amount of snowfall passing through environment
                             ! (kg/m**2/s)
 ,rain_dd(npossdd)         & ! Amount of rainfall passing through
                             ! downdraught (kg/m**2/s)
 ,snow_dd(npossdd)         & ! Amount of snowfall passing through
                             ! downdraught (kg/m**2/s)
 ,flx_strt_c(npossdd)      & ! Compressed value of flx_strt

 ,cca_c(npossdd)           & ! Compressed convective cloud amount

 ,lr_ud_ref(npossdd)         ! precipitation mixing ratio at lowest
                             ! precipitationg level of UD

type(conv_dd_wtrac_type) :: wtrac_dd(n_wtrac)
                             ! Structure containing water tracer compressed
                             ! fields needed in downdraught calculations

type(conv_ev_wtrac_type) :: wtrac_ev(n_wtrac)
                             ! Structure containing water tracer compressed
                             ! fields needed in phase change and evap of
                             ! precip below cloud base calculations

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='DD_ALL_CALL_6A'


!-----------------------------------------------------------------------
! Compression to Downdraught points  (all levels even above term level)
!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ndd=npossdd

! Allocate water tracer arrays
! Note, these require size ndd (rather than npossdd) for passing to
! downd (although ndd = npossdd here)
if (l_wtrac_conv) then
  call wtrac_alloc_conv_dd(ndd, n_wtrac, wtrac_dd)
  call wtrac_alloc_conv_ev(ndd, n_wtrac, wtrac_ev)
end if

do k=0,klev
  do i=1,ndd
    p_layer_centres_c(i,k)     = p_layer_centres(index1(i),k)
    exner_layer_centres_c(i,k) = exner_layer_centres(index1(i),k)
    p_layer_boundaries_c(i,k)  = p_layer_boundaries(index1(i),k)
  end do
end do

!----------------------------------------------------------------------
! Initailise logical arrays as .false.
!-----------------------------------------------------------------------

do i=1,npnts
  bddi(i)      = .false.
  bdd_start(i) = .false.
  bddwt_k(i)   = .false.
  bddwt_km1(i) = .false.
  bdd_on(i)    = .false.
end do

do i2=1,ndd
  bddi(index1(i2)) = .true.    ! all points with downdraught
end do

!----------------------------------------------------------------------
! Calculate initial downdraught mass flux
! Done on all npnts but uses test on bddi (klev not used).
!-----------------------------------------------------------------------

call flx_init_6a(npnts, klev, iccb, icct, flx, flx_dd_k, bddi, flx_strt)

!-----------------------------------------------------------------------
! Compress all input arrays for the downdraught calculation down to
! those points where a downdraught is possible when the convection
! terminates in the column.
!-----------------------------------------------------------------------
! Initialise various arrays before level loop

do i=1,ndd
  flx_dd_k_c(i) = flx_dd_k(index1(i))
  flx_strt_c(i) = flx_strt(index1(i))
  pstar_c(i)    = pstar(index1(i))
  recip_pstar_c(i)=recip_pstar(index1(i))
  iccb_c(i) = iccb(index1(i))
  bdd_start_c(i) = bdd_start(index1(i))

  rain_c(i) = rain(index1(i))
  snow_c(i) = snow(index1(i))

  bddwt_k_c(i)   = bddwt_k(index1(i))
  bddwt_km1_c(i) = bddwt_km1(index1(i))
  bdd_on_c(i)    = bdd_on(index1(i))

  cca_c(i) = cca(index1(i))
  lr_ud_ref(i) = 0.0

  ! Initialise compressed downdraught indicator array to false

  bddi_c(i) = .false.

end do

if (l_wtrac_conv) then
  ! Initialise water tracer compressed fields
  do i_wt = 1, n_wtrac
    do i=1,ndd
      wtrac_dd(i_wt)%rain(i) = wtrac_e(i_wt)%rain(index1(i))
      wtrac_dd(i_wt)%snow(i) = wtrac_e(i_wt)%snow(index1(i))
    end do
  end do
end if


!-----------------------------------------------------------------------
! Main level loop working from top Downwards
!-----------------------------------------------------------------------

do k = klev,2,-1

  ! Need at this stage to reset bddi_c to true if reached level where
  ! convection terminated

  do i=1,ndd
    if (kterm(index1(i)) == k) then
      bddi_c(i) = .true.
    end if

    ! Compress arrays to those points with downdraughts possible in the column
    th_k_c(i) = thp(index1(i),k)
    q_k_c(i)  = qp(index1(i),k)
    the_k_c(i)   = the(index1(i),k)
    the_km1_c(i) = the(index1(i),k-1)
    qe_k_c(i)   = qe(index1(i),k)
    qe_km1_c(i) = qe(index1(i),k-1)
    qse_k_c(i)   = qse(index1(i),k)
    qse_km1_c(i) = qse(index1(i),k-1)
    dthbydt_k_c(i)   = dthbydt(index1(i),k)
    dthbydt_km1_c(i) = dthbydt(index1(i),k-1)
    dqbydt_k_c(i)   = dqbydt(index1(i),k)
    dqbydt_km1_c(i) = dqbydt(index1(i),k-1)
    exner_km12_c(i) = exner_layer_boundaries(index1(i),k-1)
    exner_kp12_c(i) = exner_layer_boundaries(index1(i),k)
    exner_km32_c(i) = exner_layer_boundaries(index1(i),k-2)
    precip_k_c(i) = precip(index1(i),k)
    flx_ud_k_c(i) = flx(index1(i),k)
    bwater_k_c(i) = bwater(index1(i),k)
  end do

  if (l_tracer) then   ! If run has tracers

    do ktra=1,ntra
      do i=1,ndd
        tra_k_c(i,ktra)    = trap(index1(i),k,ktra)
        trae_k_c(i,ktra)   = trae(index1(i),k,ktra)
        trae_km1_c(i,ktra) = trae(index1(i),k-1,ktra)
        dtra_k_c(i,ktra)   = dtrabydt(index1(i),k,ktra)
        dtra_km1_c(i,ktra) = dtrabydt(index1(i),k-1,ktra)
      end do
    end do

  end if

  if (l_wtrac_conv) then   ! If run has water tracers
    ! Compress water tracer arrays to those points with downdraughts
    ! possible in the column
    do i_wt = 1, n_wtrac
      do i=1,ndd
        wtrac_dd(i_wt)%q_k(i)        =  wtrac_p(i_wt)%q(index1(i),k)
        wtrac_dd(i_wt)%qe_k(i)       =  wtrac_e(i_wt)%q(index1(i),k)
        wtrac_dd(i_wt)%qe_km1(i)     =  wtrac_e(i_wt)%q(index1(i),k-1)
        wtrac_dd(i_wt)%dqbydt_k(i)   =  wtrac_e(i_wt)%dqbydt(index1(i),k)
        wtrac_ev(i_wt)%dqbydt_km1(i) =  wtrac_e(i_wt)%dqbydt(index1(i),k-1)
        wtrac_dd(i_wt)%precip_k(i)   =  wtrac_p(i_wt)%precip(index1(i),k)
      end do
    end do

  end if

  !----------------------------------------------------------------------
  ! If below convective cloud base downdraught not allowed to form
  !----------------------------------------------------------------------

  do i=1,ndd
    if (k <  iccb_c(i)) bddi_c(i)=.false.
  end do

  !-----------------------------------------------------------------------
  ! Reset en/detrainment rates for downdraught
  !-----------------------------------------------------------------------
  ! Possible problem with calculation of kmin in layer_dd when looping
  ! over levels starting well above the termination level.
  ! kmin only calculated in layer_dd if k=klev+1

  if ( k /= klev) then
    do i=1,ndd

      if (kterm(index1(i)) == k) then

        kmin(i)=klev+1      ! required to ensure amdetk calculated
                            ! correctly in layer_dd
      end if

    end do
  end if

  call layer_dd_6a (ndd,k,klev,the_k_c,the_km1_c,flx_strt_c,                   &
                    P_layer_centres_c,P_layer_boundaries_c,                    &
                    exner_layer_centres_c,                                     &
                    exner_km12_c,                                              &
                    pstar_c,pk,p_km1,delpk,delpkm1,exk,                        &
                    exkm1,amdetk,ekm14,ekm34,kmin,bddi_c,                      &
                    recip_pstar_c)

  !----------------------------------------------------------------------
  ! If level k within 150mb of surface then downdraught not allowed to
  ! form
  !----------------------------------------------------------------------

  do i=1,ndd
    if (pk(i) >  (pstar_c(i)-15000.0)) bddi_c(i)=.false.
  end do

  !-----------------------------------------------------------------------
  ! Initialise downdraught
  ! downdraught not allowed to form from cloud top layer (Klev+1)
  ! or from below cloud base
  !-----------------------------------------------------------------------

  call dd_init_6a (ndd,npossdd,n_wtrac,                                        &
                   th_k_c,q_k_c,the_k_c,qe_k_c,qse_k_c,pk,exk,                 &
                   thdd_k,qdd_k,deltd,delqd,                                   &
                   bdd_start_c,k,bddi_c,bdd_on_c,                              &
                   l_tracer,ntra,tra_k_c,                                      &
                   trae_k_c,tradd_k,deltrad,wtrac_dd)

  !-----------------------------------------------------------------------
  ! Update mask for where downdraught occurs
  !-----------------------------------------------------------------------

  do i=1,ndd
    if (bdd_start_c(i) .or. bdd_on_c(i)) bdd_on_c(i)=.true.
  end do

  !
  ! If downdraught initiated set diagnostic array
  !
  if (flg_dwn_flx) then
    do i=1,ndd
      if (bdd_start_c(i)) dd_flux(index1(i),k)=flx_dd_k(index1(i))
    end do
  end if

  nddon_tmp = 0
  do i=1,ndd
    if (bdd_on_c(i)) then
      nddon_tmp = nddon_tmp+1
    end if
  end do


  !-----------------------------------------------------------------------
  ! Call downdraught routine
  !-----------------------------------------------------------------------

  call downd_6a(ndd,npossdd,k,klev,ntra,n_wtrac,nddon_tmp, iccb_c,             &
                l_tracer,bwater_k_c,                                           &
                timestep,the_k_c,the_km1_c,                                    &
                qe_k_c,qe_km1_c,qse_km1_c,p_km1,delpk,delpkm1,exk,             &
                exkm1,deltd,delqd,amdetk,ekm14,ekm34,flx_ud_k_c,cca_c,         &
                trae_k_c,trae_km1_c,deltrad,                                   &
                bdd_start_c,bddwt_k_c,bddwt_km1_c,bdd_on_c,                    &
                thdd_k,qdd_k,dthbydt_k_c,dthbydt_km1_c,dqbydt_k_c,             &
                dqbydt_km1_c,rain_c,snow_c,precip_k_c, rain_env,snow_env,      &
                rain_dd,snow_dd,flx_dd_k_c,                                    &
                tradd_k,dtra_k_c,dtra_km1_c,wtrac_dd,wtrac_ev)

  !-----------------------------------------------------------------------
  ! Decompress/expand those variables which are to be output
  !-----------------------------------------------------------------------

  do i=1,ndd
    ! Want just downdraught tendencies
    dt_dd(index1(i),k)   =  dt_dd(index1(i),k)  + (dthbydt_k_c(i)              &
                       -dthbydt(index1(i),k))*exner_layer_centres(index1(i),k)
    dt_dd(index1(i),k-1) =  dt_dd(index1(i),k-1)+ (dthbydt_km1_c(i)            &
                    -dthbydt(index1(i),k-1))*exner_layer_centres(index1(i),k-1)
    dq_dd(index1(i),k)   =  dq_dd(index1(i),k)  + (dqbydt_k_c(i)               &
                                                   -dqbydt(index1(i),k))
    dq_dd(index1(i),k-1) =  dq_dd(index1(i),k-1)+ (dqbydt_km1_c(i)             &
                                                   -dqbydt(index1(i),k-1))
    dthbydt(index1(i),k)   = dthbydt_k_c(i)
    dthbydt(index1(i),k-1) = dthbydt_km1_c(i)
    dqbydt(index1(i),k)    = dqbydt_k_c(i)
    dqbydt(index1(i),k-1)  = dqbydt_km1_c(i)
  end do
  !
  ! Need to check that point would be selected in S.R DOWND or else
  ! not sensible to set entrainment and detrainment rates  in diagnostics
  !
  if (flg_dwn_flx) then
    do i=1,ndd
      if (bdd_on_c(i)) then
        dd_flux(index1(i),k-1) = flx_dd_k_c(i)
      end if
    end do
  end if

  if (flg_entr_dwn) then
    do i=1,ndd
      if (bdd_on_c(i)) then
        entrain_dwn(index1(i),k)=(1.0-amdetk(i))*                              &
                                 (ekm14(i)+ekm34(i)*(1.0+ekm14(i)))*           &
                                   dd_flux(index1(i),k)
      end if
    end do
  end if

  if (flg_detr_dwn) then
    do i=1,ndd
      if (bdd_on_c(i)) then
        detrain_dwn(index1(i),k)=-amdetk(i)*dd_flux(index1(i),k)
      end if
    end do
  end if

  if (k == 2) then
    do i=1,ndd
      rain(index1(i)) = rain_c(i)
      snow(index1(i)) = snow_c(i)
    end do
  end if


  do i=1, ndd
    precip(index1(i),k)    = precip_k_c(i)
    rain_3d(index1(i),k-1) = rain_3d(index1(i),k-1) + rain_dd(i) + rain_env(i)
    snow_3d(index1(i),k-1) = snow_3d(index1(i),k-1) + snow_dd(i) + snow_env(i)
  end do

  if (l_tracer) then       ! Runs with tracers

    do ktra=1,ntra
      do i=1,ndd
        dtrabydt(index1(i),k,ktra)   = dtra_k_c(i,ktra)
        dtrabydt(index1(i),k-1,ktra) = dtra_km1_c(i,ktra)
      end do
    end do

  end if

  if (l_wtrac_conv) then      ! Water tracers
    ! Decompress/expand water tracer variables which are to be output
    do i_wt = 1, n_wtrac
      do i = 1, ndd
        wtrac_e(i_wt)%dqbydt(index1(i),k)   = wtrac_dd(i_wt)%dqbydt_k(i)
        wtrac_e(i_wt)%dqbydt(index1(i),k-1) = wtrac_ev(i_wt)%dqbydt_km1(i)
        wtrac_p(i_wt)%precip(index1(i),k)   = wtrac_dd(i_wt)%precip_k(i)
        if (k == 2) then
          wtrac_e(i_wt)%rain(index1(i)) = wtrac_dd(i_wt)%rain(i)
          wtrac_e(i_wt)%snow(index1(i)) = wtrac_dd(i_wt)%snow(i)
        end if
      end do
    end do
  end if

  !----------------------------------------------------------------------
  !   End of main K loop
  !----------------------------------------------------------------------

end do    ! End of main level loop

! Deallocate water tracer arrays
if (l_wtrac_conv) then
  call wtrac_dealloc_conv_ev(n_wtrac, wtrac_ev)
  call wtrac_dealloc_conv_dd(n_wtrac, wtrac_dd)
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine dd_all_call_6a
end module dd_all_call_6a_mod
