! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module evap_bcb_nodd_all_6a_mod

use um_types, only: real_umphys

implicit none

! Description:
! Calculate convective precipitation reaching the surface
!
! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection

character(len=*), parameter, private :: ModuleName = 'EVAP_BCB_NODD_ALL_6A_MOD'
contains

subroutine evap_bcb_nodd_all_6a (npnts,n_nodd,n_wtrac,klev,kterm               &
                      ,iccb, index1, bwater                                    &
                      ,exner_layer_centres                                     &
                      ,p_layer_centres, p_layer_boundaries                     &
                      ,timestep , cca, the, qe, qse, precip                    &
                      ,dthbydt, dqbydt, rain, snow ,rain_3d, snow_3d           &
                      ,dt_dd, dq_dd, wtrac_p, wtrac_e)

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

use chg_phse_mod, only: chg_phse
use pevp_bcb_mod, only: pevp_bcb

use wtrac_conv_mod,            only: l_wtrac_conv, conv_e_wtrac_type,          &
                                     conv_p_wtrac_type, conv_ev_wtrac_type,    &
                                     wtrac_alloc_conv_ev, wtrac_dealloc_conv_ev
use wtrac_conv_store_mod,      only: conv_old_wtrac_type,                      &
                                     wtrac_alloc_conv_store2
use wtrac_precip_chg_phse_mod, only: wtrac_precip_chg_phse
use wtrac_precip_evap_mod,     only: wtrac_precip_evap

implicit none

!  Description : To calculate the convective precipitation reaching the
!               surface.
!
!     Method : the evaporation below cloud base follows that done in
!              the downdraught routine for the environmental part of the
!              column. the points which are gathered here are those
!              points which have an updraught, but no downdraught.
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards vn8.2.
!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------

integer, intent(in) ::                                                         &
  npnts                &  ! Vector length
 ,n_nodd               &  ! Compressed vector length for calculation
 ,n_wtrac              &  ! No. of water tracers
 ,klev                    ! Number of levels (may be model levels-1 or a
                          ! reduced set required here)

integer, intent(in) ::                                                         &
  kterm(npnts)         &  ! Convective cloud top layer
 ,iccb(npnts)          &  ! Cloud base level
 ,index1(npnts)           ! Index of points where downdraughts not possible

logical, intent(in) ::                                                         &
  bwater(npnts,2:klev+1)   ! Mask for points at which condensate is liquid

real(kind=real_umphys), intent(in) ::                                          &
  exner_layer_centres(npnts,0:klev+1)    & ! Exner function at layer centres
                                           ! starting at level k-1/2
 ,p_layer_centres(npnts,0:klev+1)        & ! pressure at layer centre (Pa)

 ,p_layer_boundaries(npnts,0:klev+1)       ! pressure at layer boundaries (Pa)

real(kind=real_umphys), intent(in) ::                                          &
  timestep               ! timestep

real(kind=real_umphys), intent(in) ::                                          &
  cca(npnts)           & ! 2d convective cloud amount
 ,the(npnts,klev+1)    & ! Model enviromental potential temperature (K)
 ,qe(npnts,klev+1)     & ! Model enviromental mixing ratio (kg/kg)
 ,qse(npnts,klev+1)      ! Model enviromental qsat mixing ratio (kg/kg)

real(kind=real_umphys), intent(in out) ::                                      &
  precip(npnts,klev+1)  & ! precipitation added when descending from layer
                          ! k to k-1 (kg/m**2/s)
 ,dthbydt(npnts,klev+1) & ! increment to model potential temperature (K/s)

 ,dqbydt(npnts,klev+1)  & ! increment to model mixing ratio (kg/kg/s)

 ,rain(npnts)           & ! rainfall at surface (kg/m**2/s)

 ,snow(npnts)           & ! snowfall at surface (kg/m**2/s)

 ,rain_3d(npnts,klev+1) & ! rainfall profile  (kg/m**2/s)

 ,snow_3d(npnts,klev+1) & ! snowfall profile  (kg/m**2/s)
 ,dt_dd(npnts,klev+1)   & ! dT/dt from evap below cloud base (K/s)
 ,dq_dd(npnts,klev+1)     ! dq/dt  from evap below cloud base  (kg/kg/s)

type(conv_p_wtrac_type), intent(in out) :: wtrac_p(n_wtrac)
                          ! Structure containing parcel water tracer fields

type(conv_e_wtrac_type), intent(in out) :: wtrac_e(n_wtrac)
                          ! Structure containing environment water tracer fields

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

integer ::                                                                     &
  i,k,i_wt        & ! Loop counters
 ,nsofar            ! points in level

integer ::                                                                     &
  index2(n_nodd)  & ! compress index for each level
 ,iccb_c(n_nodd)    ! Compressed cloud base level

real(kind=real_umphys) ::                                                      &
  pkm1_c(n_nodd)       & ! pressure of layer k-1 (Pa)
 ,exk_c(n_nodd)        & ! exner ratio for layer k
 ,exkm1_c(n_nodd)      & ! exner ratio for layer k-1
 ,delpkm1_c(n_nodd)      ! Pressure difference across layer K-1 (Pa)


real(kind=real_umphys) ::                                                      &
  qe_km1_c(n_nodd)     & ! Compressed parcel mixing ratio of layer K-1 (kg/kg)
 ,qse_km1_c(n_nodd)    & ! Compressed qsat environment of layer K-1 (kg/kg)
 ,the_k_c(n_nodd)      & ! Compressed parcel potential temperature of
                         ! layer k (K)
 ,the_km1_c(n_nodd)    & ! Compressed parcel potential temperature of
                         ! layer k-1 (K)
 ,dthbydt_km1_c(n_nodd)& ! Compressed increment to model potential temperature
                         ! of layer k-1 (K/s)
 ,dqbydt_km1_c(n_nodd) & ! Compressed increment to model mixing ratio of
                         ! layer k-1 (kg/kg/s)
 ,rain_c(n_nodd)       & ! Amount of rainfall passing through environment
                         ! (kg/m**2/s)
 ,snow_c(n_nodd)       & ! Amount of snowfall passing through environment
                         ! (kg/m**2/s)
 ,cca_c(n_nodd)          ! Compressed convective cloud amount

real(kind=real_umphys), allocatable :: qe_km1_c_wtrac(:,:)
                         ! Compressed water tracer mixing ratio of layer K-1
type(conv_ev_wtrac_type) :: wtrac_ev(n_wtrac)
                         ! Structure containing water tracer compressed
                         ! fields needed in phase change and evap of
                         ! precip in env or below cloud base calculations

type(conv_old_wtrac_type) :: wtrac_conv_old
                         ! Store of water values prior to phase change


integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='EVAP_BCB_NODD_ALL_6A'


!----------------------------------------------------------------------
! Loop over levels working downwards from maximum termination level + 1
!----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do k = klev+1,2,-1

  ! How many points terminated at or above this layer ?
  ! Only work on these points. For the bottom levels nsofar should
  ! equal n_nodd

  nsofar = 0
  do i=1,n_nodd
    if (kterm(index1(i))+1 >= k) then
      nsofar = nsofar + 1
      index2(nsofar) = index1(i)
    end if
  end do

  ! Compress to required points

  do i = 1, nsofar
    the_k_c(i)   = the(index2(i),k)
    the_km1_c(i) = the(index2(i),k-1)
    qe_km1_c(i)  = qe(index2(i),k-1)
    qse_km1_c(i)  = qse(index2(i),k-1)
    dthbydt_km1_c(i)  = dthbydt(index2(i),k-1)
    dqbydt_km1_c(i)   = dqbydt(index2(i),k-1)
    rain_c(i)  = rain(index2(i))
    snow_c(i)  = snow(index2(i))
    iccb_c(i)  = iccb(index2(i))
    cca_c(i)   = cca(index2(i))
    pkm1_c(i)    = p_layer_centres(index2(i),k-1)
    delpkm1_c(i) = p_layer_boundaries(index2(i),k-2) -                         &
                          p_layer_boundaries(index2(i),k-1)
    exk_c(i)   = exner_layer_centres(index2(i),k)
    exkm1_c(i) = exner_layer_centres(index2(i),k-1)
  end do

  if (l_wtrac_conv) then
    ! Allocate water tracer local arrays
    ! Note, these fields have size nsofar (not n_nodd) as is required for
    ! the routines precip_chg_phse_wtrac and pevp_bcb_wtrac

    call wtrac_alloc_conv_ev(nsofar, n_wtrac, wtrac_ev)
    allocate(qe_km1_c_wtrac(nsofar,n_wtrac))

    do i_wt = 1, n_wtrac
      do i = 1, nsofar
        qe_km1_c_wtrac(i,i_wt)       = wtrac_e(i_wt)%q(index2(i),k-1)
        wtrac_ev(i_wt)%dqbydt_km1(i) = wtrac_e(i_wt)%dqbydt(index2(i),k-1)
        ! Water tracer equivalent of rain_c
        wtrac_ev(i_wt)%rain_env(i)   = wtrac_e(i_wt)%rain(index2(i))
        ! Water tracer equivalent of snow_c
        wtrac_ev(i_wt)%snow_env(i)   = wtrac_e(i_wt)%snow(index2(i))
      end do
    end do
  end if

  do i = 1,nsofar
    if (bwater(index2(i),k)) then
      rain_c(i) = rain_c(i) + precip(index2(i),k)
    else
      snow_c(i) = snow_c(i) + precip(index2(i),k)
    end if
  end do

  if (l_wtrac_conv) then
    do i_wt = 1, n_wtrac
      do i = 1,nsofar
        if (bwater(index2(i),k)) then
          wtrac_ev(i_wt)%rain_env(i) =                                         &
             wtrac_ev(i_wt)%rain_env(i) + wtrac_p(i_wt)%precip(index2(i),k)
        else
          wtrac_ev(i_wt)%snow_env(i) =                                         &
             wtrac_ev(i_wt)%snow_env(i) + wtrac_p(i_wt)%precip(index2(i),k)
        end if
      end do
    end do
  end if

  !----------------------------------------------------------------------
  ! Carry out change of phase calculation for precipitation falling
  ! through environment
  !----------------------------------------------------------------------

  if (l_wtrac_conv) then
    ! Store values before phase change for use in water tracer calculations
    call wtrac_alloc_conv_store2(nsofar, qe_km1_c, rain_c, snow_c,             &
                                 wtrac_conv_old)
  end if

  call chg_phse (nsofar,rain_c,snow_c,dthbydt_km1_c,                           &
                 exk_c,exkm1_c,delpkm1_c,the_k_c,the_km1_c,                    &
                 timestep,cca_c,wtrac_conv_old=wtrac_conv_old)

  if (l_wtrac_conv) then
    ! Update water tracers for phase change
    call wtrac_precip_chg_phse(nsofar, n_wtrac, rain_c, snow_c,                &
                               wtrac_conv_old, wtrac_ev = wtrac_ev)
  end if

  !----------------------------------------------------------------------
  ! Reset precipitation falling through environment if downdraught
  ! terminates
  !----------------------------------------------------------------------

  call pevp_bcb (nsofar,k-1,iccb_c,the_km1_c,pkm1_c,qe_km1_c,qse_km1_c,        &
                 delpkm1_c,rain_c,snow_c,dthbydt_km1_c,                        &
                 dqbydt_km1_c,exkm1_c,timestep,cca_c)

  if (l_wtrac_conv) then
    ! Update water tracers for phase change
    call wtrac_precip_evap (nsofar, n_wtrac, delpkm1_c, qe_km1_c,              &
                            rain_c, snow_c, the_km1_c, exkm1_c, qse_km1_c,     &
                            wtrac_conv_old, qe_km1_c_wtrac,                    &
                            timestep = timestep, wtrac_ev = wtrac_ev)

    deallocate(qe_km1_c_wtrac)
  end if

  ! Expand output values back into full arrays

  do i=1,nsofar
    ! Want just downdraught increment rates
    dt_dd(index2(i),k-1) =  dt_dd(index2(i),k-1)+ (dthbydt_km1_c(i)            &
                   -dthbydt(index2(i),k-1))*exner_layer_centres(index2(i),k-1)
    dq_dd(index2(i),k-1) =  dq_dd(index2(i),k-1)+ (dqbydt_km1_c(i)             &
                                                   -dqbydt(index2(i),k-1))
    dthbydt(index2(i),k-1) = dthbydt_km1_c(i)
    dqbydt(index2(i),k-1)  = dqbydt_km1_c(i)

    ! Zero precipitation, as is (slyly) done in downd3c

    precip(index2(i),k) = 0.0
    rain(index2(i)) = rain_c(i)
    snow(index2(i)) = snow_c(i)
  end do

  if (l_wtrac_conv) then
    do i_wt = 1, n_wtrac
      do i=1, nsofar
        wtrac_e(i_wt)%dqbydt(index2(i),k-1) = wtrac_ev(i_wt)%dqbydt_km1(i)
        wtrac_p(i_wt)%precip(index2(i),k)   = 0.0
        wtrac_e(i_wt)%rain(index2(i))       = wtrac_ev(i_wt)%rain_env(i)
        wtrac_e(i_wt)%snow(index2(i))       = wtrac_ev(i_wt)%snow_env(i)
      end do
    end do
  end if

  ! Capture 3d rain/snow profiles
  do i=1,nsofar
    rain_3d(index2(i),k-1) = rain_3d(index2(i),k-1) + rain_c(i)
    snow_3d(index2(i),k-1) = snow_3d(index2(i),k-1) + snow_c(i)
  end do

  ! Deallocate water tracer local arrays
  if (l_wtrac_conv) then
    call wtrac_dealloc_conv_ev(n_wtrac, wtrac_ev)
  end if

end do      ! main loop over levels

!----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine evap_bcb_nodd_all_6a
end module evap_bcb_nodd_all_6a_mod
