! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module evap_bcb_nodd_6a_mod

use um_types, only: real_umphys

implicit none

! Description:
! Calculate convective precipitation reaching the surface
!
! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection

character(len=*), parameter, private :: ModuleName = 'EVAP_BCB_NODD_6A_MOD'
contains

subroutine evap_bcb_nodd_6a(npnts, n_nodd, n_wtrac, kct, iccb, index1,         &
                         bwater, timestep,                                     &
                         p_layer_centres, p_layer_boundaries,                  &
                         exner_layer_centres,                                  &
                         the, qe, qse, cca,                                    &
                         precip, dthbydt, dqbydt,                              &
                         rain, snow, rain_3d, snow_3d,                         &
                         dt_dd, dq_dd, wtrac_p, wtrac_e)

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

! Description : To calculate the convective precipitation reaching the
!               surface.
!
!     Method : The evaporation below cloud base follows that done in
!              the downdraught routine for the environmental part of the
!              column. The points which are gathered here are those
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
 ,n_wtrac              &  ! Number of water tracers
 ,kct                     ! Convective cloud top layer

integer, intent(in) ::                                                         &
  iccb(npnts)       & ! Cloud base level
 ,index1(npnts)          ! Index of downdraughts not possible

logical, intent(in) ::                                                         &
  bwater(npnts,2:kct+1)  ! mask for those points at which condensate is
                         ! water in layer k

real(kind=real_umphys), intent(in) ::                                          &
  timestep               ! Timestep

real(kind=real_umphys), intent(in) ::                                          &
  p_layer_centres(npnts,0:kct+1)      & ! pressure at layer centre (Pa)

 ,p_layer_boundaries(npnts,0:kct+1)   & ! pressure at layer boundaries (Pa)

 ,exner_layer_centres(npnts,0:kct+1)  & ! Exner function at layer centres
                                           ! starting at level k-1/2
 ,the(npnts,kct+1)                    & ! Model environmental potential
                                           ! temperature (K)
 ,qe(npnts,kct+1)                     & ! Environmental mixing ratio (kg/kg)
 ,qse(npnts,kct+1)                    & ! Environmental qsat (kg/kg)
 ,cca(npnts)                            ! Convective cloud amount

real(kind=real_umphys), intent(in out) ::                                      &
  precip(npnts,kct+1)     & ! Precipitation added when descending from layer
                            ! K to K-1(kg/m**2/s)
 ,dthbydt(npnts,kct+1) & ! in Increment to model potential temperature (K/s)
                            ! out Updated increment to model potential
                            !     temperature (K/s)
 ,dqbydt(npnts,kct+1)    ! in Increment to model mixing ratio (kg/kg/s)
                            ! out Updated increment to model
                            !     mixing ratio (kg/kg/s)

real(kind=real_umphys), intent(in out) ::                                      &
  rain(npnts)          & ! Rainfall at surface (kg/m**2/s)
 ,snow(npnts)          & ! Snowfall at surface (kg/m**2/s)
 ,rain_3d(npnts,kct+1) & ! Rainfall flux (kg/m**2/s)
 ,snow_3d(npnts,kct+1) & ! Snowfall flux  (kg/m**2/s)
 ,dt_dd(npnts,kct+1)   & ! dt/dt  from evap below cloud base  (K/s)
 ,dq_dd(npnts,kct+1)     ! dq/dt  from evap below cloud base  (kg/kg/s)

type(conv_p_wtrac_type), intent(in out) :: wtrac_p(n_wtrac)
                                        ! Structure containing parcel
                                        ! water tracer fields

type(conv_e_wtrac_type), intent(in out) :: wtrac_e(n_wtrac)
                                        ! Structure containing environment
                                        ! water tracer fields

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

integer ::                                                                     &
  i,k,i_wt          ! Loop counters

integer ::                                                                     &
  iccb_c(n_nodd)    ! Compressed cloud base level

real(kind=real_umphys) ::                                                      &
  pkm1_c(n_nodd)       & ! Pressure of layer k-1 (Pa)
 ,exk_c(n_nodd)        & ! Exner ratio for layer k
 ,exkm1_c(n_nodd)      & ! Exner ratio for layer k-1
 ,delpkm1_c(n_nodd)      ! Pressure difference across layer k-1 (Pa)

real(kind=real_umphys) ::                                                      &
  the_k_c(n_nodd)      & ! Compressed parcel potential temperature of layer k
                         !  (K)
 ,qe_km1_c(n_nodd)     & ! Compressed parcel mixing ratio of layer k-1 (kg/kg)

 ,qse_km1_c(n_nodd)    & ! Compressed qsat environment of layer k-1 (kg/kg)
 ,the_km1_c(n_nodd)    & ! Compressed parcel potential temperature of layer
                         !  k-1 (K)
 ,dthbydt_km1_c(n_nodd)& ! Compressed increment to model potential temperature
                         ! of layer k-1 (K/s)
 ,dqbydt_km1_c(n_nodd)   ! Compressed increment to model  mixing ratio
                         ! of layer k-1 (kg/kg/s)

real(kind=real_umphys) ::                                                      &
  rain_c(n_nodd)       & ! Amount of rainfall passing through environment
                         ! (kg/m**2/s)
 ,snow_c(n_nodd)       & ! Amount of ssnowfall passing through environment
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

character(len=*), parameter :: RoutineName='EVAP_BCB_NODD_6A'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Allocate water tracer arrays
if (l_wtrac_conv) then
  call wtrac_alloc_conv_ev(n_nodd, n_wtrac, wtrac_ev)
  allocate(qe_km1_c_wtrac(n_nodd,n_wtrac))
end if

!-----------------------------------------------------------------------
! Loop over levels moving downwards from cloud top
!-----------------------------------------------------------------------
do k = kct+1,2,-1

  do i=1,n_nodd
    the_k_c(i)   = the(index1(i),k)
    the_km1_c(i) = the(index1(i),k-1)
    qe_km1_c(i)  = qe(index1(i),k-1)
    qse_km1_c(i)  = qse(index1(i),k-1)
    dthbydt_km1_c(i) = dthbydt(index1(i),k-1)
    dqbydt_km1_c(i)  = dqbydt(index1(i),k-1)
  end do
  if (k == kct+1) then
    do i=1,n_nodd
      rain_c(i) = 0.0
      snow_c(i) = 0.0
      iccb_c(i) = iccb(index1(i))
      cca_c(i) = cca(index1(i))
    end do
  end if
  if (k == kct+1) then
    do i=1,n_nodd
      exk_c(i) = exner_layer_centres(index1(i),k)
    end do
  else
    do i=1,n_nodd
      exk_c(i) = exkm1_c(i)
    end do
  end if
  do i=1,n_nodd
    pkm1_c(i)    = p_layer_centres(index1(i),k-1)
    delpkm1_c(i) = p_layer_boundaries(index1(i),k-2) -                         &
                                         p_layer_boundaries(index1(i),k-1)
    exkm1_c(i)   = exner_layer_centres(index1(i),k-1)
  end do

  if (l_wtrac_conv) then
    do i_wt = 1, n_wtrac
      do i = 1, n_nodd
        qe_km1_c_wtrac(i,i_wt)       = wtrac_e(i_wt)%q(index1(i),k-1)
        wtrac_ev(i_wt)%dqbydt_km1(i) = wtrac_e(i_wt)%dqbydt(index1(i),k-1)
      end do
    end do
    if (k == kct+1) then
      do i_wt = 1, n_wtrac
        do i = 1, n_nodd
          wtrac_ev(i_wt)%rain_env(i) = 0.0   ! Equivalent of rain_c
          wtrac_ev(i_wt)%snow_env(i) = 0.0   ! Equivalent of snow_c
        end do
      end do
    end if
  end if

  do i=1,n_nodd
    if (bwater(index1(i),k)) then
      rain_c(i) = rain_c(i) + precip(index1(i),k)
    else
      snow_c(i) = snow_c(i) + precip(index1(i),k)
    end if
  end do

  if (l_wtrac_conv) then
    do i_wt = 1, n_wtrac
      do i=1,n_nodd
        if (bwater(index1(i),k)) then
          wtrac_ev(i_wt)%rain_env(i) =                                         &
             wtrac_ev(i_wt)%rain_env(i) + wtrac_p(i_wt)%precip(index1(i),k)
        else
          wtrac_ev(i_wt)%snow_env(i) =                                         &
             wtrac_ev(i_wt)%snow_env(i) + wtrac_p(i_wt)%precip(index1(i),k)
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
    call wtrac_alloc_conv_store2(n_nodd, qe_km1_c, rain_c, snow_c,             &
                                 wtrac_conv_old)
  end if

  call chg_phse (n_nodd,rain_c,snow_c,dthbydt_km1_c,                           &
                 exk_c,exkm1_c,delpkm1_c,the_k_c,the_km1_c,                    &
                 timestep,cca_c,wtrac_conv_old=wtrac_conv_old)

  if (l_wtrac_conv) then
    ! Update water tracers for phase change
    call wtrac_precip_chg_phse(n_nodd, n_wtrac, rain_c, snow_c,                &
                               wtrac_conv_old, wtrac_ev = wtrac_ev)
  end if

  !----------------------------------------------------------------------
  ! Reset precipitation falling through environment if downdraught
  ! terminates
  !----------------------------------------------------------------------

  call pevp_bcb (n_nodd,k-1,iccb_c,the_km1_c,pkm1_c,qe_km1_c,qse_km1_c,        &
                 delpkm1_c,rain_c,snow_c,dthbydt_km1_c,                        &
                 dqbydt_km1_c,exkm1_c,timestep,cca_c)

  if (l_wtrac_conv) then
    ! Update water tracers for phase change
    call wtrac_precip_evap (n_nodd, n_wtrac, delpkm1_c, qe_km1_c,              &
                            rain_c, snow_c, the_km1_c, exkm1_c, qse_km1_c,     &
                            wtrac_conv_old, qe_km1_c_wtrac,                    &
                            timestep = timestep, wtrac_ev = wtrac_ev)
  end if

  do i=1,n_nodd
    ! Want just downdraught increments
    dt_dd(index1(i),k-1) =  dt_dd(index1(i),k-1)+ (dthbydt_km1_c(i)            &
                 -dthbydt(index1(i),k-1))*exner_layer_centres(index1(i),k-1)
    dq_dd(index1(i),k-1) =  dq_dd(index1(i),k-1)+ (dqbydt_km1_c(i)             &
                                                 -dqbydt(index1(i),k-1))
    dthbydt(index1(i),k-1) = dthbydt_km1_c(i)
    dqbydt(index1(i),k-1)  = dqbydt_km1_c(i)

    ! Zero precipitation, as is (slyly) done in downd3c
    precip(index1(i),k) = 0.0
  end do

  if (l_wtrac_conv) then
    do i_wt = 1, n_wtrac
      do i=1, n_nodd
        wtrac_e(i_wt)%dqbydt(index1(i),k-1) = wtrac_ev(i_wt)%dqbydt_km1(i)
        wtrac_p(i_wt)%precip(index1(i),k)   = 0.0
      end do
    end do
  end if

  if (k == 2) then
    do i=1,n_nodd
      rain(index1(i)) = rain(index1(i)) + rain_c(i)
      snow(index1(i)) = snow(index1(i)) + snow_c(i)
    end do
  end if

  if (l_wtrac_conv) then
    if (k == 2) then
      do i_wt = 1, n_wtrac
        do i=1,n_nodd
          wtrac_e(i_wt)%rain(index1(i)) =                                      &
                 wtrac_e(i_wt)%rain(index1(i)) + wtrac_ev(i_wt)%rain_env(i)
          wtrac_e(i_wt)%snow(index1(i)) =                                      &
                 wtrac_e(i_wt)%snow(index1(i)) + wtrac_ev(i_wt)%snow_env(i)
        end do
      end do
    end if
  end if

  do i=1,n_nodd
    rain_3d(index1(i), k-1) = rain_3d(index1(i), k-1) + rain_c(i)
    snow_3d(index1(i), k-1) = snow_3d(index1(i), k-1) + snow_c(i)
  end do

end do      !  MAIN LOOP OVER LEVELS

! Deallocate water tracer local arrays
if (l_wtrac_conv) then
  deallocate(qe_km1_c_wtrac)
  call wtrac_dealloc_conv_ev(n_wtrac, wtrac_ev)
end if

!----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine evap_bcb_nodd_6a
end module evap_bcb_nodd_6a_mod
