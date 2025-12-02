! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************


!   Compress expand around Emanuel Downdraught code

module eman_cex_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'EMAN_CEX_MOD'
contains

subroutine eman_cex(npnts,nmid,kmax_term,nlev,trlev,ntra                       &
,                      kterm,l_mid,l_tracer                                    &
,                      exner_layer_centres,exner_layer_boundaries              &
,                      p,ph,timestep,scale_factor                              &
,                      th,q,qse,tracer,precip                                  &
,                      dthbydt, dqbydt, dtrabydt                               &
,                      rain, snow ,dwn_flux, dt_dd, dq_dd                      &
                   )

use cv_run_mod, only:                                                          &
  l_snow_rain

! Subroutines
use eman_dd_rev_mod, only: eman_dd_rev

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

use eman_dd_mod, only: eman_dd
implicit none

!  Description : Interface to Emanuel downdraught routine
!
!  Method :
!   First compress to points with mid-level convection.
!   Then call  Emanuel Down draughts.
!   Expand results back to full grid.
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
  npnts                & ! No. of deep convection points
, nmid                 & ! No. of mid level convection points
, nlev                 & ! No. of model layers
, trlev                & ! No. of tracer levels
, ntra                 & ! No. of tracer fields
, kmax_term            & ! highest level reached by convection
, kterm(npnts)           ! level reached by convection


logical, intent(in) ::                                                         &
  l_tracer             & ! true in tracers present
, l_mid(npnts)           ! true  if mid level convection

real(kind=real_umphys), intent(in)    ::                                       &
  exner_layer_centres(npnts,0:nlev)    & ! Exner
, exner_layer_boundaries(npnts,0:nlev) & ! Exner at half level above
, p(npnts,0:nlev)                      & ! Pressure  (Pa)
, ph(npnts,0:nlev)                     & ! Pressure at half level (Pa)
, timestep                             & ! timestep for model physics in seconds
, scale_factor(npnts)                  & ! factor to scale sigd by
, qse(npnts,nlev)                      & ! Saturation specific humidity of cloud
                                         ! environment (kg/kg)
, q (npnts,nlev)                       & ! specific humidity of cloud
                                         ! environment(kg/kg)
, th (npnts,nlev)                      & ! theta of cloud environment(K)
, precip(npnts,nlev)                   & ! Precip from updraught G-R scheme
                                         ! (kg/m2/s) not units for wdtrain
, tracer(npnts,trlev,ntra)               !  Tracer on model levels  (kg/kg)


real(kind=real_umphys), intent(in out)    ::                                   &
  rain(npnts)             & ! rainfall at surface (kg/m**2/s)
, snow(npnts)               ! snowfall at surface (kg/m**2/s)

! increments
real(kind=real_umphys), intent(in out)    ::                                   &
  dqbydt(npnts,nlev)        & ! increments to q (kg/kg/s)
, dthbydt(npnts,nlev)       & ! increments to potential temperature(K/s)
, dtrabydt(npnts,nlev,ntra)   ! increments to model tracers(kg/kg/s)

real(kind=real_umphys), intent(in out)    ::                                   &
  dt_dd(npnts,nlev)         & ! dT/dt from convective DD (K/s)
 ,dq_dd(npnts,nlev)           ! dq/dt from convective DD (kg/kg/s)

! Arguments with intent out:

real(kind=real_umphys), intent(out)    ::                                      &
  dwn_flux(npnts,nlev)       ! Downdraught mass flux (Pa/s) diagnostic

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
! compressed version of most input arrays

integer ::                                                                     &
  kterm_c(nmid)      ! compressed termination level

real(kind=real_umphys) ::                                                      &
  rain_c(nmid)                            & ! compressed rain
, snow_c(nmid)                                                                 &
, scale_factor_c(nmid)                                                         &
, dwn_flux_c(nmid,nlev)                                                        &
, dqbydt_c(nmid,nlev)                                                          &
, dthbydt_c(nmid,nlev)                                                         &
, dtrabydt_c(nmid,nlev,ntra)                                                   &
, exner_layer_centres_c(nmid,0:nlev)                                           &
, exner_layer_boundaries_c(nmid,0:nlev)                                        &
, p_c(nmid,0:nlev)                                                             &
, ph_c(nmid,0:nlev)                                                            &
, qse_c(nmid,nlev)                                                             &
, q_c(nmid,nlev)                                                               &
, th_c(nmid,nlev)                                                              &
, dt_dd_c(nmid,nlev)                                                           &
, dq_dd_c(nmid,nlev)                                                           &
, precip_c(nmid,nlev)                                                          &
, tracer_c(nmid,trlev,ntra)

integer ::                                                                     &
 dmi(nmid)      ! index of convecting points

integer ::                                                                     &
 i ,j ,k        ! loop counters

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='EMAN_CEX'

! ----------------------------------------------------------------------
! Compress input arrays
! ----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise output arrays
do k = 1, nlev
  do i = 1, nmid
    dt_dd_c(i,k)    = 0.0
    dq_dd_c(i,k)    = 0.0
  end do
end do

j=0
do i = 1,npnts
  if (l_mid(i)) then
    j = j+1
    dmi(j) = i
  end if
end do

do j = 1, nmid
  kterm_c(j) = kterm(dmi(j))
  rain_c(j)  = rain(dmi(j))
  snow_c(j)  = snow(dmi(j))
  scale_factor_c(j)  = scale_factor(dmi(j))

end do

do k = 0, nlev
  do j = 1, nmid
    exner_layer_centres_c(j,k)   = exner_layer_centres(dmi(j),k)
    exner_layer_boundaries_c(j,k)=exner_layer_boundaries(dmi(j),k)
    p_c(j,k)  = p(dmi(j),k)
    ph_c(j,k) = ph(dmi(j),k)
  end do
end do

do k = 1, nlev
  do j = 1, nmid
    q_c(j,k)   = q(dmi(j),k)
    th_c(j,k)  = th(dmi(j),k)
    qse_c(j,k) = qse(dmi(j),k)

    precip_c(j,k) = precip(dmi(j),k)

    ! Problems as taking scaling from last i.e. highest convection in the
    ! column. If this has failed but below there is some convection then need
    ! to ensure scale_factor is not 0.0
    if (precip_c(j,k) > 0.0 .and. scale_factor_c(j) == 0.0 ) then
      scale_factor_c(j) = 1.0
    end if

    dqbydt_c(j,k)  = dqbydt(dmi(j),k)
    dthbydt_c(j,k) = dthbydt(dmi(j),k)

  end do
end do




if (l_tracer) then
  do i = 1, ntra
    do k = 1, trlev
      do j = 1, nmid
        tracer_c(j,k,i) = tracer(dmi(j),k,i)
      end do
    end do
    do k = 1, nlev
      do j = 1, nmid
        dtrabydt_c(j,k,i) = dtrabydt(dmi(j),k,i)
      end do
    end do
  end do
end if

! ----------------------------------------------------------------------
! Call Emanuel scheme for just mid level convecting points
! ----------------------------------------------------------------------
if (l_snow_rain) then
  call eman_dd_rev(nmid,kmax_term,nlev,trlev,ntra,kterm_c,l_tracer             &
,                      exner_layer_centres_c                                   &
,                      p_c, ph_c, scale_factor_c                               &
,                      th_c, q_c, qse_c, tracer_c, precip_c                    &
,                      dthbydt_c, dqbydt_c, dtrabydt_c                         &
,                      rain_c, snow_c ,dwn_flux_c,dt_dd_c,dq_dd_c  )

else
  call eman_dd (nmid,kmax_term,nlev,trlev,ntra,kterm_c,l_tracer                &
,               exner_layer_centres_c                                          &
,               p_c,ph_c,timestep,th_c,q_c,qse_c,tracer_c                      &
,               precip_c,dthbydt_c,dqbydt_c,dtrabydt_c                         &
,               rain_c, snow_c ,dwn_flux_c,dt_dd_c,dq_dd_c )
end if

! ----------------------------------------------------------------------
! Expand results
! ----------------------------------------------------------------------
do i = 1,nmid
  rain(dmi(i)) = rain_c(i)
  snow(dmi(i)) = snow_c(i)
end do
do k = 1, nlev
  do i = 1, nmid
    dwn_flux(dmi(i),k) = dwn_flux_c(i,k)
    dqbydt(dmi(i),k)   = dqbydt_c(i,k)
    dthbydt(dmi(i),k)  = dthbydt_c(i,k)
    dt_dd(dmi(i),k)    = dt_dd_c(i,k)
    dq_dd(dmi(i),k)    = dq_dd_c(i,k)
  end do
end do

if (l_tracer) then
  do j = 1, ntra
    do k = 1, nlev
      do i = 1, nmid
        dtrabydt(dmi(i),k,j) = dtrabydt_c(i,k,j)
      end do
    end do
  end do
end if
! ----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine eman_cex

end module eman_cex_mod
