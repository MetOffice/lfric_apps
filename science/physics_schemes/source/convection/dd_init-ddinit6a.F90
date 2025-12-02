! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module dd_init_6a_mod

use um_types, only: real_umphys

implicit none

! Description: Routine to initialise the downdraught
!
! Method: UM documentation paper 27
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards vn8.2.

character(len=*), parameter, private :: ModuleName = 'DD_INIT_6A_MOD'

contains

subroutine dd_init_6a(npnts, np_full, n_wtrac                                  &
                  ,th_ud_k, q_ud_k, the_k, qe_k, qse_k, pk, exk                &
                  ,thdd_k, qdd_k, deltd, delqd                                 &
                  ,bdd_start, k, bddi, bdd_on                                  &
                  ,l_tracer                                                    &
                  ,ntra, tra_ud_k, trae_k, tradd_k, deltrad, wtrac_dd)


use planet_constants_mod, only: c_virtual
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

use satcal_mod, only: satcal
use wtrac_conv_mod,       only: l_wtrac_conv, conv_dd_wtrac_type
use wtrac_calc_ratio_mod, only: wtrac_calc_ratio_fn

implicit none

!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------

integer, intent(in) ::                                                         &
  npnts                &  ! Vector length
 ,np_full              &  ! Full vector length
 ,ntra                 &  ! Number of tracers
 ,n_wtrac              &  ! Number of water tracers
 ,k                       ! Present model layer

logical, intent(in) ::                                                         &
  l_tracer                ! Switch for tracers

logical, intent(in) ::                                                         &
  bddi(npnts)          & ! Mask for points where downdraught may initiate
 ,bdd_on(npnts)          ! mask for those points where downdraught is on

real(kind=real_umphys), intent(in) ::                                          &
  th_ud_k(npnts)     & ! Parcel potential temperature of updraught, layer k (K)
 ,q_ud_k(npnts)      & ! Parcel mixing ratio of updraught, layer k (kg/kg)
 ,the_k(npnts)       & ! Potential temperature of environment in layer k (K)
 ,qe_k(npnts)        & ! Mixing ratio of environment in layer k (kg/kg)
 ,qse_k(npnts)       & ! Mixing ratio of environment qsat in layer k (kg/kg)
 ,pk(npnts)          & ! Pressure of layer k (Pa)
 ,exk(npnts)           ! Exner ratio layer k

real(kind=real_umphys), intent(in) ::                                          &
  trae_k(np_full,ntra)  & ! Tracer content of environment in layer k  (kg/kg)
 ,tra_ud_k(np_full,ntra)  ! Parcel tracer content of updraught, layer k (kg/kg)

logical, intent(in out) ::                                                     &
  bdd_start(npnts)           ! input mask for those points where DD may start ?
                             ! output set to true if DD started on level k

real(kind=real_umphys), intent(out) ::                                         &
  thdd_k(npnts)          & ! Downdraught potential temperature of layer k (K)
 ,qdd_k(npnts)           & ! Downdraught mixing ratio of layer k (kg/kg)
 ,tradd_k(np_full,ntra)  & ! Downdraught tracer content of layer k (kg/kg)
 ,deltd(npnts)           & ! Cooling necessary to achieve saturation
 ,delqd(npnts)           & ! Moistening necessary to achieve saturation
 ,deltrad(np_full,ntra)    ! Depletion of environment tracer due to formation
                           ! of Downdraught

type(conv_dd_wtrac_type), intent(in out) :: wtrac_dd(n_wtrac)
                                         ! Structure containing the water
                                         ! tracer downdraught fields

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

integer ::                                                                     &
  i, ktra, i_wt         ! Loop counters


real(kind=real_umphys) ::                                                      &
  th_mean(npnts)      & ! Mean potential temperature used in calculation of
                        ! saturated downdraught potential temperature  in
                        ! layer k (K)
 ,q_mean(npnts)       & ! Mean mixing ratio used in calculation of
                        ! saturated downdraught mixing ratio in layer k (kg/kg)
 ,t_mean(npnts)       & ! Mean temperature used in calculation of
                        ! saturated downdraught potential temperature in
                        ! layer k (kg/kg)
 ,thdds(npnts)        & ! Saturated downdraught potential temperature in
                        ! layer k (K)
 ,qdds(npnts)         & ! Saturated downdraught mixing ratio in
                        ! layer k (kg/kg)
 ,buoy(npnts)           ! Buoyancy of parcel in layer k

real(kind=real_umphys) ::                                                      &
  thdd_v            &  ! Virtual potential temperature of parcel in layer k (K)
 ,the_v                ! Virtual potential temperature of environment in
                       ! layer k (K)
real(kind=real_umphys) ::                                                      &
  q_mean_wtrac      &  ! Mean mixing ratio of water tracer used in calc of
                       ! saturated downdraught mixing ratio in layer k (kg/kg)
 ,ratio_wtrac          ! Ratio of mean water tracer mix ratio to mean
                       ! normal water mix ratio in layer k

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='DD_INIT_6A'


!-----------------------------------------------------------------------
! Calculate mean temperature, mixing ratio, U, V  and tracer.
!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do i=1,npnts
  th_mean(i) = (the_k(i)+th_ud_k(i))*0.5
  q_mean(i)  = (qe_k(i)+q_ud_k(i))*0.5
  t_mean(i)  = th_mean(i)*exk(i)
end do

!-----------------------------------------------------------------------
! Calculate saturated downdraught potential temperature for layer k
!-----------------------------------------------------------------------

call satcal(npnts,th_mean,pk,exk,q_mean,the_k,qse_k,qdds,thdds)

!-----------------------------------------------------------------------
! Is saturated parcel negatively buoyant compared to environment?
!-----------------------------------------------------------------------

do i=1,npnts
  if (.not. bdd_on(i) .and. bddi(i) ) then
    thdd_v = thdds(i)*(1.0+c_virtual*qdds(i))
    the_v  = the_k(i)*(1.0+c_virtual*qe_k(i))
    buoy(i) = thdd_v - the_v

    if (buoy(i)  <   0.5 ) then

      !-----------------------------------------------------------------------
      ! Initiate downdraught
      !-----------------------------------------------------------------------

      thdd_k(i) = thdds(i)
      qdd_k(i) = qdds(i)
      bdd_start(i) = .true.

      !-----------------------------------------------------------------------
      ! Calculate cooling and moistening to achieve saturation
      !-----------------------------------------------------------------------

      deltd(i) = thdds(i)-the_k(i)
      delqd(i) = qdds(i)-qe_k(i)
    end if
  end if
end do


if (l_tracer .and. k >= 4) then

  do ktra=1,ntra
    do i=1,npnts
      if (.not. bdd_on(i) .and. bddi(i)) then
        if (buoy(i) <  0.5) then
          tradd_k(i,ktra) = (trae_k(i,ktra)+tra_ud_k(i,ktra))*0.5
          deltrad(i,ktra) = tradd_k(i,ktra)-trae_k(i,ktra)
        end if
      end if
    end do
  end do

end if

if (l_wtrac_conv) then
  do i_wt = 1, n_wtrac
    do i=1,npnts
      if (.not. bdd_on(i) .and. bddi(i) ) then
        if (buoy(i)  <   0.5 ) then

          ! Calculate ratio of water tracer to water for the 'mean' values
          q_mean_wtrac = (wtrac_dd(i_wt)%qe_k(i)+wtrac_dd(i_wt)%q_k(i)) * 0.5
          ratio_wtrac = wtrac_calc_ratio_fn(i_wt,q_mean_wtrac,q_mean(i))

          wtrac_dd(i_wt)%qdd_k(i) = ratio_wtrac * qdds(i)
          wtrac_dd(i_wt)%delqd(i) = (ratio_wtrac * qdds(i))                    &
                                       - wtrac_dd(i_wt)%qe_k(i)

        end if
      end if
    end do
  end do
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine dd_init_6a

end module dd_init_6a_mod
