! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module dd_env_6a_mod

use um_types, only: real_umphys

implicit none

! Description: Downdraught routine
!              Calculate the effect of the downdraught on the
!              large-scale atmosphere
!
! Method: UM documentation paper 27
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards vn8.2.

character(len=*), parameter, private :: ModuleName = 'DD_ENV_6A_MOD'

contains

subroutine dd_env_6a(npnts, np_full, ntra, n_wtrac                             &
                  ,l_tracer, b_dd_end, bdd_start, bdd_on                       &
                  ,thdd_k, thdd_km1, qdd_k, qdd_km1, the_k, the_km1            &
                  ,qe_k, qe_km1, flx_dd_k, flx_dd_km1, delpk, delpkm1          &
                  ,deltd, delqd, amdetk, ekm14                                 &
                  ,tradd_k, tradd_km1, trae_k, trae_km1, deltrad               &
                  ,qdd_km1_wtrac                                               &
                  ,dthbydt_k, dthbydt_km1, dqbydt_k, dqbydt_km1                &
                  ,dtrabydt_k, dtrabydt_km1, wtrac_dd2)

use wtrac_conv_mod, only: l_wtrac_conv, conv_dd2_wtrac_type

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------
integer, intent(in) ::                                                         &
  npnts                &  ! Vector length
 ,np_full              &  ! Full vector length
 ,ntra                 &  ! Number of tracers
 ,n_wtrac                 ! Number of water tracers

logical, intent(in) ::                                                         &
  l_tracer                ! Switch for tracers

logical, intent(in) ::                                                         &
  b_dd_end(npnts)      & ! Mask for those points where downdraught is
                         ! terminating
 ,bdd_start(npnts)     & ! Mask for those points where downdraught is
                         ! starting
 ,bdd_on(npnts)          ! Mask for those points where downdraught is on

real(kind=real_umphys), intent(in) ::                                          &
  thdd_k(npnts)     & ! Potential temperature of downdraught in layer k (K)
 ,thdd_km1(npnts)   & ! Potential temperature of downdraught in layer k-1 (K)
 ,qdd_k(npnts)      & ! Mixing ratio of downdraught in layer k (kg/kg)
 ,qdd_km1(npnts)    & ! Mixing ratio of downdraught in layer k-1 (kg/kg)
 ,the_k(npnts)      & ! Potential temperature of environment in layer k (K)
 ,the_km1(npnts)    & ! Potential temperature of environment in layer k-1 (K)
 ,qe_k(npnts)       & ! Mixing ratio of environment in layer k  (kg/kg)
 ,qe_km1(npnts)     & ! Mixing ratio of environment in layer k-1 (kg/kg)
 ,flx_dd_k(npnts)   & ! Downdraught mass flux of layer k (Pa/s)
 ,flx_dd_km1(npnts) & ! Downdraught mass flux of layer k-1 (Pa/s)
 ,delpk(npnts)      & ! Change in pressure across layer k (Pa)
 ,delpkm1(npnts)    & ! Change in pressure across layer k-1  (Pa)
 ,deltd(npnts)      & ! Cooling necessary to achieve saturation (K)
 ,delqd(npnts)      & ! moistening necessary to achieve saturation (kg/kg)
 ,amdetk(npnts)     & ! Mixing detrainment rate
 ,ekm14(npnts)        ! Exner ratio at layer k-1/4

real(kind=real_umphys), intent(in) ::                                          &
  tradd_k(np_full,ntra)   & ! Downdraught tracer content of layer k (kg/kg)
 ,tradd_km1(npnts,ntra)   & ! Downdraught tracer content of layer k-1 (kg/kg)
 ,trae_k(np_full,ntra)    & ! Environment tracer content of layer k (kg/kg)
 ,trae_km1(np_full,ntra)  & ! Environment tracer content of layer k-1(kg/kg)
 ,deltrad(npnts,ntra)       ! Depletion of environment tracer due to
                            ! downdraught formation (kg/kg)

real(kind=real_umphys), intent(in) ::                                          &
  qdd_km1_wtrac(npnts,n_wtrac) ! Water tracer mixing ratio of downdraught
                               ! in layer k-1 (kg/kg)

real(kind=real_umphys), intent(in out) ::                                      &
  dthbydt_k(npnts)     & ! In  Increment to potential temperature of layer k
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
                               ! (kg/kg/s)
 ,dtrabydt_km1(np_full,ntra)   ! In  Increment to tracer content of layer k-1
                               ! Out Updated increment tracer content layer k-1
                               ! (kg/kg/s)

type(conv_dd2_wtrac_type), intent(in out) :: wtrac_dd2(n_wtrac)
                             ! Structure containing 2nd compression water
                             ! tracer arrays used in downdraught calculations

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

integer ::                                                                     &
  i, ktra, i_wt    ! Loop counters

real(kind=real_umphys) ::                                                      &
  tempry           ! Used in calculations of the effect on the environment

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='DD_ENV_6A'

!-----------------------------------------------------------------------
! Calculate the effect on the environment in layer k
!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do i=1,npnts
  if (bdd_on(i)) then

    !-----------------------------------------------------------------------
    ! Subtract the energy used to form the downdraught
    !-----------------------------------------------------------------------

    tempry = flx_dd_k(i)/delpk(i)
    if (bdd_start(i)) then
      dthbydt_k(i) = dthbydt_k(i)-tempry*deltd(i)
      dqbydt_k(i)  = dqbydt_k(i) -tempry*delqd(i)
    end if

    !-----------------------------------------------------------------------
    ! Effect of convection and downdarught upon potential temperature of
    ! layer k
    !-----------------------------------------------------------------------

    dthbydt_k(i) = dthbydt_k(i) + tempry * (                                   &

                        ! compensating subsidence term
                  (1.0+ekm14(i)) * (1.0-amdetk(i)) * (the_km1(i)-the_k(i))     &

                        ! Mixing detrainment term
                  +  amdetk(i)* (thdd_k(i)-the_k(i)) )

    !-----------------------------------------------------------------------
    ! Effect of convection and downdraught upon mixing ratio of
    ! layer k
    !-----------------------------------------------------------------------

    dqbydt_k(i) = dqbydt_k(i) + tempry * (                                     &

                        ! compensating subsidence term
                  (1.0+ekm14(i)) * (1.0-amdetk(i)) * (qe_km1(i)-qe_k(i))       &

                        ! Mixing detrainment term
                  +  amdetk(i)* (qdd_k(i)-qe_k(i)) )


    !-----------------------------------------------------------------------
    ! Terminal detrainment and subsidence in terminal layer
    !-----------------------------------------------------------------------

    if (b_dd_end(i)) then
      tempry         = flx_dd_km1(i)/delpkm1(i)

      dthbydt_km1(i) = dthbydt_km1(i)+tempry* (thdd_km1(i)-the_km1(i))

      dqbydt_km1(i)  = dqbydt_km1(i)+tempry*(qdd_km1(i)-qe_km1(i))
    end if

  end if    ! test on bdd_on
end do      ! loop over npnts

!-----------------------------------------------------------------------
! Effect of convection and downdraught upon tracer content of  layer k
!-----------------------------------------------------------------------

if (l_tracer) then

  do ktra=1,ntra
    do i=1,npnts
      if (bdd_on(i)) then

        tempry = flx_dd_k(i)/delpk(i)
        if (bdd_start(i)) then
          dtrabydt_k(i,ktra) = dtrabydt_k(i,ktra)-tempry*deltrad(i,ktra)
        end if
        dtrabydt_k(i,ktra) = dtrabydt_k(i,ktra) + tempry * (                   &

                        ! compensating subsidence term
         (1.0+ekm14(i)) * (1.0-amdetk(i)) * (trae_km1(i,ktra)-trae_k(i,ktra))  &

                        ! Mixing detrainment term
            + amdetk(i)* (tradd_k(i,ktra)-trae_k(i,ktra))  )

        !--------------------------------------------------------------------
        ! Terminal Detrainment of tracer
        !--------------------------------------------------------------------

        if (b_dd_end(i)) then
          tempry = flx_dd_km1(i)/delpkm1(i)
          dtrabydt_km1(i,ktra)=dtrabydt_km1(i,ktra)+tempry*                    &
                                   (tradd_km1(i,ktra)-trae_km1(i,ktra))
        end if

      end if
    end do
  end do

end if      ! l_tracer

!--------------------------------------------------------------------------
! Effect of convection and downdraught upon WATER tracer content of layer k
!--------------------------------------------------------------------------

if (l_wtrac_conv) then
  do i_wt = 1, n_wtrac
    do i=1,npnts
      if (bdd_on(i)) then

        !----------------------------------------------------------------
        ! Subtract the water tracer content of the moisture used to form
        !  the downdraught
        !----------------------------------------------------------------

        tempry = flx_dd_k(i)/delpk(i)
        if (bdd_start(i)) then
          wtrac_dd2(i_wt)%dqbydt_k(i)  =  wtrac_dd2(i_wt)%dqbydt_k(i)          &
                                       - tempry*wtrac_dd2(i_wt)%delqd(i)
        end if

        !----------------------------------------------------------------
        ! Effect of convection and downdraught upon mixing ratio of
        ! layer k
        !----------------------------------------------------------------

        wtrac_dd2(i_wt)%dqbydt_k(i) = wtrac_dd2(i_wt)%dqbydt_k(i)              &
           + tempry * (                                                        &

                 ! compensating subsidence term
           (1.0+ekm14(i)) * (1.0-amdetk(i)) *                                  &
           (wtrac_dd2(i_wt)%qe_km1(i)-wtrac_dd2(i_wt)%qe_k(i))                 &

                 ! Mixing detrainment term
           +  amdetk(i)* (wtrac_dd2(i_wt)%qdd_k(i)-wtrac_dd2(i_wt)%qe_k(i)) )


        !----------------------------------------------------------------
        ! Terminal detrainment and subsidence in terminal layer
        !----------------------------------------------------------------

        if (b_dd_end(i)) then
          tempry         = flx_dd_km1(i)/delpkm1(i)
          wtrac_dd2(i_wt)%dqbydt_km1(i) = wtrac_dd2(i_wt)%dqbydt_km1(i)        &
               + tempry*(qdd_km1_wtrac(i,i_wt)-wtrac_dd2(i_wt)%qe_km1(i))
        end if

      end if   ! test on bdd_on
    end do     ! loop over npnts
  end do       ! loop over water tracers
end if         ! l_wtrac_conv

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine dd_env_6a

end module dd_env_6a_mod
