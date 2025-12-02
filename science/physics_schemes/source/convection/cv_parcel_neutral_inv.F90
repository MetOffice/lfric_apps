! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate level of neutral buoyancy, CAPE and CIN for ascent
!
module cv_parcel_neutral_inv_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'CV_PARCEL_NEUTRAL_INV_MOD'
contains

subroutine cv_parcel_neutral_inv(nunstable,                                    &
               nlcl_c,k_plume,                                                 &
               z_lcl_c, thv_pert, z_full_c, z_half_c, exner_theta_levels_c,    &
               buoyancy, t_dens_env, dqsatdz, denv_bydz, dpar_bydz,            &
               zh_c, zh2,                                                      &
               k_max,k_neutral,k_inv, kshmin, shmin,                           &
               max_buoy,dt_dens_parc_t,ql_ad_c,delthvu_c,cape_c,cin_c,         &
               dt_dens_parc_t2,dt_dens_parc_tmin,ql_ad2,                       &
               delthvu2)


use planet_constants_mod, only: g => g_bl
use bl_option_mod, only: zero, one
use parkind1, only: jpim, jprb       !DrHook
use yomhook,  only: lhook, dr_hook   !DrHook
use nlsizes_namelist_mod,  only: model_levels

implicit none
!
! Description:
!   Calculate the level of neutral buoyancy for a parcel ascent.
!   Also calcuate the CAPE and CIN of the ascent
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

integer, intent(in) ::                                                         &
  nunstable             ! Number of parcel ascents

integer, intent(in) ::                                                         &
  nlcl_c(nunstable)    & ! Level number of LCL
 ,k_plume(nunstable)     ! start level for surface-driven plume

real(kind=r_bl), intent(in)    ::                                              &
  z_lcl_c(nunstable)     & ! LCL height rounded to nearest model level
 ,thv_pert(nunstable)      ! threshold thv of parcel  (K)

real(kind=r_bl), intent(in)    ::                                              &
  z_full_c(nunstable, model_levels)            & ! Height theta lev (m)
 ,z_half_c(nunstable, model_levels)            & ! Height uv lev    (m)
 ,exner_theta_levels_c(nunstable, model_levels)& ! Exner on theta lev
 ,buoyancy(nunstable, model_levels)            & ! undilute parcel buoyancy (K)
 ,t_dens_env(nunstable, model_levels)          & ! Density potential temperature
                                                 ! of environment. (K)
 ,dqsatdz(nunstable, model_levels)             & ! dqsat/dz along adiabat
 ,denv_bydz(nunstable, model_levels)           & ! Gradient of density potential
                                                 ! temp in the environment.
 ,dpar_bydz(nunstable, model_levels)             ! Gradient of density potential
                                                 ! temperature of the parcel.


real(kind=r_bl), intent(in out)    ::                                          &
  zh_c(nunstable)         & ! BL depth compressed
 ,zh2(nunstable)            !  2nd copy

integer, intent(out) ::                                                        &
  k_max(nunstable)      & ! level of max parcel buoyancy
 ,k_neutral(nunstable)  & ! level of neutral parcel buoyancy
 ,k_inv(nunstable)      & ! level from inversion testing
 ,kshmin(nunstable)       ! Position of buoyancy minimum above
                          ! topbl (used for shallow Cu diag)level

logical,  intent(out) ::                                                       &
  shmin(nunstable)         ! Flag for finding min in parcel buoyancy below
                           ! 3km (for shallow Cu)
real(kind=r_bl), intent(out)    ::                                             &
  max_buoy(nunstable)      ! Maximum buoyancy

real(kind=r_bl), intent(out)   ::                                              &
  dt_dens_parc_T(nunstable)    & ! t_dens_parc-t_dens_env at ntpar
 ,dt_dens_parc_Tmin(nunstable) & ! t_dens_parc-t_dens_env at kshmin
 ,cape_c(nunstable)            & ! CAPE from undilute parcel ascent (m2/s2)
 ,cin_c(nunstable)             & ! CIN from undilute parcel ascent (m2/s2)
 ,delthvu_c(nunstable)         & ! Integral of undilute parcel buoyancy
                                 ! over convective cloud layer (for convection)
 ,ql_ad_c(nunstable)           & ! adiabatic liquid water content at inversion
                                 ! or cloud top (kg/kg)
 ,dt_dens_parc_T2(nunstable)   & ! t_dens_parc-t_dens_env at ntpar
 ,delthvu2(nunstable)          & ! 2nd copy delthuv
 ,ql_ad2(nunstable)              ! ql_ad 2nd copy

!-------------------------------------------------------------------------
! Local variables

integer :: ii, k            !Loop counters


real(kind=r_bl)    ::                                                          &
  inc      & ! CAPE in layer
 ,dz         ! layer thickness

real(kind=r_bl)    ::                                                          &
  dtv_min(nunstable)      & ! min Tv of parcel in cld layer 1 virtual
                            ! temperature (K).
 ,dtv_min2(nunstable)       ! 2nd copy min TV of parcel

logical ::                                                                     &
  topbl(nunstable)        & ! Flag set when top of boundary layer is reached.
 ,topprof(nunstable)      & ! Flag set when top of ascent is reached.
 ,above_lcl(nunstable)      ! Flag set when parcel above LCL.

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CV_PARCEL_NEUTRAL_INV'

!-------------------------------------------------------------------------
! 1.0 Start of subroutine code: perform the calculation.

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Initialise output arrays

do ii=1,nunstable
  shmin(ii)   = .false.
  topbl(ii)   = .false.
  topprof(ii) = .false.

  kshmin(ii)     = 1
  k_max(ii)      = 1
  k_neutral(ii)  = 1
  k_inv(ii)      = 1

  dtv_min(ii)  = zero
  dtv_min2(ii) = zero
  delthvu2(ii) = zero
  max_buoy(ii)     = zero
  delthvu_c(ii) = zero
  ql_ad2(ii)   = zero
  ql_ad_c(ii)  = zero
  CAPE_c(ii)  = zero
  CIN_c(ii)  = zero
end do

do  k = 2,model_levels

  do ii=1,nunstable


    ! Set flag to true when level below is at least one level above the lcl
    ! and above the lcl transition zone
    ! Code implies ABOVE_LCL at NLCL+3 or greater.

    if (k-1 >  nlcl_c(ii)+1                                                    &
                     .and. z_full_c(ii,k-1) >  1.1_r_bl*z_lcl_c(ii)) then
      above_lcl(ii)=.true.
    else
      above_lcl(ii)=.false.
    end if

    !-----------------------------------------------------------------------
    ! Level of neutral buoyancy (LNB) & maximum buoyancy level below this
    !-----------------------------------------------------------------------
    ! Not reached LNB continue testing

    if ( .not. topprof(ii) .and. k >  k_plume(ii) ) then

      if (buoyancy(ii,k) >  max_buoy(ii)) then
        max_buoy(ii) = buoyancy(ii,k)
        k_max(ii)    = k
      end if

      ! Is parcel still buoyant ?

      if ( (buoyancy(ii,k)  <=  - thv_pert(ii))                                &
        !                      or reached top of model
                 .or. (k  >   model_levels-1)  ) then

        k_neutral(ii) = k-1
        topprof(ii) = .true.
        zh_c(ii) = z_half_c(ii,k)

        ! Buoyancy at last buoyant level

        Dt_dens_parc_T(ii) = buoyancy(ii,k-1)

        if ( delthvu_c(ii)  >   zero) then
          ! compensate for any negative buoyancy of parcel in cloud layer
          delthvu_c(ii) = delthvu_c(ii) - dtv_min(ii) *                        &
                                    ( z_half_c(ii,k) - z_lcl_c(ii) )
        end if
      end if
    end if

    !-----------------------------------------------------------------------
    ! Tests applied once found top of parcel ascent.
    ! Aim - to establish if the ascent has an inversion above the top
    !       i.e. the ascent may indicate shallow /congestus convection.
    ! Sets indicator shmin = .true. if conditions met and stops testing.
    !
    ! Conditions are ;
    ! either  denv/dz(k) < dpar/dz(k)
    !   or    par_svl(k-1) -env_svl(k-1) <= 0.0
    !
    !-----------------------------------------------------------------------

    if ( topbl(ii) .and. ( denv_bydz(ii,k)  <   dpar_bydz(ii,k)                &
                    .or. buoyancy(ii,k-1)  <=  zero )                          &
                                  .and. .not. shmin(ii) ) then
      shmin(ii) = .true.
      dt_dens_parc_tmin(ii) = buoyancy(ii,k-1)
      kshmin(ii) = k-1
    end if

    !-----------------------------------------------------------------------
    ! Tests applied to find parcel top
    !-----------------------------------------------------------------------

    if ( .not. topbl(ii) .and. k  >   k_plume(ii) .and.                        &
        (  ( buoyancy(ii,k) <=  - thv_pert(ii)) .or.                           &

      !           plume non buoyant

          (above_lcl(ii) .and. (denv_bydz(ii,k) >  1.25_r_bl*dpar_bydz(ii,k))) &

      !           or environmental virtual temperature gradient
      !           signIficantly larger than parcel gradient
      !           above lIfting condensation level

                 .or. (k  >   model_levels-1)                                  &
      !                      or reached top of model
               )                                                               &
               ) then

      topbl(ii) = .true.
      zh2(ii) = z_half_c(ii,k)
      k_inv(ii) = k-1

      dt_dens_parc_T2(ii) = buoyancy(ii,k-1)
      if ( delthvu2(ii)  >   zero) then
        ! compensate for any negative buoyancy of parcel in cloud layer
        delthvu2(ii) = delthvu2(ii) - dtv_min2(ii) *                           &
                                  ( z_half_c(ii,k) - z_lcl_c(ii) )
      end if
    end if          ! test on .not.topbl

    !-----------------------------------------------------------------------
    ! While doing parcel ascent
    ! (a) find minimum buoyancy
    ! (b) integrate CAPE over the ascent
    !-----------------------------------------------------------------------

    if (k > nlcl_c(ii) .and. k < model_levels ) then

      dz = z_half_c(ii,k+1) - z_half_c(ii,k)
      inc = g *  buoyancy(ii,k) * dz/t_dens_env(ii,k)

      !----------------------------------------------------------
      ! If not reached an inversion or level of neutral buoyancy
      !----------------------------------------------------------

      if (.not. topbl(ii)) then


        ! adiabatic liquid water content at cloud top or inversion
        !                            = -dqsat/dz * zcld

        ql_ad2(ii) = -one* dqsatdz(ii,k)                                       &
                        *(z_half_c(ii,k+1) - z_half_c(ii,nlcl_c(ii)+1))

        dtv_min2(ii) = min( dtv_min2(ii),                                      &
                           buoyancy(ii,k)/exner_theta_levels_c(ii,k) )

        delthvu2(ii) = delthvu2(ii) + buoyancy(ii,k)*dz                        &
                                           /exner_theta_levels_c(ii,k)
      end if

      !----------------------------------------------------------
      ! If not reached level of neutral buoyancy (top of ascent)
      !----------------------------------------------------------

      if (.not. topprof(ii)) then

        ! Note only calculating CIN and CAPE from ascents reaching
        ! level of neutral buoyancy. This may not always correspond
        ! to the diagnosed top for the convection scheme.

        if (inc <  zero) then
          CIN_c(ii)  = CIN_c(ii) + inc
        else      ! CAPE holds only postive part
          CAPE_c(ii) = CAPE_c(ii) + inc
        end if
        ! adiabatic liquid water content at cloud top
        !                            = -dqsat/dz * zcld

        ql_ad_c(ii) = -one* dqsatdz(ii,k)                                      &
                     *(z_half_c(ii,k+1) - z_half_c(ii,nlcl_c(ii)+1))

        dtv_min(ii) = min( dtv_min(ii),                                        &
                         buoyancy(ii,k)/exner_theta_levels_c(ii,k)  )

        delthvu_c(ii) = delthvu_c(ii) + buoyancy(ii,k)* dz                     &
                                          /exner_theta_levels_c(ii,k)

      end if    ! test on topprof

    end if

  end do   ! ii loop

end do     ! level loop


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine cv_parcel_neutral_inv
end module cv_parcel_neutral_inv_mod
