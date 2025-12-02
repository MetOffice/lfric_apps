! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate level of neutral buoyancy for a dilute and an undilute ascent.
!
module cv_parcel_neutral_dil_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'CV_PARCEL_NEUTRAL_DIL_MOD'
contains

subroutine cv_parcel_neutral_dil(nunstable,                                    &
               nlcl_c,k_plume, l_dilute,                                       &
               z_lcl_c, thv_pert, z_full_c, z_half_c, exner_theta_levels_c,    &
               buoyancy, buoyancy_dil, t_dens_env, dqsatdz,                    &
               zh_c,                                                           &
               k_max,k_max_dil,k_neutral,k_neutral_dil,                        &
               max_buoy,max_buoy_dil,ql_ad_c,delthvu_c,                        &
               cape_c,cin_c)

use planet_constants_mod, only: g => g_bl
use bl_option_mod, only: zero, one
use parkind1, only: jpim, jprb       !DrHook
use yomhook,  only: lhook, dr_hook   !DrHook
use nlsizes_namelist_mod,  only: model_levels

implicit none
!
! Description:
!   Calculate the level of neutral buoyancy for a dilute and an undilute parcel
!   ascent.
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
  nunstable              ! Number of parcel ascents

integer, intent(in) ::                                                         &
  nlcl_c(nunstable)    & ! Level number of LCL
 ,k_plume(nunstable)     ! start level for surface-driven plume

logical, intent(in) ::                                                         &
  l_dilute               ! true if dilute ascent as well as undilutre ascent

real(kind=r_bl), intent(in)    ::                                              &
  z_lcl_c(nunstable)     & ! LCL height rounded to nearest model level
 ,thv_pert(nunstable)      ! threshold thv of parcel  (K)

real(kind=r_bl), intent(in)    ::                                              &
  z_full_c(nunstable, model_levels)            & ! Height theta lev (m)
 ,z_half_c(nunstable, model_levels)            & ! Height uv lev    (m)
 ,exner_theta_levels_c(nunstable, model_levels)& ! Exner on theta lev
 ,buoyancy(nunstable, model_levels)            & ! undilute parcel buoyancy (K)
 ,buoyancy_dil(nunstable, model_levels)        & ! dilute parcel buoyancy (K)
 ,t_dens_env(nunstable, model_levels)          & ! Density potential temperature
                                                 ! of environment. (K)
 ,dqsatdz(nunstable, model_levels)               ! dqsat/dz along adiabat


real(kind=r_bl), intent(in out)    ::                                          &
  zh_c(nunstable)           ! BL depth compressed

integer, intent(out) ::                                                        &
  k_max(nunstable)          & ! level of max parcel buoyancy
 ,k_max_dil(nunstable)      & ! level of max parcel buoyancy dilute
 ,k_neutral(nunstable)      & ! level of neutral parcel buoyancy for undilute
 ,k_neutral_dil(nunstable)    ! level of neutral parcel buoyancy for dilute

real(kind=r_bl), intent(out)    ::                                             &
  max_buoy(nunstable)       & ! Maximum buoyancy undilute ascent
 ,max_buoy_dil(nunstable)     ! Maximum buoyancy dilute ascent

real(kind=r_bl), intent(out)   ::                                              &
  cape_c(nunstable)         & ! CAPE from undilute parcel ascent (m2/s2)
 ,cin_c(nunstable)          & ! CIN from undilute parcel ascent (m2/s2)
 ,delthvu_c(nunstable)      & ! Integral of undilute parcel buoyancy
                              ! over convective cloud layer (for convection)
 ,ql_ad_c(nunstable)          ! adiabatic liquid water content at inversion
                              ! or cloud top (kg/kg)

!-------------------------------------------------------------------------
! Local variables

integer :: ii, k            !Loop counters


real(kind=r_bl)    ::                                                          &
  inc      & ! CAPE in layer
 ,dz       & ! layer thickness
 ,factor     ! multiplies thv_pert

real(kind=r_bl)    ::                                                          &
  dtv_min(nunstable)      & ! min Tv of parcel in cld layer 1 virtual
                            ! temperature (K).
 ,dtv_min_dil(nunstable)    ! min Tv of parcel in cld layer dilute ascent but
                            ! undilute value
! UNUSED AT present - left in case required in future
!real    ::                &
! ,cape2_c(nunstable)      & ! undilute CAPE from dilute parcel ascent (m2/s2)
! ,cin2_c(nunstable)         ! undilute CIN from dilute parcel ascent (m2/s2)

logical ::                                                                     &
  topprof(nunstable)      & ! Flag set when top of ascent is reached.
 ,topprof_dil(nunstable)    ! Flag set when top of ascent is reached.



integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CV_PARCEL_NEUTRAL_DIL'

!-------------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

if (l_dilute) then

!$OMP  PARALLEL do DEFAULT(none) private(ii, k, factor, dz, inc)               &
!$OMP  SHARED(nunstable,k_max,k_neutral,max_buoy,ql_ad_c,delthvu_c,            &
!$OMP         CAPE_c,CIN_c,dtv_min,topprof,k_max_dil,k_neutral_dil,            &
!$OMP         max_buoy_dil,dtv_min_dil,topprof_dil,model_levels,               &
!$OMP         nlcl_c,k_plume,z_lcl_c,thv_pert,  z_full_c,z_half_c,             &
!$OMP         exner_theta_levels_c,buoyancy, buoyancy_dil,t_dens_env,          &
!$OMP         dqsatdz,zh_c, g) SCHEDULE(STATIC)
  ! The main index of the loop is "ii" so that it can be
  ! parallelised using OpenMP.
  do ii=1,nunstable

    !-----------------------------------------------------------------------
    ! Initialise output arrays

    topprof(ii)     = .false.

    k_max(ii)          = 1
    k_neutral(ii)      = 1

    dtv_min(ii)   = zero
    max_buoy(ii)  = zero
    delthvu_c(ii) = zero
    ql_ad_c(ii)   = zero
    CAPE_c(ii)    = zero
    CIN_c(ii)     = zero

    ! initialise more arrays

    topprof_dil(ii) = .false.

    k_max_dil(ii)      = 1
    k_neutral_dil(ii)  = 1

    dtv_min_dil(ii)   = zero
    max_buoy_dil(ii)  = zero
    !    CAPE2_c(ii)    = 0.0
    !    CIN2_c(ii)     = 0.0

    do  k = 2,model_levels

      !-----------------------------------------------------------------------
      ! Only perform tests if parcel ascent If unstable
      !-----------------------------------------------------------------------
      ! No flag for above_lcl required. Reduce thv_pert by a factor dependent
      ! on height relative to LCL.

      if (k-1 >  nlcl_c(ii)+1                                                  &
                        .and. z_full_c(ii,k-1) >  1.1_r_bl*z_lcl_c(ii)) then

        if ((z_full_c(ii,k)-z_lcl_c(ii)) >   1000.0_r_bl) then
          ! set to zero if z-zlcl >1000?
          factor =zero
        else
          ! decrease thv_pert by exp(-(z-zlcl)/1000.)
          factor = exp( (z_lcl_c(ii)-z_full_c(ii,k))*1.0e-3_r_bl)
        end if

      else
        factor= one
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

        if ( (buoyancy(ii,k)  <=  - thv_pert(ii))                              &
          !                      or reached top of model
                   .or. (k  >   model_levels-1)  ) then

          k_neutral(ii) = k-1
          topprof(ii) = .true.

        end if
      end if
      !-----------------------------------------------------------------------
      ! Level of neutral buoyancy (LNB) & maximum buoyancy level below this
      ! Dilute plume
      !-----------------------------------------------------------------------
      ! Not reached LNB continue testing

      if ( .not. topprof_dil(ii) .and. k >  k_plume(ii) ) then
        if (buoyancy_dil(ii,k) >  max_buoy(ii)) then
          max_buoy_dil(ii) = buoyancy_dil(ii,k)
          k_max_dil(ii)    = k
        end if

        ! Is parcel still buoyant ?

        if ( (buoyancy_dil(ii,k)  <=  - thv_pert(ii)*factor)                   &
          !                      or reached top of model
                   .or. (k  >   model_levels-1)  ) then

          k_neutral_dil(ii) = k-1
          topprof_dil(ii) = .true.
          zh_c(ii) = z_half_c(ii,k)

          if ( delthvu_c(ii)  >   zero) then
            ! Note delthvu_c is an undilute value calculated to top of dilute
            ! plume.
            ! compensate for any negative buoyancy of parcel in cloud layer
            delthvu_c(ii) = delthvu_c(ii) - dtv_min_dil(ii) *                  &
                                      ( z_half_c(ii,k) - z_lcl_C(ii) )
          end if
        end if
      end if

      !-----------------------------------------------------------------------
      ! While doing parcel ascent
      ! (a) find minimum buoyancy
      ! (b) integrate CAPE over the ascent
      !-----------------------------------------------------------------------

      if (k > nlcl_c(ii) .and. k < model_levels ) then

        dz = z_half_c(ii,k+1) - z_half_c(ii,k)
        inc = g *  buoyancy(ii,k) * dz/t_dens_env(ii,k)


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

          dtv_min(ii) = min( dtv_min(ii),                                      &
                           buoyancy(ii,k)/exner_theta_levels_c(ii,k)  )

        end if    ! test on topprof

        ! Not reached top of ascent - dilute ascent
        if (.not. topprof_dil(ii)) then
          ! require undilute buoyancy here as calculating delthvu as undilute
          dtv_min_dil(ii) = min( dtv_min_dil(ii),                              &
                            buoyancy(ii,k)/exner_theta_levels_c(ii,k) )

          dz = z_half_c(ii,k+1) - z_half_c(ii,k)
          ! undilute value
          delthvu_c(ii) = delthvu_c(ii) + buoyancy(ii,k)* dz                   &
                                            /exner_theta_levels_c(ii,k)

          ! calculation of CIN and CAPE from profiles - undilute values
          ! not USED at present so commenting out - may want in future?
          !          inc =g * buoyancy(ii,k) * dz /t_dens_env(ii,k)

          !          if (inc <  0.0) then
          !            CIN2_c(ii)  = CIN2_c(ii) + inc
          !          else      ! CAPE holds only positive part
          !            CAPE2_c(ii) = CAPE2_c(ii) + inc
          !          end if

                    ! ql_ad = -dqsat/dz * zcld

          ql_ad_c(ii) = -one* dqsatdz(ii,k)                                    &
                          *(z_half_c(ii,k+1) - z_half_c(ii,nlcl_c(ii)+1))
        end if    ! test on topprof

      end if

    end do     ! level loop

  end do   ! ii loop
!$OMP end PARALLEL do

  !     write(6,*) ' k neutral ',k_neutral(i),k_neutral_dil(i)
  !     write(6,*) ' k max buoy ',k_max(i),k_max_dil(i)
  !     write(6,*) 'max buoy ',max_buoy(i),max_buoy_dil(i)
else     ! undilute ascent only

!$OMP  PARALLEL do DEFAULT(none) private(ii, k, dz, inc)                       &
!$OMP  SHARED(nunstable,k_max,k_neutral,max_buoy,ql_ad_c,delthvu_c,            &
!$OMP         CAPE_c,CIN_c,dtv_min,topprof,model_levels,                       &
!$OMP         nlcl_c,k_plume,z_lcl_c,thv_pert,  z_full_c,z_half_c,             &
!$OMP         exner_theta_levels_c,buoyancy, t_dens_env,                       &
!$OMP         dqsatdz,zh_c, g) SCHEDULE(STATIC)
  ! The main index of the loop is "ii" so that it can be
  ! parallelised using OpenMP.
  do ii=1,nunstable

    !-----------------------------------------------------------------------
    ! Initialise output arrays

    topprof(ii)     = .false.

    k_max(ii)          = 1
    k_neutral(ii)      = 1

    dtv_min(ii)   = zero
    max_buoy(ii)  = zero
    delthvu_c(ii) = zero
    ql_ad_c(ii)   = zero
    CAPE_c(ii)    = zero
    CIN_c(ii)     = zero

    do  k = 2,model_levels

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

        if ( (buoyancy(ii,k)  <=  - thv_pert(ii))                              &
          !                      or reached top of model
                   .or. (k  >   model_levels-1)  ) then

          k_neutral(ii) = k-1
          topprof(ii) = .true.
          zh_c(ii) = z_half_c(ii,k)

          if ( delthvu_c(ii)  >   zero) then
            ! compensate for any negative buoyancy of parcel in cloud layer
            delthvu_c(ii) = delthvu_c(ii) - dtv_min(ii) *                      &
                                      ( z_half_c(ii,k) - z_lcl_c(ii) )
          end if
        end if
      end if

      !-----------------------------------------------------------------------
      ! While doing parcel ascent
      ! (a) find minimum buoyancy
      ! (b) integrate CAPE over the ascent
      !-----------------------------------------------------------------------

      if (k > nlcl_c(ii) .and. k < model_levels ) then

        dz = z_half_c(ii,k+1) - z_half_c(ii,k)
        inc = g *  buoyancy(ii,k) * dz/t_dens_env(ii,k)


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

          ql_ad_c(ii) = -one* dqsatdz(ii,k)                                    &
                       *(z_half_c(ii,k+1) - z_half_c(ii,nlcl_c(ii)+1))

          dtv_min(ii) = min( dtv_min(ii),                                      &
                           buoyancy(ii,k)/exner_theta_levels_c(ii,k)  )

          delthvu_c(ii) = delthvu_c(ii) + buoyancy(ii,k)* dz                   &
                                            /exner_theta_levels_c(ii,k)

        end if    ! test on topprof

      end if

    end do     ! level loop
  end do   ! ii loop
!$OMP end PARALLEL do

end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine cv_parcel_neutral_dil
end module cv_parcel_neutral_dil_mod
