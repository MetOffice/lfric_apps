! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

module calc_3d_cca_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'CALC_3D_CCA_MOD'
contains

!  Subroutine CALC_3D_CCA: Calculates a conv. cld amt on theta model levels.
!
!  Subroutine Interface:

subroutine calc_3d_cca                                                         &
  ( np_field, npnts, nlev, n_cca_lev, nbl, cld_base, cld_top, p_lyr_bnds       &
  , frz_lev, cca_2d, cca_3d, z_theta, z_rho )

use cv_param_mod, only:                                                        &
    deep_dp, anv_pressure, anv_height, anv_model_levels,                       &
    anv_limited_pressure_depth

use cv_run_mod, only:                                                          &
    l_cloud_deep, tower_factor, anvil_factor, anv_opt, l_ccrad

use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook
implicit none
!
! Description:
! ------------
! Calculates a 3D convective cloud amount (i.e. on theta model levels) from
! the 2D convective cloud amount array according to parameters specified in
! the gui/namelist and the position of cloud base, cloud top and freezing level.
!
! Method:
! -------
! The 2D convective cloud amount is expanded into the vertical by applying it
! between the cloud base and top with the additional constraints that:
!
!         (ia)  If the cloud base is in the boundary layer (Original)
!         (ib)  If the cloud base is below the freezing level (CCRad)
!         (ii)  Cloud top is above the freezing level and
!         (iii) The cloud is more than 500mb deep
!
! Then the cloud below the freezing level will be multiplied by TOWER_FACTOR,
! and the cloud above the freezing level will be linearly
! (model level/height/pressure(default)) increased to cloud top where it will
! be equal to the 2D fraction * ANVIL_FACTOR.
!
! NOTE: ***** The above method needs to be rewritten if these mods are********
!       ***** Implemented*****************************************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!

!-----------------------------------------------------------------------------
!   Scalar arguments with intent(in):
!-----------------------------------------------------------------------------

integer, intent(in)  :: npnts        ! Number of points
integer, intent(in)  :: np_field     ! Full data length
integer, intent(in)  :: nlev         ! Number of levels
integer, intent(in)  :: n_cca_lev    ! Number of cca levels
integer, intent(in)  :: nbl          ! Number of Boundary layer levels

!-----------------------------------------------------------------------------
!   Array  arguments with intent(in):
!-----------------------------------------------------------------------------

real(kind=real_umphys),    intent(in)  :: z_theta    (np_field,   nlev)
                                                      ! z (th layer centres)
real(kind=real_umphys),    intent(in)  :: z_rho      (np_field,   nlev)
                                                      ! z (rh level  bounds)
real(kind=real_umphys),    intent(in)  :: p_lyr_bnds (np_field, 0:nlev)
                                                ! Pressure on layer
                                                ! boundaries (rho levels-1)

real(kind=real_umphys),    intent(in)  :: cca_2d   (npnts)
                                         ! 2D convective cloud amount
integer, intent(in)  :: cld_top  (npnts) ! Conv. cloud top  (theta level)
integer, intent(in)  :: cld_base (npnts) ! Conv. cloud base (theta level)
integer, intent(in)  :: frz_lev  (npnts) ! Freezing level   (theta level)




!-----------------------------------------------------------------------------
!   Array  arguments with intent(out):
!-----------------------------------------------------------------------------

real(kind=real_umphys),    intent(out) :: cca_3d(np_field, n_cca_lev)
                                         ! Convective cloud amount on
                                         ! model levels (theta levels)


! Local variables:
! -----------------
integer             :: i, k             ! Loop counters


integer :: anv_lev     ! Base level of 'anvil'
real(kind=real_umphys)    :: anv_dep     ! Anvil depth in model levels
real(kind=real_umphys)    :: anv_p_dep   ! Anvil depth in pressure
real(kind=real_umphys)    :: anv_z_dep   ! Anvil depth in metres
real(kind=real_umphys)    :: anv_p_base  ! Anvil base pressure (rho-level)

real(kind=real_umphys)    :: p_cld_base
                       ! Pressure at lowest  cloud layer BOUNDARY
real(kind=real_umphys)    :: p_cld_top
                       ! Pressure at highest cloud layer BOUNDARY



logical :: cbct_crit (npnts) ! .true. if cloud base/top are sensible
logical :: dep_crit  (npnts) ! .true. if depth criteria met
logical :: base_crit (npnts) ! .true. if cloud base criteria met
logical :: anv_on    (npnts) ! .true. if all anvil criteria met
integer :: tp_of_lp  (npnts) ! Index of cloud top, required so that CCRad
                             ! correction can be reverted

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CALC_3D_CCA'

!-----------------------------------------------------------------------------
! Code Statements
!-----------------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!---------------------------------------------------------------------------
! 1.0 Initialise local arrays
!---------------------------------------------------------------------------

do i=1, npnts
  cbct_crit(i)  = .false.
  dep_crit(i)   = .false.
  base_crit(i)  = .false.
  anv_on(i)     = .false.
  tp_of_lp(i)   = 0
end do



!---------------------------------------------------------------------------
! 2.0 Option changes/bugfixes specified by CCRad
!---------------------------------------------------------------------------
if (l_ccrad) then

  do i=1, npnts
    ! Check for sensible Cloud base/top
    cbct_crit(i)  = ((cld_top(i)  >= cld_base(i)) .and.                        &
                     (cld_top(i)  /= 0)           .and.                        &
                     (cld_base(i) /= 0))

    ! Check for cloud base/top above/below freezing level
    base_crit(i)  = ((cld_base(i) < frz_lev(i))   .and.                        &
                     (cld_top(i)  > frz_lev(i)))

    ! Change index of cloud top (Bug fix)
    tp_of_lp(i)   = cld_top(i)
  end do

else
  ! Original test criteria
  do i=1, npnts
    cbct_crit(i) = .true.
    base_crit(i) = ((cld_base(i) < nbl)           .and.                        &
                    (cld_top(i)  > frz_lev(i)))
    tp_of_lp(i)  = cld_top(i)-1
  end do

end if


! Locate grid points where depth criteria is satisfied
if (l_cloud_deep) then
  do i=1, npnts
    p_cld_base  = p_lyr_bnds(i,max(cld_base(i)-1,0))
    p_cld_top   = p_lyr_bnds(i,cld_top(i))
    dep_crit(i) = (p_cld_base - p_cld_top) >= deep_dp
  end do
else
  do i=1,npnts
    dep_crit(i) = .true.
  end do
end if



!---------------------------------------------------------------------------
! 3.0 Locate grid points where all anvil criteria are satisfied
!---------------------------------------------------------------------------
do i=1, npnts
  if (base_crit(i) .and. dep_crit(i)) then
    anv_on(i) = .true.
  end if
end do


!---------------------------------------------------------------------------
! 4.0 Apply CCA Profiles
!---------------------------------------------------------------------------
do i=1, npnts
  if ( cbct_crit(i)       .and.                                                &
      (cca_2d(i) > 0.0) ) then

    if (anv_on(i)) then

      !---------------------------------------------------------------------
      ! 4.1a Cloud satisfies anvil criteria: Apply Anvil
      !---------------------------------------------------------------------
      select case(anv_opt)
      case (anv_height)
        !-----------------------------------------------------------------
        ! CCA increases with height from freezing level to cloud-top
        !-----------------------------------------------------------------
        anv_lev   = max(cld_base(i), frz_lev(i))
        anv_z_dep = z_rho(i,cld_top(i)) - z_rho(i,anv_lev-1)

        do k=anv_lev, cld_top(i)

          cca_3d(i,k) =                                                        &
                  (anvil_factor - tower_factor)*cca_2d(i)                      &
                * (z_theta(i,k) - z_rho(i,anv_lev-1)) / anv_z_dep              &
                + (cca_2d(i) * tower_factor)

          if (cca_3d(i,k) >= 1.0) then
            cca_3d(i,k) = 0.99
          end if

        end do


      case (anv_model_levels)
        !-----------------------------------------------------------------
        ! CCA increases with model level from freezing level to cloud-top:
        ! (original code)
        !-----------------------------------------------------------------
        anv_lev = max(cld_base(i), frz_lev(i))
        anv_dep = cld_top(i) - anv_lev

        do k=anv_lev, tp_of_lp(i)

          cca_3d(i,k) =                                                        &
                  (anvil_factor - tower_factor) * cca_2d(i)                    &
                * (k - anv_lev + 1) / anv_dep                                  &
                + (cca_2d(i) * tower_factor)

          if (cca_3d(i,k)  >=  1.0) then
            cca_3d(i,k) = 0.99
          end if

        end do


      case (anv_pressure)
        !-----------------------------------------------------------------
        ! CCA increases with pressure from freezing level to cloud-top
        !-----------------------------------------------------------------
        anv_lev   = max(cld_base(i), frz_lev(i))
        anv_p_dep = p_lyr_bnds(i,anv_lev-1) - p_lyr_bnds(i,cld_top(i))

        do k=anv_lev, tp_of_lp(i)

          cca_3d(i,k) =                                                        &
                (anvil_factor - tower_factor)*cca_2d(i)                        &
              * (p_lyr_bnds(i,anv_lev-1) -  p_lyr_bnds(i,k)) / anv_p_dep       &
              + (cca_2d(i) * tower_factor)

          if (cca_3d(i,k) >= 1.0) then
            cca_3d(i,k) = 0.99
          end if

        end do

      case (anv_limited_pressure_depth)
        !-----------------------------------------------------------------
        ! Pressure based, but limit anvil depth to 5000.0 pa
        !-----------------------------------------------------------------
        anv_p_base = p_lyr_bnds(i,cld_top(i)) + 5000.0

        ! Ensure that Anvil is at least 2 levels deep, so only loop to
        ! layer below cloud top so that it will be 2 levels deep even
        ! if k is set at top of loop
        do k=1, tp_of_lp(i)-1
          if (p_lyr_bnds(i,k) > anv_p_base) then
            anv_lev = k
          end if
        end do

        anv_p_dep = p_lyr_bnds(i,anv_lev-1) - p_lyr_bnds(i,cld_top(i))

        do k=anv_lev, tp_of_lp(i)

          cca_3d(i,k) =                                                        &
                (anvil_factor - tower_factor)*cca_2d(i)                        &
              * (p_lyr_bnds(i,anv_lev-1) -  p_lyr_bnds(i,k)) / anv_p_dep       &
              + (cca_2d(i) * tower_factor)

          if (cca_3d(i,k) >= 1.0) then
            cca_3d(i,k) = 0.99
          end if

        end do

      end select

      !---------------------------------------------------------------------
      ! 4.1b Cloud satisfies anvil criteria: Apply Tower below anvil base
      !---------------------------------------------------------------------
      do k=cld_base(i), anv_lev-1
        cca_3d(i,k) = tower_factor * cca_2d(i)
      end do

    else

      ! Anvil criteria not met
      do k=cld_base(i), tp_of_lp(i)
        cca_3d(i,k) = cca_2d(i)
      end do

    end if      ! End test on anvil criteria

    !-----------------------------------------------------------------------
    ! Finally check there is no cloud below cloud base or above cloud top:
    ! (original code)
    !-----------------------------------------------------------------------
    if (.not. l_ccrad) then
      do k=1, (cld_base(i)-1)
        cca_3d(i,k) = 0.0
      end do

      if (cld_top(i) > 0) then
        do k=cld_top(i), n_cca_lev
          cca_3d(i,k) = 0.0
        end do
      end if

    end if      ! l_ccrad
  end if      ! cca_2d > 0 and sensible ccb/cct
end do      ! loop over npnts


!
!=============================================================================
!  End of anvil calculation
!=============================================================================
!
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine calc_3d_cca
end module calc_3d_cca_mod
