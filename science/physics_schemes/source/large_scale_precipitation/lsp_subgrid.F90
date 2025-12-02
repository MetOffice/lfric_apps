! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme.
!
!  Contains the following subroutines:
!
!  lsp_subgrid (Subbgrid-scale set ups and checks)
!  lsp_qclear  (Calculates mean water vapour of clear-sky portion of gridbox)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

module lsp_subgrid_mod

implicit none


! NOTE: This interface is also called from outside the lsp scheme by:
!       atmosphere/UKCA/ukca_main1-ukca_main1.F90
!       atmosphere/GLOMAP_CLIM/prepare_fields_for_radaer_mod.F90
interface lsp_qclear
  module procedure lsp_qclear_64b, lsp_qclear_32b
end interface

character(len=*), parameter, private :: ModuleName='LSP_SUBGRID_MOD'

contains

subroutine lsp_subgrid(                                                        &
  points,                                                                      &
                                          ! Number of points
  q, qcf_cry, qcf_agg, qcftot, t,                                              &
                                          ! Water contents and temp
  qsl, qs,                                                                     &
                                          ! Saturated water contents
  snow_cry, snow_agg,                                                          &
                                          ! Fall-fluxes of ice from above
  cry_nofall, agg_nofall,                                                      &
                                          ! Fractions of ice not falling out
  q_ice, q_clear, q_ice_1, q_ice_2,                                            &
                                          ! Local vapour contents
  area_liq,area_mix,area_ice,area_clear,                                       &
                                          ! Cloud frac partitions
  area_ice_1, area_ice_2,                                                      &
                                          ! Subdivision of area_ice
  areamix_over_cfliq,                                                          &
                                          ! area_mix/cfliq
  rain_liq,rain_mix,rain_ice,rain_clear,                                       &
                                          ! Rain overlap partitions
  cftot, cfliq, cfice, cficei,                                                 &
                                          ! Cloud fractions for
                                          ! partition calculations
  cf, cff, rainfrac, rainfraci, precfrac_k, rain_new,                          &
                                          ! Cloud and rain fractions
                                          ! for updating
  dqprec_liq, dqprec_mix, dqprec_ice, dqprec_clear, dqprec_new,                &
                                          ! Partition qrain/graup increments
  lsrcp,                                                                       &
                                          ! Latent heat of sublim./cp
  rhcpt,                                                                       &
                                          ! RH crit values
  wtrac_mp_cpr, wtrac_mp_cpr_old)
                                          ! Water tracers

!Use in reals in lsprec precision, both microphysics related and general atmos
use lsprec_mod, only: qcfmin, ice_width, zerodegc,                             &
                      zero, half, one, two

use mphys_inputs_mod, only: l_mcr_precfrac, i_update_precfrac, i_homog_areas

! Temporary bug-fixes
use science_fixes_mod, only: l_fix_mcr_frac_ice

! Cloud modules- logicals and integers
use cloud_inputs_mod,  only: l_subgrid_qv

! Water tracers
use free_tracers_inputs_mod, only: l_wtrac, n_wtrac
use wtrac_mphys_mod,         only: mp_cpr_wtrac_type, mp_cpr_old_wtrac_type

! Use in kind for large scale precip, used for compressed variables passed down
! from here
use um_types,             only: real_lsprec

! Dr Hook Modules
use yomhook,           only: lhook, dr_hook
use parkind1,          only: jprb, jpim

implicit none

! Purpose:
!   Perform the subgrid-scale setting up calculations

! Method:
!   Parametrizes the width of the vapour distribution in the part
!   of the gridbox which does not have liquid water present.
!   Calculates the overlaps within each gridbox between  the cloud
!   fraction prognostics and rainfraction diagnostic.
!

! Description of Code:
!   Fortran95.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.

! The subgrid calculations are a necessary step to calculating the
! subsequent deposition and sublimation transfers and for setting up
! the partition information that is used by the subsequent transfers.

! Subroutine Arguments

integer, intent(in) ::                                                         &
  points            ! Number of points to calculate

real (kind=real_lsprec), intent(in) ::                                         &
  qs(points),                                                                  &
                        ! Saturated mixing ratio wrt ice
  qsl(points),                                                                 &
                        ! Saturated mixing ratio wrt liquid
  cfliq(points),                                                               &
                        ! Fraction of gridbox with liquid cloud
  lsrcp,                                                                       &
                        ! Latent heat of sublimation
                        ! / heat capacity of air / K
  rhcpt(points)     ! RH crit values

real (kind=real_lsprec), intent(in out) :: rainfrac(points)
real (kind=real_lsprec), intent(in out) :: rainfraci(points)
                        ! Fraction of gridbox containing rain, and inverse
real (kind=real_lsprec), intent(in) :: precfrac_k(points)
                        ! Prognostic precipitation fraction

! Area-fraction in-which any rain (and optionally graupel) is produced
! in air that didn't already contain any rain / graupel.
real (kind=real_lsprec), intent(out) :: rain_new(points)

! Increments to rain (and optionally graupel) within each of the
! 4 sub-grid partitions rain_liq, rain_mix, rain_ice, rain_clear.
! Used to update the prognostic precipitation fraction.
real (kind=real_lsprec), intent(out) :: dqprec_liq(points)
real (kind=real_lsprec), intent(out) :: dqprec_mix(points)
real (kind=real_lsprec), intent(out) :: dqprec_ice(points)
real (kind=real_lsprec), intent(out) :: dqprec_clear(points)
! Additional category for increments due to creation of new precip
! outside the initial precip region
real (kind=real_lsprec), intent(out) :: dqprec_new(points)

real (kind=real_lsprec), intent(in out) ::                                     &
  q(points),                                                                   &
                        ! Vapour content / kg kg-1
  qcf_cry(points),                                                             &
                        ! Ice crystal content / kg kg-1
  qcf_agg(points),                                                             &
                        ! Ice aggregate content / kg kg-1
  qcftot(points),                                                              &
                        ! Total ice content before advection / kg kg-1
  t(points),                                                                   &
                        ! Temperature / K
  cf(points),                                                                  &
                        ! Current cloud fraction
  cff(points)       ! Current ice cloud fraction

real (kind=real_lsprec), intent(in) ::                                         &
  snow_cry(points),                                                            &
                        ! Ice crystal fall-flux from above / kg m-2 s-1
  snow_agg(points),                                                            &
                        ! Ice aggregate fall-flux from above / kg m-2 s-1
  cry_nofall(points),                                                          &
                        ! Fraction of ice aggregate mass not falling out
  agg_nofall(points)
                        ! Fraction of ice crystal mass not falling out

real (kind=real_lsprec), intent(out) ::                                        &
  q_clear(points),                                                             &
                        ! Local vapour in clear-sky region / kg kg-1
  q_ice(points),                                                               &
                        ! Local vapour in ice-only region  / kg kg-1
  q_ice_1(points),                                                             &
                        ! Local vapour in ice-only regions that are:
  q_ice_2(points),                                                             &
                        !   1, depositing; and 2, subliming.
  cftot(points),                                                               &
                        ! Modified cloud fraction for partition calc.
  cfice(points),                                                               &
                        ! Modified ice cloud frac. for partition calc.
  cficei(points),                                                              &
                        ! 1/cfice
  area_liq(points),                                                            &
                        ! Frac of gridbox with liquid cloud but no ice
  area_mix(points),                                                            &
                        ! Frac of gridbox with liquid and ice cloud
  area_ice(points),                                                            &
                        ! Frac of gridbox with ice cloud but no liquid
  area_clear(points),                                                          &
                        ! Frac of gridbox with no cloud
  area_ice_1(points),                                                          &
                        ! Frac of gridbox where ice-only cloud is:
  area_ice_2(points),                                                          &
                        !  1, depositing; and 2, subliming.
  areamix_over_cfliq(points),                                                  &
                        ! area_mix/cfliq
  rain_liq(points),                                                            &
                        ! Frac of gbox with rain and liquid but no ice
  rain_mix(points),                                                            &
                        ! Frac of gbox with rain and liquid and ice
  rain_ice(points),                                                            &
                        ! Frac of gbox with rain and ice but no liquid
  rain_clear(points)! Frac of gbox with rain but no condensate

! Water tracer fields
type(mp_cpr_wtrac_type), intent(in out)     :: wtrac_mp_cpr(n_wtrac)
type(mp_cpr_old_wtrac_type), intent(in out) :: wtrac_mp_cpr_old

! Local Variables

integer ::                                                                     &
  i, i_wt               ! Loop counters

real (kind=real_lsprec) ::                                                     &
  temp7,                                                                       &
                        ! Temporary in width of PDF calculation
  width             ! Full width of vapour distribution in ice and
                        ! clear sky.

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_SUBGRID'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! If using prognostic precipitation fraction; reset rainfrac consistent
! with the latest value of the prognostic (which has been updated
! by sedimentation).  The sub-partitioning of the rain-fraction into
! bits that overlap with different cloud regions
! (rain_liq, rain_ice, rain_mix, rain_clear)
! will now represent sub-partitions of the prognostic precip fraction.
if ( l_mcr_precfrac ) then
  do i = 1, points
    ! Impose min limit on rainfrac, to avoid rare div-by-zero
    rainfrac(i) = max( precfrac_k(i), 0.001_real_lsprec )
    rainfraci(i) = one / rainfrac(i)
    rain_new(i) = zero
  end do
  if ( i_update_precfrac == i_homog_areas ) then
    ! If updating the precip fraction based on assumed homogeneous areas,
    ! we need to store the precip mass increments within each area partition.
    ! Initialise precip mass increments from each partition to zero:
    do i = 1, points
      dqprec_liq(i) = zero
      dqprec_mix(i) = zero
      dqprec_ice(i) = zero
      dqprec_clear(i) = zero
      dqprec_new(i) = zero
    end do
  end if
end if

do i = 1, points

      !-----------------------------------------------
      ! Check that ice cloud fraction is sensible.
      !-----------------------------------------------
      ! Difference between the way PC2 and non-PC2 code operates
      ! is kept here in order to be tracable across model versions.
      ! However, perhaps the code ought to be the same.
        ! 0.001 is to avoid divide by zero problems

  cfice(i)  = max( cff(i), 0.001_real_lsprec )
  cficei(i) = one/cfice(i)
  cftot(i)  = cf(i)
  cftot(i)  = min( max( cftot(i), cfice(i) ),( cfice(i) + cfliq(i)) )

    ! -----------------------------------------------
    ! Calculate overlaps of liquid, ice and rain fractions
    ! -----------------------------------------------
  area_liq(i) = max(cftot(i)-cfice(i),zero)
  area_mix(i) = max(cfice(i)+cfliq(i)-cftot(i),zero)

  if (cfliq(i) /= zero) then
    areamix_over_cfliq(i) = area_mix(i)/cfliq(i)
  end if

  if (cfice(i) == cftot(i)) then
    area_mix(i)           = cfliq(i)
    areamix_over_cfliq(i) = one
  end if

      ! Remove tiny mixed phase areas
  if (area_mix(i) < 0.0001_real_lsprec) then
    area_mix(i)           = zero
    areamix_over_cfliq(i) = zero
  end if

  if ( l_mcr_precfrac ) then
    ! Ensure we still have cfliq = area_liq + area_mix
    ! after the above checks; this is assumed (and required)
    ! by the prognostic precip fraction code.
    area_liq(i) = cfliq(i) - area_mix(i)
  end if

  area_ice(i) = max(cftot(i)-cfliq(i),zero)
  area_clear(i) = max(one-cftot(i),zero)

  ! Assume maximal overlap of rain with liq, then mix, then ice.
  rain_liq(i) = max(min(area_liq(i),rainfrac(i)),zero)
  rain_mix(i) = max(min(area_mix(i),rainfrac(i)-rain_liq(i)),zero)
  rain_ice(i) =                                                                &
    max(min(area_ice(i),rainfrac(i)-rain_liq(i)-rain_mix(i)),zero)

  ! Assign any remaining rain fraction to the clear-sky partition
  rain_clear(i) =                                                              &
    max(rainfrac(i)-rain_liq(i)-rain_mix(i)-rain_ice(i),zero)

end do

if ( l_fix_mcr_frac_ice ) then
  ! The copy of ice-cloud fraction limited to be above 0.001 above
  ! is used to calculate the increments to cff in various parts
  ! of the code. But this causes problems because if the actual
  ! ice cloud fraction is smaller than 0.001, the increments can
  ! remove nearly all the ice cloud fraction whilst leaving a
  ! significant ice mass. The implied in-cloud ice qcf/cff then
  ! becomes stupidly large.
  ! Under this bug-fix switch, avoid this by resetting the
  ! prognostic ice fraction consistent with the copies
  ! used to compute the increments
  do i = 1, points
    if ( qcf_agg(i) > zero  .or. qcf_cry(i) > zero .or.                        &
         snow_agg(i) > zero .or. snow_cry(i) > zero ) then
      ! Only reset where either there is already ice mass present
      ! or ice is falling in from above
      cff(i) = cfice(i)
      cf(i)  = cftot(i)
    end if
  end do
end if

call lsp_qclear(                                                               &
!     Input fields
     q, qs, qsl, qcftot, cfliq, cftot, rhcpt,                                  &
!     Output field
     q_clear,                                                                  &
!     Array dimensions
     points)

do i = 1, points
  width   = one
  if (cfliq(i)  <   one) then

    if (area_ice(i)  >   zero) then
      if (l_subgrid_qv) then
        width = two *(one-rhcpt(i))*qsl(i)                                     &
            *max((one-half*qcftot(i)/(ice_width * qsl(i))), 0.001_real_lsprec)

        ! The full width cannot be greater than 2q because otherwise
        ! part of the gridbox would have negative q. Also ensure that
        ! the full width is not zero (possible if rhcpt is 1).
        ! 0.001 is to avoid divide by zero problems
        width = min(width, max(two*q(i),0.001_real_lsprec*qs(i)))

        ! If there is no partitioning of qv between clear-sky and
        ! ice-cloud, then use the same value of qv everywhere outside
        ! of the liquid-cloud
        q_ice(i) = (q(i)-cfliq(i)*qsl(i)-area_clear(i)*q_clear(i))             &
                   / area_ice(i)
      else
        q_ice(i) = q_clear(i)
      end if !l_subgrid_qv
    else
      q_ice(i) = zero               ! q_ice is a dummy value here
    end if  ! area_ice gt 0

  else ! cf_liq lt 1

        ! -----------------------------------------------
        ! Specify dummy values for q_clear and q_ice
        ! -----------------------------------------------
    q_clear(i) = zero
    q_ice(i)   = zero

  end if ! cf_liq lt 1

      ! -------------------------------------------------
      ! Remove any small amount of ice to be tidy.
      ! -------------------------------------------------
      ! If QCF is less than QCFMIN and is not growing by deposition
      ! (assumed to be given by RHCPT) then evaporate it.
  if ((qcf_cry(i)*cry_nofall(i) + qcf_agg(i)*agg_nofall(i)) <  qcfmin) then
    if (t(i) >  zerodegc .or.                                                  &
       (q_ice(i)  <=  qs(i) .and. area_mix(i)  <=  zero)                       &
       .or. (qcf_cry(i)+qcf_agg(i)) <  zero) then
      q(i) = q(i) +qcf_cry(i)+qcf_agg(i)
      t(i) = t(i) - lsrcp * (qcf_cry(i)+qcf_agg(i))
      qcf_cry(i)=zero
      qcf_agg(i)=zero
      ! Update water tracers consistently
      ! (No isotopic fractionation as qcf is fully removed)
      if (l_wtrac) then
        do i_wt = 1, n_wtrac
          wtrac_mp_cpr(i_wt)%q(i)   = wtrac_mp_cpr(i_wt)%q(i)                  &
                                      + wtrac_mp_cpr(i_wt)%qcf(i)
          wtrac_mp_cpr(i_wt)%qcf(i) = zero
        end do
      end if
    end if ! T gt 0 etc.
  end if ! qcf_cry+qcf_agg lt qcfmin

      ! -------------------------------------------------
      ! First estimate of partition sizes for ice sublimation
      ! and deposition and vapour contents within these partitions
      ! -------------------------------------------------
  if (q_ice(i)  >   qs(i)) then
        ! First estimate is to use a deposition process
    area_ice_1(i) = area_ice(i)
    area_ice_2(i) = zero
    q_ice_1(i)    = q_ice(i)
    q_ice_2(i)    = qs(i)       ! Dummy value

  else ! q_ice gt qs
        ! First estimate is to use a sublimation process
    area_ice_1(i) = zero
    area_ice_2(i) = area_ice(i)
    q_ice_1(i)    = qs(i)       ! Dummy value
    q_ice_2(i)    = q_ice(i)

  end if ! q_ice gt qs


  ! -----------------------------------------------------
  ! If the water vapor is being partitioned between sublimating
  ! and deposition ice-regions, then calculate the detailed
  ! partition. Else there is no further
  ! subgrid partitioning of qv and the first-estimate(above) is used
  ! -----------------------------------------------------
  if ( l_subgrid_qv ) then
      ! -------------------------------------------------
      ! Detailed estimate of partition sizes for ice sublimation
      ! and deposition and vapour contents within these partitions
      ! -------------------------------------------------
    if (area_ice(i)  >   zero) then
      ! Temp7 is the estimate of the proportion of the gridbox
      ! which contains ice and has local q > than qs (wrt ice)
      temp7 = half*area_ice(i) + (q_ice(i)-qs(i)) / width

      if (temp7 >  zero .and. temp7 <  area_ice(i)) then
          ! Calculate sizes of regions and q in each region
          ! These overwrite previous estimates
        area_ice_1(i) = temp7
        area_ice_2(i) = area_ice(i) - area_ice_1(i)
        q_ice_1(i) = qs(i) + half * area_ice_1(i) * width
        q_ice_2(i) = qs(i) - half * area_ice_2(i) * width
      end if ! temp7 gt 0 etc.

    end if ! area_ice gt 0

  end if   ! l_subgrid_qv=.true.


end do  ! Points

! If water tracers, store current values of q and qcf
if (l_wtrac) then
  do i = 1, points
    wtrac_mp_cpr_old%q(i)   = q(i)
    wtrac_mp_cpr_old%qcf(i) = qcf_agg(i)
    wtrac_mp_cpr_old%t(i)   = t(i)
  end do
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_subgrid

!Note regarding variable precision:
!This routine is used in various places beyond the lsp scheme which may require
!either 64 or 32 bit calculation. Therefore, we cannot use the real_lsprec
!parameter as seen generally in the lsp_* modules. Instead, we create an
!interface with both 32 and 64 bit versions available.

subroutine lsp_qclear_64b(                                                     &
!     Input fields
     q, qsmr, qsmr_wat, qcf, cloud_liq_frac, cloud_frac, rhcrit,               &
!     Output field
     q_clear,                                                                  &
!     Array dimensions
     npnts)

  ! Cloud modules
use cloud_inputs_mod, only: l_subgrid_qv
use cloud_inputs_mod, only: ice_width
use um_types,         only: real_64

! Dr Hook Modules
use yomhook,           only: lhook, dr_hook
use parkind1,          only: jprb, jpim

implicit none
integer, parameter :: prec = real_64
character(len=*), parameter :: RoutineName='LSP_QCLEAR_64B'
#include "lsp_subgrid_lsp_qclear.h"
return
end subroutine lsp_qclear_64b

subroutine lsp_qclear_32b(                                                     &
!     Input fields
     q, qsmr, qsmr_wat, qcf, cloud_liq_frac, cloud_frac, rhcrit,               &
!     Output field
     q_clear,                                                                  &
!     Array dimensions
     npnts)

  ! Cloud modules
use cloud_inputs_mod, only: l_subgrid_qv
use cloud_inputs_mod, only: ice_width
use um_types,         only: real_32

! Dr Hook Modules
use yomhook,           only: lhook, dr_hook
use parkind1,          only: jprb, jpim

implicit none
integer, parameter :: prec = real_32
character(len=*), parameter :: RoutineName='LSP_QCLEAR_32B'
#include "lsp_subgrid_lsp_qclear.h"
return
end subroutine lsp_qclear_32b

end module lsp_subgrid_mod

