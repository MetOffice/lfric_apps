! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Autoconversion of graupel.
! Subroutine Interface:
module lsp_graup_autoc_mod

implicit none

character(len=*), parameter, private :: ModuleName='LSP_GRAUP_AUTOC_MOD'

contains

subroutine lsp_graup_autoc(                                                    &
  points, timestep,                                                            &
                                          ! Number of points and tstep
  qrain, qcf_agg, qgraup, t, rho,                                              &
                                          ! Water contents and temp
  psacw, psdep, ptransfer,                                                     &
                                          ! Mass transfer diagnostic
  agg_nofall,                                                                  &
                                          ! Fraction of qcf_agg not falling out
  one_over_tsi,                                                                &
                                          ! 1/(timestep*iterations)
  area_mix, area_mix_orog, l_enh_rime,                                         &
                                          ! Orographic enhancement of riming
  dqprec_mix, dqprec_ice, dqprec_new,                                          &
                                          ! Precip increment in mixed-phase
  rain_mix, rain_ice, rain_new, precfrac_k                                     &
                                          ! Area fractions
  )

!Use in reals in lsprec precision, both microphysics related and general atmos
use lsprec_mod, only: auto_graup_qcf_thresh, auto_graup_t_thresh,              &
                      auto_graup_coeff, zerodegc,                              &
                      zero, one, small_number

use mphys_inputs_mod, only: graupel_option, gr_field_psd,                      &
                            l_mcr_precfrac, l_subgrid_graupel_frac,            &
                            i_update_precfrac, i_homog_areas, i_sg_correl,     &
                            l_orogrime

use lsp_combine_precfrac_mod, only: lsp_combine_precfrac

! Use in kind for large scale precip, used for compressed variables passed down
! from here
use um_types,             only: real_lsprec

! Dr Hook Modules
use yomhook,             only: lhook, dr_hook
use parkind1,            only: jprb, jpim

implicit none

! Purpose:
!   Update cloud prognostics as a result of the autoconversion of
!   snow aggregates to graupel.

!  Method:
!   This is a source term for graupel when snow growth is dominated
!   by riming liquid water cloud and is the only term that can create
!   graupel where there was no graupel before.
!   The conversion rate is proportional to how much the riming rate
!   of aggregates (PSACW) exceeds the rate of growth due to vapour
!   deposition (PSDEP) and collection/accretion of ice crystals (PSACI).
!   The coefficient reduces the conversion rate as the riming
!   aggregates will not immediately increase their density to that of
!   graupel.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.


! Subroutine Arguments

integer, intent(in) ::                                                         &
  points
                        ! Number of points to process

real (kind=real_lsprec), intent(in) ::                                         &
  timestep,                                                                    &
                        ! Timestep / s
  one_over_tsi,                                                                &
                        ! 1/(timestep*iterations)
  t(points),                                                                   &
                        ! Temperature / K
  rho(points),                                                                 &
                        ! Air density / kg kg s-1
  psacw(points),                                                               &
                        ! Riming rate of snow aggs. / kg kg-1 s-1
  psdep(points),                                                               &
                        ! Deposition rate of snow aggs. / kg kg-1 s-1
  agg_nofall(points),                                                          &
                        ! Fraction of qcf_agg not falling out
  qrain(points)
                        ! Rain water content / kg kg-1

real (kind=real_lsprec), intent(in out) ::                                     &
  qgraup(points),                                                              &
                        ! Graupel content / kg kg-1
  qcf_agg(points),                                                             &
                        ! Ice water content of snow aggs. / kg kg-1
  ptransfer(points)
                        ! Autoconversion rate / kg kg-1 s-1

! Area fraction of mixed-phase cloud
real (kind=real_lsprec), intent(in) :: area_mix(points)
! Mixed-phase cloud area as diagnosed in the orographic riming enhancement
! (accounts for orographically-generated cloud that doesn't exist in the model)
real (kind=real_lsprec), intent(in) :: area_mix_orog(points)
! Flag for points where orographic riming enhancement has occured
logical, intent(in) :: l_enh_rime(points)

! Precip-mass increment within overlap of rain/graupel and mixed-phase cloud
real (kind=real_lsprec), intent(in out) :: dqprec_mix(points)
! Precip-mass increment within overlap of rain/graupel and ice-only cloud
real (kind=real_lsprec), intent(in out) :: dqprec_ice(points)
! Precip mass increment within area not yet containing rain or graupel
real (kind=real_lsprec), intent(in out) :: dqprec_new(points)

! Area of overlap of existing rain/graupel with mixed-phase cloud
real (kind=real_lsprec), intent(in) :: rain_mix(points)
! Area of overlap of existing rain/graupel with ice-only cloud
real (kind=real_lsprec), intent(in) :: rain_ice(points)
! Area where rain/graupel is produced in air that didn't have any already
real (kind=real_lsprec), intent(in out) :: rain_new(points)
! Prognostic precipitation fraction
real (kind=real_lsprec), intent(in out) :: precfrac_k(points)


! Local Variables

integer :: i


real (kind=real_lsprec) ::                                                     &
  dqi(points)
                        ! Transfer amount from ice to snow / kg kg-1

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_GRAUP_AUTOC'

real(kind=real_lsprec) :: ff       ! Thompson et al (2008) factor
real(kind=real_lsprec) :: ff_ratio ! Ratio of riming to deposition
                                   ! (pscaw / psdep)

! Area fraction where graupel created in air  not already containing graupel
real(kind=real_lsprec) :: area_new
! Fraction of increment occuring in air not already containing graupel
real(kind=real_lsprec) :: frac
! Work variable storing 1 / area_mix
real(kind=real_lsprec) :: tmp
! Area fraction where graupel created
real(kind=real_lsprec) :: area_inc(points)
! Precip content (rain+graupel) before the increment is applied
real(kind=real_lsprec) :: qprec(points)

! Flag for points where autoconversion has been calculated
logical :: l_calc(points)

! Parameters used in calculation of ff from Thompson et al (2008):
real(kind=real_lsprec), parameter :: ff_lower  = 0.05_real_lsprec
! Lower limit of ff

real(kind=real_lsprec), parameter :: ff_upper  = 0.75_real_lsprec
! Upper limit of ff

real(kind=real_lsprec), parameter :: ff_coeff  = 0.028_real_lsprec
! Coefficient used in calculation of ff

real(kind=real_lsprec), parameter :: ff_thresh = 5.0_real_lsprec
! Threshold value which ff_ratio must exceed to start graupel autoconversion

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
do i = 1, points

  if (rho(i)*qcf_agg(i)*agg_nofall(i)  >   auto_graup_qcf_thresh               &
      .and. t(i)  <  (zerodegc+auto_graup_t_thresh)) then
    l_calc(i) = .true.

        !-----------------------------------------------
        ! Calculate conversion rate / kg kg-1 s-1
        !-----------------------------------------------
    if ( graupel_option == gr_field_psd ) then
      ! Follow Thompson et al. (2008) method for generating
      ! graupel autoconversion

      ff_ratio = psacw(i)/psdep(i)

      if ( ff_ratio >= ff_thresh ) then
        ! Calculate the initial value of ff based on Thompson et al. (2008).
        ff = (ff_ratio - ff_thresh) * ff_coeff + ff_lower

        ! Apply the upper and lower limits to ff
        ff = max(ff_lower, min(ff_upper,  ff))

        ! Now modify the value of dqi (the transfer rate between snow and
        ! graupel) based on the values of riming, deposition and the ff
        ! factor.
        dqi(i) = ff * max(zero, psacw(i)-psdep(i))

      else ! ratio < ff_thresh

        ! Set dqi to zero as no graupel autoconversion will take place here.
        dqi(i) = zero

      end if ! ratio >= ff_thresh

    else ! graupel_option not gr_field_psd

      dqi(i) = auto_graup_coeff * max(zero, psacw(i)-psdep(i))

    end if ! test on graupel option

    dqi(i) = min(dqi(i)*timestep, qcf_agg(i))

        !-----------------------------------------------
        ! Store process rate / kg kg-1 s-1
        !-----------------------------------------------
    ptransfer(i) = ptransfer(i) + dqi(i) * one_over_tsi

        !-----------------------------------------------
        ! Recalculate water contents
        !-----------------------------------------------
    qgraup(i)  = qgraup(i)  + dqi(i)
    qcf_agg(i) = qcf_agg(i) - dqi(i)

        ! No cloud fraction updates

  else

    l_calc(i) = .false.
    dqi(i) = zero

  end if !  qcf_agg > threshold etc.

end do  ! Points


if ( l_mcr_precfrac .and. l_subgrid_graupel_frac ) then
  if ( i_update_precfrac == i_homog_areas ) then
    ! If using prognostic sub-grid fraction for rain and graupel,
    ! need to update the precip mass increments in each sub-grid
    ! partition.
    ! Assuming that all autoconversion of ice-cloud to graupel occurs
    ! within the mixed-phase cloud partition, where riming occurs.

    if ( l_orogrime ) then
      ! If orographic enhancement of riming is in use, then the area in
      ! which riming occurs maybe much larger than area_mix.  Need to
      ! account for the orographic mixed-phase cloud area

      do i = 1, points
        if ( l_calc(i) ) then
          if ( l_enh_rime(i) ) then
            ! If orographically enhanced riming occured at the current point

            tmp = one / max( area_mix_orog(i), small_number )

            ! Assign the increment between sub-partitions of the rain-fraction,
            ! assuming maximal overlap with first mixed-phase cloud and
            ! ice-only cloud (the 2 partitions where rain/graupel overlaps
            ! with ice), and finally create new precip area outside
            ! the existing rainfrac if it doesn't fit within rain_mix+rain_ice.
            frac = area_mix_orog(i)
            if ( frac > zero ) then
              dqprec_mix(i) = dqprec_mix(i) + dqi(i)*min(rain_mix(i),frac)*tmp
              frac = frac - rain_mix(i)
            end if
            if ( frac > zero ) then
              dqprec_ice(i) = dqprec_ice(i) + dqi(i)*min(rain_ice(i),frac)*tmp
              frac = frac - rain_ice(i)
            end if
            if ( frac > zero ) then
              dqprec_new(i) = dqprec_new(i) + dqi(i)*frac*tmp
              ! (this should only be updated where dqi(i)>0; this is a bug...)
              rain_new(i) = max( rain_new(i), frac )
            end if

          else
            ! No orographic enhancement at this point; area of riming is
            ! just area_mix (some of which will be rain_mix, and some of
            ! which will be new precip area in the mixed-phase cloud)

            tmp = one / max( area_mix(i), small_number)

            frac = area_mix(i)
            if ( frac > zero ) then
              dqprec_mix(i) = dqprec_mix(i) + dqi(i)*min(rain_mix(i),frac)*tmp
              frac = frac - rain_mix(i)
            end if
            if ( frac > zero ) then
              dqprec_new(i) = dqprec_new(i) + dqi(i)*frac*tmp
              ! (this should only be updated where dqi(i)>0; this is a bug...)
              rain_new(i) = max( rain_new(i), frac )
            end if

          end if  ! ( l_enh_rime(i) )
        end if  ! ( l_calc(i) )
      end do

    else  ! ( l_orogrime )
      ! Not using orographic riming enhancement, so all of the riming
      ! always occurs within area_mix.

      do i = 1, points
        if ( l_calc(i) ) then

          ! If part of mixed-phase area doesn't already contain rain/graupel,
          ! need to assign part of the precip mass increment to "new" area...

          ! Compute area of newly-created graupel
          area_new = max( zero, area_mix(i) - rain_mix(i) )

          ! Compute fraction of increment applied to the "new" area vs existing
          frac = area_new / max( area_mix(i), small_number )

          ! divvy the increment up between the two in proportion to area
          dqprec_new(i) = dqprec_new(i) +      frac  * dqi(i)
          dqprec_mix(i) = dqprec_mix(i) + (one-frac) * dqi(i)

          ! Update area of new precip passed back up
          ! (this should only be updated where dqi(i)>0; this is a bug...)
          rain_new(i) = max( rain_new(i), area_new )

        end if  ! ( l_calc(i) )
      end do

    end if  ! ( l_orogrime )

  else if ( i_update_precfrac == i_sg_correl ) then

    ! Calculate rain + graupel content before the increment was applied
    ! (ignoring negative values)
    do i = 1, points
      qprec(i) = max( qgraup(i) - dqi(i), zero )
      dqi(i) = max( qgraup(i), zero ) - qprec(i)
      qprec(i) = qprec(i) + max( qrain(i), zero )
    end do
    ! Set area fraction of the increment
    if ( l_orogrime ) then
      do i = 1, points
        if ( l_enh_rime(i) ) then
          area_inc(i) = area_mix_orog(i)
        else
          area_inc(i) = area_mix(i)
        end if
      end do
    else
      do i = 1, points
        area_inc(i) = area_mix(i)
      end do
    end if

    ! Calculate combined fraction of the existing precip and increment
    ! (precfrac_k is updated)
    call lsp_combine_precfrac( points,                                         &
                               qprec, dqi, precfrac_k, area_inc )

  end if  ! ( i_update_precfrac )
end if  ! ( l_mcr_precfrac .and. l_subgrid_graupel_frac )


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_graup_autoc
end module lsp_graup_autoc_mod
