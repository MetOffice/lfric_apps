! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_env_region_tq_nb_mod

implicit none

contains


! Subroutine to parameterise the temperature and water-vapour content
! of the 4 sub-grid regions of the grid-box:
! - clear-sky         (dry)
! - liquid-only cloud (liq)
! - mixed-phase cloud (mph)
! - ice and rain only (icr)
! Each region may differ from the grid-mean, such that the grid-mean
! T, qv are retrieved by averaging over all 4 regions:
! f_dry       + f_liq       + f_mph       + f_icr       = 1
! f_dry T_dry + f_liq T_liq + f_mph T_mph + f_icr T_icr = T_env
! f_dry q_dry + f_liq q_liq + f_mph q_mph + f_icr q_icr = q_env
!
! This routine assumes that:
! a) All 4 regions have equal virtual temperature, such that they are
!    neutrally buoyant.
! b) The liquid-cloud (and mixed-phase) regions are saturated
!    w.r.t. to liquid water.
! c) The subsaturated ice and rain region contains more vapour than the
!    mean, in proportion to its water-loading.
subroutine calc_env_region_tq_nb( n_points, n_points_super,                    &
                                  temperature, q_vap, cloudfracs,              &
                                  qc_tot, q_cond_loc,                          &
                                  qsat_liq, dqsatdt_liq, qsat_ice, dqsatdt_ice,&
                                  dtv_dt, dtv_dqv, dtv_dqc,                    &
                                  frac_r, temperature_r, q_vap_r )

use comorph_constants_mod, only: real_cvprec, zero, one,                       &
                                 min_delta, sqrt_min_delta, sqrt_min_float,    &
                                 n_cond_species, n_cond_species_liq,           &
                                 cond_params, i_sg_homog, i_sg_frac_liq,       &
                                 i_sg_frac_ice, i_sg_frac_prec,                &
                                 melt_temp
use subregion_mod, only: n_regions, i_dry, i_liq, i_mph, i_icr
use cloudfracs_type_mod, only: n_cloudfracs, i_frac_liq

use set_qsat_mod, only: set_qsat_liq

implicit none

! Number of points
integer, intent(in) :: n_points

! Size of input super-arrays (maybe larger than n_points,
! as arrays are dimensioned with the max size they could need
! on any level and reused for all levels)
integer, intent(in) :: n_points_super

! Grid-mean temperature and water-vapour mixing-ratio
real(kind=real_cvprec), intent(in) :: temperature(n_points)
real(kind=real_cvprec), intent(in) :: q_vap(n_points)

! Super-array containing:
! frac_liq    : Total liquid cloud fraction (including mixed-phase)
! frac_ice    : Total ice cloud fraction (including mixed-phase)
! frac_bulk   : Bulk cloud fraction
! frac_precip : Rain / graupel fraction
real(kind=real_cvprec), intent(in) :: cloudfracs                               &
                                      ( n_points_super, n_cloudfracs )

! Grid-mean total condensed water
real(kind=real_cvprec), intent(in) :: qc_tot(n_points)

! Local condensate mixing ratios within the regions in which
! each species is allowed to be non-zero
real(kind=real_cvprec), intent(in) :: q_cond_loc ( n_points, n_cond_species )

! Saturation water-vapour mixing-ratio w.r.t. liquid at the
! grid-mean temperature, and its gradient with temperature
real(kind=real_cvprec), intent(in) :: qsat_liq(n_points)
real(kind=real_cvprec), intent(in) :: dqsatdt_liq(n_points)
! Values w.r.t. ice (only valid where below freezing)
real(kind=real_cvprec), intent(in) :: qsat_ice(n_points)
real(kind=real_cvprec), intent(in) :: dqsatdt_ice(n_points)

! Derivatives of virtual temperature w.r.t. temperature,
! water-vapour mixing-ratio, and condensed water mixing-ratio
real(kind=real_cvprec), intent(in) :: dtv_dt(n_points)
real(kind=real_cvprec), intent(in) :: dtv_dqv(n_points)
real(kind=real_cvprec), intent(in) :: dtv_dqc(n_points)

! Fractional area of each sub-region
real(kind=real_cvprec), intent(in) :: frac_r ( n_points, n_regions )

! Temperature and water-vapour in each sub-region
real(kind=real_cvprec), intent(in out) :: temperature_r ( n_points, n_regions )
! Water vapour mixing ratio in each sub-region
real(kind=real_cvprec), intent(in out) :: q_vap_r       ( n_points, n_regions )


! Re-occuring ratio dTv/dT / dqsat/dT
real(kind=real_cvprec) :: dtv_dqsat_liq(n_points)
real(kind=real_cvprec) :: dtv_dqsat_ice(n_points)

! Local total condensed water within each region
real(kind=real_cvprec) :: qc_tot_loc(n_points,n_regions)
! Local total condensed water outside the liquid-cloud
real(kind=real_cvprec) :: qc_tot_noliq

! Local total ice mixing-ratio within the mph,icr regions
real(kind=real_cvprec) :: q_ice_loc(n_points)

! Target supersaturation for cloudy region (usually 0)
real(kind=real_cvprec) :: supersat(n_points)

! Temperature and vapour in the non-liquid-cloud region
! (dry and icr regions combined)
real(kind=real_cvprec) :: temperature_noliq(n_points)
real(kind=real_cvprec) :: q_vap_noliq(n_points)

! Total-condensate in the icr region (excess relative to dry region)
real(kind=real_cvprec) :: qc_icr

! Temporary work variable
real(kind=real_cvprec) :: tmp

! Fraction of condensate that is ice in icr region
real(kind=real_cvprec) :: ice_frac_icr

! Max or min limit applied to q_vap for safety
real(kind=real_cvprec) :: q_lim

! Indices of sub-set of points to do certain calculations at
integer :: nc
integer :: index_ic(n_points)

! Loop counters
integer :: ic, ic2, i_cond, i_region


!------------------------------------------------------------------------------
! 1) Set work arrays needed for the calculations
!------------------------------------------------------------------------------

! Initialise total-condensed-water in each region to zero
do i_region = 1, n_regions
  do ic = 1, n_points
    qc_tot_loc(ic,i_region) = zero
  end do
end do

! Find the total condensed water in each region
! (needed for the water-loading term to calculate neutral buoyancy)
do i_cond = 1, n_cond_species
  ! Which sub-grid region does this species belong to?
  select case ( cond_params(i_cond)%pt % i_sg )
  case (i_sg_homog)
    ! Species contributes to total-condensate in all 4 regions
    do ic = 1, n_points
      qc_tot_loc(ic,i_dry) = qc_tot_loc(ic,i_dry) + q_cond_loc(ic,i_cond)
      qc_tot_loc(ic,i_liq) = qc_tot_loc(ic,i_liq) + q_cond_loc(ic,i_cond)
      qc_tot_loc(ic,i_mph) = qc_tot_loc(ic,i_mph) + q_cond_loc(ic,i_cond)
      qc_tot_loc(ic,i_icr) = qc_tot_loc(ic,i_icr) + q_cond_loc(ic,i_cond)
    end do
  case (i_sg_frac_liq)
    ! Species contributes to total-condensate in liq and mph regions
    do ic = 1, n_points
      qc_tot_loc(ic,i_liq) = qc_tot_loc(ic,i_liq) + q_cond_loc(ic,i_cond)
      qc_tot_loc(ic,i_mph) = qc_tot_loc(ic,i_mph) + q_cond_loc(ic,i_cond)
    end do
  case (i_sg_frac_ice)
    ! Species contributes to total-condensate in mph and icr regions
    do ic = 1, n_points
      qc_tot_loc(ic,i_mph) = qc_tot_loc(ic,i_mph) + q_cond_loc(ic,i_cond)
      qc_tot_loc(ic,i_icr) = qc_tot_loc(ic,i_icr) + q_cond_loc(ic,i_cond)
    end do
  case (i_sg_frac_prec)
    ! Species contributes to total-condensate in liq, mph and icr regions
    do ic = 1, n_points
      qc_tot_loc(ic,i_liq) = qc_tot_loc(ic,i_liq) + q_cond_loc(ic,i_cond)
      qc_tot_loc(ic,i_mph) = qc_tot_loc(ic,i_mph) + q_cond_loc(ic,i_cond)
      qc_tot_loc(ic,i_icr) = qc_tot_loc(ic,i_icr) + q_cond_loc(ic,i_cond)
    end do
  end select
end do  ! i_cond = 1, n_cond_species

! Set frequently used ratio dTv/dT / dqsat/dT
do ic = 1, n_points
  dtv_dqsat_liq(ic) = dtv_dt(ic) / dqsatdt_liq(ic)
end do
if ( any( temperature < melt_temp ) ) then
  do ic = 1, n_points
    dtv_dqsat_ice(ic) = dtv_dt(ic) / dqsatdt_ice(ic)
  end do
end if


!------------------------------------------------------------------------------
! 2) Set liquid-cloud q_vap to saturated and neutrally-buoyant
!------------------------------------------------------------------------------

! Set target supersaturation for the liquid cloud region.
! This is normally zero (i.e. saturation).  However,
! in the event that the grid-mean state is supersaturated,
! we avoid making the cloud-free region more supersaturated
! than the cloudy region by resetting to the grid-mean
! supersaturation.
do ic = 1, n_points
  supersat(ic) = max( zero, q_vap(ic) - qsat_liq(ic) )
end do

! Find points with non-zero fraction for liq region
nc = 0
do ic = 1, n_points
  if ( frac_r(ic,i_liq) > min_delta ) then
    nc = nc + 1
    index_ic(nc) = ic
  end if
end do

if ( nc > 0 ) then
  ! If any liq points
  ! Calculate properties of the liquid cloud only region...

  do ic2 = 1, nc
    ic = index_ic(ic2)
    ! Calculate water-vapour mixing ratio; we have:
    !
    ! Must have saturation:
    ! qv_loc = ( q_sat(T_mean) + dqsat/dT (T_loc - T_mean) ) + supersat
    !
    ! Must have specified buoyancy:
    ! (T_loc - T_mean) dTv/dT + (qv_loc - qv_mean) dTv/dqv
    !                         + (qc_loc - qc_mean) dTv/dqc = 0
    !
    ! Combining to eliminate (T_loc - T_mean) and rearranging,
    ! we have:
    ! ( qv_loc - q_sat(T_mean) - supersat ) dTv/dT / dqsat/dT
    !                         + (qv_loc - qv_mean) dTv/dqv
    !                         + (qc_loc - qc_mean) dTv/dqc = 0
    ! => qv_loc ( dTv/dT / dqsat/dT + dTv/dqv )
    !  - ( q_sat(T_mean) + supersat ) dTv/dT / dqsat/dT
    !  - qv_mean dTv/dqv + (qc_loc - qc_mean) dTv/dqc = 0
    ! => qv_loc = ( ( q_sat(T_mean) + supersat ) dTv/dT / dqsat/dT
    !             + qv_mean dTv/dqv
    !             - (qc_loc - qc_mean) dTv/dqc
    !             ) / ( dTv/dT / dqsat/dT + dTv/dqv )
    q_vap_r(ic,i_liq) =                                                        &
         ( ( qsat_liq(ic) + supersat(ic) ) * dtv_dqsat_liq(ic)                 &
         + q_vap(ic) * dtv_dqv(ic)                                             &
         - ( qc_tot_loc(ic,i_liq) - qc_tot(ic) ) * dtv_dqc(ic) )               &
       / ( dtv_dqsat_liq(ic) + dtv_dqv(ic) )
  end do

end if  ! ( nc > 0 ) THEN


!------------------------------------------------------------------------------
! 3) Repeat the above for mixed-phase cloud (may have different water-loading)
!------------------------------------------------------------------------------

! Find points with non-zero fraction for mph region
nc = 0
do ic = 1, n_points
  if ( frac_r(ic,i_mph) > min_delta ) then
    nc = nc + 1
    index_ic(nc) = ic
  end if
end do

if ( nc > 0 ) then
  ! If any mph points
  ! Calculate properties of the mixed-phase region...

  do ic2 = 1, nc
    ic = index_ic(ic2)
    ! Calculate q_vap as above
    q_vap_r(ic,i_mph) =                                                        &
         ( ( qsat_liq(ic) + supersat(ic) ) * dtv_dqsat_liq(ic)                 &
         + q_vap(ic) * dtv_dqv(ic)                                             &
         - ( qc_tot_loc(ic,i_mph) - qc_tot(ic) ) * dtv_dqc(ic) )               &
       / ( dtv_dqsat_liq(ic) + dtv_dqv(ic) )
  end do

end if  ! ( nc > 0 )


!------------------------------------------------------------------------------
! 4) Find mean properties of air outside the liquid cloud,
!    so-as to get the right grid-mean
!------------------------------------------------------------------------------

! Initialise non-liquid-cloud region T,q to grid-mean
do ic = 1, n_points
  temperature_noliq(ic) = temperature(ic)
  q_vap_noliq(ic)       = q_vap(ic)
end do

! Find points where there is partial liquid cloud cover
nc = 0
do ic = 1, n_points
  if ( ( frac_r(ic,i_liq)>min_delta .or. frac_r(ic,i_mph)>min_delta ) .and.    &
       ( frac_r(ic,i_dry)>min_delta .or. frac_r(ic,i_icr)>min_delta ) ) then
    nc = nc + 1
    index_ic(nc) = ic
  end if
end do

! If any points
if ( nc > 0 ) then

  do ic2 = 1, nc
    ic = index_ic(ic2)
    ! We have:
    ! f_liq q_liq + f_mph q_mph + f_noliq q_noliq = q_mean
    ! Therefore:
    ! q_noliq = ( q_mean - f_liq q_liq - f_mph q_mph ) / f_noliq
    !
    ! To improve numerical precision, we rearrange this to:
    ! q_noliq = q_mean - ( f_liq(q_liq-q_mean)
    !                    + f_mph(q_mph-q_mean) ) / f_noliq
    !
    ! (Note: sub-partition of the noliq area's properties
    !  calculated here into the icr and dry regions will be
    !  done subsequently...)
    !
    ! (note frac_liq = f_liq + f_mph)
    ! Find q_vap outside the liquid cloud
    q_vap_noliq(ic) = q_vap(ic)                                                &
          - ( frac_r(ic,i_liq) * ( q_vap_r(ic,i_liq) - q_vap(ic) )             &
            + frac_r(ic,i_mph) * ( q_vap_r(ic,i_mph) - q_vap(ic) ) )           &
          / ( one - cloudfracs(ic,i_frac_liq) )
  end do

  ! Safety checks...
  do ic2 = 1, nc
    ic = index_ic(ic2)
    ! If making the liquid cloud region saturated has implied
    ! negative vapour content in the non-liquid cloud region:
    if ( q_vap_noliq(ic) < zero ) then
      ! Fix problem by reducing q_vap in the liquid-cloud regions
      ! and recalculating.  We want to rescale the moisture in
      ! the liq and mph regions such that:
      ! 0 = qv - fac ( f_liq qv_liq + f_mph qv_mph )
      ! => fac = qv / ( f_liq qv_liq + f_mph qv_mph )
      tmp = q_vap(ic) / ( frac_r(ic,i_liq) * q_vap_r(ic,i_liq)                 &
                        + frac_r(ic,i_mph) * q_vap_r(ic,i_mph) )
      ! Rescale liq region
      if ( frac_r(ic,i_liq) > min_delta ) then
        q_vap_r(ic,i_liq) = q_vap_r(ic,i_liq) * tmp
      end if
      ! Rescale mph region
      if ( frac_r(ic,i_mph) > min_delta ) then
        q_vap_r(ic,i_mph) = q_vap_r(ic,i_mph) * tmp
      end if
      ! Reset non-liquid cloud region q_vap to zero
      q_vap_noliq(ic) = zero
    end if
  end do

  ! Find temperature of the non-liquid-cloud region, so-as to have
  ! same virtual temperature as the grid-mean
  do ic2 = 1, nc
    ic = index_ic(ic2)
    ! Find local total-condensed-water in the no-liquid cloud region
    qc_tot_noliq = qc_tot(ic)                                                  &
          - ( frac_r(ic,i_liq) * ( qc_tot_loc(ic,i_liq) - qc_tot(ic) )         &
            + frac_r(ic,i_mph) * ( qc_tot_loc(ic,i_mph) - qc_tot(ic) ) )       &
          / ( one - cloudfracs(ic,i_frac_liq) )
    ! Set T for specified Tv
    temperature_noliq(ic) = temperature(ic)                                    &
     - ( ( q_vap_noliq(ic) - q_vap(ic) ) * dtv_dqv(ic)                         &
       + ( qc_tot_noliq - qc_tot(ic) ) * dtv_dqc(ic) ) / dtv_dt(ic)
  end do

end if  ! ( nc > 0 )


!------------------------------------------------------------------------------
! 5) Parameterise the difference in water-vapour mixing ratio
!    between the dry and icr regions...
!------------------------------------------------------------------------------

! Find total ice in the icr region
do ic = 1, n_points
  q_ice_loc(ic) = zero
end do
do i_cond = n_cond_species_liq+1, n_cond_species
  ! Sum local values for all ice species
  ! (all ice species are assumed to be present in the icr region).
  do ic = 1, n_points
    q_ice_loc(ic) = q_ice_loc(ic) + q_cond_loc(ic,i_cond)
  end do
end do

do ic = 1, n_points

  ! If the icr and dry regions both have non-zero fraction...
  if ( frac_r(ic,i_icr) > min_delta .and.                                      &
       frac_r(ic,i_dry) > min_delta ) then

    ! Initialise vapour in the icr region to
    ! liquid-saturated value with specified buoyancy
    q_vap_r(ic,i_icr) =                                                        &
         ( qsat_liq(ic) * dtv_dqsat_liq(ic)                                    &
         + q_vap(ic) * dtv_dqv(ic)                                             &
         - ( qc_tot_loc(ic,i_icr) - qc_tot(ic) ) * dtv_dqc(ic) )               &
       / ( dtv_dqsat_liq(ic) + dtv_dqv(ic) )

    ! If below freezing and ice is present:
    ! (note: need both the grid-mean T and local T to be below
    !  0oC for the linearisation of qsat_ice to work correctly)
    if ( temperature_noliq(ic) < melt_temp                                     &
         .and. temperature(ic) < melt_temp                                     &
         .and. q_ice_loc(ic) > zero ) then
      ! Initialise vapour in the icr region to specified buoyancy with
      ! ice-saturated value (weighted by fraction of condensate
      ! that is ice)
      ice_frac_icr = q_ice_loc(ic) / qc_tot_loc(ic,i_icr)
      ! ice_frac_icr stores fraction of that which is ice
      q_vap_r(ic,i_icr) = (one-ice_frac_icr) * q_vap_r(ic,i_icr)               &
                        +      ice_frac_icr  * (                               &
           ( qsat_ice(ic) * dtv_dqsat_ice(ic)                                  &
           + q_vap(ic) * dtv_dqv(ic)                                           &
           - ( qc_tot_loc(ic,i_icr) - qc_tot(ic) ) * dtv_dqc(ic) )             &
         / ( dtv_dqsat_ice(ic) + dtv_dqv(ic) )   )
    end if

    ! Find excess of icr region condensate above the dry region condensate
    qc_icr = qc_tot_loc(ic,i_icr) - qc_tot_loc(ic,i_dry)

    ! If mean q_vap in the no-liquid-cloud region
    ! (stored in qv_noliq) is greater than the saturated value
    ! (stored in qv_icr)
    if ( q_vap_noliq(ic) > q_vap_r(ic,i_icr) ) then
      ! Where the no-liquid-cloud environment is supersaturated,
      ! assume the condensate has grown at the expense of the
      ! vapour, so that qv_icr + qc_icr = qv_dry.
      ! Now:
      ! f_icr qv_icr + f_dry qv_dry = (f_icr + f_dry) qv_noliq
      ! Substituting:
      ! qv_dry = qv_icr + qc_icr
      ! we get:
      ! f_icr qv_icr + f_dry (qv_icr + qc_icr)
      !              = (f_icr + f_dry) qv_noliq
      ! =>
      ! qv_icr (f_icr + f_dry) = -f_dry qc_icr
      !                        + (f_icr + f_dry) qv_noliq
      ! =>
      ! qv_icr = qv_noliq - f_dry/(f_icr + f_dry) qc_icr
      q_vap_r(ic,i_icr) = max( q_vap_r(ic,i_icr),                              &
           q_vap_noliq(ic) - qc_icr * frac_r(ic,i_dry)                         &
                   / ( frac_r(ic,i_icr) + frac_r(ic,i_dry) ) )
      ! Note safety check to stop q_vap going below its
      ! initialised saturated value.
    else
      ! Where the no-liquid-cloud environment is subsaturated,
      ! assume:  qv_icr = qv_dry + qc_icr
      ! i.e. the icr region is moister than the dry region
      ! by an amount equal to its condensed water content.
      ! There is no particular justification for this precise
      ! formula, but it should have a sensible general effect
      ! of making the icr region's moisture excess increase
      ! when there is more rain / ice in it.
      !
      ! We end up with the same formula as above, but with one
      ! sign change:
      q_vap_r(ic,i_icr) = min( q_vap_r(ic,i_icr),                              &
           q_vap_noliq(ic) + qc_icr * frac_r(ic,i_dry)                         &
                   / ( frac_r(ic,i_icr) + frac_r(ic,i_dry) ) )
      ! Note: now the safety check stops q_vap going above its
      ! initialised saturated value.
    end if

    ! Apply bounds to qv_icr...

    ! Upper bound is liquid saturation
    q_lim = ( ( qsat_liq(ic) + supersat(ic) ) * dtv_dqsat_liq(ic)              &
            + q_vap(ic) * dtv_dqv(ic)                                          &
            - ( qc_tot_loc(ic,i_icr) - qc_tot(ic) ) * dtv_dqc(ic) )            &
          / ( dtv_dqsat_liq(ic) + dtv_dqv(ic) )
    q_vap_r(ic,i_icr) = min( q_vap_r(ic,i_icr), q_lim )

    ! Upper bound for dry region is also liquid saturation
    q_lim = ( ( qsat_liq(ic) + supersat(ic) ) * dtv_dqsat_liq(ic)              &
            + q_vap(ic) * dtv_dqv(ic)                                          &
            - ( qc_tot_loc(ic,i_dry) - qc_tot(ic) ) * dtv_dqc(ic) )            &
          / ( dtv_dqsat_liq(ic) + dtv_dqv(ic) )
    ! This translates into a lower bound for the icr region:
    ! f_dry qv_dry + f_icr qv_icr = (f_dry + f_icr) qv_noliq
    ! => qv_icr = ( (f_dry + f_icr) qv_noliq - f_dry qv_dry )
    !           / f_icr
    ! => qv_icr = qv_noliq + f_dry/f_icr ( qv_noliq - qv_dry )
    q_lim = q_vap_noliq(ic)                                                    &
          + ( frac_r(ic,i_dry) / frac_r(ic,i_icr) )                            &
          * ( q_vap_noliq(ic) - q_lim )
    q_vap_r(ic,i_icr) = max( q_vap_r(ic,i_icr), q_lim )

    ! Finally, we must have qv_icr > 0
    q_vap_r(ic,i_icr) = max( q_vap_r(ic,i_icr), zero )

    ! Compute dry region q_vap to get right grid-mean
    ! f_dry qv_dry + f_icr qv_icr = (f_dry + f_icr) qv_noliq
    ! => qv_dry = ( (f_dry + f_icr) qv_noliq - f_icr qv_icr )
    !           / f_dry
    ! => qv_dry = qv_noliq + f_icr/f_dry ( qv_noliq - qv_icr )
    q_vap_r(ic,i_dry) = q_vap_noliq(ic)                                        &
         + ( frac_r(ic,i_icr) / frac_r(ic,i_dry) )                             &
         * ( q_vap_noliq(ic) - q_vap_r(ic,i_icr) )

    ! If q_vap_dry is negative, reset to zero and recalculate
    ! q_vap_icr to get right grid-mean
    if ( q_vap_r(ic,i_dry) < zero ) then
      ! We will now have:
      ! f_dry * 0 + f_icr qv_icr = (f_dry+f_icr) qv_noliq
      ! => qv_icr = ( 1 + f_dry/f_icr ) qv_noliq
      q_vap_r(ic,i_icr) = ( one + ( frac_r(ic,i_dry) / frac_r(ic,i_icr) ) )    &
                        * q_vap_noliq(ic)
      q_vap_r(ic,i_dry) = zero
    end if

  else  ! ( frac_r(ic,i_icr) > min_delta .AND.                                 &
        ! ( frac_r(ic,i_dry) > min_delta )

    if ( frac_r(ic,i_icr) > min_delta ) then
      ! If icr region has non-zero fraction but dry region is
      ! zero, just set icr region equal to mean outside the
      ! liquid cloud region
      q_vap_r(ic,i_icr) = q_vap_noliq(ic)
    end if

    if ( frac_r(ic,i_dry) > min_delta ) then
      ! If dry region has non-zero fraction but icr region is
      ! zero, just set dry region equal to mean outside the
      ! liquid cloud region
      q_vap_r(ic,i_dry) = q_vap_noliq(ic)
    end if

  end if
end do


!------------------------------------------------------------------------------
! 6) Set the temperatures of all 4 regions so-as to yield neutral buoyancy
!------------------------------------------------------------------------------
! Note: this should also yield region temperatures which average
! to the correct grid-mean temperature, because we've assumed
! Tv is a linear function of q_vap and q_cond, which should
! average up correctly.

do i_region = 1, n_regions
  ! Find points with non-zero fraction for current region
  nc = 0
  do ic = 1, n_points
    if ( frac_r(ic,i_region) > min_delta ) then
      nc = nc + 1
      index_ic(nc) = ic
    end if
  end do
  if ( nc > 0 ) then
    ! If any points...
    do ic2 = 1, nc
      ic = index_ic(ic2)

      ! Calculate T so-as to yield the correct buoyancy
      temperature_r(ic,i_region) = temperature(ic)                             &
      - ( ( q_vap_r(ic,i_region) - q_vap(ic) ) * dtv_dqv(ic)                   &
        + ( qc_tot_loc(ic,i_region) - qc_tot(ic) ) * dtv_dqc(ic) ) / dtv_dt(ic)

    end do
  end if
end do


return
end subroutine calc_env_region_tq_nb

end module calc_env_region_tq_nb_mod
