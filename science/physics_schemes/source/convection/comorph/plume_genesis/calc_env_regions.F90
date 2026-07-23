! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_env_regions_mod

implicit none

contains


! Subroutine to calculate the properties of the clear, cloudy
! and rainy sub-regions of the grid-box.  The 4 sub-regions calculated are:
! - clear-sky         (dry)
! - liquid-only cloud (liq)
! - mixed-phase cloud (mph)
! - ice and rain only (icr)
!
! The local values of T, q_vap and condensed water mixing ratios
! are calculated in each region so-as to satisfy some chosen
! assumptions, and the conservation laws applying to the
! fractional area f of each region:
! f_dry       + f_liq       + f_mph       + f_icr       = 1
! f_dry T_dry + f_liq T_liq + f_mph T_mph + f_icr T_icr = T_env
! f_dry q_dry + f_liq q_liq + f_mph q_mph + f_icr q_icr = q_env
!
! This routine also outputs non-turbulent parcel temperature and moisture
! perturbations, since their calculation requires some of the same
! factors as are used to compute the sub-region properties...

subroutine calc_env_regions( n_points, n_points_super,                         &
                             cmpr, k, call_string,                             &
                             pressure, cloudfracs,                             &
                             temperature, q_vap, q_cond,                       &
                             frac_r, temperature_r, q_vap_r,                   &
                             q_cond_loc,                                       &
                             genesis_diags, diags_super,                       &
                             delta_temp_neut, delta_qvap_neut )

use comorph_constants_mod, only: real_cvprec, zero, one,                       &
                                 min_delta, sqrt_min_float, n_cond_species,    &
                                 cond_params, i_sg_homog, i_sg_frac_liq,       &
                                 i_sg_frac_ice, i_sg_frac_prec,                &
                                 par_gen_rhpert, melt_temp, R_dry, R_vap,      &
                                 name_length, i_check_bad_values_cmpr,         &
                                 i_check_bad_none
use subregion_mod, only: n_regions, i_dry, i_liq, i_mph, i_icr,                &
                         region_names
use cloudfracs_type_mod, only: n_cloudfracs,                                   &
                               i_frac_liq, i_frac_ice,                         &
                               i_frac_bulk, i_frac_precip
use cmpr_type_mod, only: cmpr_type
use check_bad_values_mod, only: check_bad_values_cmpr
use fields_type_mod, only: field_names, n_fields,                              &
                           i_temperature, i_q_vap, i_qc_first
use genesis_diags_type_mod, only: genesis_diags_type

use calc_virt_temp_mod, only: calc_virt_temp
use set_qsat_mod, only: set_qsat_liq, set_qsat_ice
use set_dqsatdt_mod, only: set_dqsatdt_liq, set_dqsatdt_ice
use calc_env_region_tq_nb_mod, only: calc_env_region_tq_nb
use set_region_cond_fields_mod, only: set_region_cond_fields
use fields_diags_type_mod, only: fields_diags_copy

implicit none

! Number of points
integer, intent(in) :: n_points

! Size of input super-arrays (maybe larger than n_points,
! as arrays are dimensioned with the max size they could need
! on any level and reused for all levels)
integer, intent(in) :: n_points_super

! Stuff used to make error messages more informative:
! Structure storing compression indices of the current points
type(cmpr_type), intent(in) :: cmpr
! Current model-level index
integer, intent(in) :: k
! String identifying where this routine has been called
character(len=name_length), intent(in) :: call_string

! Ambient pressure (assumed homogeneous)
real(kind=real_cvprec), intent(in) :: pressure(n_points)

! Super-array containing:
! frac_liq    : Total liquid cloud fraction (including mixed-phase)
! frac_ice    : Total ice cloud fraction (including mixed-phase)
! frac_bulk   : Bulk cloud fraction
! frac_precip : Rain / graupel fraction
real(kind=real_cvprec), intent(in) :: cloudfracs                               &
                                ( n_points_super, n_cloudfracs )

! Grid-mean temperature
real(kind=real_cvprec), intent(in) :: temperature(n_points)
! Grid-mean water-vapour mixing-ratio
real(kind=real_cvprec), intent(in) :: q_vap(n_points)
! Super-array containing grid-means of condensed water species
real(kind=real_cvprec), intent(in) :: q_cond                                   &
                              ( n_points_super, n_cond_species )

! Fractional area of each sub-region
real(kind=real_cvprec), intent(out) :: frac_r                                  &
                                       ( n_points, n_regions )
! Temperature in each sub-region
real(kind=real_cvprec), intent(out) :: temperature_r                           &
                                       ( n_points, n_regions )
! Water vapour mixing ratio in each sub-region
real(kind=real_cvprec), intent(out) :: q_vap_r                                 &
                                       ( n_points, n_regions )

! Local condensate mixing ratios within the regions in which
! each species is allowed to be non-zero
real(kind=real_cvprec), intent(out) :: q_cond_loc                              &
                                     ( n_points, n_cond_species )

! Structure storing diagnostics and associated meta-data
type(genesis_diags_type), intent(in out) :: genesis_diags
! Super-array storing diagnostics to be output
real(kind=real_cvprec), intent(in out) :: diags_super                          &
                          ( n_points_super, genesis_diags % n_diags )

! Neutrally buoyant perturbations to T,q needed for a fractional increase of RH
real(kind=real_cvprec), intent(out) :: delta_temp_neut(n_points)
real(kind=real_cvprec), intent(out) :: delta_qvap_neut(n_points)

! Saturation water-vapour mixing-ratio w.r.t. liquid at the
! grid-mean temperature, and its gradient with temperature
real(kind=real_cvprec) :: qsat_liq(n_points)
real(kind=real_cvprec) :: dqsatdt_liq(n_points)
! Values w.r.t. ice (only valid where below freezing)
real(kind=real_cvprec) :: qsat_ice(n_points)
real(kind=real_cvprec) :: dqsatdt_ice(n_points)

! Virtual temperature (assumed to be equal in all regions)
real(kind=real_cvprec) :: virt_temp(n_points)

! Derivatives of virtual temperature w.r.t. temperature,
! water-vapour mixing-ratio, and condensed water mixing-ratio
real(kind=real_cvprec) :: dtv_dt(n_points)
real(kind=real_cvprec) :: dtv_dqv(n_points)
real(kind=real_cvprec) :: dtv_dqc(n_points)

! Grid-mean total condensed water
real(kind=real_cvprec) :: qc_tot(n_points)

! Temporary work variable
real(kind=real_cvprec) :: tmp

! Fields super-array only used for copying diagnostics if requested
real(kind=real_cvprec), allocatable :: fields_reg(:,:)

! Flag passed into fields_diags_copy (always false  here)
logical, parameter :: l_conserved_form_false = .false.

! Indices of sub-set of points to do certain calculations at
integer :: nc
integer :: index_ic(n_points)

! Description of where we are in the code, for error messages
character(len=name_length) :: where_string
character(len=name_length) :: field_name
logical :: field_positive

! Loop counter
integer :: ic, ic2, i_cond, i_region, i_field, i_super


!------------------------------------------------------------------------------
! 1) Set work arrays needed for the calculations
!------------------------------------------------------------------------------

! Calculate saturation water-vapour mixing-ratio and its
! gradient with temperature
call set_qsat_liq( n_points, temperature, pressure, qsat_liq )
call set_dqsatdt_liq( n_points, temperature, qsat_liq, dqsatdt_liq )
if ( any( temperature < melt_temp ) ) then
  call set_qsat_ice( n_points, temperature, pressure, qsat_ice )
  call set_dqsatdt_ice( n_points, temperature, qsat_ice, dqsatdt_ice )
end if

! Calculate grid-mean virtual temperature
call calc_virt_temp( n_points, n_points_super,                                 &
                     temperature, q_vap, q_cond, virt_temp )

! Calculate grid-mean total condensed water
do ic = 1, n_points
  qc_tot(ic) = zero
end do
do i_cond = 1, n_cond_species
  do ic = 1, n_points
    qc_tot(ic) = qc_tot(ic) + q_cond(ic,i_cond)
  end do
end do

! Calculate gradients of virtual temperature with T, q_vap, qc
do ic = 1, n_points
  ! Store 1 + q_tot
  tmp = one + q_vap(ic) + qc_tot(ic)
  ! dTv/dT = ( 1 + Rv/Rd q_vap ) / ( 1 + q_tot ) = Tv/T
  dtv_dt(ic) = virt_temp(ic) / temperature(ic)
  ! dTv/dqv = ( (Rv/Rd) T - Tv ) / ( 1 + q_tot )
  dtv_dqv(ic) = ( (R_vap/R_dry) * temperature(ic) - virt_temp(ic) ) / tmp
  ! dTv/dqc = -Tv / ( one + q_tot(ic) )
  dtv_dqc(ic) = -virt_temp(ic) / tmp
end do


!------------------------------------------------------------------------------
! 2) Find fractional areas of the 4 partitions...
!------------------------------------------------------------------------------

do ic = 1, n_points
  ! Mixed-phase cloud fraction
  ! = overlap between cf_liq and cf_ice
  ! = cf_liq + cf_ice - cf_bulk
  ! To avoid cf_mph spuriously falling to zero due to rounding-error,
  ! write this as either:
  !   cf_mph = cf_liq - ( cf_bulk - cf_ice )
  ! or
  !   cf_mph = cf_ice - ( cf_bulk - cf_liq )
  ! with the smaller of cf_liq, cf_ice outside the brackets.
  ! This ensures that if the larger one equals cf_bulk, it correctly cancels
  ! it leaving cf_mph equal to the smaller one.  Otherwise the smaller one
  ! can get wiped-out by rounding-error when added to the larger one,
  ! wrongly leaving zero residual when cf_bulk is subtracted.
  ! This formula also avoids cf_mph > cf_liq or cf_ice due to rounding error.
  frac_r(ic,i_mph) = max(                                                      &
        min( cloudfracs(ic,i_frac_liq), cloudfracs(ic,i_frac_ice) )            &
    - ( cloudfracs(ic,i_frac_bulk)                                             &
      - max( cloudfracs(ic,i_frac_liq), cloudfracs(ic,i_frac_ice) ) ),         &
                          zero )
  ! Still need to limit to avoid negative values due to rounding-error.

  ! Liquid-only cloud fraction
  ! = part of liquid cloud that is not mixed-phase
  ! = cf_liq - f_mph
  frac_r(ic,i_liq) = cloudfracs(ic,i_frac_liq) - frac_r(ic,i_mph)

  ! Ice / rain but no liquid cloud fraction
  ! = area containing ice or rain but no liquid cloud
  ! = MAX( cf_ice - cf_mph, cf_precip - cf_liq )
  frac_r(ic,i_icr) = max(                                                      &
      cloudfracs(ic,i_frac_ice) - frac_r(ic,i_mph),                            &
      cloudfracs(ic,i_frac_precip) - cloudfracs(ic,i_frac_liq)                 &
                        )

  ! Clear-sky fraction = any space remaining
  frac_r(ic,i_dry) = one - ( frac_r(ic,i_liq) + frac_r(ic,i_mph)               &
                           + frac_r(ic,i_icr) )
end do

do ic = 1, n_points
  ! Remove spurious not-quite-zero values of frac_dry due to rounding-errors
  ! (frac_dry should be exactly 0 if either cloud or precip fraction is 1)
  if ( .not. ( cloudfracs(ic,i_frac_bulk) < one .and.                          &
               cloudfracs(ic,i_frac_precip) < one ) ) then
    frac_r(ic,i_dry) = zero
  end if
end do


!------------------------------------------------------------------------------
! 3) Find local condensate mixing-rataios within the partitions...
!------------------------------------------------------------------------------

! Loop over condensed water species
do i_cond = 1, n_cond_species
  ! Which sub-grid region does this species belong to?
  select case ( cond_params(i_cond)%pt % i_sg )
  case (i_sg_homog)
    ! Species is homogeneously distributed over the whole grid-box
    do ic = 1, n_points
      ! Local value equals grid-mean
      q_cond_loc(ic,i_cond) = q_cond(ic,i_cond)
    end do
  case (i_sg_frac_liq)
    ! Species is only present within the liquid-cloud fraction
    do ic = 1, n_points
      ! This corresponds to both the liq and mph regions
      q_cond_loc(ic,i_cond) = q_cond(ic,i_cond)                                &
        / max( frac_r(ic,i_liq) + frac_r(ic,i_mph), sqrt_min_float )
    end do
  case (i_sg_frac_ice)
    ! Species is only present within the ice-cloud fraction
    do ic = 1, n_points
      ! This corresponds to both the mph and icr regions
      q_cond_loc(ic,i_cond) = q_cond(ic,i_cond)                                &
        / max( frac_r(ic,i_mph) + frac_r(ic,i_icr), sqrt_min_float )
    end do
  case (i_sg_frac_prec)
    ! Species is only present within the precipitation fraction
    do ic = 1, n_points
      ! This corresponds to the liq, mph and icr regions
      ! (precip assumed to be maximally overlapped with cloud)
      q_cond_loc(ic,i_cond) = q_cond(ic,i_cond)                                &
        / max( frac_r(ic,i_liq) + frac_r(ic,i_mph)                             &
                                + frac_r(ic,i_icr), sqrt_min_float )
    end do
  end select
end do  ! i_cond = 1, n_cond_species


!------------------------------------------------------------------------------
! 4) Parameterise T, q_vap within each region...
!------------------------------------------------------------------------------

! Initialise output T,q of regions to grid-box means
do i_region = 1, n_regions
  do ic = 1, n_points
    temperature_r(ic,i_region) = temperature(ic)
    q_vap_r(ic,i_region)       = q_vap(ic)
  end do
end do

! Assume the 4 sub-regions are neutrally buoyant (equal virtual temperature)
! TEMPORARY COMPILER DIRECTIVE TO FORCE CCE-FAST-DEBUG TO PRESERVE KGO:
!DIR$ INLINE
call calc_env_region_tq_nb( n_points, n_points_super,                          &
                            temperature, q_vap, cloudfracs,                    &
                            qc_tot, q_cond_loc,                                &
                            qsat_liq, dqsatdt_liq, qsat_ice, dqsatdt_ice,      &
                            dtv_dt, dtv_dqv, dtv_dqc,                          &
                            frac_r, temperature_r, q_vap_r )
!DIR$ RESETINLINE


!------------------------------------------------------------------------------
! 5) Calculate neutrally-buoyant T,q perturbations
!------------------------------------------------------------------------------

do ic = 1, n_points

  ! We want T and q perturbations such that RH is increased
  ! by a fixed fractional perturbation RHpert:
  ! (q + dq) / qsat(T+dT) = q/qsat(T) (1 + RHpert)
  !
  ! Linearising qsat:
  ! => q + dq = ( qsat(T) + dqsat/dT dT ) q/qsat(T) (1 + RHpert)
  ! => dq = q RHpert + dqsat/dT dT q/qsat(T) (1 + RHpert)
  !
  ! Now for neutral buoyancy:
  ! dTv/dT dT + dTv/dq dq = 0
  !
  ! Combine to eliminate dT:
  ! dq = q RHpert - dqsat/dT q/qsat(T) (1 + RHpert) ( dTv/dq / dTv/dT ) dq
  ! => dq = q RHpert
  !       / ( 1 + dqsat/dT q/qsat(T) (1 + RHpert) ( dTv/dq / dTv/dT ) )

  ! i.e. the q perturbation scales with the grid-mean q:

  delta_qvap_neut(ic) = q_vap(ic) * par_gen_rhpert                             &
      / ( one + (dtv_dqv(ic)/dtv_dt(ic)) * dqsatdt_liq(ic)                     &
              * (q_vap(ic)/qsat_liq(ic)) * (one + par_gen_rhpert) )

  delta_temp_neut(ic) = -delta_qvap_neut(ic)                                   &
                      * (dtv_dqv(ic)/dtv_dt(ic))

end do


!------------------------------------------------------------------------------
! 6) Copy outputs to diagnostics if requested
!------------------------------------------------------------------------------

do i_region = 1, n_regions
  if ( genesis_diags % subregion_diags(i_region) % n_diags > 0 ) then
    ! If any diags requested for this region
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

      ! If area fraction requested, copy into super-array
      if ( genesis_diags % subregion_diags(i_region) % frac % flag ) then
        i_super = genesis_diags % subregion_diags(i_region) % frac % i_super
        do ic2 = 1, nc
          ic = index_ic(ic2)
          diags_super(ic,i_super) = frac_r(ic,i_region)
        end do
      end if

      ! If any field diagnostics requested for this region
      if ( genesis_diags % subregion_diags(i_region)%fields%n_diags > 0 ) then

        ! Initialise a fields super-array to compute fields diags
        if ( .not. allocated(fields_reg) ) then
          allocate( fields_reg(n_points,n_fields) )
        end if
        do i_field = 1, n_fields
          do ic = 1, n_points
            fields_reg(ic,i_field) = zero
          end do
        end do

        ! Copy fields from the current region into a fields super-array
        do ic2 = 1, nc
          ic = index_ic(ic2)
          fields_reg(ic2,i_temperature) = temperature_r(ic,i_region)
          fields_reg(ic2,i_q_vap)       = q_vap_r(ic,i_region)
        end do
        ! Set local condensed water species mixing ratios and
        ! cloud-fractions based on current region index:
        call set_region_cond_fields( n_points, nc, index_ic,                   &
                                     n_points, i_region,                       &
                                     q_cond_loc, cloudfracs(:,i_frac_ice),     &
                                     frac_r, fields_reg )
        ! Expand fields onto full list for input to fields_diags_copy
        do ic2 = nc, 1, -1
          ic = index_ic(ic2)
          if ( ic > ic2 ) then
            do i_field = 1, n_fields
              fields_reg(ic,i_field) = fields_reg(ic2,i_field)
              fields_reg(ic2,i_field) = zero
            end do
          end if
        end do
        ! Copy to diagnostics super-array
        call fields_diags_copy(                                                &
               n_points, n_points, n_points_super,                             &
               genesis_diags % n_diags, n_fields,                              &
               l_conserved_form_false,                                         &
               genesis_diags % subregion_diags(i_region) % fields,             &
               fields_reg, virt_temp, pressure,                                &
               diags_super )

      end if ! ( genesis_diags % subregion_diags(i_region)%fields%n_diags > 0 )

    end if  ! ( nc > 0 ) THEN
  end if  ! ( genesis_diags % subregion_diags(i_region) % n_diags > 0 )
end do  ! i_region = 1, n_regions


! Check outputs for bad values (NaN, Inf, negatives, etc)
if ( i_check_bad_values_cmpr > i_check_bad_none ) then

  where_string = "calc_env_regions call for " // trim(adjustl(call_string))
  field_positive = .true.

  do i_region = 1, n_regions
    ! Check area fractions
    field_name = trim(adjustl(region_names(i_region))) // "_frac"
    call check_bad_values_cmpr( cmpr, k,                                       &
                                frac_r(:,i_region),                            &
                                where_string, field_name,                      &
                                field_positive )
    ! Check temperatures
    field_name = trim(adjustl(region_names(i_region))) // "_temperature"
    call check_bad_values_cmpr( cmpr, k,                                       &
                                temperature_r(:,i_region),                     &
                                where_string, field_name,                      &
                                field_positive )
    ! Check water-vapour mixing-ratio
    field_name = trim(adjustl(region_names(i_region))) // "_q_vap"
    call check_bad_values_cmpr( cmpr, k,                                       &
                                q_vap_r(:,i_region),                           &
                                where_string, field_name,                      &
                                field_positive )
  end do

  ! Check condensed water species
  do i_cond = 1, n_cond_species
    field_name = "loc_" // trim(adjustl(field_names(i_qc_first-1+i_cond)))
    call check_bad_values_cmpr( cmpr, k,                                       &
                                q_cond_loc(:,i_cond),                          &
                                where_string, field_name,                      &
                                field_positive )
  end do

end if  ! ( i_check_bad_values_cmpr > i_check_bad_none )


return
end subroutine calc_env_regions

end module calc_env_regions_mod
