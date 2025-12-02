! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_subregion_diags_mod

implicit none

contains


! Subroutine to calculate diagnostics of the properties of the
! clear, cloudy, and rainy / icy sub-regions of the grid-box,
! as assumed in the convective initiation mass-source
! calculation.  This is basically just a compression
! wrapper around the routine calc_env_partitions, which
! is called from init_mass_moist_frac.
subroutine calc_subregion_diags( grid, fields, cloudfracs,                     &
                                 subregion_diags )

use comorph_constants_mod, only: real_cvprec, zero,                            &
                     nx_full, ny_full, k_bot_conv, k_top_conv,                 &
                     l_cv_cloudfrac,                                           &
                     n_cond_species, name_length
use cmpr_type_mod, only: cmpr_type, cmpr_alloc, cmpr_dealloc
use grid_type_mod, only: grid_type
use fields_type_mod, only: fields_type, n_fields,                              &
                           i_temperature, i_q_vap,                             &
                           i_qc_first, i_qc_last,                              &
                           i_cf_liq, field_positive
use cloudfracs_type_mod, only: cloudfracs_type, n_cloudfracs,                  &
                               i_frac_liq, i_frac_ice,                         &
                               i_frac_bulk, i_frac_precip
use subregion_mod, only: n_regions
use subregion_diags_type_mod, only: subregion_diags_type
use fields_diags_type_mod, only: fields_diags_copy
use calc_virt_temp_mod, only: calc_virt_temp
use calc_env_partitions_mod, only: calc_env_partitions
use set_region_cond_fields_mod, only: set_region_cond_fields
use compress_mod, only: compress
use decompress_mod, only: decompress
use force_cloudfrac_consistency_mod, only:                                     &
                                 force_cloudfrac_consistency

implicit none

! Structure containing pointers to model grid fields
type(grid_type), intent(in) :: grid

! Structure containing pointers to the primary fields
type(fields_type), intent(in) :: fields

! Structure containing pointers to diagnostic cloud fractions,
! in the event that cloud fractions are not included as
! prognostic primary fields
type(cloudfracs_type), intent(in) :: cloudfracs

! Structures containing meta-data and pointers to the diagnostic
! arrays to be calculated (separate structure for each region)
type(subregion_diags_type), intent(in out) :: subregion_diags                  &
                                             (n_regions)


! Structure storing compression / decompression indices
type(cmpr_type) :: cmpr

! Compressed fields required:
! Pressure
real(kind=real_cvprec) :: pressure( nx_full*ny_full )
! Compressed environment virtual temperature
real(kind=real_cvprec) :: virt_temp( nx_full*ny_full )
! Primary fields
real(kind=real_cvprec) :: fields_cmpr                                          &
                    ( nx_full*ny_full, n_fields )
! Cloud fractions
real(kind=real_cvprec) :: cloudfracs_cmpr                                      &
                    ( nx_full*ny_full, n_cloudfracs )

! Separate region fields output by calc_env_partitions:
! Fractional area of each region
real(kind=real_cvprec) :: frac_r                                               &
                    ( nx_full*ny_full, n_regions )
! Temperature of each region
real(kind=real_cvprec) :: temperature_r                                        &
                    ( nx_full*ny_full, n_regions )
! Water-vapour mixing-ratio of each region
real(kind=real_cvprec) :: q_vap_r                                              &
                    ( nx_full*ny_full, n_regions )
! Local values of condensed-water mixing-ratios
real(kind=real_cvprec) :: q_cond_loc                                           &
                    ( nx_full*ny_full, n_cond_species )

! Array for storing any requested fields diagnostics
real(kind=real_cvprec), allocatable :: diags_super(:,:)

! Redundant compression indices needed for input to
! set_region_cond_fields
integer :: index_ic( nx_full*ny_full )

! Flag input to fields_diags_copy
logical :: l_conserved_form

! Lower and upper bounds of arrays
integer :: lb(3), ub(3)

character(len=name_length) :: call_string

! Loop counters
integer :: i, j, k, ic, i_field, i_frac,                                       &
           i_diag, i_region, i_super


! Setup compression indices for all points on a slice
cmpr % n_points = nx_full*ny_full
call cmpr_alloc( cmpr, cmpr % n_points )
do j = 1, ny_full
  do i = 1, nx_full
    ic = nx_full*(j-1) + i
    cmpr % index_i(ic) = i
    cmpr % index_j(ic) = j
  end do
end do

do ic = 1, cmpr % n_points
  index_ic(ic) = ic
end do


! Loop over levels
!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP SHARED(  k_bot_conv, k_top_conv, cmpr, index_ic, l_cv_cloudfrac,         &
!$OMP          grid, fields, cloudfracs, n_fields, field_positive,             &
!$OMP          i_cf_liq, i_temperature, i_q_vap, i_qc_first, i_qc_last,        &
!$OMP          subregion_diags )                                               &
!$OMP private( k, lb, ub, i_field, i_frac, i_region, i_diag, i_super, ic,      &
!$OMP          pressure, fields_cmpr, cloudfracs_cmpr,                         &
!$OMP          virt_temp, call_string, frac_r, temperature_r, q_vap_r,         &
!$OMP          q_cond_loc, diags_super, l_conserved_form )
do k = k_bot_conv, k_top_conv

  ! Compress pressure
  lb = lbound( grid % pressure_full )
  ub = ubound( grid % pressure_full )
  call compress( cmpr, lb(1:2), ub(1:2),                                       &
                 grid % pressure_full(:,:,k), pressure )

  ! Compress all primary fields
  do i_field = 1, n_fields
    lb = lbound( fields % list(i_field)%pt )
    ub = ubound( fields % list(i_field)%pt )
    call compress( cmpr, lb(1:2), ub(1:2),                                     &
                   fields % list(i_field)%pt(:,:,k), fields_cmpr(:,i_field) )
    if ( field_positive(i_field) ) then
      ! Remove any spurious negative values from input data
      do ic = 1, cmpr % n_points
        fields_cmpr(ic,i_field) = max( fields_cmpr(ic,i_field), zero )
      end do
    end if
  end do

  ! Compress cloud fractions
  if ( l_cv_cloudfrac ) then
    ! If cloud fractions are included as primary fields,
    ! just copy them from there
    do i_frac = i_frac_liq, i_frac_bulk
      i_field = i_cf_liq-1 + i_frac
      do ic = 1, cmpr % n_points
        cloudfracs_cmpr(ic,i_frac) = fields_cmpr(ic,i_field)
      end do
    end do
  else
    ! Cloud fractions are not primary fields; compress from
    ! separate structure
    lb = lbound( cloudfracs % frac_liq )
    ub = ubound( cloudfracs % frac_liq )
    call compress( cmpr, lb(1:2), ub(1:2),                                     &
                   cloudfracs % frac_liq(:,:,k),                               &
                   cloudfracs_cmpr(:,i_frac_liq) )
    lb = lbound( cloudfracs % frac_ice )
    ub = ubound( cloudfracs % frac_ice )
    call compress( cmpr, lb(1:2), ub(1:2),                                     &
                   cloudfracs % frac_ice(:,:,k),                               &
                   cloudfracs_cmpr(:,i_frac_ice) )
    lb = lbound( cloudfracs % frac_bulk )
    ub = ubound( cloudfracs % frac_bulk )
    call compress( cmpr, lb(1:2), ub(1:2),                                     &
                   cloudfracs % frac_bulk(:,:,k),                              &
                   cloudfracs_cmpr(:,i_frac_bulk) )
  end if

  ! Rounding errors when converting the cloud-fractions
  ! to 32-bit in compress can cause them to become
  ! slightly inconsistent; correct if needed:
  call force_cloudfrac_consistency( cmpr%n_points, cmpr%n_points,              &
                                    cloudfracs_cmpr(:,i_frac_liq:i_frac_bulk) )

  ! Compress rain fraction (always a separate diagnostic field)
  lb = lbound( cloudfracs % frac_precip )
  ub = ubound( cloudfracs % frac_precip )
  call compress( cmpr, lb(1:2), ub(1:2),                                       &
                 cloudfracs % frac_precip(:,:,k),                              &
                 cloudfracs_cmpr(:,i_frac_precip) )
  ! Remove spurious negatives
  do ic = 1, cmpr % n_points
    cloudfracs_cmpr(ic,i_frac_precip)                                          &
      = max( cloudfracs_cmpr(ic,i_frac_precip), zero )
  end do

  ! Compute grid-mean virtual temperature
  call calc_virt_temp( cmpr%n_points, cmpr%n_points,                           &
                       fields_cmpr(:,i_temperature),                           &
                       fields_cmpr(:,i_q_vap),                                 &
                       fields_cmpr(:,i_qc_first:i_qc_last),                    &
                       virt_temp )


  ! Calculate env sub-region properties
  call_string = "diagnostics"
  call calc_env_partitions( cmpr%n_points, cmpr%n_points,                      &
                            cmpr, k, call_string,                              &
                            pressure, cloudfracs_cmpr,                         &
                            fields_cmpr(:,i_temperature),                      &
                            fields_cmpr(:,i_q_vap),                            &
                            fields_cmpr(:,i_qc_first:i_qc_last),               &
                            frac_r, temperature_r, q_vap_r,                    &
                            q_cond_loc )


  ! Loop over regions
  do i_region = 1, n_regions
    ! If any diags requested for this region
    if ( subregion_diags(i_region) % n_diags > 0 ) then

      ! Allocate diagnostics super-array
      allocate( diags_super( cmpr % n_points,                                  &
                             subregion_diags(i_region) % n_diags ) )

      ! If area fraction requested, copy into super-array
      if ( subregion_diags(i_region) % frac % flag ) then
        i_super = subregion_diags(i_region) % frac % i_super
        do ic = 1, cmpr % n_points
          diags_super(ic,i_super) = frac_r(ic,i_region)
        end do
      end if

      ! If any field diagnostics requested for this region
      if ( subregion_diags(i_region)%fields % n_diags > 0 ) then

        ! Copy fields from the current region back into the
        ! fields super-array
        do ic = 1, cmpr % n_points
          fields_cmpr(ic,i_temperature)                                        &
                                    = temperature_r(ic,i_region)
          fields_cmpr(ic,i_q_vap)   = q_vap_r(ic,i_region)
        end do
        ! Set local condensed water species mixing ratios and
        ! cloud-fractions based on current region index:
        call set_region_cond_fields( cmpr%n_points,cmpr%n_points,              &
                                     index_ic, cmpr%n_points,                  &
                                     i_region, q_cond_loc,                     &
                                     cloudfracs_cmpr(:,i_frac_ice),            &
                                     frac_r, fields_cmpr )

        ! Copy requested diags into the super-array,
        ! and calculate any derived field diagnostics
        l_conserved_form = .false.
        call fields_diags_copy( cmpr%n_points, cmpr%n_points, cmpr%n_points,   &
                                subregion_diags(i_region)%n_diags, n_fields,   &
                                l_conserved_form,                              &
                                subregion_diags(i_region)%fields,              &
                                fields_cmpr,                                   &
                                virt_temp, pressure,                           &
                                diags_super )

      end if  !( subregion_diags(i_region)%fields % n_diags > 0 )

      ! Decompress diags into the output 3-D arrays
      do i_diag = 1, subregion_diags(i_region) % n_diags
        i_super = subregion_diags(i_region) % list(i_diag)%pt % i_super
        lb = lbound( subregion_diags(i_region) % list(i_diag)%pt % field_3d )
        ub = ubound( subregion_diags(i_region) % list(i_diag)%pt % field_3d )
        call decompress( cmpr, diags_super(:,i_super), lb(1:2), ub(1:2),       &
                         subregion_diags(i_region) % list(i_diag)%pt           &
                         % field_3d(:,:,k) )
      end do

      ! Deallocate diags super-array
      deallocate( diags_super )

    end if  ! ( subregion_diags(i_region) % n_diags > 0 )
  end do  ! i_region = 1, n_regions


end do  ! k = k_bot_conv, k_top_conv
!$OMP end PARALLEL do


call cmpr_dealloc( cmpr )


return
end subroutine calc_subregion_diags


end module calc_subregion_diags_mod
