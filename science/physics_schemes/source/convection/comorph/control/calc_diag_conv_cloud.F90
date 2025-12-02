! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_diag_conv_cloud_mod

implicit none

contains


!----------------------------------------------------------------
! Subroutine sums the diagnosed convective cloud contributions
! over all convective drafts / types / layers to calculate
! the overall convective cloud fraction and grid-mean convective
! cloud liquid and ice mixing-ratios.
!----------------------------------------------------------------
subroutine calc_diag_conv_cloud( n_updraft_layers,                             &
                                 n_dndraft_layers,                             &
                                 n_points_super, n_fields_tot,                 &
                                 cmpr_any, k,                                  &
                                 ij_first, ij_last, index_ic,                  &
                                 updraft_res_source,                           &
                                 updraft_fallback_res_source,                  &
                                 dndraft_res_source,                           &
                                 dndraft_fallback_res_source,                  &
                                 lb_p, ub_p, pressure,                         &
                                 fields_k, cloudfracs )

use comorph_constants_mod, only: real_cvprec, real_hmprec, zero,               &
                     k_bot_conv, k_top_conv,                                   &
                     n_updraft_types, n_dndraft_types,                         &
                     l_updraft_fallback, l_dndraft_fallback,                   &
                     max_sigma, i_convcloud, i_convcloud_mph

use cmpr_type_mod, only: cmpr_type
use res_source_type_mod, only: res_source_type
use cloudfracs_type_mod, only: cloudfracs_type,                                &
                               n_convcloud, i_frac_liq_conv,                   &
                               i_frac_ice_conv, i_frac_bulk_conv,              &
                               i_q_cl_conv, i_q_cf_conv
use fields_type_mod, only: i_temperature,                                      &
                           i_q_vap, i_q_cl, i_q_cf
use add_conv_cloud_mod, only: add_conv_cloud
use set_qsat_mod, only: set_qsat_liq
use compress_mod, only: compress
use decompress_mod, only: decompress

implicit none

! Highest number of distinct updraft and downdraft layers
! occuring on the current segment
integer, intent(in) :: n_updraft_layers
integer, intent(in) :: n_dndraft_layers

! Number of points in the compressed fields super-array
integer, intent(in) :: n_points_super

! Number of fields in the compressed fields super-array
integer, intent(in) :: n_fields_tot

! Structure storing compression list indices for the points where
! any kind of convection is occuring
type(cmpr_type), intent(in) :: cmpr_any( k_bot_conv:k_top_conv )

! Current model-level index
integer, intent(in) :: k

! Position indices ( nx*(j-1) + i ) of the first and last
! convecting point on the current segment
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last

! cmpr_any compression indices scattered into a common
! array-space for referencing from the individual convection
! type compression lists
integer, intent(in) :: index_ic( ij_first : ij_last,                           &
                                 k_bot_conv-1 : k_top_conv+1 )

! Arrays of structures containing resolved-scale source terms
! in separate compression lists for updrafts & downdrafts
! and their fall-back flows
type(res_source_type), intent(in) :: updraft_res_source (:,:,:)
type(res_source_type), intent(in) :: updraft_fallback_res_source               &
                                                        (:,:,:)
type(res_source_type), intent(in) :: dndraft_res_source (:,:,:)
type(res_source_type), intent(in) :: dndraft_fallback_res_source               &
                                                        (:,:,:)
! These are allocatable arrays, and we don't know for certain
! what shape they will have been allocated to.
! If n_updraft_layers > 0, shape is
! ( n_updraft_types, n_updraft_layers, k_bot_conv:k_top_conv )
! But if n_updraft_layers = 0, shape is
! ( 1, 1, k_bot_conv:k_top_conv )

! Full 3-D array of pressure
! (and lower and upper bounds in case it has halos)
integer, intent(in) :: lb_p(3), ub_p(3)
real(kind=real_hmprec), intent(in) :: pressure                                 &
                          ( lb_p(1):ub_p(1), lb_p(2):ub_p(2), lb_p(3):ub_p(3) )

! Compressed super-array containing primary fields on the
! current level, updated with increments from convective
! entrainment, detrainment, subsidence etc.
real(kind=real_cvprec), intent(in) :: fields_k                                 &
                                ( n_points_super, n_fields_tot )

! Structure containing pointers to the output full 3-D
! convective cloud fields.
type(cloudfracs_type), intent(in out) :: cloudfracs


! Super-array containing sum of convective cloud variables over
! all drafts, types and layers on the current level
real(kind=real_cvprec) :: conv_cloud                                           &
                          ( cmpr_any(k) % n_points, n_convcloud )

! Pressure at convecting points on current level
real(kind=real_cvprec) :: pressure_k( cmpr_any(k) % n_points )
! qsat w.r.t. liquid water
real(kind=real_cvprec) :: qsat_liq( cmpr_any(k) % n_points )
! Work variables for safety-check to avoid impying little or no
! water left outside the convective clouds
real(kind=real_cvprec) :: q_tot_mean, q_tot_conv, tmp

! Array lower and upper bounds
integer :: lb(3), ub(3)

! Limit on fraction of grid-mean water which is allowed to be
! assumed to be held inside the convective clouds
real(kind=real_cvprec), parameter :: max_water_frac                            &
                                     = 0.8_real_cvprec

! Loop counters
integer :: ic, i_field


! Initialise cloud properties to zero
do i_field = 1, n_convcloud
  do ic = 1, cmpr_any(k) % n_points
    conv_cloud(ic,i_field) = zero
  end do
end do

! Add on the contributions from each draft...

if ( n_updraft_types > 0 .and. n_updraft_layers > 0 ) then
  ! Contribution from primary updrafts
  call add_conv_cloud( cmpr_any(k) % n_points,                                 &
                       n_updraft_types, n_updraft_layers,                      &
                       ij_first, ij_last, index_ic(:,k),                       &
                       updraft_res_source(:,:,k),                              &
                       conv_cloud )
  ! Contribution from updraft fall-backs
  if ( l_updraft_fallback ) then
    call add_conv_cloud( cmpr_any(k) % n_points,                               &
                         n_updraft_types, n_updraft_layers,                    &
                         ij_first, ij_last, index_ic(:,k),                     &
                         updraft_fallback_res_source(:,:,k),                   &
                         conv_cloud )
  end if
end if

if ( n_dndraft_types > 0 .and. n_dndraft_layers > 0 ) then
  ! Contribution from primary downdrafts
  call add_conv_cloud( cmpr_any(k) % n_points,                                 &
                       n_dndraft_types, n_dndraft_layers,                      &
                       ij_first, ij_last, index_ic(:,k),                       &
                       dndraft_res_source(:,:,k),                              &
                       conv_cloud )
  ! Contribution from downdraft fall-backs
  if ( l_dndraft_fallback ) then
    call add_conv_cloud( cmpr_any(k) % n_points,                               &
                         n_dndraft_types, n_dndraft_layers,                    &
                         ij_first, ij_last, index_ic(:,k),                     &
                         dndraft_fallback_res_source(:,:,k),                   &
                         conv_cloud )
  end if
end if

! Safety-check to avoid convective cloud-fraction > 1/2
do ic = 1, cmpr_any(k) % n_points
  ! If convective bulk cloud-fraction exceeds max, scale down
  ! all the convective cloud variables consistently
  if ( conv_cloud(ic,i_frac_bulk_conv) > max_sigma ) then
    tmp = max_sigma / conv_cloud(ic,i_frac_bulk_conv)
    ! Scale down all convective cloud fields consistently
    do i_field = 1, n_convcloud
      conv_cloud(ic,i_field) = conv_cloud(ic,i_field) * tmp
    end do
  end if
end do

! Safety-check to avoid implying that the convective cloud
! contains more water than exists in the grid-box mean fields:

! Compress pressure
call compress( cmpr_any(k), lb_p(1:2), ub_p(1:2),                              &
               pressure(:,:,k), pressure_k )
! Calculate qsat_liq
call set_qsat_liq( cmpr_any(k) % n_points,                                     &
                   fields_k(:,i_temperature), pressure_k, qsat_liq )

! Calculate approximate water-content of the convective cloud
! (taking approximation that cloud has same temperature as
! the environment), and compare with the grid-mean
! water content...
! (Note: only dealing with vapour + liquid and ice cloud.
!  Also neglecting increase in T and qsat due to condensation
!  when adding the cloud).
do ic = 1, cmpr_any(k) % n_points

  ! Calculate grid-mean vapour + cloud water
  ! (with check to exclude negative values which can
  !  occasionally occur in the input moisture fields).
  q_tot_mean = max( fields_k(ic,i_q_vap)                                       &
                  + fields_k(ic,i_q_cl) + fields_k(ic,i_q_cf),                 &
                    zero )

  q_tot_conv = conv_cloud(ic,i_q_cl_conv)                                      &
  ! Add water-vapour contained in the convective
  ! liquid cloud fraction, assuming it is saturated.
             + conv_cloud(ic,i_frac_liq_conv) * qsat_liq(ic)                   &
  ! Add water-vapour contained in any convective
  ! ice-only cloud fraction, assuming same q_vap as grid-mean
             + ( conv_cloud(ic,i_frac_bulk_conv)                               &
               - conv_cloud(ic,i_frac_liq_conv) )                              &
               * max( fields_k(ic,i_q_vap), zero )

  ! Also include ice content if it is output
  if ( i_convcloud == i_convcloud_mph ) then
    q_tot_conv = q_tot_conv + conv_cloud(ic,i_q_cf_conv)
  end if

  ! Don't allow the implied water content of the convective
  ! cloud to exceed a certain fraction of the grid-mean water
  if ( q_tot_conv > max_water_frac * q_tot_mean ) then
    tmp = ( max_water_frac * q_tot_mean ) / q_tot_conv
    ! Scale down all convective cloud fields consistently
    do i_field = 1, n_convcloud
      conv_cloud(ic,i_field) = conv_cloud(ic,i_field) * tmp
    end do
  end if

end do


! Scatter convective cloud properties back to full 3-D arrays
lb = lbound( cloudfracs % frac_liq_conv )
ub = ubound( cloudfracs % frac_liq_conv )
call decompress( cmpr_any(k), conv_cloud(:,i_frac_liq_conv),                   &
                 lb(1:2), ub(1:2), cloudfracs % frac_liq_conv(:,:,k) )
lb = lbound( cloudfracs % frac_bulk_conv )
ub = ubound( cloudfracs % frac_bulk_conv )
call decompress( cmpr_any(k), conv_cloud(:,i_frac_bulk_conv),                  &
                 lb(1:2), ub(1:2), cloudfracs % frac_bulk_conv(:,:,k) )
lb = lbound( cloudfracs % q_cl_conv )
ub = ubound( cloudfracs % q_cl_conv )
call decompress( cmpr_any(k), conv_cloud(:,i_q_cl_conv),                       &
                 lb(1:2), ub(1:2), cloudfracs % q_cl_conv(:,:,k) )
if ( i_convcloud == i_convcloud_mph ) then
  lb = lbound( cloudfracs % frac_ice_conv )
  ub = ubound( cloudfracs % frac_ice_conv )
  call decompress( cmpr_any(k), conv_cloud(:,i_frac_ice_conv),                 &
                   lb(1:2), ub(1:2), cloudfracs % frac_ice_conv(:,:,k) )
  lb = lbound( cloudfracs % q_cf_conv )
  ub = ubound( cloudfracs % q_cf_conv )
  call decompress( cmpr_any(k), conv_cloud(:,i_q_cf_conv),                     &
                   lb(1:2), ub(1:2), cloudfracs % q_cf_conv(:,:,k) )
end if


return
end subroutine calc_diag_conv_cloud


!----------------------------------------------------------------
! Subroutine to do the convective cloud calculations at points
! where there is no convection this timestep.
! Currently just zeros the fields, but if we implement
! prognostic convective cloud, it will need to add increments
! due to evaporating any pre-existing convective cloud.
!----------------------------------------------------------------
subroutine zero_diag_conv_cloud( ij_0, ij_1, k, cmpr_any,                      &
                                 cloudfracs )

use comorph_constants_mod, only: k_bot_conv, k_top_conv, nx_full,              &
                     i_convcloud, i_convcloud_mph
use cmpr_type_mod, only: cmpr_type, cmpr_alloc, cmpr_dealloc
use cloudfracs_type_mod, only: cloudfracs_type
use init_zero_mod, only: init_zero_cmpr

implicit none

! Position of the first and last point within the current
! segment's domain, including non-convecting points
integer, intent(in) :: ij_0
integer, intent(in) :: ij_1

! Current model-level index
integer, intent(in) :: k

! Structure storing compression list indices for the points where
! any kind of convection is occuring
type(cmpr_type), intent(in) :: cmpr_any( k_bot_conv:k_top_conv )

! Structure containing pointers to the convective cloud fields
type(cloudfracs_type), intent(in out) :: cloudfracs

! Mask of points on the current segment where there is
! not any convection
logical :: l_mask(ij_0:ij_1)

! Points where calculations are needed
type(cmpr_type) :: cmpr_none

! Array bounds
integer :: lb(3), ub(3)

! Loop counters
integer :: i, j, ic, ij


! Set number of non-convecting points and allocate compression
! list for those points:
cmpr_none % n_points = ij_1 - ij_0 + 1 - cmpr_any(k) % n_points
call cmpr_alloc( cmpr_none, cmpr_none % n_points )

! If any convecting points on this level:
if ( cmpr_any(k) % n_points > 0 ) then

  ! Initialise mask to true everywhere
  do ij = ij_0, ij_1
    l_mask(ij) = .true.
  end do

  ! Reset to false at points where there is convection on the
  ! current model-level
  do ic = 1, cmpr_any(k) % n_points
    i = cmpr_any(k) % index_i(ic)
    j = cmpr_any(k) % index_j(ic)
    ij = nx_full*(j-1)+i
    l_mask(ij) = .false.
  end do

  ! Find list of points where no convection
  cmpr_none % n_points = 0
  do ij = ij_0, ij_1
    if ( l_mask(ij) ) then
      cmpr_none % n_points = cmpr_none % n_points + 1
      j = 1 + (ij-1) / nx_full
      i = ij - nx_full*(j-1)
      cmpr_none % index_i(cmpr_none % n_points) = i
      cmpr_none % index_j(cmpr_none % n_points) = j
    end if
  end do

else
  ! No convecting points; set compression list to all points on
  ! the current segment.

  do ij = ij_0, ij_1
    j = 1 + (ij-1) / nx_full
    i = ij - nx_full*(j-1)
    ic = ij - ij_0 + 1
    cmpr_none % index_i(ic) = i
    cmpr_none % index_j(ic) = j
  end do

end if

! Set convective cloud fields to zero at these points
lb = lbound( cloudfracs % frac_liq_conv )
ub = ubound( cloudfracs % frac_liq_conv )
call init_zero_cmpr( cmpr_none, lb(1:2), ub(1:2),                              &
                     cloudfracs % frac_liq_conv(:,:,k) )
lb = lbound( cloudfracs % frac_bulk_conv )
ub = ubound( cloudfracs % frac_bulk_conv )
call init_zero_cmpr( cmpr_none, lb(1:2), ub(1:2),                              &
                     cloudfracs % frac_bulk_conv(:,:,k) )
lb = lbound( cloudfracs % q_cl_conv )
ub = ubound( cloudfracs % q_cl_conv )
call init_zero_cmpr( cmpr_none, lb(1:2), ub(1:2),                              &
                     cloudfracs % q_cl_conv(:,:,k) )

if ( i_convcloud == i_convcloud_mph ) then
  lb = lbound( cloudfracs % frac_ice_conv )
  ub = ubound( cloudfracs % frac_ice_conv )
  call init_zero_cmpr( cmpr_none, lb(1:2), ub(1:2),                            &
                       cloudfracs % frac_ice_conv(:,:,k) )
  lb = lbound( cloudfracs % q_cf_conv )
  ub = ubound( cloudfracs % q_cf_conv )
  call init_zero_cmpr( cmpr_none, lb(1:2), ub(1:2),                            &
                       cloudfracs % q_cf_conv(:,:,k) )
end if

call cmpr_dealloc( cmpr_none )


return
end subroutine zero_diag_conv_cloud



end module calc_diag_conv_cloud_mod
