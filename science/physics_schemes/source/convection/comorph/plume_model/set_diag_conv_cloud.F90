! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module set_diag_conv_cloud_mod

implicit none

contains

! Subroutine to set the diagnosed convective cloud properties
! at the current level, based on the convective updraft or
! downdraft parcel properties
subroutine set_diag_conv_cloud( n_points, n_points_env,                        &
                                n_points_res,                                  &
                                n_points_prev, n_points_next,                  &
                                n_fields_tot,                                  &
                                max_buoy_heights, n_buoy_vars,                 &
                                grid_prev, grid_next, frac_level_step,         &
                                prev_virt_temp, next_virt_temp,                &
                                prev_massflux_d, next_massflux_d,              &
                                prev_radius, next_radius,                      &
                                par_prev_mean, par_next_mean,                  &
                                buoyancy_super, i_next, i_sat,                 &
                                convcloud_super )

use comorph_constants_mod, only: real_cvprec, zero, half, one,                 &
                                 gravity, wind_w_fac, cf_conv_fac,             &
                                 l_cf_conv_buoy, wind_w_buoy_fac, w_min,       &
                                 l_cv_cloudfrac
use fields_type_mod, only: n_fields, i_temperature,                            &
                           i_q_vap, i_qc_first, i_qc_last,                     &
                           i_q_cl, i_q_cf, i_cf_liq, i_cf_bulk
use cloudfracs_type_mod, only: n_convcloud, i_frac_liq, i_frac_bulk
use buoyancy_mod, only: i_prev, i_mean_buoy
use grid_type_mod, only: n_grid, i_height, i_pressure

use calc_q_tot_mod, only: calc_q_tot
use calc_rho_dry_mod, only: calc_rho_dry
use interp_diag_conv_cloud_mod, only: interp_diag_conv_cloud
use set_par_cloudfrac_mod, only: set_par_cloudfrac

implicit none

! Number of points
integer, intent(in) :: n_points

! Size of grid field super-arrays
integer, intent(in) :: n_points_env
! Size of output convective cloud super-arrays
integer, intent(in) :: n_points_res

! Size of primary fields super-arrays
integer, intent(in) :: n_points_prev
integer, intent(in) :: n_points_next
integer, intent(in) :: n_fields_tot

! Dimensions of the buoyancy super-array
integer, intent(in) :: max_buoy_heights
integer, intent(in) :: n_buoy_vars

! Heights and pressures of model-levels:
! Previous model-level interface
real(kind=real_cvprec), intent(in) :: grid_prev                                &
                                      ( n_points_env, n_grid )
! Next model-level interface
real(kind=real_cvprec), intent(in) :: grid_next                                &
                                      ( n_points_env, n_grid )
! Fraction of full-level-step performed by this call
real(kind=real_cvprec), intent(in) :: frac_level_step(n_points)

! Previous and next environment virtual temperatures
real(kind=real_cvprec), intent(in) :: prev_virt_temp(n_points)
real(kind=real_cvprec), intent(in) :: next_virt_temp(n_points)

! Parcel mass-flux at previous and next model-level interfaces
real(kind=real_cvprec), intent(in) :: prev_massflux_d(n_points)
real(kind=real_cvprec), intent(in) :: next_massflux_d(n_points)

! Parcel radius at previous and next model-level interfaces
real(kind=real_cvprec), intent(in) :: prev_radius(n_points)
real(kind=real_cvprec), intent(in) :: next_radius(n_points)

! In-parcel mean primary fields at previous and next
real(kind=real_cvprec), intent(in) :: par_prev_mean                            &
                             ( n_points_prev, i_temperature:n_fields )
real(kind=real_cvprec), intent(in) :: par_next_mean                            &
                             ( n_points_next, n_fields_tot )

! Super-array storing the buoyancies of the parcel mean and core
! at up to 4 sub-level heights within the current level-step:
! a) Previous model-level interface
! b) Next model-level interface
! c) Height where parcel mean properties first hit saturation
! d) Height where parcel core first hit saturation
real(kind=real_cvprec), intent(in) :: buoyancy_super                           &
                     ( n_points, n_buoy_vars, max_buoy_heights )

! Address of next model-level interface in buoyancy_super
integer, intent(in) :: i_next(n_points)
! Address of saturation height in buoyancy_super
integer, intent(in) :: i_sat(n_points)

! Super-array containing diagnosed convective cloud properties
real(kind=real_cvprec), intent(in out) :: convcloud_super                      &
                                          ( n_points_res, n_convcloud )

! Dry density
real(kind=real_cvprec) :: rho_dry(n_points)

! Total-water mixing ratio
real(kind=real_cvprec) :: q_tot(n_points)

! Mean parcel buoyancy over the current model-level
real(kind=real_cvprec) :: mean_buoy(n_points)

! Convective fractions at prev and next
real(kind=real_cvprec) :: prev_cf_conv(n_points)
real(kind=real_cvprec) :: next_cf_conv(n_points)

! Convective grid-mean liq and ice mixing ratios at prev and next
real(kind=real_cvprec) :: prev_q_cl_conv(n_points)
real(kind=real_cvprec) :: next_q_cl_conv(n_points)
real(kind=real_cvprec) :: prev_q_cf_conv(n_points)
real(kind=real_cvprec) :: next_q_cf_conv(n_points)

! In-parcel cloud-fractions at prev and next (separate arrays for
! calculating them in the case where cloud-fractions aren't included
! in the primary fields super-arrays)
real(kind=real_cvprec), allocatable :: par_prev_cloudfracs(:,:)
real(kind=real_cvprec), allocatable :: par_next_cloudfracs(:,:)

! Height where parcel reaches saturation (if between prev and next)
real(kind=real_cvprec) :: sat_height(n_points)

! Super-array for convective cloud-fields from only this call
real(kind=real_cvprec) :: convcloud_step ( n_points, n_convcloud )

! Loop counter
integer :: ic, i_lev, i_field


! Calculate the convective fraction irrespective of whether
! or not it contains any cloud (cf_conv), at prev and next...

! Now the dry-mass flux  Md = rho_dry w' sigma
!
! => sigma = Md / ( rho_dry w' )

if ( l_cf_conv_buoy ) then

  ! Assuming that w' = fac * sqrt( buoyancy * radius )
  !
  ! => sigma = Md / ( rho_dry fac * sqrt( buoyancy * radius ) )

  ! To reduce noisiness where a very small buoyancy occurs
  ! on just one sub-level, take the vertical mean of the buoyancy
  ! profile between prev and next...
  do ic = 1, n_points
    mean_buoy(ic) = zero
  end do
  do ic = 1, n_points
    do i_lev = i_prev, i_next(ic) - 1
      mean_buoy(ic) = mean_buoy(ic)                                            &
        + half * ( buoyancy_super(ic,i_mean_buoy,i_lev)                        &
                 + buoyancy_super(ic,i_mean_buoy,i_lev+1) )                    &
               * ( buoyancy_super(ic,i_height,i_lev+1)                         &
                 - buoyancy_super(ic,i_height,i_lev) )
    end do
  end do
  do ic = 1, n_points
    mean_buoy(ic) = mean_buoy(ic)                                              &
                  / ( buoyancy_super(ic,i_height,i_next(ic))                   &
                    - buoyancy_super(ic,i_height,i_prev) )
  end do

  ! Calc dry-density of parcel at previous model-level interface
  call calc_rho_dry( n_points,                                                 &
                     par_prev_mean(:,i_temperature),                           &
                     par_prev_mean(:,i_q_vap),                                 &
                     grid_prev(:,i_pressure),                                  &
                     rho_dry )
  ! Calculate sigma
  do ic = 1, n_points
    prev_cf_conv(ic) = prev_massflux_d(ic)                                     &
      / ( rho_dry(ic) * wind_w_buoy_fac                                        &
        * sqrt( max( ( gravity / prev_virt_temp(ic) )                          &
                   * mean_buoy(ic)                                             &
                   * prev_radius(ic), w_min*w_min ) ) )
  end do

  ! Same again to cf_conv at next model-level interface
  call calc_rho_dry( n_points,                                                 &
                     par_next_mean(:,i_temperature),                           &
                     par_next_mean(:,i_q_vap),                                 &
                     grid_next(:,i_pressure),                                  &
                     rho_dry )
  do ic = 1, n_points
    next_cf_conv(ic) = next_massflux_d(ic)                                     &
      / ( rho_dry(ic) * wind_w_buoy_fac                                        &
        * sqrt( max( ( gravity / next_virt_temp(ic) )                          &
                   * mean_buoy(ic)                                             &
                   * next_radius(ic), w_min*w_min ) ) )
  end do

else  ! ( l_cf_conv_buoy )

  ! For now we have prescribed  w' = w_fac / rho_wet
  !
  ! => sigma = Md/w_fac rho_wet/rho_dry
  !          = Md (1+q_tot) / w_fac
  !          = Mw / w_fac

  ! Calc parcel total-water at previous model-level interface
  call calc_q_tot( n_points, n_points_prev,                                    &
                   par_prev_mean(:,i_q_vap),                                   &
                   par_prev_mean(:,i_qc_first:i_qc_last), q_tot )
  ! Calculate sigma
  do ic = 1, n_points
    prev_cf_conv(ic) = prev_massflux_d(ic) * ( one + q_tot(ic) ) / wind_w_fac
  end do

  ! Same again to cf_conv at next model-level interface
  call calc_q_tot( n_points, n_points_next,                                    &
                   par_next_mean(:,i_q_vap),                                   &
                   par_next_mean(:,i_qc_first:i_qc_last), q_tot )
  do ic = 1, n_points
    next_cf_conv(ic) = next_massflux_d(ic) * ( one + q_tot(ic) ) / wind_w_fac
  end do

end if  ! ( l_cf_conv_buoy )

! Compute grid-mean convective liquid and ice mixing ratios
do ic = 1, n_points
  prev_q_cl_conv(ic) = prev_cf_conv(ic)*par_prev_mean(ic,i_q_cl)
  prev_q_cf_conv(ic) = prev_cf_conv(ic)*par_prev_mean(ic,i_q_cf)
  next_q_cl_conv(ic) = next_cf_conv(ic)*par_next_mean(ic,i_q_cl)
  next_q_cf_conv(ic) = next_cf_conv(ic)*par_next_mean(ic,i_q_cf)
end do

! Inflate the fraction without inflating the grid-mean
! condensate fields.  This is a temporary approximation to
! acount for in-plume w-variation; the slower-rising parts of
! the plume have greater area but lower condensate and so
! bring down the horizontal mean.
! When the in-parcel w-equation / PDF is fully implemented,
! we should be able to constrain this more accurately.
do ic = 1, n_points
  prev_cf_conv(ic) = prev_cf_conv(ic) * cf_conv_fac
  next_cf_conv(ic) = next_cf_conv(ic) * cf_conv_fac
end do

! Extract saturation height where applicable
do ic = 1, n_points
  if ( i_sat(ic) > 0 ) then
    sat_height(ic) = buoyancy_super(ic,i_height,i_sat(ic))
  else
    sat_height(ic) = zero
  end if
end do

! Interpolate the convective cloud condensate and fraction from
! the previous and next interfaces (where the parcel properties are
! defined) to the intervening full level...

if ( l_cv_cloudfrac ) then
  ! CoMorph is carrying cloud-fractions within the parcel;
  ! Pass in the in-parcel cloud-fractions from the parcel mean super-arrays

  call interp_diag_conv_cloud( n_points, n_points,                             &
                               n_points_prev, n_points_next,                   &
                               grid_prev(:,i_height),                          &
                               grid_next(:,i_height), sat_height,              &
                               prev_cf_conv, next_cf_conv,                     &
                               prev_q_cl_conv, next_q_cl_conv,                 &
                               prev_q_cf_conv, next_q_cf_conv,                 &
                               par_prev_mean(:,i_cf_liq:i_cf_bulk),            &
                               par_next_mean(:,i_cf_liq:i_cf_bulk),            &
                               convcloud_step )
  ! Note: i_cf_liq etc are the addresses of the cloud-fractions in the
  ! primary fields super-array; when cloud-fractions are not treated
  ! as primary fields (l_cv_cloudfrac=.false.), these are not set.

else
  ! CoMorph not carrying cloud-fractions within the parcel;
  ! Diagnose in-parcel cloud-fractions in separate arrays...

  ! Note: using the indices i_frac_liq etc, which are set even when
  ! cloud-fractions are not treated as primary fields, unlike i_cf_liq etc.
  allocate( par_prev_cloudfracs( n_points, i_frac_liq:i_frac_bulk ) )
  allocate( par_next_cloudfracs( n_points, i_frac_liq:i_frac_bulk ) )

  ! Call routine to diagnose in-parcel cloud-fractions based on q_cl, q_cf
  call set_par_cloudfrac( n_points, n_points,                                  &
                          par_prev_mean(:,i_q_cl), par_prev_mean(:,i_q_cf),    &
                          par_prev_cloudfracs )
  call set_par_cloudfrac( n_points, n_points,                                  &
                          par_next_mean(:,i_q_cl), par_next_mean(:,i_q_cf),    &
                          par_next_cloudfracs )

  ! Use diagnosed in-parcel cloud-fractions to interpolate to full-level
  call interp_diag_conv_cloud( n_points, n_points,                             &
                               n_points, n_points,                             &
                               grid_prev(:,i_height),                          &
                               grid_next(:,i_height), sat_height,              &
                               prev_cf_conv, next_cf_conv,                     &
                               prev_q_cl_conv, next_q_cl_conv,                 &
                               prev_q_cf_conv, next_q_cf_conv,                 &
                               par_prev_cloudfracs,                            &
                               par_next_cloudfracs,                            &
                               convcloud_step )

  deallocate( par_next_cloudfracs )
  deallocate( par_prev_cloudfracs )

end if  ! ( l_cv_cloudfrac )

! We have now calculated the convective cloud-fields from the current
! level-step performed by this call; now add these onto the total,
! which may be summed over sub-level-steps.  We need to weight by
! the fraction of the full level-step to get correct vertical sums
! of the convcloud fields on full-levels at the end.
do i_field = 1, n_convcloud
  do ic = 1, n_points
    convcloud_super(ic,i_field) = convcloud_super(ic,i_field)                  &
                                + convcloud_step(ic,i_field)                   &
                                  * frac_level_step(ic)
  end do
end do


return
end subroutine set_diag_conv_cloud


end module set_diag_conv_cloud_mod
