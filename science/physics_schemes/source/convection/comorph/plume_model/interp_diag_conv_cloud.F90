! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module interp_diag_conv_cloud_mod

implicit none

contains

! Subroutine to interpolate the diagnosed convective cloud
! properties from the previous and next model-level interfaces
! (where the parcel properties are held) onto the intervening
! full model-level.
subroutine interp_diag_conv_cloud( n_points, n_points_res,                     &
                                   n_points_prev, n_points_next,               &
                                   prev_height, next_height, sat_height,       &
                                   prev_cf_conv, next_cf_conv,                 &
                                   prev_q_cl_conv, next_q_cl_conv,             &
                                   prev_q_cf_conv, next_q_cf_conv,             &
                                   par_prev_cloudfracs, par_next_cloudfracs,   &
                                   convcloud )

use comorph_constants_mod, only: real_cvprec, zero, half,                      &
                                 i_convcloud, i_convcloud_mph
use cloudfracs_type_mod, only: i_frac_liq, i_frac_ice, i_frac_bulk,            &
                               n_convcloud,                                    &
                               i_frac_liq_conv, i_frac_ice_conv,               &
                               i_frac_bulk_conv,                               &
                               i_q_cl_conv, i_q_cf_conv

implicit none

! Number of points
integer, intent(in) :: n_points

! Size of primary fields super-arrays
integer, intent(in) :: n_points_res
integer, intent(in) :: n_points_prev
integer, intent(in) :: n_points_next

! Heights of the previous and next model-level interfaces,
! and the saturation-point
real(kind=real_cvprec), intent(in) :: prev_height(n_points)
real(kind=real_cvprec), intent(in) :: next_height(n_points)
real(kind=real_cvprec), intent(in) :: sat_height(n_points)

! Convective area fraction at prev and next
real(kind=real_cvprec), intent(in) :: prev_cf_conv(n_points)
real(kind=real_cvprec), intent(in) :: next_cf_conv(n_points)

! Grid-mean convective liquid and ice cloud mixing-ratios at prev and next
real(kind=real_cvprec), intent(in) :: prev_q_cl_conv(n_points)
real(kind=real_cvprec), intent(in) :: next_q_cl_conv(n_points)
real(kind=real_cvprec), intent(in) :: prev_q_cf_conv(n_points)
real(kind=real_cvprec), intent(in) :: next_q_cf_conv(n_points)

! In-parcel cloud-fractions at prev and next
real(kind=real_cvprec), intent(in) :: par_prev_cloudfracs                      &
                        ( n_points_prev, i_frac_liq:i_frac_bulk )
real(kind=real_cvprec), intent(in) :: par_next_cloudfracs                      &
                        ( n_points_next, i_frac_liq:i_frac_bulk )

! Super-array containing diagnosed convective cloud properties
real(kind=real_cvprec), intent(out) :: convcloud                               &
                                       ( n_points_res, n_convcloud )

! Weight for vertical interpolation
real(kind=real_cvprec) :: interp
! Convective ice-cloud fraction, for check on the bulk cloud fraction
real(kind=real_cvprec) :: frac_ice_tmp

! Loop counter
integer :: ic


! Find mean convective cloud properties between prev and next
! Note:grid-mean fraction of convective cloud = convective fraction
! * in-convection cloud-fraction
do ic = 1, n_points

  convcloud(ic,i_frac_liq_conv)                                                &
    = half * prev_cf_conv(ic) * par_prev_cloudfracs(ic,i_frac_liq)             &
    + half * next_cf_conv(ic) * par_next_cloudfracs(ic,i_frac_liq)

  convcloud(ic,i_frac_bulk_conv)                                               &
    = half * prev_cf_conv(ic) * par_prev_cloudfracs(ic,i_frac_bulk)            &
    + half * next_cf_conv(ic) * par_next_cloudfracs(ic,i_frac_bulk)

  convcloud(ic,i_q_cl_conv)                                                    &
    = half * prev_q_cl_conv(ic)                                                &
    + half * next_q_cl_conv(ic)

end do

if ( i_convcloud == i_convcloud_mph ) then
  ! If outputting convective ice cloud properties as well,
  ! repeat the above for the ice cloud properties...

  do ic = 1, n_points

    convcloud(ic,i_frac_ice_conv)                                              &
      = half * prev_cf_conv(ic) * par_prev_cloudfracs(ic,i_frac_ice)           &
      + half * next_cf_conv(ic) * par_next_cloudfracs(ic,i_frac_ice)

    convcloud(ic,i_q_cf_conv)                                                  &
      = half * prev_q_cf_conv(ic)                                              &
      + half * next_q_cf_conv(ic)

  end do

end if  ! ( i_convcloud == i_convcloud_mph )

! If at the saturation height, reset interp for liquid cloud to
! volume-average fraction of the model-level containing liquid cloud
do ic = 1, n_points
  if ( sat_height(ic) > zero ) then

    if ( par_prev_cloudfracs(ic,i_frac_liq) == zero .and.                      &
         par_next_cloudfracs(ic,i_frac_liq) > zero ) then

      interp = ( next_height(ic) - sat_height(ic) )                            &
             / ( next_height(ic) - prev_height(ic) )

      convcloud(ic,i_frac_liq_conv) = interp * next_cf_conv(ic)                &
                                      * par_next_cloudfracs(ic,i_frac_liq)
      convcloud(ic,i_q_cl_conv)     = interp * next_q_cl_conv(ic)

    else if ( par_prev_cloudfracs(ic,i_frac_liq) > zero .and.                  &
              par_next_cloudfracs(ic,i_frac_liq) == zero ) then

      interp = ( sat_height(ic) - prev_height(ic) )                            &
             / ( next_height(ic) - prev_height(ic) )

      convcloud(ic,i_frac_liq_conv) = interp * prev_cf_conv(ic)                &
                                      * par_prev_cloudfracs(ic,i_frac_liq)
      convcloud(ic,i_q_cl_conv)     = interp * prev_q_cl_conv(ic)

    end if

    ! Adjust bulk convective cloud fraction consistent with the new
    ! liquid convective cloud fraction:
    ! Find convective ice cloud fraction
    frac_ice_tmp                                                               &
      = half * prev_cf_conv(ic) * par_prev_cloudfracs(ic,i_frac_ice)           &
      + half * next_cf_conv(ic) * par_next_cloudfracs(ic,i_frac_ice)
    ! We must have max(frac_liq,frac_ice) <= frac_bulk <= frac_liq+frac_ice
    convcloud(ic,i_frac_bulk_conv)                                             &
      = max( min( convcloud(ic,i_frac_bulk_conv),                              &
                  convcloud(ic,i_frac_liq_conv) + frac_ice_tmp ),              &
             max( convcloud(ic,i_frac_liq_conv),  frac_ice_tmp ) )

  end if
end do


return
end subroutine interp_diag_conv_cloud

end module interp_diag_conv_cloud_mod
