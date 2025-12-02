! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_fields_next_mod

implicit none

contains

! Subroutine to calculate environment Tl,qt,u,v,w interpolated to
! the next model-level interface.
subroutine calc_fields_next( n_points, n_points_super,                         &
                             grid_k, grid_next, grid_kpdk,                     &
                             fields_k, fields_kpdk, tl_k, qt_k,                &
                             tl_next, qt_next, winds_next )

use comorph_constants_mod, only: real_cvprec, one
use grid_type_mod, only: n_grid, i_height, i_pressure
use fields_type_mod, only: i_wind_u, i_wind_w,                                 &
                           i_temperature, i_q_vap,                             &
                           i_q_cl, i_qc_first, i_qc_last

use dry_adiabat_mod, only: dry_adiabat
use set_cp_tot_mod, only: set_cp_tot
use lat_heat_mod, only: lat_heat_incr, i_phase_change_evp


implicit none

! Number of points
integer, intent(in) :: n_points

! Size of input super-arrays (maybe larger than n_points,
! as arrays are dimensioned with the max size they could need
! on any level and reused for all levels)
integer, intent(in) :: n_points_super

! Height and pressure at k, next model-level interface, next full level
real(kind=real_cvprec), intent(in) :: grid_k                                   &
                                      ( n_points_super, n_grid )
real(kind=real_cvprec), intent(in) :: grid_next                                &
                                      ( n_points_super, n_grid )
real(kind=real_cvprec), intent(in) :: grid_kpdk                                &
                                      ( n_points_super, n_grid )

! Primary model-fields at level k and next full level
real(kind=real_cvprec), intent(in) :: fields_k                                 &
                                ( n_points_super, i_wind_u:i_qc_last )
real(kind=real_cvprec), intent(in) :: fields_kpdk                              &
                                ( n_points_super, i_wind_u:i_qc_last )

! Liquid water temperature and q_vap + q_cl at level k
real(kind=real_cvprec), intent(in) :: tl_k(n_points)
real(kind=real_cvprec), intent(in) :: qt_k(n_points)

! Liquid water temperature and q_vap + q_cl at next model-level interface
real(kind=real_cvprec), intent(out) :: tl_next(n_points)
real(kind=real_cvprec), intent(out) :: qt_next(n_points)

! Winds at next model-level interface
real(kind=real_cvprec), intent(out) :: winds_next                              &
                                       ( n_points, i_wind_u:i_wind_w )


! Factor for subsiding next level temperature to level k
real(kind=real_cvprec) :: exner_ratio(n_points)
! Total heat capacity for phase-change
real(kind=real_cvprec) :: cp_tot(n_points)
! Vertical interpolation weight
real(kind=real_cvprec) :: interp(n_points)

! Loop counters
integer :: ic, i_field


! Calculate temperature the air at the next full level (k+dk)
! would have at level k, and the vapour + liquid-cloud
call dry_adiabat( n_points, n_points_super,                                    &
                  grid_kpdk(:,i_pressure),                                     &
                  grid_k(:,i_pressure),                                        &
                  fields_kpdk(:,i_q_vap),                                      &
                  fields_kpdk(:,i_qc_first:i_qc_last),                         &
                  exner_ratio )
do ic = 1, n_points
  tl_next(ic) = fields_kpdk(ic,i_temperature) * exner_ratio(ic)
  qt_next(ic) = fields_kpdk(ic,i_q_vap) + fields_kpdk(ic,i_q_cl)
end do

! Convert temperature to liquid-water temperature
call set_cp_tot( n_points, n_points_super,                                     &
                 fields_kpdk(:,i_q_vap),                                       &
                 fields_kpdk(:,i_qc_first:i_qc_last),                          &
                 cp_tot )
call lat_heat_incr( n_points, n_points, i_phase_change_evp,                    &
                    cp_tot, tl_next,                                           &
                    dq=fields_kpdk(:,i_q_cl) )

! Store interpolation weight to get to next model-level interface
do ic = 1, n_points
  interp(ic) = ( grid_next(ic,i_height) - grid_k(ic,i_height) )                &
             / ( grid_kpdk(ic,i_height) - grid_k(ic,i_height) )
end do

! Interpolate Tl and qt to the next model-level interface
! (tl_next, qt_next currently store values at the next full-level)
do ic = 1, n_points
  tl_next(ic) = (one-interp(ic)) * tl_k(ic)                                    &
              +      interp(ic)  * tl_next(ic)
  qt_next(ic) = (one-interp(ic)) * qt_k(ic)                                    &
              +      interp(ic)  * qt_next(ic)
end do

! Interpolate winds onto the next model-level interface
do i_field = i_wind_u, i_wind_w
  do ic = 1, n_points
    winds_next(ic,i_field) = (one-interp(ic)) * fields_k(ic,i_field)           &
                           +      interp(ic)  * fields_kpdk(ic,i_field)
  end do
end do


return
end subroutine calc_fields_next

end module calc_fields_next_mod
