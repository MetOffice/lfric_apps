! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_turb_parcel_mod

implicit none

contains

! Subroutine to calculate parcel initial perturbations to winds and
! scalars, on full-levels (involves some fiddly interpolation from
! the fluxes defined on half-levels).
subroutine calc_turb_parcel( n_points, n_points_super, cmpr_init, k,           &
                             grid_km1, grid_kmh, grid_k,                       &
                             grid_kph, grid_kp1,                               &
                             fields_km1, fields_k, fields_kp1,                 &
                             turb_kmh, turb_kph, turb_pert_k )

use comorph_constants_mod, only: real_cvprec, zero, one, min_float,            &
                                 k_bot_conv, k_top_conv,                       &
                                 l_turb_par_gen, par_gen_qpert,                &
                                 name_length,                                  &
                                 i_check_bad_values_cmpr, i_check_bad_none
use cmpr_type_mod, only: cmpr_type
use grid_type_mod, only: n_grid, i_height
use fields_type_mod, only: i_wind_u, i_wind_v, i_wind_w,                       &
                           i_temperature, i_q_vap,                             &
                           i_q_cl, i_qc_first, i_qc_last,                      &
                           field_names
use turb_type_mod, only: n_turb
use set_cp_tot_mod, only: set_cp_tot
use lat_heat_mod, only: lat_heat_incr, i_phase_change_evp

use calc_fields_next_mod, only: calc_fields_next
use calc_turb_perts_mod, only: calc_turb_perts
use check_bad_values_mod, only: check_bad_values_cmpr

implicit none

! Number of points
integer, intent(in) :: n_points

! Size of input super-arrays (maybe larger than n_points,
! as arrays are dimensioned with the max size they could need
! on any level and reused for all levels)
integer, intent(in) :: n_points_super

! Indices of points being processed, for making error messages more informative
type(cmpr_type), intent(in) :: cmpr_init
integer, intent(in) :: k

! Height and pressure at k, model-level interfaces, and
! neighbouring full levels
real(kind=real_cvprec), intent(in) :: grid_km1                                 &
                                      ( n_points_super, n_grid )
real(kind=real_cvprec), intent(in) :: grid_kmh                                 &
                                      ( n_points_super, n_grid )
real(kind=real_cvprec), intent(in) :: grid_k                                   &
                                      ( n_points_super, n_grid )
real(kind=real_cvprec), intent(in) :: grid_kph                                 &
                                      ( n_points_super, n_grid )
real(kind=real_cvprec), intent(in) :: grid_kp1                                 &
                                      ( n_points_super, n_grid )

! Primary model-fields at level k and neighbouring full levels
real(kind=real_cvprec), intent(in) :: fields_km1                               &
                                ( n_points_super, i_wind_u:i_qc_last )
real(kind=real_cvprec), intent(in) :: fields_k                                 &
                                ( n_points_super, i_wind_u:i_qc_last )
real(kind=real_cvprec), intent(in) :: fields_kp1                               &
                                ( n_points_super, i_wind_u:i_qc_last )

! Super-arrays containing turbulence fields at the
! model-level interfaces
real(kind=real_cvprec), intent(in) :: turb_kmh                                 &
                                      ( n_points_super, n_turb )
real(kind=real_cvprec), intent(in) :: turb_kph                                 &
                                      ( n_points_super, n_turb )

! Turbulence-based perturbations to u,v,w,Tl,qt, at level k
real(kind=real_cvprec), intent(out) :: turb_pert_k                             &
                                       ( n_points, i_wind_u:i_q_vap )

! Toal heat capacity at k
real(kind=real_cvprec) :: cp_tot(n_points)
! Tl,qt at level k
real(kind=real_cvprec) :: tl_k(n_points)
real(kind=real_cvprec) :: qt_k(n_points)

! u,v,w,Tl,qt interpolated to the model-level interfaces
real(kind=real_cvprec) :: tl_kph(n_points)
real(kind=real_cvprec) :: qt_kph(n_points)
real(kind=real_cvprec) :: winds_kph(n_points,i_wind_u:i_wind_w)
real(kind=real_cvprec) :: tl_kmh(n_points)
real(kind=real_cvprec) :: qt_kmh(n_points)
real(kind=real_cvprec) :: winds_kmh(n_points,i_wind_u:i_wind_w)

! Turbulence-based perturbations on their native half-levels
real(kind=real_cvprec) :: turb_pert_kph                                        &
                          ( n_points, i_wind_u:i_q_vap )
real(kind=real_cvprec) :: turb_pert_kmh                                        &
                          ( n_points, i_wind_u:i_q_vap )

! Weight for interpolating the turbulence-based perturbations
real(kind=real_cvprec) :: interp

! String containing info to print in error messages
character(len=name_length) :: call_string
character(len=name_length) :: field_name
logical :: l_positive

! Loop counters
integer :: ic, i_field


if ( l_turb_par_gen ) then
  ! If using turbulence-based parcel perturbations...


  ! Calculate turbulence-based perturbations at the model-level interfaces
  call calc_turb_perts( n_points, n_points_super,                              &
                        turb_kph, turb_pert_kph )
  call calc_turb_perts( n_points, n_points_super,                              &
                        turb_kmh, turb_pert_kmh )

  if ( k > k_bot_conv .and. k < k_top_conv ) then
    ! If not at the bottom or top model-level
    ! (need fields at k+1 and k-1 for the interpolations we're about to do!)

    ! Calculate Tl,qt at k
    do ic = 1, n_points
      tl_k(ic) = fields_k(ic,i_temperature)
      qt_k(ic) = fields_k(ic,i_q_vap) + fields_k(ic,i_q_cl)
    end do
    call set_cp_tot( n_points, n_points_super, fields_k(:,i_q_vap),            &
                     fields_k(:,i_qc_first:i_qc_last),                         &
                     cp_tot )
    call lat_heat_incr( n_points, n_points, i_phase_change_evp,                &
                        cp_tot, tl_k,                                          &
                        dq=fields_k(:,i_q_cl) )

    ! Calculate u,v,w,Tl,qt interpolated to the model-level interfaces
    call calc_fields_next( n_points, n_points_super,                           &
                           grid_k, grid_kph, grid_kp1,                         &
                           fields_k, fields_kp1, tl_k, qt_k,                   &
                           tl_kph, qt_kph, winds_kph )
    call calc_fields_next( n_points, n_points_super,                           &
                           grid_k, grid_kmh, grid_km1,                         &
                           fields_k, fields_km1, tl_k, qt_k,                   &
                           tl_kmh, qt_kmh, winds_kmh )

    ! Interpolate turbulence-based perturbations onto level k.
    ! The method is to interpolate using a weight based on how close
    ! the level k value of a field is to its values at the neighbouring
    ! model-level interfaces.  e.g. if Tl_k is very similar to Tl_k-1/2,
    ! then the Tl perturbation is weighted towards its value at k-1/2.
    !
    ! This is to avoid nasty things happening at the tops of mixed-layers.
    ! If pert_kph has a much larger value, due to an inversion jump between
    ! k and k+1, but k itself is within the mixed-layer, then we shouldn't
    ! add the large inversion perturbation at k.

    ! Interpolations for Tl and qt
    do ic = 1, n_points
      interp = abs(tl_k(ic) - tl_kmh(ic))                                      &
        / max( abs(tl_k(ic) - tl_kmh(ic)) + abs(tl_kph(ic) - tl_k(ic)),        &
               min_float )
      turb_pert_k(ic,i_temperature)                                            &
        = (one-interp) * turb_pert_kmh(ic,i_temperature)                       &
        +      interp  * turb_pert_kph(ic,i_temperature)
      interp = abs(qt_k(ic) - qt_kmh(ic))                                      &
        / max( abs(qt_k(ic) - qt_kmh(ic)) + abs(qt_kph(ic) - qt_k(ic)),        &
               min_float )
      turb_pert_k(ic,i_q_vap)                                                  &
        = (one-interp) * turb_pert_kmh(ic,i_q_vap)                             &
        +      interp  * turb_pert_kph(ic,i_q_vap)
    end do

    ! Interpolations for winds
    do i_field = i_wind_u, i_wind_v
      do ic = 1, n_points
        interp = abs(fields_k(ic,i_field) - winds_kmh(ic,i_field))             &
          / max( abs(fields_k(ic,i_field) - winds_kmh(ic,i_field))             &
               + abs(winds_kph(ic,i_field) - fields_k(ic,i_field)),            &
                 min_float )
        turb_pert_k(ic,i_field)                                                &
          = (one-interp) * turb_pert_kmh(ic,i_field)                           &
          +      interp  * turb_pert_kph(ic,i_field)
      end do
    end do

    ! The above interpolation approach doesn't make sense for the
    ! w perturbation, as w' is not driven by the vertical gradient in
    ! grid-mean w, like for the other variables.
    ! Just use linear interpolation in height for the w perturbation
    do ic = 1, n_points
      interp = ( grid_k(ic,i_height) - grid_kmh(ic,i_height) )                 &
             / ( grid_kph(ic,i_height) - grid_kmh(ic,i_height) )
      turb_pert_k(ic,i_wind_w)                                                 &
        = (one-interp) * turb_pert_kmh(ic,i_wind_w)                            &
        +      interp  * turb_pert_kph(ic,i_wind_w)
    end do

  else  ! ( k > k_bot_conv .AND. k < k_top_conv )
    ! Special case of the top and bottom model-levels, where
    ! we don't have the fields at the levels above and below to
    ! use in the interpolation.
    ! Just interpolate all the perturbations linearly in height space.

    do i_field = i_wind_u, i_q_vap
      do ic = 1, n_points
        interp = ( grid_k(ic,i_height) - grid_kmh(ic,i_height) )               &
               / ( grid_kph(ic,i_height) - grid_kmh(ic,i_height) )
        turb_pert_k(ic,i_field)                                                &
          = (one-interp) * turb_pert_kmh(ic,i_field)                           &
          +      interp  * turb_pert_kph(ic,i_field)
      end do
    end do

  end if  ! ( k > k_bot_conv .AND. k < k_top_conv )


else  ! ( l_turb_par_gen )
  ! If NOT using turbulence-based parcel perturbations...


  ! Set u,v,w,T perturbations to zero
  do i_field = i_wind_u, i_temperature
    do ic = 1, n_points
      turb_pert_k(ic,i_field) = zero
    end do
  end do

  ! Set fixed fractional q perturbation
  do ic = 1, n_points
    turb_pert_k(ic,i_q_vap) = par_gen_qpert * fields_k(ic,i_q_vap)
  end do


end if  ! ( l_turb_par_gen )


if ( i_check_bad_values_cmpr > i_check_bad_none ) then
  ! Check outputs for bad values (NaN, Inf etc).

  call_string = "On output from calc_turb_parcel"

  ! Parcel initial perturbations
  do i_field = i_wind_u, i_q_vap
    ! Perturbations can be +ive or -ive, except w' which must be positive
    if ( i_field==i_wind_w) then
      l_positive = .true.
    else
      l_positive = .false.
    end if
    field_name = "turb_pert_k_" // trim(adjustl(field_names(i_field)))
    call check_bad_values_cmpr( cmpr_init, k, turb_pert_k(:,i_field),          &
                                call_string, field_name, l_positive )
  end do

end if  ! ( i_check_bad_values_cmpr > i_check_bad_none )


return
end subroutine calc_turb_parcel

end module calc_turb_parcel_mod
