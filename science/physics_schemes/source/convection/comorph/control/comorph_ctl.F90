! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module comorph_ctl_mod

implicit none

contains


! Top-level subroutine for the CoMorph convection scheme.
subroutine comorph_ctl( l_tracer, n_segments,                                  &
                        grid, turb, cloudfracs,                                &
                        fields_n, fields_np1,                                  &
                        comorph_diags )

! Get array sizes etc.
use comorph_constants_mod, only: real_hmprec, l_init_constants, name_length,   &
                                 nx_full, ny_full,                             &
                                 k_bot_conv, k_top_conv, k_top_init,           &
                                 n_tracers,                                    &
                                 l_cv_rain, l_cv_cf, l_cv_snow, l_cv_graup,    &
                                 l_turb_par_gen, i_check_bad_values_3d,        &
                                 i_check_turb_consistent, i_check_bad_none

use fields_type_mod, only: fields_type, l_init_fields_type_mod,                &
                           fields_set_addresses,                               &
                           fields_list_make, fields_list_clear,                &
                           field_names, field_positive, n_fields
use grid_type_mod, only: grid_type, grid_check_bad_values
use turb_type_mod, only: turb_type,                                            &
                         turb_list_make, turb_list_clear,                      &
                         n_turb, turb_names, turb_positive,                    &
                         turb_check_consistent
use cloudfracs_type_mod, only: cloudfracs_type,                                &
                               l_init_cloudfracs_type_mod,                     &
                               cloudfracs_set_addresses,                       &
                               cloudfracs_check_bad_values
use fields_2d_mod, only: l_init_fields_2d_mod, fields_2d_set_addresses
use comorph_diags_type_mod, only: comorph_diags_type,                          &
                                  comorph_diags_assign,                        &
                                  comorph_diags_dealloc

use set_dependent_constants_mod, only: set_dependent_constants
use check_bad_values_mod, only: check_bad_values_3d
use calc_layer_mass_mod, only: calc_layer_mass
use copy_field_mod, only: copy_field_3d
use calc_turb_diags_mod, only: calc_turb_diags
use calc_subregion_diags_mod, only: calc_subregion_diags
use calc_virt_temp_mod, only: calc_virt_temp_3d
use init_test_mod, only: init_test
use comorph_main_mod, only: comorph_main

implicit none


!----------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------

! Flag for whether to compute convective transport of
! passive tracers
logical, intent(in) :: l_tracer
! If this flag is set to true, the list of pointers to the
! tracers in the fields_np1 structure passed into comorph_ctl
! must be correctly allocated and assigned.
! Tracers at start-of-timestep (fields_n) are not required.

! Number of segments for segmented call
! (used to reduce the memory overhead, and/or for
!  shared-memory parallelisation using OMP)
integer, intent(in) :: n_segments
! This might want to change on different iterations of the
! solver; e.g. a smaller segment size might be desirable
! on the last iteration when tracer variables are included,
! giving a higher memory footprint per grid-column.

! Structure containing pointers to model grid and related
! fields; contains model-level heights, pressures and
! dry-density
type(grid_type), intent(in) :: grid

! Structure containing pointers to turbulence fields passed
! in from the boundary-layer scheme
type(turb_type), intent(in out) :: turb

! Structure containing cloud fractions, in the event that they
! are diagnostic rather than prognostic and so aren't included
! in the primary fields structure.
! Also contains the diagnosed convective cloud fraction and
! water content, which are inout:
! in: convective cloud at start of timestep (from previous step).
! out: updated convective cloud from this call.
type(cloudfracs_type), intent(in out) :: cloudfracs

! Structure containing pointers to input start-of-timestep fields
type(fields_type), intent(in out) :: fields_n
! (needs intent(inout) so that the contained fields list
!  can be setup).

! Structure containing pointers to the _np1 fields:
! in: updated with any increments added before convection
! out: updated by convection as well
type(fields_type), intent(in out) :: fields_np1

! Note: all the input primary fields need to be on the same
! levels.  In the UM, u,v are on different levels (rho-levels)
! to the rest of the primary fields (theta-levels).
! They will need to be interpolated before input and output to
! this routine.

! Structure containing pointers to output diagnostics and
! miscellaneous outputs required elsewhere in the model
type(comorph_diags_type), intent(in out) :: comorph_diags
! (holds outputs, but intent inout to retain the assigned status
!  of the contained pointers on input)

!----------------------------------------------------------------
! Local variables...
!----------------------------------------------------------------

! The following fields are calculated directly from the input
! fields on the host model's grid, hence declared with
! the same precision as the host model:

! Dry-mass on each model-level, per unit surface area
real(kind=real_hmprec) :: layer_mass( nx_full, ny_full,                        &
                                      k_bot_conv:k_top_conv )
! This is rho_dry * grid-cell volume per unit surface area.
! For a Cartesian grid, this is just rho_dry * delta z
! But for spherical coordinates, we also need to account for the
! grid area getting wider with height.

! Environment virtual temperature profile (full 3-D array):
! Start-of-timestep
real(kind=real_hmprec) :: virt_temp_n( nx_full, ny_full,                       &
                                       k_bot_conv:k_top_conv )
! Latest fields
real(kind=real_hmprec) :: virt_temp_np1( nx_full, ny_full,                     &
                                         k_bot_conv:k_top_conv )

! 3-D mask of points where convective initiation mass-sources
! for updrafts or downdrafts might be possible
logical :: l_init_poss( nx_full, ny_full,                                      &
                        k_bot_conv:k_top_init )


! Number of points on each segment
integer :: seg_n_points(n_segments)
! ij-index of the last point on each segment
integer :: seg_ij_last(0:n_segments)


! Local flag input to fields_list_make to indicate no
! tracers exist in the start-of-timestep fields structure
logical, parameter :: l_tracer_false = .false.

! Total number of fields to act upon, to pass into conv_sweep_ctl
integer :: n_fields_tot

! Flag for input to comorph_diags_assign
logical :: l_count_diags

! Lower and upper bounds of arrays
integer :: lb_rho(3), ub_rho(3)
integer :: lb_rs(2), ub_rs(2)
integer :: lb_h(3), ub_h(3)
integer :: lb_t(3), ub_t(3)
integer :: lb_v(3), ub_v(3)
integer :: lb_l(3), ub_l(3)
integer :: lb_r(3), ub_r(3)
integer :: lb_f(3), ub_f(3)
integer :: lb_s(3), ub_s(3)
integer :: lb_g(3), ub_g(3)
integer :: lb_1(3), ub_1(3)
integer :: lb_2(3), ub_2(3)

! String indicating where in the code bad-value checks are done
character(len=name_length) :: where_string

! Loop counters
integer :: i_seg, i_field


!----------------------------------------------------------------
! 1) Initialisations
!----------------------------------------------------------------

! Set the values of constants that depend on other constants
! (if not already set)
if ( .not. l_init_constants )  call set_dependent_constants()

! Set super-array addresses for primary fields
! (if not already set)
if ( .not. l_init_fields_type_mod )  call fields_set_addresses()

! Setup lists of pointers to the primary fields
call fields_list_make( l_tracer_false, fields_n )
call fields_list_make( l_tracer, fields_np1 )

! Set super-array addresses for convective cloud fields
if ( .not. l_init_cloudfracs_type_mod )  call cloudfracs_set_addresses()

! Set super-array addresses for 2D fields
if ( .not. l_init_fields_2d_mod )  call fields_2d_set_addresses()

! Setup list of pointers to turbulence fields
if ( l_turb_par_gen  )  call turb_list_make( turb )

! Set number of fields including tracers; (an array size
! required by conv_sweep_ctl)
n_fields_tot = n_fields
if ( l_tracer )  n_fields_tot = n_fields_tot + n_tracers


! Check the fields to ensure no bad values on input
if ( i_check_bad_values_3d > i_check_bad_none ) then

  ! Check start-of-timestep fields:
  where_string = "On input to CoMorph: start-of-timestep fields:"
  do i_field = 1, n_fields
    lb_1 = lbound( fields_n % list(i_field)%pt )
    ub_1 = ubound( fields_n % list(i_field)%pt )
    call check_bad_values_3d( lb_1, ub_1, fields_n%list(i_field)%pt,           &
                              where_string, field_names(i_field),              &
                              field_positive(i_field) )
  end do
  ! Check latest fields:
  where_string = "On input to CoMorph: latest fields:"
  do i_field = 1, n_fields_tot
    lb_1 = lbound( fields_np1 % list(i_field)%pt )
    ub_1 = ubound( fields_np1 % list(i_field)%pt )
    call check_bad_values_3d( lb_1, ub_1, fields_np1%list(i_field)%pt,         &
                              where_string, field_names(i_field),              &
                              field_positive(i_field) )
  end do

  ! Check input grid fields
  where_string = "On input to CoMorph: model-grid fields:"
  call grid_check_bad_values( grid, where_string )

  ! Check cloud and rain fractions
  where_string = "On input to CoMorph: cloud-fraction fields:"
  call cloudfracs_check_bad_values( cloudfracs, where_string )

  ! Also check turbulence fields if used
  if ( l_turb_par_gen ) then
    where_string = "On input to CoMorph: turbulence fields:"
    do i_field = 1, n_turb
      lb_1 = lbound( turb % list(i_field)%pt )
      ub_1 = ubound( turb % list(i_field)%pt )
      call check_bad_values_3d( lb_1, ub_1, turb%list(i_field)%pt,             &
                                where_string,turb_names(i_field),              &
                                turb_positive(i_field),                        &
                                l_half=.true., l_init=.true. )
    end do
  end if

end if  ! ( i_check_bad_values_3d > i_check_bad_none )

if ( l_turb_par_gen .and. i_check_turb_consistent > i_check_bad_none ) then
  ! Check consistency between the different turbulence fields
  call turb_check_consistent( turb )
end if

! Initialise comorph diagnostics...
! Note that the "diagnostics" system needs to be active
! even when not on the last solver outer loop iteration,
! because some "diagnostic" fields might be required by other
! parts of the model, not just for diagnostic output.
! Also note that the diagnostic switches maybe different
! on different calls to comorph_ctl within the same run
! (e.g. if a diag is only requested on certain timesteps).
! Therefore the following code to count the number of required
! fields must be called every time.

! Count the number of active diagnostics
l_count_diags = .true.
call comorph_diags_assign( l_count_diags, comorph_diags )

! Now that the required number of diags is known...
! (stored in comorph_diags % n_diags)

! If any diagnostics requested
if ( comorph_diags % n_diags > 0 ) then

  ! Allocate list of pointers to active diagnostics and
  ! point them at all the required arrays
  l_count_diags = .false.
  call comorph_diags_assign( l_count_diags, comorph_diags )
  ! This routine also initialises the diagnostics to zero.

end if

! Calculate the layer masses at all points, for conservation
! (note need to pass in the bounds of the full 3-D arrays,
!  in case they have halos).
lb_rho = lbound(grid%rho_dry)
ub_rho = ubound(grid%rho_dry)
lb_rs = lbound(grid%r_surf)
ub_rs = ubound(grid%r_surf)
lb_h = lbound(grid%height_half)
ub_h = ubound(grid%height_half)
call calc_layer_mass( lb_rho, ub_rho, grid%rho_dry, lb_rs, ub_rs, grid%r_surf, &
                      lb_h, ub_h, grid%height_half,                            &
                      layer_mass )

! Initialise array bounds for optional water-species
lb_r = [1,1,1]
ub_r = [1,1,1]
lb_f = [1,1,1]
ub_f = [1,1,1]
lb_s = [1,1,1]
ub_s = [1,1,1]
lb_g = [1,1,1]
ub_g = [1,1,1]

! Find array bounds needed for virtual temperature calculation routine
lb_t = lbound(fields_n % temperature)
ub_t = ubound(fields_n % temperature)
lb_v = lbound(fields_n % q_vap)
ub_v = ubound(fields_n % q_vap)
lb_l = lbound(fields_n % q_cl)
ub_l = ubound(fields_n % q_cl)
if ( l_cv_rain ) then
  lb_r = lbound(fields_n % q_rain)
  ub_r = ubound(fields_n % q_rain)
end if
if ( l_cv_cf ) then
  lb_f = lbound(fields_n % q_cf)
  ub_f = ubound(fields_n % q_cf)
end if
if ( l_cv_snow ) then
  lb_s = lbound(fields_n % q_snow)
  ub_s = ubound(fields_n % q_snow)
end if
if ( l_cv_graup ) then
  lb_g = lbound(fields_n % q_graup)
  ub_g = ubound(fields_n % q_graup)
end if
! Calculate start-of-timestep virtual temperature profile
call calc_virt_temp_3d( lb_t, ub_t, fields_n % temperature,                    &
                        lb_v, ub_v, fields_n % q_vap,                          &
                        lb_l, ub_l, fields_n % q_cl,                           &
                        lb_r, ub_r, fields_n % q_rain,                         &
                        lb_f, ub_f, fields_n % q_cf,                           &
                        lb_s, ub_s, fields_n % q_snow,                         &
                        lb_g, ub_g, fields_n % q_graup,                        &
                        virt_temp_n )

! Find array bounds for the latest fields
lb_t = lbound(fields_np1 % temperature)
ub_t = ubound(fields_np1 % temperature)
lb_v = lbound(fields_np1 % q_vap)
ub_v = ubound(fields_np1 % q_vap)
lb_l = lbound(fields_np1 % q_cl)
ub_l = ubound(fields_np1 % q_cl)
if ( l_cv_rain ) then
  lb_r = lbound(fields_np1 % q_rain)
  ub_r = ubound(fields_np1 % q_rain)
end if
if ( l_cv_cf ) then
  lb_f = lbound(fields_np1 % q_cf)
  ub_f = ubound(fields_np1 % q_cf)
end if
if ( l_cv_snow ) then
  lb_s = lbound(fields_np1 % q_snow)
  ub_s = ubound(fields_np1 % q_snow)
end if
if ( l_cv_graup ) then
  lb_g = lbound(fields_np1 % q_graup)
  ub_g = ubound(fields_np1 % q_graup)
end if
! Calculate latest virtual temperature profile
call calc_virt_temp_3d( lb_t, ub_t, fields_np1 % temperature,                  &
                        lb_v, ub_v, fields_np1 % q_vap,                        &
                        lb_l, ub_l, fields_np1 % q_cl,                         &
                        lb_r, ub_r, fields_np1 % q_rain,                       &
                        lb_f, ub_f, fields_np1 % q_cf,                         &
                        lb_s, ub_s, fields_np1 % q_snow,                       &
                        lb_g, ub_g, fields_np1 % q_graup,                      &
                        virt_temp_np1 )


!----------------------------------------------------------------
! 2) Compute any diagnostics not calculated inside the convection
!----------------------------------------------------------------

! If layer-mass requested as a diagnostic, copy into diag array
if ( comorph_diags % layer_mass % flag ) then
  lb_1 = lbound( layer_mass )
  ub_1 = ubound( layer_mass )
  lb_2 = lbound( comorph_diags%layer_mass%field_3d )
  ub_2 = ubound( comorph_diags%layer_mass%field_3d )
  call copy_field_3d( lb_1, ub_1, layer_mass,                                  &
                      lb_2, ub_2, comorph_diags%layer_mass%field_3d )
end if

! If turbulence-based parcel perturbations requested as
! diagnostics, calculate them here (these are calculated
! in conv_genesis_ctl, but only at initiation mass-source
! points; we might want the diags at all points).
if ( l_turb_par_gen .and.                                                      &
  ! Diags only available if turbulence-based properties are on.
     ( comorph_diags % turb_fields_pert % n_diags > 0 .or.                     &
       comorph_diags % turb_radius % flag ) ) then
  call calc_turb_diags( turb, comorph_diags )
end if

! If diagnostics for separate dry, cloudy, icy/rainy sub-grid
! regions are requested, calculate them here
if ( l_cv_rain .and.                                                           &
     comorph_diags % n_diags_subregions > 0 ) then
  call calc_subregion_diags( grid, fields_np1, cloudfracs,                     &
                             comorph_diags % subregion_diags )
end if


!----------------------------------------------------------------
! 3) Find mask of points where convective initiation is possible
!----------------------------------------------------------------

! Currently this sets the mask to true at any points with
! negative dry N^2 (dry statically unstable), or with some
! condensate of any type (will be used later to compute moist
! N^2 and see if that is negative).
call init_test( grid, fields_n, virt_temp_n, l_init_poss )
! l_init_poss stores 3-D mask of points to test for initiation.


!----------------------------------------------------------------
! 4) Setup segments for OMP-parallelised loop over points...
!----------------------------------------------------------------

! The convective parcel ascents and descents are performed
! sequentially in the vertical, so for optimal use of
! vectorisation, multiple columns need to be processed
! together, with calculations done on lists of points on
! the same level.  However, processing all the convecting
! points at a given level together like this would require
! a lot of memory for all the work arrays.  We therefore
! divide the horizontal domain into "segments", and process
! a sub-set of the columns together in each segment.
! Further, multiple segments can be processed in parallel
! using shared-memory OMP, to greatly speed-up the code.
! Also, oranges are really tasty, and they're segmented,
! so there you go.

! So, we need to divvy-up the grid-columns among the segments....

! Fictitious segment number 0, so that we can always
! get ij_first(i_seg) = ij_last(i_seg-1) + 1
seg_ij_last(0) = 0

! Set the number of points on each segment
do i_seg = 1, n_segments
  ! This formula divvies up the points among the segments
  ! as evenly as possible:
  ! Points on this segment = number of not yet assigned points
  !                        / number of not yet assigned segments
  seg_n_points(i_seg) = ( nx_full*ny_full                                      &
                        - seg_ij_last(i_seg-1) )                               &
                      / ( n_segments - (i_seg-1) )
  seg_ij_last(i_seg) = seg_ij_last(i_seg-1) + seg_n_points(i_seg)
end do


! BEGIN THE GREAT PARALLEL LOOP OVER SEGMENTS...
! Note that often the OMP threads will be very poorly
! load-balanced (in terms of both memory and CPU), since
! convection tends to be highly unevenly distributed in space.
! Use of dynamic scheduling should help with this.

!$OMP  PARALLEL DEFAULT(none)                                                  &
!$OMP  SHARED(  n_segments, seg_n_points, seg_ij_last,                         &
!$OMP           n_fields_tot, l_tracer,                                        &
!$OMP           grid, turb, cloudfracs, fields_np1,                            &
!$OMP           layer_mass, virt_temp_n, virt_temp_np1,                        &
!$OMP           l_init_poss, comorph_diags )                                   &
!$OMP  private( i_seg )
!$OMP do SCHEDULE(DYNAMIC)
do i_seg = 1, n_segments

  ! Call main comorph routine; calls various other routines to diagnose
  ! initiation mass-sources, perform updraft and downdraft parcel sweeps,
  ! and update primary fields with convective increments
  ! (all only at points on the current segment).
  call comorph_main( seg_n_points(i_seg),                                      &
                     seg_ij_last(i_seg-1)+1, seg_ij_last(i_seg),               &
                     n_fields_tot, l_tracer,                                   &
                     grid, turb, cloudfracs, fields_np1,                       &
                     layer_mass, virt_temp_n, virt_temp_np1,                   &
                     l_init_poss, comorph_diags )

end do  ! i_seg = 1, n_segments
!$OMP end do NOWAIT
!$OMP end PARALLEL


! Check the fields to ensure no bad values on output
if ( i_check_bad_values_3d > i_check_bad_none ) then

  ! Check latest fields:
  where_string = "On output from CoMorph: latest fields:"
  do i_field = 1, n_fields_tot
    lb_1 = lbound( fields_np1 % list(i_field)%pt )
    ub_1 = ubound( fields_np1 % list(i_field)%pt )
    call check_bad_values_3d( lb_1, ub_1, fields_np1%list(i_field)%pt,         &
                              where_string, field_names(i_field),              &
                              field_positive(i_field) )
  end do

  ! Check cloud and rain fractions (the convective cloud fields
  ! will have been updated).
  where_string = "On output from CoMorph: cloud-fraction fields:"
  call cloudfracs_check_bad_values( cloudfracs, where_string )

end if


! Deallocate miscellaneous allocatables in the diag structures
call comorph_diags_dealloc( comorph_diags )

! Clear list of pointers to turbulence fields
if ( l_turb_par_gen )  call turb_list_clear( turb )

! Clear lists of pointers to the primary fields
call fields_list_clear( fields_np1 )
call fields_list_clear( fields_n )


return
end subroutine comorph_ctl


end module comorph_ctl_mod
