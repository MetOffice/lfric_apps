! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfric2um_initialise_um_mod

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only: int64, real64
use, intrinsic :: iso_c_binding, only: c_bool

! Shumlib modules
use f_shum_file_mod, only: shum_file_type
use f_shum_field_mod, only: shum_field_type
use f_shum_ff_status_mod, only: shum_ff_status_type

! lfricinputs modules
use lfricinp_check_shumlib_status_mod, only: shumlib

use lfricinp_um_parameters_mod, only:                                          &
    um_imdi, um_rmdi, rh_deltaEW, rh_deltaNS,                                  &
    rh_baselat, rh_baselong, rh_polelat, rh_polelong, rh_model_top,            &
    ih_row_length, ih_rows, ih_model_levels, ih_wet_levels, ih_soilT_levels,   &
    ih_cloud_levels, ih_tracer_levels, ih_boundary_levels, ih_N_types,         &
    ih_height_gen, ih_1_c_rho_level, ih_land_points, ih_ozone_levels,          &
    ih_soilQ_levels, ih_convect_levels, ldc_eta_theta, ldc_eta_rho,            &
    ldc_Zsea_theta, ldc_C_theta, ldc_Zsea_rho, ldc_C_rho

use lfricinp_grid_type_mod, only: lfricinp_grid_type

! lfric modules
use mesh_mod, only: mesh_type
use log_mod,  only: log_event, log_scratch_space,                              &
                    LOG_LEVEL_INFO, LOG_LEVEL_ERROR

implicit none

type(shum_file_type), save, public :: um_output_file

private
public :: lfric2um_initialise_um

contains

!------------------------------------------------------------------

subroutine lfric2um_initialise_um()
! Description:
!  Create UM file object using shumlib Populates metadata using
!  information from UM grid object
use lfricinp_um_grid_mod, only: um_grid
use lfricinp_lfric_driver_mod, only: mesh
use lfricinp_stashmaster_mod, only: lfricinp_read_stashmaster
use lfric2um_namelists_mod, only: lfric2um_config

implicit none

character(len=*), parameter :: routinename='lfric2um_initialise_um'
integer(kind=int64) :: num_lookup
! Create UM file
! Future work todo, calculate lookups based on requested fields, for
! now use default 4096
num_lookup = 4096
call shumlib(routinename//'::open_file', &
     um_output_file%open_file(lfric2um_config%output_filename, &
     num_lookup=num_lookup, overwrite=.true._C_BOOL))

call lfric2um_set_fixed_length_header(um_output_file, um_grid)

call lfric2um_set_integer_constants(um_output_file, um_grid, mesh)

call lfric2um_set_real_constants(um_output_file, um_grid, mesh)

call lfric2um_set_level_dep_constants(um_output_file, mesh)

! Write output header
call shumlib(routinename//'::write_header', &
         um_output_file%write_header())

end subroutine lfric2um_initialise_um

!------------------------------------------------------------------

subroutine lfric2um_set_fixed_length_header(um_output_file, um_grid)

use f_shum_fieldsfile_mod, only: f_shum_fixed_length_header_len

use f_shum_fixed_length_header_indices_mod, only: data_set_format_version, &
                                                  sub_model,               &
                                                  vert_coord_type,         &
                                                  horiz_grid_type,         &
                                                  dataset_type,            &
                                                  run_identifier,          &
                                                  calendar,                &
                                                  projection_number,       &
                                                  model_version,           &
                                                  grid_staggering,         &
                                                  sub_model,               &
                                                  t1_year,                 &
                                                  t1_month,                &
                                                  t1_day,                  &
                                                  t1_hour,                 &
                                                  t1_minute,               &
                                                  t1_second,               &
                                                  t2_year,                 &
                                                  t2_month,                &
                                                  t2_day,                  &
                                                  t2_hour,                 &
                                                  t2_minute,               &
                                                  t2_second,               &
                                                  t3_year,                 &
                                                  t3_month,                &
                                                  t3_day,                  &
                                                  t3_hour,                 &
                                                  t3_minute,               &
                                                  t3_second

use lfric2um_namelists_mod, only: lfric2um_config

implicit none

type(shum_file_type), intent(in out) :: um_output_file
type(lfricinp_grid_type), intent(in) :: um_grid

integer(kind=int64) :: fixed_length_header(            &
                           f_shum_fixed_length_header_len) = um_imdi
character(len=*), parameter :: routinename = 'lfric2um_set_fixed_length_header'

fixed_length_header(data_set_format_version) = 20
fixed_length_header(sub_model) = 1 ! Atmosphere
fixed_length_header(vert_coord_type) = 5 ! Charney Phillips
fixed_length_header(horiz_grid_type) = um_grid % horiz_grid_type
fixed_length_header(dataset_type) = 1 ! Dump
fixed_length_header(run_identifier) = 0 ! Undefined
fixed_length_header(grid_staggering) = 6 ! Only Arakawa C - ENDGame supported
fixed_length_header(calendar) = 1 ! Gregorian
fixed_length_header(projection_number) = 1 ! Lat/lon grid
fixed_length_header(model_version) = lfric2um_config%um_version_int

! Data time of forecast, set equal to validity time
fixed_length_header(t1_year) = lfric2um_config%dump_validity_time(1)
fixed_length_header(t1_month) = lfric2um_config%dump_validity_time(2)
fixed_length_header(t1_day) = lfric2um_config%dump_validity_time(3)
fixed_length_header(t1_hour) = lfric2um_config%dump_validity_time(4)
fixed_length_header(t1_minute) = lfric2um_config%dump_validity_time(5)
fixed_length_header(t1_second) = lfric2um_config%dump_validity_time(6)
! Validity time of data
fixed_length_header(t2_year) = lfric2um_config%dump_validity_time(1)
fixed_length_header(t2_month) = lfric2um_config%dump_validity_time(2)
fixed_length_header(t2_day) = lfric2um_config%dump_validity_time(3)
fixed_length_header(t2_hour) = lfric2um_config%dump_validity_time(4)
fixed_length_header(t2_minute) = lfric2um_config%dump_validity_time(5)
fixed_length_header(t2_second) = lfric2um_config%dump_validity_time(6)

! This is normally the time the file was produced, don't have access
! to that information right now (presumably can get it from system)
fixed_length_header(t3_year) = 0
fixed_length_header(t3_month) = 0
fixed_length_header(t3_day) = 0
fixed_length_header(t3_hour) = 0
fixed_length_header(t3_minute) = 0
fixed_length_header(t3_second) = 0

call shumlib(routinename//'::set_fixed_length_header', &
     um_output_file% set_fixed_length_header(fixed_length_header))

end subroutine lfric2um_set_fixed_length_header
!------------------------------------------------------------------

subroutine lfric2um_set_integer_constants(um_output_file, um_grid, lfric_mesh)

implicit none

type(shum_file_type), intent(in out) :: um_output_file
type(lfricinp_grid_type), intent(in) :: um_grid
type(mesh_type), pointer, intent(in) :: lfric_mesh
character(len=*), parameter :: routinename = 'lfric2um_set_integer_constants'
! Dumps and fieldsfiles have int consts length 46
integer(kind=int64) :: int_constants(46) = um_imdi

int_constants(ih_row_length) = um_grid%num_p_points_x
int_constants(ih_rows) = um_grid%num_p_points_y
int_constants(ih_model_levels) = lfric_mesh%get_nlayers()
int_constants(ih_wet_levels) = int_constants(ih_model_levels)
! hardcode for now, no info on soil in lfric dump
int_constants(ih_soilT_levels) = 4
int_constants(ih_soilQ_levels) = int_constants(ih_soilT_levels)
int_constants(ih_cloud_levels) = lfric_mesh%get_nlayers()
int_constants(ih_tracer_levels) = lfric_mesh%get_nlayers()
! hardcode for now, doesn't exist as a concept in lfric yet
int_constants(ih_boundary_levels) = 29
int_constants(ih_height_gen) = 2 ! Smooth height generation method
! hardcode for now - doesn't exist in lfric yet
int_constants(ih_1_c_rho_level) = 30
!  hardcode for now - aquaplanet
int_constants(ih_land_points) = 0
int_constants(ih_ozone_levels) = int_constants(ih_model_levels)
int_constants(ih_convect_levels) = 0

call shumlib(routinename//'::set_integer_constants', &
     um_output_file%set_integer_constants(int_constants))

end subroutine lfric2um_set_integer_constants

!------------------------------------------------------------------

subroutine lfric2um_set_real_constants(um_output_file, um_grid, lfric_mesh)

implicit none

type(shum_file_type), intent(in out) :: um_output_file
type(lfricinp_grid_type), intent(in) :: um_grid
type(mesh_type), pointer, intent(in) :: lfric_mesh
character(len=*), parameter :: routinename = 'lfric2um_set_real_constants'
! Dumps and fieldsfiles have real consts length 38
real(kind=real64) :: real_constants(38) = um_rmdi

real_constants(rh_deltaEW) = um_grid%spacing_x
real_constants(rh_deltaNS) = um_grid%spacing_y
real_constants(rh_baselat) = um_grid%grid_origin_y
real_constants(rh_baselong) = um_grid%grid_origin_x
real_constants(rh_polelat) = um_grid%pole_lat
real_constants(rh_polelong) = um_grid%pole_long
real_constants(rh_model_top) = lfric_mesh%get_domain_top()

call shumlib(routinename//'::set_real_constants', &
     um_output_file%set_real_constants(real_constants))

end subroutine lfric2um_set_real_constants

!------------------------------------------------------------------

subroutine lfric2um_set_level_dep_constants(um_output_file, lfric_mesh)

implicit none

type(shum_file_type), intent(in out) :: um_output_file
type(mesh_type), pointer, intent(in) :: lfric_mesh
character(len=*), parameter :: routinename = 'lfric2um_set_level_dep_constants'
real(kind=real64), allocatable :: level_dep_constants(:,:)
real(kind=real64) :: model_height
integer(kind=int64) :: num_levels
integer(kind=int64) :: k
integer(kind=int64) :: level_num_of_first_constant_rho

! Get num_levels
call shumlib(routinename //'::get_integer_constants_by_index',       &
      um_output_file%get_integer_constants_by_index(ih_model_levels, &
      num_levels))
! Get height at top of model
call shumlib(routinename //'::get_real_constants_by_index',          &
     um_output_file%get_real_constants_by_index(rh_model_top,        &
     model_height))
! Get level number of first constant rho level
call shumlib(routinename//'::get_integer_constants_by_index',        &
     um_output_file%get_integer_constants_by_index(ih_1_c_rho_level, &
     level_num_of_first_constant_rho))

allocate(level_dep_constants(num_levels + 1, 8))
level_dep_constants(:,:) = um_rmdi

! Set eta values for theta and rho levels
call lfric_mesh%get_eta(level_dep_constants(:, ldc_eta_theta))
do k = 1, num_levels
  ! Rho levels are midpoint between theta
  level_dep_constants(k, ldc_eta_rho) = &
      (level_dep_constants(k, ldc_eta_theta) +  &
       level_dep_constants(k+1, ldc_eta_theta))  / 2.0_real64
end do
! No soil level information yet - leave as rmdi
! No RHCrit information yet - leave as rmdi

! Calculate Zsea and C Values - formulae taken from UMDP F03 Appendix A
! model level definition.
! Zsea(k) = eta_value(k) * Height_at_top_of_model
! C(k) = (1 - ( eta_value(k) / eta_value(k_of_first_constant_rho_level)) ) ^2
do k = 1, num_levels + 1
  level_dep_constants(k, ldc_Zsea_theta) = &
       level_dep_constants(k, ldc_eta_theta) * model_height
  level_dep_constants(k, ldc_C_theta) = ( 1 - &
    (level_dep_constants(k, ldc_eta_theta)  / &
     level_dep_constants(level_num_of_first_constant_rho, ldc_eta_theta))) ** 2
  if (k > level_num_of_first_constant_rho) then
    level_dep_constants(k, ldc_C_theta) = 0.0_real64
  end if
end do
do k = 1, num_levels
  level_dep_constants(k, ldc_Zsea_rho) = &
       level_dep_constants(k, ldc_eta_rho) * model_height
  level_dep_constants(k, ldc_C_rho) = ( 1 - &
    (level_dep_constants(k, ldc_eta_rho)  / &
     level_dep_constants(level_num_of_first_constant_rho, ldc_eta_rho))) ** 2
  if (k >= level_num_of_first_constant_rho) then
    level_dep_constants(k, ldc_C_rho) = 0.0_real64
  end if
end do

call shumlib(routinename//'::set_level_dependent_constants', &
     um_output_file%set_level_dependent_constants(level_dep_constants))

end subroutine lfric2um_set_level_dep_constants

end module lfric2um_initialise_um_mod
