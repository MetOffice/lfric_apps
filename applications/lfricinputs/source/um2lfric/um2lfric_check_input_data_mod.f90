! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module um2lfric_check_input_data_mod

implicit none

private
public :: um2lfric_check_input_data
contains

subroutine um2lfric_check_input_data(um_input_file)
! Description:  Perform checks on the input data to ensure
!               that it can be suitably handled by UM2LFRic
!

! Intrinsic modules
use, intrinsic :: iso_fortran_env,          only: int64, real64

! LFRic modules
use log_mod,                                only: log_event,       &
                                                  LOG_LEVEL_ERROR, &
                                                  log_scratch_space
use extrusion_config_mod,                   only: number_of_layers, &
                                                  domain_height
! Shumlib modules
use f_shum_file_mod,                        only: shum_file_type
use f_shum_fixed_length_header_indices_mod, only: horiz_grid_type, &
                                                  grid_staggering, &
                                                  dataset_type
! lfricinp modules
use lfricinp_check_shumlib_status_mod,      only: shumlib

implicit none

type(shum_file_type), intent(INOUT) :: um_input_file

! Local variables
! Parameters for accessing UM header information
! Dataset type indicator values - fixed header
integer(kind=int64), parameter :: inst_dump = 1
integer(kind=int64), parameter :: fieldsfile = 3
! Horizontal grid indicator values- fixed header
integer(kind=int64), parameter :: global_grid = 0, lam_no_wrap = 3,            &
                                  rotated_grid = 100
! Grid staggering indicator values - fixed header
integer(kind=int64), parameter :: arakawa_C_endgame = 6
integer(kind=int64), parameter :: arakawa_C_nd = 3

! Values needed from dump headers
integer(kind=int64) ::                                                         &
               um_file_type, dump_stagger, dump_grid_type, dump_num_levels
real(kind=real64)   :: dump_model_top

! Indices of required values in integer constants
integer(kind=int64), parameter :: ih_num_levels = 8
integer(kind=int64), parameter :: rh_height_model_top = 16

! Error reporting
character(len=*), parameter :: routinename='um2lfric_check_input_data'

! Tolerance used for floating point comparisons
real(kind=real64) :: tolerance

tolerance  = tiny(1.0_real64)

! Check that the input file is a UM dump or a Fieldsfile - we can't run with
! any other data format
call shumlib(routinename//'::get_fixed_length_header_by_index',                &
     um_input_file % get_fixed_length_header_by_index(                         &
     dataset_type, um_file_type))
if ((um_file_type /= inst_dump) .and. (um_file_type /= fieldsfile)) then
  write(log_scratch_space, "(3(A,I0))" )                                       &
       "Input file is not a UM dump or fieldsfile, dataset type "              &
        // "found was: ", um_file_type, " but input file should have: ",       &
        inst_dump, " or ", fieldsfile
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end if

! Get number of levels from integer constants
call shumlib(routinename//'::get_integer_constants_by_index',                  &
     um_input_file % get_integer_constants_by_index(                           &
     ih_num_levels, dump_num_levels))
! Check that the number of layers in LFRic namelist matches the number
! of levels in the UM input file. Interpolation is not supported.
if (dump_num_levels /= number_of_layers) then
   write(log_scratch_space, '(2(A,I0))')                                       &
        "Mismatch between number of levels in UM "                             &
        // " dump: ", dump_num_levels, " and number"                           &
        // " of layers in LFRic mesh: ", number_of_layers
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end if

! Get model top from real constants
call shumlib(routinename//'::get_real_constants_by_index',                     &
     um_input_file % get_real_constants_by_index(                              &
     rh_height_model_top, dump_model_top))
! Check that height of the top of the model is the same.
if (abs(dump_model_top - domain_height) > tolerance ) then
   write(log_scratch_space, '(2(A,F0.12))')                                    &
        "Mismatch between top of model in UM "                                 &
        // " dump: ", dump_model_top,                                          &
        " and top of model in LFRic namelist: ", domain_height
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end if

! Get horizontal grid type indicator
call shumlib(routinename//'::get_fixed_length_header_by_index',                &
     um_input_file % get_fixed_length_header_by_index(                         &
     horiz_grid_type, dump_grid_type))
! Check that we have a global grid or LAM with no wrapping (rotated or not)
if (mod(dump_grid_type, rotated_grid) /= global_grid .and.                     &
    mod(dump_grid_type, rotated_grid) /= lam_no_wrap                           &
   ) then
   write(log_scratch_space, '(A,I0)')                                          &
     "Unsupported horiz grid type. Fixed header(4) = ", dump_grid_type
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end if

! Get grid staggering indicator
call shumlib(routinename//'::get_fixed_length_header_by_index',                &
     um_input_file % get_fixed_length_header_by_index(                         &
     grid_staggering, dump_stagger))
if (dump_stagger /= arakawa_C_endgame) then
  if (dump_stagger == arakawa_C_nd) then
    call log_event("New Dynamics Arakawa C grid not supported",                &
         LOG_LEVEL_ERROR)
  else
    write(log_scratch_space, '(A,I0)')                                         &
      "Unsupported grid stagger. Fixed header(9) = ", dump_stagger
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if
end if

end subroutine um2lfric_check_input_data

end module um2lfric_check_input_data_mod
