! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE um2lfric_check_input_data_mod

IMPLICIT NONE

PRIVATE
PUBLIC :: um2lfric_check_input_data
CONTAINS

SUBROUTINE um2lfric_check_input_data(um_input_file)
! Description:  Perform checks on the input data to ensure
!               that it can be suitably handled by UM2LFRic
!

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY: int64, real64

! LFRic modules
USE log_mod,       ONLY: log_event, LOG_LEVEL_ERROR, log_scratch_space
USE extrusion_config_mod, ONLY: number_of_layers, domain_top
! Shumlib modules
USE f_shum_file_mod,  ONLY: shum_file_type
USE f_shum_fixed_length_header_indices_mod, ONLY:                              &
    horiz_grid_type, grid_staggering, dataset_type
! lfricinp modules
USE lfricinp_check_shumlib_status_mod, ONLY: shumlib

IMPLICIT NONE

TYPE(shum_file_type), INTENT(INOUT) :: um_input_file

! Local variables
! Parameters for accessing UM header information
! Dataset type indicator values - fixed header
INTEGER(KIND=int64), PARAMETER :: inst_dump = 1
INTEGER(KIND=int64), PARAMETER :: fieldsfile = 3
! Horizontal grid indicator values- fixed header
INTEGER(KIND=int64), PARAMETER :: global_grid = 0, lam_no_wrap = 3,            &
                                  rotated_grid = 100
! Grid staggering indicator values - fixed header
INTEGER(KIND=int64), PARAMETER :: arakawa_C_endgame = 6
INTEGER(KIND=int64), PARAMETER :: arakawa_C_nd = 3

! Values needed from dump headers
INTEGER(KIND=int64) ::                                                         &
               um_file_type, dump_stagger, dump_grid_type, dump_num_levels
REAL(KIND=real64)   :: dump_model_top

! Indices of required values in integer constants
INTEGER(KIND=int64), PARAMETER :: ih_num_levels = 8
INTEGER(KIND=int64), PARAMETER :: rh_height_model_top = 16

! Error reporting
CHARACTER(LEN=*), PARAMETER :: routinename='um2lfric_check_input_data'

! Tolerance used for floating point comparisons
REAL(KIND=real64) :: tolerance

tolerance  = TINY(1.0_real64)

! Check that the input file is a UM dump or a Fieldsfile - we can't run with
! any other data format
CALL shumlib(routinename//'::get_fixed_length_header_by_index',                &
     um_input_file % get_fixed_length_header_by_index(                         &
     dataset_type, um_file_type))
IF ((um_file_type /= inst_dump) .AND. (um_file_type /= fieldsfile)) THEN
  WRITE(log_scratch_space, "(3(A,I0))" )                                       &
       "Input file is not a UM dump or fieldsfile, dataset type "              &
        // "found was: ", um_file_type, " but input file should have: ",       &
        inst_dump, " or ", fieldsfile
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

! Get number of levels from integer constants
CALL shumlib(routinename//'::get_integer_constants_by_index',                  &
     um_input_file % get_integer_constants_by_index(                           &
     ih_num_levels, dump_num_levels))
! Check that the number of layers in LFRic namelist matches the number
! of levels in the UM input file. Interpolation is not supported.
IF (dump_num_levels /= number_of_layers) THEN
   WRITE(log_scratch_space, '(2(A,I0))')                                       &
        "Mismatch between number of levels in UM "                             &
        // " dump: ", dump_num_levels, " and number"                           &
        // " of layers in LFRic mesh: ", number_of_layers
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

! Get model top from real constants
CALL shumlib(routinename//'::get_real_constants_by_index',                     &
     um_input_file % get_real_constants_by_index(                              &
     rh_height_model_top, dump_model_top))
! Check that height of the top of the model is the same.
IF (ABS(dump_model_top - domain_top) > tolerance ) THEN
   WRITE(log_scratch_space, '(2(A,F0.12))')                                    &
        "Mismatch between top of model in UM "                                 &
        // " dump: ", dump_model_top,                                          &
        " and top of model in LFRic namelist: ", domain_top
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

! Get horizontal grid type indicator
CALL shumlib(routinename//'::get_fixed_length_header_by_index',                &
     um_input_file % get_fixed_length_header_by_index(                         &
     horiz_grid_type, dump_grid_type))
! Check that we have a global grid or LAM with no wrapping (rotated or not)
IF (mod(dump_grid_type, rotated_grid) /= global_grid .AND.                     &
    mod(dump_grid_type, rotated_grid) /= lam_no_wrap                           &
   ) THEN
   WRITE(log_scratch_space, '(A,I0)')                                          &
     "Unsupported horiz grid type. Fixed header(4) = ", dump_grid_type
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

! Get grid staggering indicator
CALL shumlib(routinename//'::get_fixed_length_header_by_index',                &
     um_input_file % get_fixed_length_header_by_index(                         &
     grid_staggering, dump_stagger))
IF (dump_stagger /= arakawa_C_endgame) THEN
  IF (dump_stagger == arakawa_C_nd) THEN
    CALL log_event("New Dynamics Arakawa C grid not supported",                &
         LOG_LEVEL_ERROR)
  ELSE
    WRITE(log_scratch_space, '(A,I0)')                                         &
      "Unsupported grid stagger. Fixed header(9) = ", dump_stagger
    CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
  END IF
END IF

END SUBROUTINE um2lfric_check_input_data

END MODULE um2lfric_check_input_data_mod
