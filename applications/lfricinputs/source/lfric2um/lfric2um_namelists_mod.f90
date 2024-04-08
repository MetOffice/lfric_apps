! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!> @brief   Module containing lfric2um configuration TYPE
!> @details Hold information on input and output files and fields to be
!!          regridded, including proceedures to read this information from
!!          the relevant namelist.
MODULE lfric2um_namelists_mod

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : int64, real64

! lfricinputs modules
USE lfricinp_um_parameters_mod,     ONLY: fnamelen, um_imdi

IMPLICIT NONE
PRIVATE

PUBLIC :: lfric2um_config, required_lfric_namelists

TYPE :: config
  CHARACTER(LEN=fnamelen) :: output_filename = 'unset'
  CHARACTER(LEN=fnamelen) :: target_grid_namelist = 'unset'
  CHARACTER(LEN=fnamelen) :: stashmaster_file = 'unset'
  CHARACTER(LEN=fnamelen) :: weights_file_face_centre_to_p_bilinear = 'unset'
  CHARACTER(LEN=fnamelen) :: weights_file_face_centre_to_p_neareststod = 'unset'
  CHARACTER(LEN=fnamelen) :: weights_file_face_centre_to_u_bilinear = 'unset'
  CHARACTER(LEN=fnamelen) :: weights_file_face_centre_to_v_bilinear = 'unset'
  INTEGER(KIND=int64), ALLOCATABLE ::  stash_list(:)
  INTEGER(KIND=int64) :: um_version_int = um_imdi
  INTEGER(KIND=int64) :: dump_validity_time(6) = um_imdi
  INTEGER :: num_fields

  INTEGER :: status = -1
  CHARACTER(LEN=512) :: message = 'No namelist read'
  INTEGER :: unit_number
CONTAINS

  PROCEDURE :: load_namelists

END TYPE config

! Input namelist configuration
TYPE(config) :: lfric2um_config

! Namelist filenames read from command line
CHARACTER(LEN=fnamelen), PUBLIC :: lfric2um_nl_fname

INTEGER(KIND=int64), PARAMETER :: max_stash_list = 999

CHARACTER(*), PARAMETER  :: required_lfric_namelists(6) = ['logging         ', &
                                                           'finite_element  ', &
                                                           'base_mesh       ', &
                                                           'planet          ', &
                                                           'extrusion       ', &
                                                           'io              ']

CONTAINS

SUBROUTINE load_namelists(self)

! Descriptions:
!  Reads in lfric2um namelists. Performs checking on namelist values.
!  Populates UM grid object and prints out grid diagnostics.

! LFRic modules
USE log_mod,          ONLY: log_event, LOG_LEVEL_ERROR, LOG_LEVEL_INFO
USE constants_mod,    ONLY: imdi, rmdi

! lfricinputs modules
USE lfricinp_grid_namelist_mod, ONLY: grid, lambda_origin_targ,                &
     phi_origin_targ, phi_pole, lambda_pole, delta_lambda_targ,                &
     delta_phi_targ, points_lambda_targ, points_phi_targ,                      &
     igrid_targ, rotated
USE lfricinp_unit_handler_mod, ONLY: get_free_unit

! lfric2um modules
USE lfricinp_um_grid_mod, ONLY: um_grid

IMPLICIT NONE
CLASS(config) :: self


! Local variables
INTEGER :: i_stash

! Namelist variables
CHARACTER(LEN=fnamelen) :: output_filename = 'unset'
CHARACTER(LEN=fnamelen) :: stashmaster_file = 'unset'
CHARACTER(LEN=fnamelen) :: target_grid_namelist = 'unset'
CHARACTER(LEN=fnamelen) :: weights_file_face_centre_to_p_bilinear = 'unset'
CHARACTER(LEN=fnamelen) :: weights_file_face_centre_to_p_neareststod = 'unset'
CHARACTER(LEN=fnamelen) :: weights_file_face_centre_to_u_bilinear = 'unset'
CHARACTER(LEN=fnamelen) :: weights_file_face_centre_to_v_bilinear = 'unset'
INTEGER(KIND=int64) ::  stash_list(max_stash_list)
INTEGER(KIND=int64) :: um_version_int = um_imdi
INTEGER(KIND=int64) :: dump_validity_time(6) = um_imdi

NAMELIST /configure_lfric2um/ output_filename,                                 &
                              target_grid_namelist,                            &
                              stashmaster_file,                                &
                              weights_file_face_centre_to_p_bilinear,          &
                              weights_file_face_centre_to_p_neareststod,       &
                              weights_file_face_centre_to_u_bilinear,          &
                              weights_file_face_centre_to_v_bilinear,          &
                              stash_list,                                      &
                              um_version_int,                                  &
                              dump_validity_time

stash_list(:) = um_imdi

self%status = 0
self%message = 'Reading namelist from ' // TRIM(lfric2um_nl_fname)

CALL get_free_unit(self%unit_number)

OPEN(UNIT=self%unit_number, FILE=lfric2um_nl_fname, IOSTAT=self%status,                    &
                            IOMSG=self%message)
IF (self%status /= 0) CALL log_event(self%message, LOG_LEVEL_ERROR)

READ(self%unit_number, NML=configure_lfric2um, IOSTAT=self%status,             &
                       IOMSG=self%message)
IF (self%status /= 0) CALL log_event(self%message, LOG_LEVEL_ERROR)

IF (TRIM(output_filename) == 'unset') THEN
  self%status = 1
  self%message='Target filename is unset'
  CALL log_event(self%message, LOG_LEVEL_ERROR)
END IF

IF (TRIM(stashmaster_file) == 'unset') THEN
  self%status = 1
  self%message='Stashmaster filename is unset'
  CALL log_event(self%message, LOG_LEVEL_ERROR)
END IF

IF (TRIM(weights_file_face_centre_to_p_bilinear) == 'unset') THEN
  self%status = 1
  self%message='weights_file_face_centre_to_p_bilinear filename is unset'
  CALL log_event(self%message, LOG_LEVEL_ERROR)
END IF

IF (TRIM(weights_file_face_centre_to_p_neareststod) == 'unset') THEN
  self%status = 1
  self%message='weights_file_face_centre_to_p_neareststod filename is unset'
  CALL log_event(self%message, LOG_LEVEL_ERROR)
END IF

IF (TRIM(weights_file_face_centre_to_u_bilinear) == 'unset') THEN
  self%status = 1
  self%message='weights_file_face_centre_to_u_bilinear filename is unset'
  CALL log_event(self%message, LOG_LEVEL_ERROR)
END IF

IF (TRIM(weights_file_face_centre_to_v_bilinear) == 'unset') THEN
  self%status = 1
  self%message='weights_file_face_centre_to_v_bilinear filename is unset'
  CALL log_event(self%message, LOG_LEVEL_ERROR)
END IF

IF (TRIM(target_grid_namelist) == 'unset') THEN
  self%status = 1
  self%message='Target grid namelist is unset'
  CALL log_event(self%message, LOG_LEVEL_ERROR)
END IF

CLOSE(self%unit_number)

! Now read grid namelist to define UM grid
OPEN(UNIT=self%unit_number, FILE=target_grid_namelist, IOSTAT=self%status,     &
                            IOMSG=self%message)
IF (self%status /= 0) CALL log_event(self%message, LOG_LEVEL_ERROR)
READ(self%unit_number, NML=grid, IOSTAT=self%status,                           &
                       IOMSG=self%message)
IF (self%status /= 0) CALL log_event(self%message, LOG_LEVEL_ERROR)

! Load namelist variables into objects
self%output_filename = output_filename
self%target_grid_namelist = target_grid_namelist
self%stashmaster_file = stashmaster_file
self%weights_file_face_centre_to_p_bilinear =                                  &
                                 weights_file_face_centre_to_p_bilinear
self%weights_file_face_centre_to_p_neareststod =                               &
                                 weights_file_face_centre_to_p_neareststod
self%weights_file_face_centre_to_u_bilinear =                                  &
                                 weights_file_face_centre_to_u_bilinear
self%weights_file_face_centre_to_v_bilinear =                                  &
                                 weights_file_face_centre_to_v_bilinear
self%um_version_int = um_version_int
self%dump_validity_time(:) = dump_validity_time(:)

self%num_fields=0
! Count how many fields have been requested
DO i_stash = 1, max_stash_list
  IF (stash_list(i_stash) == um_imdi) THEN
    EXIT
  ELSE
    self%num_fields = self%num_fields + 1
  END IF
END DO

IF (self%num_fields <= 0) THEN
  CALL log_event('No fields selected in stash_list namelist variable',         &
       LOG_LEVEL_ERROR)
END IF

! Can now allocate type variable
ALLOCATE(self%stash_list(self%num_fields))
self%stash_list(:) = stash_list(1:self%num_fields)

CALL um_grid%set_grid_coords(                                                  &
     grid_staggering = igrid_targ ,                                            &
     num_p_points_x = points_lambda_targ,                                      &
     num_p_points_y = points_phi_targ,                                         &
     grid_spacing_x = delta_lambda_targ,                                       &
     grid_spacing_y = delta_phi_targ,                                          &
     ! namelist provides p grid values and routine
     ! wants base grid origin, so apply offset
     grid_origin_x = lambda_origin_targ - (0.5 * delta_lambda_targ),           &
     grid_origin_y = phi_origin_targ - (0.5 * delta_phi_targ))

CALL um_grid%print_grid_coords()

self%status = 0
self%message = 'Successfully read namelists from ' // TRIM(lfric2um_nl_fname)
CALL log_event(self%message, LOG_LEVEL_INFO)
CLOSE(self%unit_number)

END SUBROUTINE load_namelists

END MODULE lfric2um_namelists_mod
