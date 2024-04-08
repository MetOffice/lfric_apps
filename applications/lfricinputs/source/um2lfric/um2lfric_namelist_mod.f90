! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!> @brief   Module containing um2lfric configuration TYPE
!> @details Hold information on input and output files and fields to be
!!          regridded, including proceedures to read this information from
!!          the relevant namelist.
MODULE um2lfric_namelist_mod

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : int64, real64

! UM2LFRic modules
USE lfricinp_um_parameters_mod,     ONLY: fnamelen, um_imdi

IMPLICIT NONE
PRIVATE
PUBLIC :: um2lfric_config, required_lfric_namelists

TYPE :: config
  CHARACTER(LEN=fnamelen) :: um_file = 'unset'
  CHARACTER(LEN=fnamelen) :: stashmaster_file = 'unset'
  CHARACTER(LEN=fnamelen) :: weights_file_p_to_face_centre_bilinear = 'unset'
  CHARACTER(LEN=fnamelen) :: weights_file_p_to_face_centre_neareststod = 'unset'
  CHARACTER(LEN=fnamelen) :: weights_file_u_to_face_centre_bilinear = 'unset'
  CHARACTER(LEN=fnamelen) :: weights_file_v_to_face_centre_bilinear = 'unset'
  INTEGER(KIND=int64), ALLOCATABLE ::  stash_list(:)
  INTEGER(KIND=int64) :: num_snow_layers = um_imdi
  INTEGER(KIND=int64) :: num_surface_types = um_imdi

  INTEGER :: num_fields

  INTEGER :: status = -1
  CHARACTER(LEN=512) :: message = 'No namelist read'
  INTEGER :: unit_number

CONTAINS

  PROCEDURE :: load_namelist

END TYPE config

INTEGER(KIND=int64), PARAMETER :: max_stash_list = 999

! Input namelist configuration
TYPE(config) :: um2lfric_config

! Namelist filenames read from command line
CHARACTER(LEN=fnamelen), PUBLIC :: um2lfric_nl_fname

CHARACTER(*), PARAMETER :: required_lfric_namelists(6) =  &
    ['logging             ', &
     'finite_element      ', &
     'base_mesh           ', &
     'planet              ', &
     'extrusion           ', &
     'io                  ']

CONTAINS

SUBROUTINE load_namelist(self)

  ! lfricinp modules
  USE lfricinp_unit_handler_mod, ONLY: get_free_unit
  USE lfricinp_ancils_mod, ONLY: lfricinp_l_land_area_fraction => &
                                 l_land_area_fraction

  ! LFRic modules
  USE log_mod,          ONLY: log_event, LOG_LEVEL_ERROR

  IMPLICIT NONE
  CLASS(config) :: self

  ! Local variables
  INTEGER :: i_stash

  ! Namelist variables
  CHARACTER(LEN=fnamelen) :: um_file = 'unset'
  CHARACTER(LEN=fnamelen) :: stashmaster_file = 'unset'
  CHARACTER(LEN=fnamelen) :: weights_file_p_to_face_centre_bilinear = 'unset'
  CHARACTER(LEN=fnamelen) :: weights_file_p_to_face_centre_neareststod = 'unset'
  CHARACTER(LEN=fnamelen) :: weights_file_u_to_face_centre_bilinear = 'unset'
  CHARACTER(LEN=fnamelen) :: weights_file_v_to_face_centre_bilinear = 'unset'
  INTEGER(KIND=int64)     :: stash_list(max_stash_list)
  INTEGER(KIND=int64)     :: num_snow_layers = um_imdi
  INTEGER(KIND=int64)     :: num_surface_types = um_imdi
  LOGICAL :: l_land_area_fraction = .FALSE.

  NAMELIST /configure_um2lfric/ um_file,                                       &
                                stashmaster_file,                              &
                                weights_file_p_to_face_centre_bilinear,        &
                                weights_file_p_to_face_centre_neareststod,     &
                                weights_file_u_to_face_centre_bilinear,        &
                                weights_file_v_to_face_centre_bilinear,        &
                                stash_list, num_snow_layers, num_surface_types,&
                                l_land_area_fraction

  stash_list(:) = um_imdi

  self%status = 0
  self%message = 'Reading namelist from ' // TRIM(um2lfric_nl_fname)

  CALL get_free_unit(self%unit_number)

  OPEN(UNIT=self%unit_number, FILE=um2lfric_nl_fname, IOSTAT=self%status,                  &
                              IOMSG=self%message)
  IF (self%status /= 0) CALL log_event(self%message, LOG_LEVEL_ERROR)

  READ(self%unit_number, NML=configure_um2lfric, IOSTAT=self%status,           &
                         IOMSG=self%message)
  IF (self%status /= 0) CALL log_event(self%message, LOG_LEVEL_ERROR)

  IF (TRIM(um_file) == 'unset') THEN
    self%status = 1
    self%message='UM filename is unset'
    CALL log_event(self%message, LOG_LEVEL_ERROR)
  END IF

 ! Further error checking goes here

  ! Load namelist variables into object
  self%um_file = um_file
  self%stashmaster_file = stashmaster_file
  self%weights_file_p_to_face_centre_bilinear =                                &
                                    weights_file_p_to_face_centre_bilinear
  self%weights_file_p_to_face_centre_neareststod =                             &
                                    weights_file_p_to_face_centre_neareststod
  self%weights_file_u_to_face_centre_bilinear =                                &
                                    weights_file_u_to_face_centre_bilinear
  self%weights_file_v_to_face_centre_bilinear =                                &
                                    weights_file_v_to_face_centre_bilinear
  self%num_snow_layers = num_snow_layers
  self%num_surface_types = num_surface_types
  ! Pass the um2lfric namelist variable to the lfricinputs module variable
  lfricinp_l_land_area_fraction = l_land_area_fraction

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
    CALL log_event('No fields selected in stash_list namelist variable', &
         LOG_LEVEL_ERROR)
  END IF

  ! Can now allocate type variable
  ALLOCATE(self%stash_list(self%num_fields))
  self%stash_list(:) = stash_list(1:self%num_fields)

  self%status = 0
  self%message = 'Successfully read namelist from ' // TRIM(um2lfric_nl_fname)

  CLOSE(self%unit_number)

! Broadcasting?

END SUBROUTINE load_namelist

END MODULE um2lfric_namelist_mod
