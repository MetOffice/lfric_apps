! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE um2lfric_regrid_fields_mod


IMPLICIT NONE

PRIVATE

PUBLIC :: um2lfric_regrid_fields

CONTAINS

SUBROUTINE um2lfric_regrid_fields(fctime)

! Description:
!
!  This routine performs a level by level regridding of the UM fields to
!  an intermediate regridded field array for each field, which is then
!  unpacked/copied into the corresponding lfric field data array as follows:
!
!  The LFRic field is a single 1D array containing all levels where the field
!  dof ordering in the LFRic array loops over each column before moving onto
!  the next horizontal point.
!
!  In this current UM2LFRic code we use a serial implementation of the LFRic
!  infrastructure. This means that there are no complications from the
!  MPI partitioning or from field stencils.
!
!  Simple example of a 2 layer mesh, with 4 points on each layer and both
!  the target and LFRic field indicies labelled for each point
!
!  Note that the target array indices repeat for each layer but the LFRic
!  field index refers to the full 3D field

!             Layer 1                                  Layer 2
!
!
! x  regridded_1      x  regridded_3       x  regridded_1       x  regridded_3
!    LFRIC_1             LFRIC_5              LFRIC_2              LFRIC_6
!
! x  regridded_2      x  regridded_4       x  regridded_2       x  regridded_4
!    LFRIC_3             LFRIC_7              LFRIC_4              LFRIC_8

! 2D array field containing the regridded data. First dimension corresponds
! to number of points per level and the second dimension to the number levels

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env,  ONLY: real64, int32, int64
USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_BOOL

! lfricinputs modules
USE lfricinp_check_shumlib_status_mod, ONLY: shumlib
USE lfricinp_regrid_weights_type_mod, ONLY: lfricinp_regrid_weights_type
USE lfricinp_stash_to_lfric_map_mod, ONLY: get_field_name
USE lfricinp_lfric_driver_mod, ONLY: lfric_fields
USE lfricinp_regrid_options_mod, ONLY: regrid_type

! um2lfric modules
USE um2lfric_namelist_mod,       ONLY: um2lfric_config
USE um2lfric_read_um_file_mod,  ONLY: um_input_file
USE um2lfric_regrid_weights_mod, ONLY: get_weights
USE um2lfric_post_process_fields_mod, ONLY: um2lfric_post_process_fields
USE um2lfric_apply_masked_field_adjustments_mod, ONLY:                         &
                                       um2lfric_apply_masked_field_adjustments

! shumlib modules
USE f_shum_field_mod, ONLY: shum_field_type

! lfric modules
USE field_mod,          ONLY: lfric_field_type => field_type,                  &
                              lfric_proxy_type => field_proxy_type
USE function_space_mod, ONLY: function_space_type
USE fs_continuity_mod,  ONLY: W3, Wtheta, W2H
USE log_mod,            ONLY: log_event, LOG_LEVEL_INFO, LOG_LEVEL_ERROR,      &
                              log_scratch_space

IMPLICIT NONE

! Forecast time of a field
REAL(KIND=real64), INTENT(IN) :: fctime

! Array of shumlib field objects that will be returned from UM file
TYPE(shum_field_type), ALLOCATABLE  :: um_input_fields(:)

! Intermediate target field for regridding
REAL(KIND=real64), ALLOCATABLE :: regridded_field(:,:)

! Pointers to lfric objects
TYPE(lfric_field_type), POINTER :: lfric_field => NULL()
TYPE(function_space_type), POINTER :: lfric_field_fs => NULL()
TYPE(lfricinp_regrid_weights_type), POINTER :: weights => NULL()

! LFRic field proxy
TYPE(lfric_proxy_type) :: lfric_field_proxy

! Other variables
INTEGER(KIND=int64) :: stashcode
INTEGER(KIND=int32) :: fs_type
! Iterators
INTEGER :: i_field, level, regridded_index, len_regridded_field, lfric_index
INTEGER :: errorstatus
! Number levels of UM field
INTEGER(KIND=int32) :: num_levels, num_dofs_per_level, num_lfric_levels
LOGICAL(KIND=C_BOOL) :: true_cbool

true_cbool = LOGICAL(.TRUE., KIND=C_BOOL)

!-------------------------------------------------------------------------------
! Initialise lfric fields to zero
!-------------------------------------------------------------------------------
DO i_field = 1, um2lfric_config%num_fields
  stashcode = um2lfric_config%stash_list(i_field)
  call lfric_fields % get_field(get_field_name(stashcode), lfric_field)
  lfric_field_proxy = lfric_field % get_proxy()
  lfric_field_proxy % data(:) = 0.0
  lfric_field => NULL()
END DO

!-------------------------------------------------------------------------------
! Main loop over fields for regridding
!-------------------------------------------------------------------------------
WRITE(log_scratch_space, '(A,I0,A)') 'Will process ', &
     um2lfric_config%num_fields, ' fields'
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)

DO i_field = 1, um2lfric_config%num_fields

  stashcode = um2lfric_config%stash_list(i_field)
  WRITE(log_scratch_space, '(A,I0)') 'Processing/regrid STASH code: ',         &
                                     stashcode
  CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
  WRITE(log_scratch_space,'(A,A)') 'LFRic field name: ',                       &
                                   TRIM(get_field_name(stashcode))
  CALL log_event(log_scratch_space, LOG_LEVEL_INFO)

  !-----------------------------------------------------------------------------
  ! Get UM field array and number of UM field levels
  !-----------------------------------------------------------------------------
  CALL shumlib("um2lfric::find_fields_in_file",                                &
                um_input_file%find_fields_in_file(um_input_fields,             &
                stashcode = stashcode, lbproc = 0_int64,                       &
                fctime = fctime), ignore_warning = true_cbool,                 &
                errorstatus = errorstatus)

  IF (errorstatus /= 0) THEN ! Field has not been found in dump

    WRITE(log_scratch_space,'(A,I0,A)') 'WARNING: stashcode ', stashcode,      &
                                        'not found in input dump. Data set '// &
                                        'to zero in LFRic output field.'
    CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
    IF (ALLOCATED(um_input_fields)) DEALLOCATE(um_input_fields)

  ELSE ! Field has sucessfully been found

    num_levels = SIZE(um_input_fields)

    !---------------------------------------------------------------------------
    ! Create pointers to lfric field and function space field lives on
    !---------------------------------------------------------------------------
    call lfric_fields % get_field(get_field_name(stashcode), lfric_field)
    lfric_field_fs => lfric_field % get_function_space()

    !---------------------------------------------------------------------------
    ! Allocate and set dimensions of intermediate regridded field data array
    !---------------------------------------------------------------------------
    fs_type = lfric_field % which_function_space()
    IF ((fs_type == W3) .OR.  (fs_type == W2H)) THEN
      num_lfric_levels = lfric_field_fs % get_nlayers()
    ELSE IF (fs_type == Wtheta) THEN
      num_lfric_levels = lfric_field_fs % get_nlayers() + 1
    ELSE
      WRITE(log_scratch_space,'(A)') 'Function space for field ' //            &
                                     get_field_name(stashcode) //              &
                                     ' not currently supported in regridding'
      CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
    END IF
    num_dofs_per_level = lfric_field_fs % get_undf() /                         &
                         (num_lfric_levels * lfric_field_fs % get_ndata())
    ALLOCATE(regridded_field(num_dofs_per_level, num_levels))

    !---------------------------------------------------------------------------
    ! Get required regridding weights
    !---------------------------------------------------------------------------
    weights => get_weights(stashcode)

    !---------------------------------------------------------------------------
    ! Loop for level by level regridding
    !---------------------------------------------------------------------------
    DO level = 1, num_levels

      !-------------------------------------------------------------------------
      ! Perform regridding/copying data from one grid/mesh to another
      !-------------------------------------------------------------------------
      CALL weights % regrid_src_2d_dst_1d (src=um_input_fields(level)%rdata,   &
                                           dst=regridded_field(:, level))
      !-------------------------------------------------------------------------
      ! Perform post regridding masked field adjustments, if regridding was not
      ! a simple copy of data (as is case for LAMs)
      !-------------------------------------------------------------------------
      IF (TRIM(regrid_type) /= 'lam_to_lam') THEN
        CALL um2lfric_apply_masked_field_adjustments(                          &
                                            stashcode,                         &
                                            src=um_input_fields(level)%rdata,  &
                                            dst=regridded_field(:, level))
      END IF

    END DO ! loop over levels

    ! Tidy up input field memory
    IF (ALLOCATED(um_input_fields)) DEALLOCATE(um_input_fields)

    !---------------------------------------------------------------------------
    ! Do any final post-processing to field, if required
    !---------------------------------------------------------------------------
    CALL um2lfric_post_process_fields(regridded_field, stashcode)
    ! Update number of levels for field if it changed during post-processing
    num_levels = SIZE(regridded_field, 2)

    !---------------------------------------------------------------------------
    ! Copy the regridded data to the lfric field
    !---------------------------------------------------------------------------
    lfric_field_proxy = lfric_field % get_proxy()
    IF (SIZE(regridded_field) /= SIZE(lfric_field_proxy % data)) THEN
      WRITE(log_scratch_space,'(A,I0,A,I0)')                                   &
                              'Mismatch between lfric field data array size ', &
                              SIZE(lfric_field_proxy % data), ' and ' //       &
                              'regridded array size ',                         &
                              SIZE(regridded_field)
      CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
    END IF
    WRITE(log_scratch_space,'(A)') 'Fill LFRic field object with regridded data'
    CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
    DO level = 1, num_levels
       len_regridded_field = SIZE(regridded_field, 1)
       ! Split up the data here and insert it into the LFRic field in the
       ! correct place. Loop over all points in the 2D regridded field that
       ! represents a single level. The first lfric array index will match the
       ! current level number
       lfric_index = level
       DO regridded_index = 1, len_regridded_field
         lfric_field_proxy % data(lfric_index) =                               &
                                      lfric_field_proxy % data(lfric_index) +  &
                                      regridded_field(regridded_index, level)
         ! Need to step by the total number of levels to get to the lfric
         ! array index that corresponds to the next regridded point
         lfric_index = lfric_index + num_levels
       END DO
     END DO

    DEALLOCATE(regridded_field)

    !---------------------------------------------------------------------------
    ! Nullify pointers before redefining
    !---------------------------------------------------------------------------
    weights => NULL()
    lfric_field => NULL()
    lfric_field_fs => NULL()

  END IF

END DO ! loop over stashcodes

END SUBROUTINE um2lfric_regrid_fields

END MODULE um2lfric_regrid_fields_mod
