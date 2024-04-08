! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfric2um_main_loop_mod

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY: int64, real64

! lfric2um modules
USE lfric2um_namelists_mod, ONLY: lfric2um_config
USE lfricinp_lfric_driver_mod, ONLY: lfric_fields, local_rank, comm, &
     twod_mesh
USE lfric2um_initialise_um_mod, ONLY: um_output_file
USE lfricinp_um_grid_mod, ONLY: um_grid
USE lfric2um_regrid_weights_mod, ONLY: get_weights

! lfricinp modules
USE lfricinp_stashmaster_mod, ONLY: get_stashmaster_item, levelt, &
     rho_levels, theta_levels, single_level
USE lfricinp_add_um_field_to_file_mod, ONLY: lfricinp_add_um_field_to_file
USE lfricinp_um_level_codes_mod, ONLY: lfricinp_get_num_levels
USE lfricinp_check_shumlib_status_mod, ONLY: shumlib
USE lfricinp_regrid_weights_type_mod, ONLY: lfricinp_regrid_weights_type
USE lfricinp_stash_to_lfric_map_mod, ONLY: get_field_name
USE lfricinp_gather_lfric_field_mod, ONLY: lfricinp_gather_lfric_field

! lfric modules
USE field_mod, ONLY: field_type
USE log_mod, ONLY: log_event, log_scratch_space, LOG_LEVEL_INFO

IMPLICIT NONE

PRIVATE
PUBLIC :: lfric2um_main_loop

CONTAINS

SUBROUTINE lfric2um_main_loop()
! Description:
!  Main processing field loop in lfric2um. Loops over requested fields,
!  reads field on all ranks, loops over levels, gather full level data
!  onto rank 0, perform regridding and write to disk.
IMPLICIT NONE

INTEGER(KIND=int64) :: i_stash, level, i_field
INTEGER(KIND=int64) :: stashcode, num_levels
CHARACTER(LEN=*), PARAMETER :: routinename='lfric2um_main_loop'
TYPE(field_type), POINTER :: lfric_field
TYPE(lfricinp_regrid_weights_type), POINTER :: weights
REAL(KIND=real64), ALLOCATABLE :: global_field_array(:)

! Main loop over requested stashcodes
DO i_stash = 1, lfric2um_config%num_fields
  stashcode = lfric2um_config%stash_list(i_stash)
  WRITE(log_scratch_space, '(A,I0)') "Processing stashcode ", stashcode
  CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
  num_levels = lfricinp_get_num_levels(um_output_file, stashcode)

  ! Select appropriate weights
  weights => get_weights(stashcode)

  ! Get pointer to lfric field + read
  call lfric_fields%get_field(get_field_name(stashcode), lfric_field)
  CALL lfric_field%read_field("read_"//lfric_field%get_name())

  ! Allocate space for global data, only need full field on rank 0
  IF (ALLOCATED(global_field_array)) DEALLOCATE(global_field_array)
  IF (local_rank == 0) THEN
    ALLOCATE(global_field_array(weights%num_points_src))
  ELSE
    ALLOCATE(global_field_array(1))
  END IF

  ! Loop over number of levels in field
  DO level = 1, num_levels
    global_field_array(:) = 0.0
    ! Gather local lfric fields into global array on base rank
    CALL lfricinp_gather_lfric_field( lfric_field, global_field_array, comm, &
         num_levels, level, twod_mesh )
    IF (local_rank == 0 ) THEN
      ! Create UM field in output file
      CALL lfricinp_add_um_field_to_file(um_output_file, stashcode, &
           level, um_grid)
      ! Adding a 2D UM field to the file increments num_fields by 1 each time
      ! Get the index of the last field added, which is just num_fields
      i_field = um_output_file%num_fields
      CALL weights%validate_dst(SIZE(um_output_file%fields(i_field)%rdata))
      ! Perform the regridding from lfric to um
      CALL weights%regrid_src_1d_dst_2d(global_field_array(:), &
           um_output_file%fields(i_field)%rdata(:,:))

      ! Write to file
      CALL shumlib(routinename//'::write_field',  &
           um_output_file%write_field(i_field))
      ! Unload data from field
      CALL shumlib(routinename//'::unload_field', &
           um_output_file%unload_field(i_field))
    END IF
  END DO ! end loop over levels
  ! Unload field data from memory
  CALL lfric_field%field_final()
END DO  ! end loop over field stashcodes

IF (local_rank == 0) THEN ! Only write UM file with rank 0
  ! Write output header
  CALL shumlib(routinename//'::write_header', &
           um_output_file%write_header())
  ! Close file
  CALL shumlib(routinename//'::close_file', &
           um_output_file%close_file())
END IF

END SUBROUTINE lfric2um_main_loop

END MODULE lfric2um_main_loop_mod
