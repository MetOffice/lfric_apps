! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
PROGRAM um2lfric

! lfricinputs modules
USE lfricinp_lfric_driver_mod,          ONLY: lfricinp_initialise_lfric,       &
                                              lfricinp_finalise_lfric, mesh,   &
                                              twod_mesh, lfric_fields
USE lfricinp_ancils_mod,                ONLY: lfricinp_create_ancil_fields,    &
                                              ancil_fields
USE lfricinp_create_lfric_fields_mod,   ONLY: lfricinp_create_lfric_fields
USE lfricinp_um_grid_mod,               ONLY: um_grid
USE lfricinp_datetime_mod,              ONLY: datetime
USE lfricinp_initialise_mod,            ONLY: lfricinp_initialise

! um2lfric modules
USE um2lfric_namelist_mod,              ONLY: um2lfric_nl_fname,               &
                                              um2lfric_config,                 &
                                              required_lfric_namelists
USE um2lfric_initialise_um2lfric_mod,   ONLY: um2lfric_initialise_um2lfric
USE um2lfric_regrid_weights_mod,        ONLY: um2lfric_regrid_weightsfile_ctl
USE um2lfric_init_masked_field_adjustments_mod,  ONLY:                         &
                                       um2lfric_init_masked_field_adjustments
USE um2lfric_regrid_and_output_data_mod,ONLY: um2lfric_regrid_and_output_data
USE um2lfric_check_input_data_mod,      ONLY: um2lfric_check_input_data
USE um2lfric_read_um_file_mod,          ONLY: um2lfric_close_um_file,          &
                                              um_input_file

IMPLICIT NONE

!==========================================================================
! Read inputs and initialise setup
!==========================================================================

! Read command line arguments and return details of filenames.
! Initialise common infrastructure
CALL lfricinp_initialise(um2lfric_nl_fname)

! Initialise um2lfric
CALL um2lfric_initialise_um2lfric()

! Initialise LFRic Infrastructure
CALL lfricinp_initialise_lfric(program_name_arg="um2lfric",                    &
     required_lfric_namelists = required_lfric_namelists,                      &
     start_date = datetime % first_validity_time,                              &
     time_origin = datetime % first_validity_time,                             &
     first_step = datetime % first_step,                                       &
     last_step = datetime % last_step,                                         &
     spinup_period = datetime % spinup_period,                                 &
     seconds_per_step = datetime % seconds_per_step)

!==========================================================================
! Further input and output file setup
!==========================================================================

! Check input file
CALL um2lfric_check_input_data(um_input_file)

! Initialise LFRic field collection
CALL lfricinp_create_lfric_fields( mesh, twod_mesh, lfric_fields,              &
                                   um2lfric_config%stash_list,                 &
                                   um_grid, um_input_file )

! Initialise LFRic ancils field collection
CALL lfricinp_create_ancil_fields( ancil_fields, mesh, twod_mesh )

! Read in, process and partition regridding weights
CALL um2lfric_regrid_weightsfile_ctl()

!==========================================================================
! um2lfric main loop
!==========================================================================
! Now initialise masked points that requires post regridding adjustments
CALL um2lfric_init_masked_field_adjustments()

! Perform regridding and output data
CALL um2lfric_regrid_and_output_data(datetime)

!==========================================================================
! Close files and finalise lfric infrastructure
!==========================================================================
! Unloads data from memory and closes UM input file
CALL um2lfric_close_um_file()

! Finalise YAXT, XIOS, MPI, logging
CALL lfricinp_finalise_lfric()

END PROGRAM um2lfric
