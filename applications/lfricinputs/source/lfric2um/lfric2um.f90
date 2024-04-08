! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
PROGRAM lfric2um

! lfricinputs modules
USE lfricinp_create_lfric_fields_mod,     ONLY: lfricinp_create_lfric_fields
USE lfricinp_um_grid_mod,                 ONLY: um_grid
USE lfricinp_datetime_mod,                ONLY: datetime
USE lfricinp_initialise_mod,              ONLY: lfricinp_initialise
USE lfricinp_lfric_driver_mod,            ONLY: lfricinp_initialise_lfric,     &
                                                lfricinp_finalise_lfric, mesh, &
                                                twod_mesh,                     &
                                                lfric_fields

! lfric2um modules
USE lfric2um_namelists_mod,               ONLY: lfric2um_nl_fname,             &
                                                lfric2um_config,               &
                                                required_lfric_namelists
USE lfric2um_initialise_um_mod,           ONLY: lfric2um_initialise_um,        &
                                                um_output_file
USE lfric2um_initialise_lfric2um_mod,     ONLY: lfric2um_initialise_lfric2um
USE lfric2um_main_loop_mod,               ONLY: lfric2um_main_loop


IMPLICIT NONE

!==========================================================================
! Read inputs and initialise setup
!==========================================================================

! Read command line arguments and return details of filenames.
! Initialise common infrastructure
CALL lfricinp_initialise(lfric2um_nl_fname)

! Initialise lfric2um
CALL lfric2um_initialise_lfric2um()

! Initialise LFRic Infrastructure
CALL lfricinp_initialise_lfric(program_name_arg="lfric2um",                    &
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

! Create UM output file and initialise the headers
CALL lfric2um_initialise_um()

! Create LFRic field collection based on list of stashcodes
CALL lfricinp_create_lfric_fields( mesh, twod_mesh, lfric_fields,              &
                                   lfric2um_config%stash_list, um_grid,        &
                                   um_output_file )

!==========================================================================
! lfric2um main loop
!==========================================================================
! Main loop over fields to be read, regridded and written to output dump
CALL lfric2um_main_loop()

!==========================================================================
! Close files and finalise lfric infrastructure
!==========================================================================
! Finalise LFRic infrastructure
CALL lfricinp_finalise_lfric()

END PROGRAM lfric2um
