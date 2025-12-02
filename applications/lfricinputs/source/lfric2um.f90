! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
program lfric2um

! lfricinputs modules
use lfricinp_create_lfric_fields_mod,     only: lfricinp_create_lfric_fields
use lfricinp_um_grid_mod,                 only: um_grid
use lfricinp_datetime_mod,                only: datetime
use lfricinp_initialise_mod,              only: lfricinp_initialise
use lfricinp_lfric_driver_mod,            only: lfricinp_initialise_lfric,     &
                                                lfricinp_finalise_lfric, mesh, &
                                                twod_mesh,                     &
                                                lfric_fields

! lfric2um modules
use lfric2um_namelists_mod,               only: lfric2um_nl_fname,             &
                                                lfric2um_config,               &
                                                required_lfric_namelists
use lfric2um_initialise_um_mod,           only: lfric2um_initialise_um,        &
                                                um_output_file
use lfric2um_initialise_lfric2um_mod,     only: lfric2um_initialise_lfric2um
use lfric2um_main_loop_mod,               only: lfric2um_main_loop


implicit none

!==========================================================================
! Read inputs and initialise setup
!==========================================================================

! Read command line arguments and return details of filenames.
! Initialise common infrastructure
call lfricinp_initialise(lfric2um_nl_fname)

! Initialise lfric2um
call lfric2um_initialise_lfric2um()

! Initialise LFRic Infrastructure
call lfricinp_initialise_lfric(program_name_arg="lfric2um",                    &
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
call lfric2um_initialise_um()

! Create LFRic field collection based on list of stashcodes
call lfricinp_create_lfric_fields( mesh, twod_mesh, lfric_fields,              &
                                   lfric2um_config%stash_list, um_grid,        &
                                   um_output_file )

!==========================================================================
! lfric2um main loop
!==========================================================================
! Main loop over fields to be read, regridded and written to output dump
call lfric2um_main_loop()

!==========================================================================
! Close files and finalise lfric infrastructure
!==========================================================================
! Finalise LFRic infrastructure
call lfricinp_finalise_lfric()

end program lfric2um
