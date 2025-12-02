! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
program um2lfric

! lfricinputs modules
use lfricinp_lfric_driver_mod,     only: lfricinp_initialise_lfric,     &
                                         lfricinp_finalise_lfric, mesh, &
                                         twod_mesh, lfric_fields
use lfricinp_ancils_mod,           only: lfricinp_create_ancil_fields, &
                                         ancil_fields
use lfricinp_create_lfric_fields_mod, &
                                   only: lfricinp_create_lfric_fields
use lfricinp_um_grid_mod,          only: um_grid
use lfricinp_datetime_mod,         only: datetime
use lfricinp_initialise_mod,       only: lfricinp_initialise

! um2lfric modules
use um2lfric_namelist_mod,         only: um2lfric_nl_fname, &
                                         um2lfric_config,   &
                                         required_lfric_namelists
use um2lfric_initialise_um2lfric_mod, &
                                   only: um2lfric_initialise_um2lfric
use um2lfric_regrid_weights_mod,   only: um2lfric_regrid_weightsfile_ctl
use um2lfric_init_masked_field_adjustments_mod, &
                                   only: um2lfric_init_masked_field_adjustments
use um2lfric_regrid_and_output_data_mod, &
                                   only: um2lfric_regrid_and_output_data
use um2lfric_check_input_data_mod, only: um2lfric_check_input_data
use um2lfric_read_um_file_mod,     only: um2lfric_close_um_file, &
                                         um_input_file

implicit none

!==========================================================================
! Read inputs and initialise setup
!==========================================================================

! Read command line arguments and return details of filenames.
! Initialise common infrastructure
call lfricinp_initialise(um2lfric_nl_fname)

! Initialise um2lfric
call um2lfric_initialise_um2lfric()

! Initialise LFRic Infrastructure
call lfricinp_initialise_lfric(program_name_arg="um2lfric",                    &
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
call um2lfric_check_input_data(um_input_file)

! Initialise LFRic field collection
call lfricinp_create_lfric_fields( mesh, twod_mesh, lfric_fields,              &
                                   um2lfric_config%stash_list,                 &
                                   um_grid, um_input_file )

! Initialise LFRic ancils field collection
call lfricinp_create_ancil_fields( ancil_fields, mesh, twod_mesh )

! Read in, process and partition regridding weights
call um2lfric_regrid_weightsfile_ctl()

!==========================================================================
! um2lfric main loop
!==========================================================================
! Now initialise masked points that requires post regridding adjustments
call um2lfric_init_masked_field_adjustments()

! Perform regridding and output data
call um2lfric_regrid_and_output_data(datetime)

!==========================================================================
! Close files and finalise lfric infrastructure
!==========================================================================
! Unloads data from memory and closes UM input file
call um2lfric_close_um_file()

! Finalise YAXT, XIOS, MPI, logging
call lfricinp_finalise_lfric()

end program um2lfric
