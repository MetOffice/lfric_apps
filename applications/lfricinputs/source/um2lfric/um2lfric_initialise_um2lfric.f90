! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module um2lfric_initialise_um2lfric_mod

implicit none

private

public :: um2lfric_initialise_um2lfric

contains

!> @brief Initialises um2lfric specific infrastructure
!> This includes reading namelists, and stashmaster file,
!> initialising dates and times and reading in gridding weights
subroutine um2lfric_initialise_um2lfric()

! um2lfric modules
use um2lfric_namelist_mod,              only: um2lfric_config
use um2lfric_read_um_file_mod,          only: um2lfric_read_um_file, um_input_file

! lfricinputs modules
use lfricinp_stashmaster_mod,           only: lfricinp_read_stashmaster
use lfricinp_read_um_time_data_mod,     only: lfricinp_read_um_time_data
use lfricinp_um_grid_mod,               only: lfricinp_set_grid_from_file
implicit none

! Read um2lfric configuration namelist
call um2lfric_config%load_namelist()

! Open the UM file
!call log_event('Initialising UM input file', LOG_LEVEL_INFO)
call um2lfric_read_um_file(um2lfric_config%um_file)

! Load date and time information for requested stash items from um input file
call lfricinp_read_um_time_data(um_input_file,                                 &
                                um2lfric_config%stash_list)

! Read in UM stashmaster
call lfricinp_read_stashmaster(um2lfric_config%stashmaster_file)

! Initialise the grid
call lfricinp_set_grid_from_file(um_input_file,                     &
                                 um2lfric_config%num_snow_layers,   &
                                 um2lfric_config%num_surface_types, &
                                 um2lfric_config%num_ice_cats)

end subroutine um2lfric_initialise_um2lfric

end module um2lfric_initialise_um2lfric_mod
