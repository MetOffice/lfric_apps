! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE um2lfric_initialise_um2lfric_mod

IMPLICIT NONE

PRIVATE

PUBLIC :: um2lfric_initialise_um2lfric

CONTAINS

!> @brief Initialises um2lfric specific infrastructure
!> This includes reading namelists, and stashmaster file,
!> initialising dates and times and reading in gridding weights
SUBROUTINE um2lfric_initialise_um2lfric()

! um2lfric modules
USE um2lfric_namelist_mod,              ONLY: um2lfric_config
USE um2lfric_read_um_file_mod,          ONLY: um2lfric_read_um_file, um_input_file

! lfricinputs modules
USE lfricinp_stashmaster_mod,           ONLY: lfricinp_read_stashmaster
USE lfricinp_read_um_time_data_mod,     ONLY: lfricinp_read_um_time_data
USE lfricinp_um_grid_mod,               ONLY: lfricinp_set_grid_from_file
IMPLICIT NONE

! Read um2lfric configuration namelist
CALL um2lfric_config%load_namelist()

! Open the UM file
!CALL log_event('Initialising UM input file', LOG_LEVEL_INFO)
CALL um2lfric_read_um_file(um2lfric_config%um_file)

! Load date and time information for requested stash items from um input file
CALL lfricinp_read_um_time_data(um_input_file,                                 &
                                um2lfric_config%stash_list)

! Read in UM stashmaster
CALL lfricinp_read_stashmaster(um2lfric_config%stashmaster_file)

! Initialise the grid
CALL lfricinp_set_grid_from_file(um_input_file,                     &
                                 um2lfric_config%num_snow_layers,   &
                                 um2lfric_config%num_surface_types)

END SUBROUTINE um2lfric_initialise_um2lfric

END MODULE um2lfric_initialise_um2lfric_mod
