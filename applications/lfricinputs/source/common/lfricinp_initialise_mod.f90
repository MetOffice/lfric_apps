! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!> @brief A module providing common initialisation routines.
!>
!> @details Provides a common basis for initialising um2lfric and lfric2um.

MODULE lfricinp_initialise_mod

USE lfricinp_um_parameters_mod, ONLY: fnamelen
! LFRic modules
USE log_mod,         ONLY: log_event, LOG_LEVEL_INFO

IMPLICIT NONE

PRIVATE

PUBLIC :: lfricinp_initialise

CONTAINS

!> Initialises generic lfricinputs infrastructure
!> @param [out] program_fname Filename for program specific namelists provided
!>                            on the command line
SUBROUTINE lfricinp_initialise(program_fname)
  USE lfricinp_read_command_line_args_mod, ONLY: lfricinp_read_command_line_args
  USE lfricinp_setup_io_mod,          ONLY: io_config, io_fname
  USE lfricinp_regrid_options_mod, ONLY: lfricinp_init_regrid_options
  USE lfricinp_datetime_mod, ONLY: datetime
  USE lfricinp_stash_to_lfric_map_mod, ONLY: lfricinp_init_stash_to_lfric_map
  USE lfricinp_lfric_driver_mod, ONLY: lfric_nl_fname

  IMPLICIT NONE

  CHARACTER(LEN=fnamelen), INTENT(OUT) :: program_fname

  CALL log_event('Reading command line', LOG_LEVEL_INFO)
  CALL lfricinp_read_command_line_args(program_fname, lfric_nl_fname, io_fname)

  CALL log_event('Loading IO namelist', LOG_LEVEL_INFO)
  CALL io_config%load_namelist()

  CALL log_event('Loading global regridding options', LOG_LEVEL_INFO)
  CALL lfricinp_init_regrid_options(program_fname)

  CALL log_event('Initialise stashcode to lfric field mapping', LOG_LEVEL_INFO)
  CALL lfricinp_init_stash_to_lfric_map()

  CALL log_event('Initialise datetime class', LOG_LEVEL_INFO)
  CALL datetime % initialise()


END SUBROUTINE lfricinp_initialise

END MODULE lfricinp_initialise_mod
