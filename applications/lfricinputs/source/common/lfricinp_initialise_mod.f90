! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!> @brief A module providing common initialisation routines.
!>
!> @details Provides a common basis for initialising um2lfric and lfric2um.

module lfricinp_initialise_mod

use lfricinp_um_parameters_mod, only: fnamelen
! LFRic modules
use log_mod,         only: log_event, LOG_LEVEL_INFO

implicit none

private

public :: lfricinp_initialise

contains

!> Initialises generic lfricinputs infrastructure
!> @param [out] program_fname Filename for program specific namelists provided
!>                            on the command line
subroutine lfricinp_initialise(program_fname)
  use lfricinp_read_command_line_args_mod, only: lfricinp_read_command_line_args
  use lfricinp_setup_io_mod,          only: io_config, io_fname
  use lfricinp_regrid_options_mod, only: lfricinp_init_regrid_options
  use lfricinp_datetime_mod, only: datetime
  use lfricinp_stash_to_lfric_map_mod, only: lfricinp_init_stash_to_lfric_map
  use lfricinp_lfric_driver_mod, only: lfric_nl_fname

  implicit none

  character(len=fnamelen), intent(out) :: program_fname

  call log_event('Reading command line', LOG_LEVEL_INFO)
  call lfricinp_read_command_line_args(program_fname, lfric_nl_fname, io_fname)

  call log_event('Loading IO namelist', LOG_LEVEL_INFO)
  call io_config%load_namelist()

  call log_event('Loading global regridding options', LOG_LEVEL_INFO)
  call lfricinp_init_regrid_options(program_fname)

  call log_event('Initialise stashcode to lfric field mapping', LOG_LEVEL_INFO)
  call lfricinp_init_stash_to_lfric_map()

  call log_event('Initialise datetime class', LOG_LEVEL_INFO)
  call datetime % initialise()


end subroutine lfricinp_initialise

end module lfricinp_initialise_mod
