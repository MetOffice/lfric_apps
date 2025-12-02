! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfricinp_read_command_line_args_mod

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only: int32

implicit none

private

public :: lfricinp_read_command_line_args

contains

  !> @brief   Read the command line arguments provided to the program
  !>
  !> @param[out] lfricinputs_fname String holding namelist for the Lfric Inputs
  !!                               application being run
  !> @param[out] lfric_fname       String holding namelists for the generic lfric
  !!                               LFRic infrastructure
  !> @param[out] io_fname          String holding namelists for the input and
  !!                               output files being used including ancilaries
  subroutine lfricinp_read_command_line_args(lfricinputs_fname,                &
                                           lfric_fname,                        &
                                           io_fname)

    use lfricinp_um_parameters_mod, only: fnamelen
    implicit none

    character(len=fnamelen), intent(out) :: lfricinputs_fname
    character(len=fnamelen), intent(out) :: lfric_fname
    character(len=fnamelen), intent(out) :: io_fname

    ! Other variables
    integer :: arglen
    integer(kind=int32) :: icode_32

    ! Read LFRic Inputs namelist filename from command line
    call get_command_argument(1, lfricinputs_fname, arglen, icode_32)
    call check_command_line_errors(icode_32)

    ! Read LFRic infrastructure namelist filename from command line
    call get_command_argument(2, lfric_fname, arglen, icode_32)
    call check_command_line_errors(icode_32)

    ! Read IO namelist filename from command line
    call get_command_argument(3, io_fname, arglen, icode_32)
    call check_command_line_errors(icode_32)

    return
  end subroutine lfricinp_read_command_line_args

  ! Checks the output from "get_command_argument()" and writes appropriate error
  ! messages. Either returns nothing or finishes with an error.
  subroutine check_command_line_errors(icode_32)

    use log_mod,         only: log_event, LOG_LEVEL_INFO, LOG_LEVEL_ERROR,     &
                               log_scratch_space
    implicit none

    integer(kind=int32), intent(in out):: icode_32

    select case(icode_32)
      case (0)
        continue
      case (1)
        log_scratch_space = 'No filename provided on command line'
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
      case (-1)
        log_scratch_space =                                                    &
             'The filename and path to the namelist file is too long '//       &
             'for currently compiled string declaration. Please '     //       &
             'recompile with a larger filenamelength parameter or '   //       &
             'reconsider location of the provided namelist.'
        icode_32 = 1 ! Force fatal error
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
      case DEFAULT
        log_scratch_space =                                                    &
             'Unknown error reading namelist file from command line.'
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end select

    return
  end subroutine check_command_line_errors

end module lfricinp_read_command_line_args_mod
