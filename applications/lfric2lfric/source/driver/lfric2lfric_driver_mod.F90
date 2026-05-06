!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Drives the execution of the lfric2lfric miniapp.
!>
module lfric2lfric_driver_mod

  use constants_mod,            only: str_def, i_def, l_def, r_second
  use driver_fem_mod,           only: final_fem
  use driver_io_mod,            only: final_io
  use driver_modeldb_mod,       only: modeldb_type
  use field_collection_mod,     only: field_collection_type
  use lfric_xios_action_mod,    only: advance
  use lfric_xios_context_mod,   only: lfric_xios_context_type
  use lfric_xios_read_mod,      only: read_checkpoint, &
                                      read_state
  use lfric_xios_write_mod,     only: write_checkpoint, &
                                      write_state
  use log_mod,                  only: log_scratch_space, &
                                      log_event,         &
                                      log_level_info
  use model_clock_mod,          only: model_clock_type
  use namelist_mod,             only: namelist_type
  use sci_checksum_alg_mod,     only: checksum_alg
  use xios,                     only: xios_date, xios_get_current_date, &
                                      xios_date_convert_to_string

  !------------------------------------
  ! lfric2lfric modules
  !------------------------------------
  use lfric2lfric_config_mod,         only: mode_ics, mode_lbc
  use lfric2lfric_infrastructure_mod, only: initialise_infrastructure, &
                                            context_dst, context_src,  &
                                            source_collection_name,    &
                                            target_collection_name
  use lfric2lfric_regrid_mod,         only: lfric2lfric_regrid

  implicit none

  private
  public initialise, run, finalise

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief          Sets up required state in preparation for run.
  !> @details        Calls the `initialise_infrastructure` subroutine that
  !!                 checks the configuration namelist, initialises meshes,
  !!                 extrusions, XIOS contexts and files, field collections
  !!                 and fields.
  !> @param [in,out] modeldb                 The structure holding model state
  !> @param [in,out] oasis_clock             Clock for OASIS exchanges
  subroutine initialise( modeldb, oasis_clock )

    implicit none

    type(modeldb_type),     intent(inout) :: modeldb
    type(model_clock_type), allocatable, &
                            intent(inout) :: oasis_clock

    call initialise_infrastructure( modeldb, oasis_clock )

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief    Performs regridding of a field collection
  !> @details  Populates fields in the source field collection with data
  !!           read from the source XIOS context file, regrids the field
  !!           collection to a destination mesh, switches to the destination
  !!           XIOS context and writes to an output file.
  !> @param [in,out] modeldb      The structure that holds model state
  !> @param [in,out] oasis_clock  Clock for OASIS exchanges
  subroutine run( modeldb, oasis_clock )

    implicit none

    type(modeldb_type),     intent(inout) :: modeldb
    type(model_clock_type), allocatable, &
                            intent(inout) :: oasis_clock

    ! LFRic-XIOS constants
    integer(kind=i_def), parameter :: start_timestep = 1_i_def

    ! Namelist variables
    character(len=str_def) :: start_dump_filename
    character(len=str_def) :: checkpoint_stem_name
    integer(kind=i_def)    :: mode
    integer(kind=i_def)    :: regrid_method

    ! Local parameters
    type(namelist_type), pointer :: files_nml
    type(namelist_type), pointer :: lfric2lfric_nml

    integer(kind=i_def)          :: step, time_steps
    logical(kind=l_def)          :: is_running
    type(xios_date)              :: current_date
    character(len=32)            :: current_date_str
    real(kind=r_second)          :: checkpoint_times(1)

    type(field_collection_type),   pointer :: source_fields
    type(field_collection_type),   pointer :: target_fields

    type(lfric_xios_context_type), pointer :: io_context

    ! Namelist pointers
    files_nml       => modeldb%configuration%get_namelist('files')
    lfric2lfric_nml => modeldb%configuration%get_namelist('lfric2lfric')

    ! Extract configuration variables
    call files_nml%get_value( 'start_dump_filename', start_dump_filename )
    call files_nml%get_value( 'checkpoint_stem_name', checkpoint_stem_name )
    call lfric2lfric_nml%get_value( 'mode', mode )
    call lfric2lfric_nml%get_value( 'regrid_method', regrid_method )

    ! Point to source and target field collections
    source_fields => modeldb%fields%get_field_collection(source_collection_name)
    target_fields => modeldb%fields%get_field_collection(target_collection_name)

    ! Read fields and perform the regridding
    if (mode == mode_ics) then
      call read_checkpoint(source_fields,      &
                           start_timestep,     &
                           start_dump_filename )

      call lfric2lfric_regrid(modeldb, oasis_clock, source_fields,   &
                              target_fields, regrid_method)

      ! Write output
      call modeldb%io_contexts%get_io_context(context_dst, io_context)
      call io_context%set_current()

      checkpoint_times(1) = modeldb%clock%seconds_from_steps(modeldb%clock%get_step())
      call write_checkpoint(target_fields, modeldb%values, modeldb%clock, &
                            checkpoint_stem_name, checkpoint_times)

    else if (mode == mode_lbc) then
      time_steps = modeldb%clock%get_last_step() - &
                   modeldb%clock%get_first_step() + 1

      do step=1, time_steps
        call xios_get_current_date(current_date)
        call xios_date_convert_to_string(current_date, current_date_str)
        write(log_scratch_space, '(A)') 'Regridding at xios date: ' // &
                                                   current_date_str
        call log_event( log_scratch_space, log_level_info )

        call read_state(source_fields)

        call lfric2lfric_regrid(modeldb, oasis_clock, source_fields, &
                                target_fields, regrid_method)

        is_running = modeldb%clock%tick()

        call modeldb%io_contexts%get_io_context(context_dst, io_context)
        call io_context%set_current()
        call advance(io_context, modeldb%clock)

        call write_state(target_fields, prefix='lbc_')

        call modeldb%io_contexts%get_io_context(context_src, io_context)
        call io_context%set_current()
        call advance(io_context, modeldb%clock)
      end do
    end if

    ! Write checksum
    call checksum_alg("lfric2lfric", field_collection=target_fields)

  end subroutine run

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief   Tidies up after a run.
  !>
  !> @param [in]     program_name An identifier given to the model being run
  !> @param [in,out] modeldb      The structure that holds model state
  subroutine finalise( program_name, modeldb )

    implicit none

    character(len=*),   intent(in)     :: program_name
    type(modeldb_type), intent(inout)  :: modeldb


    !--------------------------------------------------------------------------
    ! Model finalise
    !--------------------------------------------------------------------------

    call log_event( program_name//': Miniapp completed', log_level_info )

    !-------------------------------------------------------------------------
    ! Driver layer finalise
    !-------------------------------------------------------------------------

    ! Finalise IO
    call final_io( modeldb )

    call final_fem()

  end subroutine finalise

end module lfric2lfric_driver_mod
