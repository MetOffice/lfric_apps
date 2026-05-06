!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Provides a method to setup the IO and associated XIOS context for JEDI-LFRIC
!>
module jedi_lfric_io_setup_mod

  use calendar_mod,              only: calendar_type
  use constants_mod,             only: i_def
  use driver_fem_mod,            only: init_fem, final_fem
  use empty_io_context_mod,      only: empty_io_context_type
  use event_actor_mod,           only: event_actor_type
  use event_mod,                 only: event_action
  use field_mod,                 only: field_type
  use inventory_by_mesh_mod,     only: inventory_by_mesh_type
  use io_config_mod,             only: subroutine_timers
  use lfric_mpi_mod,             only: lfric_mpi_type, &
                                       lfric_comm_type
  use log_mod,                   only: log_event, log_level_error
  use linked_list_mod,           only: linked_list_type
  use mesh_mod,                  only: mesh_type
  use mesh_collection_mod,       only: mesh_collection
  use model_clock_mod,           only: model_clock_type
  use namelist_collection_mod,   only: namelist_collection_type
  use namelist_mod,              only: namelist_type

  use jedi_lfric_file_meta_mod,  only: jedi_lfric_file_meta_type
  use jedi_lfric_init_files_mod, only: jedi_lfric_init_files

#ifdef USE_XIOS
  use io_context_mod,           only: io_context_type, callback_clock_arg
  use lfric_xios_context_mod,   only: lfric_xios_context_type
  use lfric_xios_action_mod,    only: advance
#endif

  implicit none

  private

  public initialise_io

contains

  !> @brief Initialise the XIOS context and IO
  !>
  !> @param [in]    context_name  The name of the context
  !> @param [in]    mpi           The mpi communicator
  !> @param [in]    file_meta     The file meta data
  !> @param [in]    mesh_name     The name of the mesh
  !> @param [in]    calendar      The model calendar
  !> @param [inout] io_context  The LFRic context object
  !> @param [inout] model_clock   The model clock
  subroutine initialise_io( context_name, mpi, file_meta, mesh_name, &
                            calendar, io_context, model_clock )

    implicit none

    character(len=*),                       intent(in) :: context_name
    class(lfric_mpi_type),                  intent(in) :: mpi
    type(jedi_lfric_file_meta_type),        intent(in) :: file_meta(:)
    character(len=*),                       intent(in) :: mesh_name
    class(calendar_type),                   intent(in) :: calendar
    class(io_context_type), allocatable, intent(inout) :: io_context
    type(model_clock_type),              intent(inout) :: model_clock

    ! Local
    type(inventory_by_mesh_type) :: chi_inventory
    type(inventory_by_mesh_type) :: panel_id_inventory
    type(mesh_type), pointer     :: mesh
    type(field_type), pointer    :: chi(:)
    type(field_type), pointer    :: panel_id
    type(lfric_comm_type)        :: lfric_comm

    nullify(mesh, chi, panel_id)

    ! Create FEM specifics (function spaces and chi field)
    call init_fem( mesh_collection, chi_inventory, panel_id_inventory )

    ! Get coordinate fields for prime mesh
    mesh => mesh_collection%get_mesh( mesh_name )
    call chi_inventory%get_field_array( mesh, chi )
    call panel_id_inventory%get_field( mesh, panel_id )

    ! Initialise I/O context and setup file to use
    lfric_comm = mpi%get_comm()
    call init_io( context_name, lfric_comm%get_comm_mpi_val(), file_meta, &
                  calendar, io_context, chi, panel_id, model_clock )

    ! Do initial step
    if ( model_clock%is_initialisation() ) then
      select type (io_context)
      type is (lfric_xios_context_type)
        call advance(io_context, model_clock)
      end select
    end if

    call final_fem()

  end subroutine initialise_io

  ! The following method (init_io) is based on driver_io_mod init_io method.
  ! The main difference is the way that the file_list is created. Only the code
  ! required for JEDI-LFRIC is included.

  !> @brief  Initialises the model I/O and context
  !>
  !> @param[in] context_name  A string identifier for the context
  !> @param[in] communicator  The ID for the model MPI communicator
  !> @param[in] file_meta     The file meta data
  !> @param[in] calendar      The model calendar
  !> @param[in] io_context    The LFRic context object
  !> @param[in] chi           The model's coordinate fields
  !> @param[in] panel_id      The model's panel ID fields
  !> @param[in] model_clock   The model clock
  !> @param[in] before_close  Optional routine to be called before
  !>                          context closes
  subroutine init_io( context_name,  &
                      communicator,  &
                      file_meta,     &
                      calendar,      &
                      io_context,    &
                      chi,           &
                      panel_id,      &
                      model_clock,   &
                      before_close )

    implicit none

    character(*),                           intent(in) :: context_name
    integer(i_def),                         intent(in) :: communicator
    type(jedi_lfric_file_meta_type),        intent(in) :: file_meta(:)
    class(calendar_type),                   intent(in) :: calendar
    class(io_context_type), target, allocatable, intent(inout) :: io_context
    type(field_type),                       intent(in) :: chi(:)
    type(field_type),                       intent(in) :: panel_id
    type(model_clock_type),              intent(inout) :: model_clock
    procedure(callback_clock_arg), optional            :: before_close


    ! Local
    procedure(callback_clock_arg), pointer :: before_close_ptr => null()
    integer(i_def) :: rc
    type(linked_list_type), pointer :: file_list
    class(event_actor_type), pointer :: event_actor_ptr
    procedure(event_action), pointer :: context_advance
    type(lfric_comm_type)            :: lfric_comm

    ! Allocate XIOS IO context types
    if (present(before_close)) then
      before_close_ptr => before_close
    end if

    allocate( lfric_xios_context_type::io_context, stat=rc )
    if (rc /= 0) then
      call log_event( "Unable to allocate LFRic-XIOS context object", &
                      log_level_error )
    end if

    ! Select the lfric_xios_context_type
    select type(io_context)
    type is (lfric_xios_context_type)

      ! Populate list of I/O files if procedure passed through
      file_list => io_context%get_filelist()
      call jedi_lfric_init_files(file_list, file_meta)

      ! Setup the context
      call io_context%initialise( context_name )
      call lfric_comm%set_comm_mpi_val(communicator)
      call io_context%initialise_xios_context( lfric_comm,            &
                                               chi, panel_id,         &
                                               model_clock, calendar, &
                                               before_close_ptr )
      ! Attach context advancement to the model's clock
      context_advance => advance
      event_actor_ptr => io_context
      call model_clock%add_event( context_advance, event_actor_ptr )
    end select

  end subroutine init_io

end module jedi_lfric_io_setup_mod
