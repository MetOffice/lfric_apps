!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Controls the initialisation and finalisation of multifile IO
!>        for the IAU
!>
module iau_multifile_io_mod

  use base_mesh_config_mod,        only: prime_mesh_name
  use calendar_mod,                only: calendar_type
  use constants_mod,               only: str_def, i_def
  use driver_modeldb_mod,          only: modeldb_type
  use event_mod,                   only: event_action
  use event_actor_mod,             only: event_actor_type
  use field_mod,                   only: field_type
  use sci_geometric_constants_mod, only: get_chi_inventory,     &
                                         get_panel_id_inventory
  use iau_multifile_file_setup_mod,only: init_iau_inc_files
#ifdef UM_PHYSICS
  use files_config_mod,            only: iau_addinf_path, &
                                         iau_bcorr_path
  use iau_config_mod,              only: iau_ainc_multifile, &
                                         iau_use_addinf,     &
                                         iau_use_bcorr
  use iau_firstfile_io_mod,        only: iau_incs_firstfile_io
#endif
  use iau_time_control_mod,        only: calc_iau_ts_num
  use inventory_by_mesh_mod,       only: inventory_by_mesh_type
  use io_context_mod,              only: callback_clock_arg
  use lfric_xios_context_mod,      only: lfric_xios_context_type
  use linked_list_mod,             only: linked_list_type
  use lfric_xios_action_mod,       only: advance_read_only
  use mesh_mod,                    only: mesh_type
  use mesh_collection_mod,         only: mesh_collection
  use model_clock_mod,             only: model_clock_type
  use namelist_mod,                only: namelist_type
  use step_calendar_mod,           only: step_calendar_type

  implicit none

  private

  public :: init_multifile_io
  public :: setup_step_multifile_io
  public :: step_multifile_io
  public :: finalise_multifile_io
  private :: init_iau_incs_io
  private :: context_init

contains

  !> @brief Initialise the multifile IO
  !>
  !> @param[in]    io_context_name name of main context
  !> @param[inout] modeldb Modeldb object
  subroutine init_multifile_io(io_context_name, modeldb)
    implicit none

    character(*),         intent(in)  :: io_context_name
    type(modeldb_type), intent(inout) :: modeldb

#ifdef UM_PHYSICS
    if ( iau_use_addinf ) then
      call iau_incs_firstfile_io ( io_context_name,     &
                                   modeldb,             &
                                   "iau_addinf_fields", &
                                   iau_addinf_path )
      call init_iau_incs_io( modeldb,             &
                             "iau_addinf_fields", &
                             "iau_addinf_io" )
    end if
    if ( iau_ainc_multifile ) then
      call init_iau_incs_io( modeldb,     &
                            "iau_fields", &
                            "iau_ainc_io"  )
    end if
    if ( iau_use_bcorr ) then
      call iau_incs_firstfile_io ( io_context_name,    &
                                   modeldb,            &
                                   "iau_bcorr_fields", &
                                   iau_bcorr_path )
      call init_iau_incs_io( modeldb,            &
                             "iau_bcorr_fields", &
                             "iau_bcorr_io" )
    end if
#endif

  end subroutine init_multifile_io

  !> @brief Initialise the IO for the different IAU increment types
  !>
  !> @param[inout] modeldb  Modeldb object
  !> @param[in]    iau_incs Type of IAU increment
  !> @param[in]    nml_name Name of multifile namelist
  subroutine init_iau_incs_io(modeldb, iau_incs, nml_name)
    implicit none

    type(modeldb_type), intent(inout), target :: modeldb
    character(*),       intent(in)            :: iau_incs
    character(*),       intent(in)            :: nml_name

    type(lfric_xios_context_type), pointer :: io_context
    type(namelist_type), pointer           :: multifile_nml
    type(linked_list_type), pointer        :: file_list
    class( model_clock_type ), pointer     :: model_clock
    integer(i_def)     :: iau_time
    integer(i_def)     :: multifile_start_timestep
    integer(i_def)     :: multifile_stop_timestep
    integer(i_def)     :: i
    character(str_def) :: context_name
    character(str_def) :: filename
    character(str_def), pointer :: multifile_io_profiles(:)

    allocate(multifile_io_profiles, &
             source=modeldb%configuration%get_namelist_profiles(nml_name))
    model_clock => modeldb % clock

    do i=1, size(multifile_io_profiles)

      multifile_nml => modeldb%configuration%get_namelist(trim(nml_name), &
                            profile_name=trim(multifile_io_profiles(i)))
      call multifile_nml%get_value('filename', filename)
      call multifile_nml%get_value('start_time', iau_time)
      multifile_start_timestep = calc_iau_ts_num (model_clock, iau_time)
      multifile_stop_timestep = multifile_start_timestep + 1_i_def

      context_name = "multifile_context_" // trim(filename)
      call context_init(modeldb, context_name, multifile_start_timestep, &
                        multifile_stop_timestep)

      call modeldb%io_contexts%get_io_context(context_name, io_context)

      file_list => io_context%get_filelist()
      call init_iau_inc_files(file_list, modeldb, iau_incs, filename)

    end do

    deallocate(multifile_io_profiles)

  end subroutine init_iau_incs_io

  !> @brief Wrapper to step the multifile IO
  !>
  !> @param[inout] modeldb            Model database object
  subroutine setup_step_multifile_io( io_context_name, modeldb )
    implicit none

    character(*),            intent(in)    :: io_context_name ! main context
    type(modeldb_type),      intent(inout) :: modeldb

#ifdef UM_PHYSICS
    if ( iau_use_addinf ) then
      call step_multifile_io( io_context_name, &
                              modeldb,         &
                              "iau_addinf_io" )
    end if

    if ( iau_ainc_multifile ) then
      call step_multifile_io( io_context_name, &
                              modeldb,         &
                              "iau_ainc_io" )
    end if

    if ( iau_use_bcorr ) then
      call step_multifile_io( io_context_name, &
                              modeldb,         &
                              "iau_bcorr_io" )
    end if
#endif

  end subroutine setup_step_multifile_io

  !> @brief Step the multifile IO
  !>
  !> @param[in]    io_context_name Name of main context
  !> @param[inout] modeldb         Model database object
  !> @param[in]    nml_name        Name of multifile namelist
  subroutine step_multifile_io( io_context_name, modeldb, nml_name )
    implicit none

    character(*),       intent(in)    :: io_context_name ! main context name
    character(*),       intent(in)    :: nml_name
    type(modeldb_type), intent(inout) :: modeldb

    type(inventory_by_mesh_type),  pointer :: chi_inventory
    type(inventory_by_mesh_type),  pointer :: panel_id_inventory
    type(lfric_xios_context_type), pointer :: io_context
    class(event_actor_type),       pointer :: event_actor_ptr
    type(mesh_type),               pointer :: mesh
    type(field_type),              pointer :: chi(:)
    type(field_type),              pointer :: panel_id
    type(namelist_type),           pointer :: multifile_nml
    type(namelist_type),           pointer :: time_nml
    class(calendar_type), allocatable      :: tmp_calendar
    character(str_def) :: context_name
    character(str_def) :: filename
    character(str_def) :: time_origin
    character(str_def) :: time_start
    character(str_def), pointer :: multifile_io_profiles(:)
    integer(i_def)     :: i
    procedure(event_action), pointer       :: context_advance
    procedure(callback_clock_arg), pointer :: before_close

    nullify(before_close)

    chi_inventory => get_chi_inventory()
    panel_id_inventory => get_panel_id_inventory()

    allocate(multifile_io_profiles, &
             source=modeldb%configuration%get_namelist_profiles(nml_name))

    do i=1, size(multifile_io_profiles)

      multifile_nml => modeldb%configuration%get_namelist(nml_name, &
                    profile_name=trim(multifile_io_profiles(i)))
      call multifile_nml%get_value('filename', filename)

      context_name = "multifile_context_" // trim(filename)
      call modeldb%io_contexts%get_io_context(context_name, io_context)

      if (modeldb%clock%get_step() == io_context%get_stop_time()) then
        ! Finalise XIOS context
        call io_context%set_current()
        call io_context%set_active(.false.)
        call modeldb%clock%remove_event(context_name)
        call io_context%finalise_xios_context()

      elseif (modeldb%clock%get_step() == io_context%get_start_time()) then
        ! Initialise XIOS context
        mesh => mesh_collection%get_mesh(prime_mesh_name)
        call chi_inventory%get_field_array(mesh, chi)
        call panel_id_inventory%get_field(mesh, panel_id)

        time_nml => modeldb%configuration%get_namelist('time')

        call time_nml%get_value('calendar_origin', time_origin)
        call time_nml%get_value('calendar_start', time_start)

        allocate(tmp_calendar, source=step_calendar_type(time_origin, time_start))

        call io_context%initialise_xios_context( modeldb%mpi%get_comm(),      &
                                                 chi, panel_id,               &
                                                 modeldb%clock, tmp_calendar, &
                                                 before_close,                &
                                                 start_at_zero=.true. )

        ! Attach context advancement to the model's clock
        context_advance => advance_read_only
        event_actor_ptr => io_context
        call modeldb%clock%add_event( context_advance, event_actor_ptr )
        call io_context%set_active(.true.)
      end if
    end do

    call modeldb%io_contexts%get_io_context(io_context_name, io_context)
    call io_context%set_current()

    deallocate(multifile_io_profiles)

    nullify ( chi_inventory, panel_id_inventory, mesh, chi, panel_id )

  end subroutine step_multifile_io

  !> @brief Finalise the multifile IO
  !>
  !> @param[inout] modeldb Model database object
  subroutine finalise_multifile_io(modeldb)
    implicit none

    type(modeldb_type), intent(inout) :: modeldb

  end subroutine finalise_multifile_io

  !> @brief Initialise the IO context
  !>
  !> @param[inout] modeldb Model            database object
  !> @param[in]    context_name             name of context
  !> @param[in]    multifile_start_timestep timestep to start context
  !> @param[in]    multifile_stop_timestep  timestep to stop context
  subroutine context_init(modeldb, &
                          context_name, &
                          multifile_start_timestep, &
                          multifile_stop_timestep)

    implicit none

    type(modeldb_type), intent(inout) :: modeldb
    character(*), intent(in)          :: context_name
    integer(i_def), intent(in)        :: multifile_start_timestep
    integer(i_def), intent(in)        :: multifile_stop_timestep

    type(lfric_xios_context_type)     :: tmp_io_context

    call tmp_io_context%initialise( context_name,                   &
                                    start=multifile_start_timestep, &
                                    stop=multifile_stop_timestep )
    call modeldb%io_contexts%add_context(tmp_io_context)

  end subroutine context_init

end module iau_multifile_io_mod
