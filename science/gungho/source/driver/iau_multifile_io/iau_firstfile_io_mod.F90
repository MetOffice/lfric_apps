!-------------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Handles opening and closing of contexts to read the first files of
!>        IAU multifile increments during model initialisation
!>
module iau_firstfile_io_mod

  use calendar_mod,              only: calendar_type
  use constants_mod,             only: str_def, l_def
  use driver_modeldb_mod,        only: modeldb_type
  use field_collection_mod,      only: field_collection_type
  use field_mod,                 only: field_type
  use file_mod,                  only: FILE_MODE_READ
  use inventory_by_mesh_mod,     only: inventory_by_mesh_type
  use io_context_mod,            only: callback_clock_arg
  use lfric_xios_context_mod,    only: lfric_xios_context_type
  use lfric_xios_file_mod,       only: lfric_xios_file_type, &
                                       OPERATION_ONCE
  use linked_list_mod,           only: linked_list_type
  use mesh_collection_mod,       only: mesh_collection
  use mesh_mod,                  only: mesh_type
  use namelist_mod,              only: namelist_type
  use sci_geometric_constants_mod, only: get_chi_inventory, &
                                         get_panel_id_inventory
  use step_calendar_mod,         only: step_calendar_type

  implicit none

  public :: iau_incs_firstfile_io

contains
  !> @brief Opens and closes context to read the first set of
  !>        IAU multifile increments
  !> @param [in]    io_context_name Name of main context
  !> @param [inout] modeldb         The structure that holds model state
  !> @param [in]    iau_incs        Name of IAU increment field collection
  !> @param [in]    iau_incs_path   Path for IAU increment file
  !>
  subroutine iau_incs_firstfile_io ( io_context_name, modeldb, &
                                     iau_incs, iau_incs_path )

    implicit none

    character(*),            intent(in)    :: io_context_name ! main context
    character(*),            intent(in)    :: iau_incs
    character(*),            intent(in)    :: iau_incs_path
    type(modeldb_type),      intent(inout) :: modeldb

    type(mesh_type),               pointer :: mesh
    type(field_type),              pointer :: chi(:)
    type(field_type),              pointer :: panel_id
    type(inventory_by_mesh_type),  pointer :: chi_inventory
    type(inventory_by_mesh_type),  pointer :: panel_id_inventory
    type(lfric_xios_context_type)          :: tmp_io_context
    type(lfric_xios_context_type), pointer :: io_context
    type(linked_list_type),        pointer :: file_list
    type(field_collection_type),   pointer :: multifile_fields
    type(namelist_type),           pointer :: time_nml
    type(namelist_type),           pointer :: base_mesh_nml
    type(namelist_type),           pointer :: files_nml
    type(namelist_type),           pointer :: io_nml

    class(calendar_type), allocatable :: tmp_calendar

    character(str_def) :: time_origin
    character(str_def) :: time_start
    character(str_def) :: prime_mesh_name
    character(str_def) :: context_name
    character(str_def) :: iau_addinf_path
    character(str_def) :: iau_bcorr_path

    logical(l_def) :: use_xios_io

    procedure(callback_clock_arg), pointer :: before_close

    nullify(before_close)

    chi_inventory => get_chi_inventory()
    panel_id_inventory => get_panel_id_inventory()

    time_nml      => modeldb%configuration%get_namelist('time')
    base_mesh_nml => modeldb%configuration%get_namelist('base_mesh')
    files_nml     => modeldb%configuration%get_namelist('files')
    io_nml        => modeldb%configuration%get_namelist('io')

    call time_nml%get_value('calendar_origin', time_origin)
    call time_nml%get_value('calendar_start', time_start)
    call base_mesh_nml%get_value('prime_mesh_name', prime_mesh_name)
    call files_nml%get_value('iau_addinf_path', iau_addinf_path)
    call files_nml%get_value('iau_bcorr_path', iau_bcorr_path)
    call io_nml%get_value('use_xios_io', use_xios_io)

    ! get filename and set up context name for this file
    context_name = "multifile_context_" // trim(iau_incs_path)
    call tmp_io_context%initialise( context_name )

    !add context to modeldb
    call modeldb%io_contexts%add_context(tmp_io_context)

    !get context
    call modeldb%io_contexts%get_io_context(context_name, io_context)

    !set up file list
    file_list => io_context%get_filelist()
    multifile_fields  => modeldb%fields%get_field_collection(iau_incs)

    if ( use_xios_io) then

      call file_list%insert_item( &
           lfric_xios_file_type( trim(iau_incs_path), &
           xios_id = iau_incs, &
           io_mode=FILE_MODE_READ, &
           freq=1, &
           operation=OPERATION_ONCE, &
           fields_in_file=multifile_fields))

    end if ! use_xios_io

    ! Initialise XIOS context
    mesh => mesh_collection%get_mesh(prime_mesh_name)
    call chi_inventory%get_field_array(mesh, chi)
    call panel_id_inventory%get_field(mesh, panel_id)

    allocate(tmp_calendar, source=step_calendar_type(time_origin, time_start))

    call io_context%initialise_xios_context( modeldb%mpi%get_comm(),      &
                                             chi, panel_id,               &
                                             modeldb%clock, tmp_calendar, &
                                             before_close,                &
                                             start_at_zero=.true. )
    ! Finalise XIOS context
    call io_context%finalise_xios_context()

    ! Remove io context from collection
    call modeldb%io_contexts%remove_context(context_name)

    ! set context back to main context
    call modeldb%io_contexts%get_io_context(io_context_name, io_context)
    call io_context%set_current()

  end subroutine iau_incs_firstfile_io

end module iau_firstfile_io_mod
