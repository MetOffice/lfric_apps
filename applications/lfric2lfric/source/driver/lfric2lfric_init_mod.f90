!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief   Initialisation functionality for the lfric2lfric miniapp.

!> @details Handles initialisation of prognostic fields through the call to
!!          field_maker_mod.

module lfric2lfric_init_mod

  use constants_mod,              only: i_def, r_def, str_def
  use driver_modeldb_mod,         only: modeldb_type
  use field_collection_mod,       only: field_collection_type
  use lfric_xios_context_mod,     only: lfric_xios_context_type
  use log_mod,                    only: log_event, &
                                        log_level_info
  use mesh_mod,                   only: mesh_type
  use netcdf,                     only: nf90_max_name

  ! lfric2lfric mods
  use lfric2lfric_config_mod,     only: mode_ics, mode_lbc
  use lfric2lfric_field_init_mod, only: get_field_list, field_maker

  implicit none
  private
  public :: init_lfric2lfric

  contains

  !> @brief    Initialises the required fields for lfric2lfric miniapp.
  !> @details  Calls out to field_maker to initialise fields on given field
  !!           collection.
  !> @param [in,out]   modeldb                Holds model state
  !> @param [in]       context_src            The name of the XIOS context that
  !!                                          will hold the source file
  !> @param [in]       context_dst            The name of the XIOS context that
  !!                                          will hold the file to be written
  !> @param [in]       start_dump_filename    File to get field names from
  !> @param [in]       mode                   Process ics or lbcs
  !> @param [in]       origin_collection_name Holds the origin fields
  !> @param [in]       origin_mesh            Mesh to initialise 3D fields
  !> @param [in]       origin_twod_mesh       Mesh to initialise 2D fields
  !> @param [in]       target_collection_name Holds target fields
  !> @param [in]       target_mesh            Mesh for target 3D fields
  !> @param [in]       target_twod_mesh       Mesh for target 2D fields
  subroutine init_lfric2lfric( modeldb, context_src, context_dst, &
                               start_dump_filename, mode,         &
                               origin_collection_name,            &
                               origin_mesh, origin_twod_mesh,     &
                               target_collection_name,            &
                               target_mesh, target_twod_mesh  )

    implicit none

    type(modeldb_type), intent(inout)       :: modeldb
    character(len=*),   intent(in)          :: context_src
    character(len=*),   intent(in)          :: context_dst
    character(len=*),   intent(in)          :: start_dump_filename
    integer(i_def),     intent(in)          :: mode
    character(len=*),   intent(in)          :: origin_collection_name
    type(mesh_type),    intent(in), pointer :: origin_mesh
    type(mesh_type),    intent(in), pointer :: origin_twod_mesh
    ! Optionals
    character(len=*),   intent(in)          :: target_collection_name
    type(mesh_type),    intent(in), pointer :: target_mesh
    type(mesh_type),    intent(in), pointer :: target_twod_mesh

    ! For field creation and storage
    type(field_collection_type), pointer :: field_collection

    ! For get_field_list returns
    integer(kind=i_def)                 :: num_fields
    character(len=str_def), allocatable :: config_list(:)
    character(len=nf90_max_name)        :: prefix

    ! Source context pointer and temporary context for setup
    type(lfric_xios_context_type), pointer :: io_context

    ! Looping variable
    integer(kind=i_def) :: i

    call log_event( 'lfric2lfric: Initialising miniapp ...', log_level_info )

    if (mode == mode_ics) then
      prefix = 'restart_'
    else if (mode == mode_lbc) then
      prefix = ''
    end if

    ! Get field names from filename and validate presence in iodef.xml
    call get_field_list( num_fields, config_list, start_dump_filename, prefix )

    !--------------------------------------------------------------------------
    ! Initialise Source Fields
    !--------------------------------------------------------------------------
    ! Initialise our field collection
    call modeldb%fields%add_empty_field_collection(origin_collection_name)
    field_collection => modeldb%fields%get_field_collection(origin_collection_name)

    ! Now need to loop over length of config_list make field for each
    do i = 1, num_fields
      call field_maker( field_collection, &
                        config_list(i),   &
                        origin_mesh,      &
                        origin_twod_mesh, &
                        prefix )
    end do

    !--------------------------------------------------------------------------
    ! Initialise Target Fields
    !--------------------------------------------------------------------------
    call modeldb%fields%add_empty_field_collection(target_collection_name)
    field_collection => &
                    modeldb%fields%get_field_collection(target_collection_name)

    call modeldb%io_contexts%get_io_context(context_dst, io_context)
    call io_context%set_current()

    if (mode == mode_ics) then
      prefix = 'checkpoint_'
    else if (mode == mode_lbc) then
      prefix = 'lbc_'
    end if

    do i = 1, num_fields
      call field_maker( field_collection, &
                        config_list(i),   &
                        target_mesh,      &
                        target_twod_mesh, &
                        prefix )
    end do

    call modeldb%io_contexts%get_io_context(context_src, io_context)
    call io_context%set_current()

    ! Now finished with config_list, deallocate
    deallocate(config_list)

    call log_event( 'lfric2lfric: Miniapp initialised', log_level_info )

  end subroutine init_lfric2lfric

end module lfric2lfric_init_mod
