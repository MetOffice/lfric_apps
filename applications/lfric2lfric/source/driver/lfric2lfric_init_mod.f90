!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief   Initialisation functionality for the lfric2lfric miniapp.

!> @details Handles initialisation of prognostic fields through the call to
!!          field_maker_mod.

module lfric2lfric_init_mod

  use constants_mod,              only : i_def, r_def, str_def
  use driver_modeldb_mod,         only : modeldb_type
  use field_collection_mod,       only : field_collection_type
  use log_mod,                    only : log_event, &
                                         log_level_info
  use mesh_mod,                   only : mesh_type

  ! lfric2lfric mods
  use lfric2lfric_field_init_mod, only : get_field_list, field_maker

  implicit none
  private
  public :: init_lfric2lfric

  contains

  !> @brief    Initialises the required fields for lfric2lfric miniapp.
  !> @details  Calls out to field_maker to initialise fields on given field
  !!           collection.
  !> @param [in,out]   modeldb                Holds model state
  !> @param [in]       start_dump_filename    File to get field names from
  !> @param [in]       origin_collection_name Holds the origin fields
  !> @param [in]       origin_mesh            Mesh to initialise 3D fields
  !> @param [in]       origin_twod_mesh       Mesh to initialise 2D fields
  !> @param [in]       target_collection_name Holds target fields
  !> @param [in]       target_mesh            Mesh for target 3D fields
  !> @param [in]       target_twod_mesh       Mesh for target 2D fields
  subroutine init_lfric2lfric( modeldb, start_dump_filename,  &
                               origin_collection_name,        &
                               origin_mesh, origin_twod_mesh, &
                               target_collection_name,        &
                               target_mesh, target_twod_mesh  )

    implicit none

    type(modeldb_type), intent(inout)       :: modeldb
    character(len=*),   intent(in)          :: start_dump_filename
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

    ! Looping variable
    integer(kind=i_def) :: i

    call log_event( 'lfric2lfric: Initialising miniapp ...', log_level_info )

    ! Get field names from filename and validate presence in iodef.xml
    call get_field_list( num_fields, config_list, start_dump_filename )

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
                        origin_twod_mesh  )
    end do

    !--------------------------------------------------------------------------
    ! Initialise Target Fields
    !--------------------------------------------------------------------------
    call modeldb%fields%add_empty_field_collection(target_collection_name)
    field_collection => &
                    modeldb%fields%get_field_collection(target_collection_name)

    do i = 1, num_fields
      call field_maker( field_collection, &
                        config_list(i),   &
                        target_mesh,      &
                        target_twod_mesh  )
    end do

    ! Now finished with config_list, deallocate
    deallocate(config_list)

    call log_event( 'lfric2lfric: Miniapp initialised', log_level_info )

  end subroutine init_lfric2lfric

end module lfric2lfric_init_mod
