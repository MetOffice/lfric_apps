!-------------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Create fields from XIOS metadata

module field_from_metadata_mod

  use log_mod,                              only: log_event,         &
                                                  log_scratch_space, &
                                                  log_level_error
  use constants_mod,                        only: i_def, l_def, str_def
  use space_from_metadata_mod,              only: space_from_metadata
  use mesh_mod,                             only: mesh_type
  use field_mod,                            only: field_type
  use integer_field_mod,                    only: integer_field_type
  use function_space_mod,                   only: function_space_type
  use empty_data_mod,                       only: empty_real_data,    &
                                                  empty_integer_data
  implicit none

  private

  interface init_field_from_metadata
    module procedure init_real_field_from_metadata
    module procedure init_integer_field_from_metadata
  end interface init_field_from_metadata

  public :: init_field_from_metadata

contains

  !> @brief Derive field status for prognostic and diagnostic fields
  !> @param[in, optional] empty            Create empty field?
  !> @param[in, optional] diag_status      Diagnostic enabling status
  !> @return                               Field status
  function combined_status(empty, diag_status) result(status)
    implicit none

    logical(l_def),         intent(in) :: empty
    character(*), optional, intent(in) :: diag_status

    character(:), allocatable :: status

    if (present(diag_status)) then
      status = trim(diag_status)
    else
      if (empty) then
        status = 'Empty'
      else
        status = 'Progno'
      end if
    end if
  end function combined_status

  !> @brief Initialise a real-valued field from XIOS metadata.
  !> @param[out]          field            Field to initialise
  !> @param[in]           xios_id          XIOS id of field
  !> @param[in, optional] empty            Create empty field?
  !> @param[in, optional] mesh_3d          Override derived 3D mesh
  !> @param[in, optional] mesh_2d          Override derived 2D mesh
  !> @param[in, optional] force_mesh       Override derived mesh
  !> @param[in, optional] force_rad_levels Override derived radiation levels
  !> @param[in, optional] force_order_h    Override derived order in horizontal
  !> @param[in, optional] force_order_v    Override derived order in vertical
  !> @param[in, optional] force_ndata      Override derived ndata
  !> @param[in, optional] diag_status      Diagnostic enabling status
  subroutine init_real_field_from_metadata(field, xios_id, empty,           &
                                           mesh_3d, mesh_2d, force_mesh,    &
                                           force_rad_levels, force_order_h, &
                                           force_order_v, force_ndata,      &
                                           diag_status)

    implicit none

    type(field_type),                   intent(out) :: field
    character(*),                       intent(in)  :: xios_id
    logical(l_def), optional,           intent(in)  :: empty
    type(mesh_type), pointer, optional, intent(in)  :: mesh_3d
    type(mesh_type), pointer, optional, intent(in)  :: mesh_2d
    type(mesh_type), pointer, optional, intent(in)  :: force_mesh
    integer(kind=i_def), optional,      intent(in)  :: force_rad_levels
    integer(kind=i_def), optional,      intent(in)  :: force_order_h
    integer(kind=i_def), optional,      intent(in)  :: force_order_v
    integer(kind=i_def), optional,      intent(in)  :: force_ndata
    character(*), optional,             intent(in)  :: diag_status

    character(str_def), parameter :: routine_name = 'init_real_field_from_metadata'

    logical(l_def)                      :: make_empty
    character(:), allocatable           :: status
    type(function_space_type),  pointer :: vector_space => null()

    make_empty = .false.
    if (present(empty)) make_empty = empty

    status = trim(combined_status(make_empty, diag_status))

    vector_space => space_from_metadata(xios_id, status,                 &
                                        mesh_3d, mesh_2d, force_mesh,    &
                                        force_rad_levels, force_order_h, &
                                        force_order_v, force_ndata)

    if (make_empty) then
      call field%initialise(                                                  &
        vector_space = vector_space,                                          &
        override_data = empty_real_data,                                      &
        name = xios_id)
    else
      call field%initialise(                                                  &
        vector_space = vector_space,                                          &
        name = xios_id)
    end if

    ! check post-condition
    if (field%get_name() /= xios_id) then
      write(log_scratch_space,*) trim(routine_name) // ': ' //                &
       'name mismatch: ' // trim(field%get_name()) //' vs '// trim(xios_id)
       call log_event(log_scratch_space, log_level_error)
    end if

    ! paranoia
    nullify(vector_space)
  end subroutine init_real_field_from_metadata

  !> @brief Initialise an integer-valued field from XIOS metadata.
  !> @param[out]          field            Field to initialise
  !> @param[in]           xios_id          XIOS id of field
  !> @param[in, optional] empty            Create empty field?
  !> @param[in, optional] mesh_3d          Override derived 3D mesh
  !> @param[in, optional] mesh_2d          Override derived 2D mesh
  !> @param[in, optional] force_mesh       Override derived mesh
  !> @param[in, optional] force_rad_levels Override derived radiation levels
  !> @param[in, optional] force_order_h    Override derived order in horizontal
  !> @param[in, optional] force_order_v    Override derived order in vertical
  !> @param[in, optional] force_ndata      Override derived ndata
  !> @param[in, optional] diag_status      Diagnostic enabling status
  subroutine init_integer_field_from_metadata(field, xios_id, empty,           &
                                              mesh_3d, mesh_2d, force_mesh,    &
                                              force_rad_levels, force_order_h, &
                                              force_order_v, force_ndata,      &
                                              diag_status)

    implicit none

    type(integer_field_type),           intent(out) :: field
    character(*),                       intent(in)  :: xios_id
    logical(l_def), optional,           intent(in)  :: empty
    type(mesh_type), pointer, optional, intent(in)  :: mesh_3d
    type(mesh_type), pointer, optional, intent(in)  :: mesh_2d
    type(mesh_type), pointer, optional, intent(in)  :: force_mesh
    integer(kind=i_def), optional,      intent(in)  :: force_rad_levels
    integer(kind=i_def), optional,      intent(in)  :: force_order_h
    integer(kind=i_def), optional,      intent(in)  :: force_order_v
    integer(kind=i_def), optional,      intent(in)  :: force_ndata
    character(*), optional,             intent(in)  :: diag_status

    character(str_def), parameter :: routine_name = 'init_field_from_metadata'

    logical(l_def)                      :: make_empty
    character(str_def)                  :: status
    type(function_space_type),  pointer :: vector_space => null()

    make_empty = .false.
    if (present(empty)) make_empty = empty

    status = combined_status(make_empty, diag_status)

    vector_space => space_from_metadata(xios_id, status,                 &
                                        mesh_3d, mesh_2d, force_mesh,    &
                                        force_rad_levels, force_order_h, &
                                        force_order_v, force_ndata)

    if (make_empty) then
      call field%initialise(                                                  &
        vector_space = vector_space,                                          &
        override_data = empty_integer_data,                                   &
        name = xios_id)
    else
      call field%initialise(                                                  &
        vector_space = vector_space,                                          &
        name = xios_id)
    end if

    ! check post-condition
    if (field%get_name() /= xios_id) then
      write(log_scratch_space,*) trim(routine_name) // ': ' //                &
       'name mismatch: ' // trim(field%get_name()) //' vs '// trim(xios_id)
       call log_event(log_scratch_space, log_level_error)
    end if

    ! paranoia
    nullify(vector_space)
  end subroutine init_integer_field_from_metadata

end module field_from_metadata_mod
