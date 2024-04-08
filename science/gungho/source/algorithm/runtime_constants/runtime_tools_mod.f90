!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Provides some generic routines used in setting up runtime constants
!>
!> @details This module provides a few tools common to each of the sets of
!>          runtime constants groups:
!>          - Global mesh ID lists are provided here as they are needed
!>            throughout run-time by the runtime constants getter functions.
!>          - Meshes are given labels to describe what type of mesh they are.
!>            This is used to determine which constants to set up on which mesh,
!>            as many constants aren't needed on every mesh. The enumerated
!>            labels are provided here.
!>          Having these in a separate module is necessary for avoiding circular
!>          dependency trees and duplicated code.
module runtime_tools_mod

  use constants_mod,     only: i_def, str_def
  use field_r32_mod,     only: field_r32_type
  use field_r64_mod,     only: field_r64_type
  use field_mod,         only: field_type
  use r_tran_field_mod,  only: r_tran_field_type
  use integer_field_mod, only: integer_field_type
  use fs_continuity_mod, only: W0, W1, W2, W2H, W2V, W3,    &
                               W2trace, W2Htrace, W2Vtrace, &
                               W2broken, Wtheta, Wchi
  use log_mod,           only: log_event, log_scratch_space, LOG_LEVEL_ERROR
  use operator_mod,      only: operator_type

  implicit none

  private

  ! Mesh IDs
  integer(kind=i_def), allocatable :: global_mesh_id_list(:)
  integer(kind=i_def), allocatable :: hierarchical_mesh_id_list(:)

  ! Enumerated types to label what type the mesh is (randomly selected integers)
  integer(kind=i_def), parameter, public :: primary_mesh_label = 19
  integer(kind=i_def), parameter, public :: shifted_mesh_label = 964
  integer(kind=i_def), parameter, public :: double_level_mesh_label = 209
  integer(kind=i_def), parameter, public :: multigrid_mesh_label = 624
  integer(kind=i_def), parameter, public :: twod_mesh_label = 71
  integer(kind=i_def), parameter, public :: extra_mesh_label = 117

  ! Public functions to create and access the module contents
  public :: init_mesh_id_list
  public :: final_mesh_id_list
  public :: find_mesh_index
  public :: init_hierarchical_mesh_id_list
  public :: final_hierarchical_mesh_id_list
  public :: get_hierarchical_mesh_id
  public :: check_initialised_field
  public :: check_initialised_operator

  interface check_initialised_field
     procedure check_initialised_field_r32,   &
               check_initialised_field_integer, &
               check_initialised_field_r64
  end interface check_initialised_field

contains
  !> @brief Subroutine to initialise mesh ID list
  !> @param[in] mesh_id_list          List of mesh IDs
  subroutine init_mesh_id_list(mesh_id_list)

    implicit none

    integer(kind=i_def), intent(in) :: mesh_id_list(:)

    ! Internal variables
    integer(kind=i_def) :: i, number_of_meshes

    number_of_meshes = size(mesh_id_list)

    allocate(global_mesh_id_list(number_of_meshes))

    do i = 1, number_of_meshes
      global_mesh_id_list(i) = mesh_id_list(i)
    end do

  end subroutine init_mesh_id_list

  !> @brief Deallocates the mesh ID list
  subroutine final_mesh_id_list()
    implicit none

    if (allocated(global_mesh_id_list)) deallocate(global_mesh_id_list)

  end subroutine final_mesh_id_list

  !> @brief Returns the mesh index in the mesh list from the mesh_id
  !> @param[in] mesh_id
  !> @return The mesh index
  function find_mesh_index(mesh_id) result(mesh_index)
    ! TODO: This routine should be removed in #2790. It should be possible
    ! to return the appropriate constant without relying on a mesh_index

    implicit none
    integer(kind=i_def), intent(in) :: mesh_id
    integer(kind=i_def)             :: i, mesh_index

    mesh_index = 0_i_def
    do i = 1, size(global_mesh_id_list)
      if ( global_mesh_id_list(i) == mesh_id ) then
        mesh_index = i
        exit
      end if
    end do

    if ( mesh_index == 0_i_def ) then
      call log_event( "mesh_id cannot be found in mesh_id_list", LOG_LEVEL_ERROR )
    end if
  end function find_mesh_index

  !> @brief Subroutine to initialise hierarchical mesh ID list
  !> @param[in] mesh_id_list          List of mesh IDs
  subroutine init_hierarchical_mesh_id_list(mesh_id_list)

    implicit none

    integer(kind=i_def), intent(in) :: mesh_id_list(:)

    ! Internal variables
    integer(kind=i_def)             :: num_meshes, i

    num_meshes = size(mesh_id_list)

    allocate(hierarchical_mesh_id_list(num_meshes))

    do i = 1, num_meshes
      hierarchical_mesh_id_list(i) = mesh_id_list(i)
    end do

  end subroutine init_hierarchical_mesh_id_list

  !> @brief Deallocates the mesh ID list
  subroutine final_hierarchical_mesh_id_list()
    implicit none

    if (allocated(hierarchical_mesh_id_list)) deallocate(hierarchical_mesh_id_list)

  end subroutine final_hierarchical_mesh_id_list

  !> @brief Gets the mesh id of a given hierarchical level
  function get_hierarchical_mesh_id(level) result(mesh_id)
    implicit none
    integer(kind=i_def), intent(in) :: level
    integer(kind=i_def)             :: mesh_id

    mesh_id = hierarchical_mesh_id_list(level)

  end function get_hierarchical_mesh_id

  !> @brief Checks whether a field is initialised and returns an error if not
  !> @param[in] field      The field to check
  !> @param[in] field_name A name of the field to include in the error message
  !> @param[in] mesh_id    ID of the mesh
  !> @param[in] space      An optional integer representing the function space
  subroutine check_initialised_field_r32(field, field_name, mesh_id, space)
    implicit none
    type(field_r32_type),          intent(in) :: field
    character(str_def),            intent(in) :: field_name
    integer(kind=i_def),           intent(in) :: mesh_id
    integer(kind=i_def), optional, intent(in) :: space

    if (.not. field%is_initialised()) then
      if (present(space)) then
        write(log_scratch_space, '(A,A,I3,A,A)') &
        trim(field_name), ' on mesh ', mesh_id, ' not initialised for ', find_space_name(space)
      else
        write(log_scratch_space, '(A,A,I3,A)') &
        trim(field_name), ' on mesh ', mesh_id, ' not initialised'
      end if
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if

  end subroutine check_initialised_field_r32

  !> @brief Checks whether a field is initialised and returns an error if not
  !> @param[in] field      The field to check
  !> @param[in] field_name A name of the field to include in the error message
  !> @param[in] mesh_id    ID of the mesh
  !> @param[in] space      An optional integer representing the function space
  subroutine check_initialised_field_r64(field, field_name, mesh_id, space)
    implicit none
    type(field_r64_type),          intent(in) :: field
    character(str_def),            intent(in) :: field_name
    integer(kind=i_def),           intent(in) :: mesh_id
    integer(kind=i_def), optional, intent(in) :: space

    if (.not. field%is_initialised()) then
      if (present(space)) then
        write(log_scratch_space, '(A,A,I3,A,A)') &
        trim(field_name), ' on mesh ', mesh_id, ' not initialised for ', find_space_name(space)
      else
        write(log_scratch_space, '(A,A,I3,A)') &
        trim(field_name), ' on mesh ', mesh_id, ' not initialised'
      end if
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if

  end subroutine check_initialised_field_r64

  !> @brief Checks whether a field is initialised and returns an error if not
  !> @param[in] field      The field to check
  !> @param[in] field_name A name of the field to include in the error message
  !> @param[in] mesh_id    ID of the mesh
  !> @param[in] space      An optional integer representing the function space
  subroutine check_initialised_field_integer(field, field_name, mesh_id, space)
    implicit none
    type(integer_field_type),       intent(in) :: field
    character(str_def),            intent(in) :: field_name
    integer(kind=i_def),           intent(in) :: mesh_id
    integer(kind=i_def), optional, intent(in) :: space

    if (.not. field%is_initialised()) then
      if (present(space)) then
        write(log_scratch_space, '(A,A,I3,A,A)') &
        trim(field_name), ' on mesh ', mesh_id, ' not initialised for ', find_space_name(space)
      else
        write(log_scratch_space, '(A,A,I3,A)') &
        trim(field_name), ' on mesh ', mesh_id, ' not initialised'
      end if
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if

  end subroutine check_initialised_field_integer

  !> @brief Checks whether an operator is initialised and returns an error if not
  !> @param[in] operator      The operator to check
  !> @param[in] operator_name A name of the operator to include in the error message
  !> @param[in] mesh_id       ID of the mesh
  !> @param[in] space         An optional integer representing the function space
  subroutine check_initialised_operator(operator, operator_name, mesh_id, space)
    implicit none
    type(operator_type),           intent(in) :: operator
    character(str_def),            intent(in) :: operator_name
    integer(kind=i_def),           intent(in) :: mesh_id
    integer(kind=i_def), optional, intent(in) :: space

    if (.not. operator%is_initialised()) then
      if (present(space)) then
        write(log_scratch_space, '(A,A,I3,A,A)') &
        trim(operator_name), ' on mesh ', mesh_id, ' not initialised for ', find_space_name(space)
      else
        write(log_scratch_space, '(A,A,I3,A)') &
        trim(operator_name), ' on mesh ', mesh_id, ' not initialised'
      end if
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if

  end subroutine check_initialised_operator

  !> Converts the integer function space name into a human readable string
  function find_space_name(space) result(space_name)
    implicit none
    integer(kind=i_def), intent(in) :: space
    character(str_def)              :: space_name

    select case(space)
    case (W0)
      space_name = 'W0'
    case (W1)
      space_name = 'W1'
    case (W2)
      space_name = 'W2'
    case (W2V)
      space_name = 'W2V'
    case (W2H)
      space_name = 'W2H'
    case (W2broken)
      space_name = 'W2broken'
    case (W2trace)
      space_name = 'W2trace'
    case (W2Vtrace)
      space_name = 'W2Vtrace'
    case (W2Htrace)
      space_name = 'W2Htrace'
    case (W3)
      space_name = 'W3'
    case (Wtheta)
      space_name = 'Wtheta'
    case (Wchi)
      space_name = 'Wchi'
    case default
      call log_event("Space not identified", LOG_LEVEL_ERROR)
    end select

  end function find_space_name

end module runtime_tools_mod
