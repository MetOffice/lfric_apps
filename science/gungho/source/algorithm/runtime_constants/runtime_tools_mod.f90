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
  use field_parent_mod,  only: field_parent_type
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
  integer(kind=i_def), allocatable :: hierarchical_mesh_id_list(:)

  ! Public functions to create and access the module contents
  public :: init_hierarchical_mesh_id_list
  public :: final_hierarchical_mesh_id_list
  public :: get_hierarchical_mesh_id
  public :: check_initialised_field
  public :: check_initialised_operator

contains

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
  subroutine check_initialised_field(field, field_name, mesh_id, space)
    implicit none
    class(field_parent_type),      intent(in) :: field
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

  end subroutine check_initialised_field

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
