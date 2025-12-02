!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
module test_helpers_mod

  use, intrinsic :: iso_fortran_env, only : real32, real64, output_unit

  use field_parent_mod,   only : field_parent_type
  use field_real32_mod,   only : field_real32_type
  use field_real64_mod,   only : field_real64_type
  use function_space_mod, only : function_space_type
  use lfric_mpi_mod,      only : global_mpi
  use test_field_mod,     only : test_field_32_type, test_field_64_type

  implicit none

  private
  public :: create_test_field_flood_constant, print_test_field

  interface create_test_field_flood_constant
    procedure create_test_field_32_flood_constant
    procedure create_test_field_64_flood_constant
  end interface create_test_field_flood_constant

  interface print_test_field
    procedure print_test_field_32
    procedure print_test_field_64
  end interface print_test_field

contains

  subroutine create_test_field_32_flood_constant( fs, name,                &
                                                  lfric_field, test_field, &
                                                  value )

    implicit none

    class(function_space_type), intent(in), pointer      :: fs
    character(*),               intent(in)               :: name
    class(field_real32_type),   intent(out), &
                                     allocatable, target :: lfric_field
    class(test_field_32_type),  intent(out), allocatable :: test_field
    real(real32),               intent(in)               :: value

    class(field_real32_type), pointer :: lfric_field_ptr
    real(real32),             pointer :: test_data(:)

    allocate(lfric_field)
    call lfric_field%initialise( fs, name=name )
    lfric_field_ptr => lfric_field
    allocate(test_field, source=test_field_32_type( lfric_field_ptr ))
    test_data => test_field%get_data()
    test_data = value
    call test_field%copy_to_lfric()

  end subroutine create_test_field_32_flood_constant


  subroutine create_test_field_64_flood_constant( fs, name,                &
                                                  lfric_field, test_field, &
                                                  value )

    implicit none

    type(function_space_type), intent(in), pointer      :: fs
    character(*),              intent(in)               :: name
    type(field_real64_type),   intent(out), &
                                    allocatable, target :: lfric_field
    type(test_field_64_type),  intent(out), allocatable :: test_field
    real(real64),              intent(in)               :: value

    type(field_real64_type), pointer :: lfric_field_ptr
    real(real64),            pointer :: test_data(:)

    allocate(lfric_field)
    call lfric_field%initialise( fs, name=name )
    lfric_field_ptr => lfric_field
    allocate(test_field, source=test_field_64_type( lfric_field_ptr ))
    test_data => test_field%get_data()
    test_data = value
    call test_field%copy_to_lfric()

  end subroutine create_test_field_64_flood_constant


  subroutine print_test_field_32( test_field )

    implicit none

    type(test_field_32_type), intent(inout) :: test_field

    class(field_parent_type),  pointer :: lfric_field
    type(function_space_type), pointer :: fs
    character(255)                     :: format_str
    real(real32),              pointer :: test_data(:)

    lfric_field => test_field%get_lfric_field_ptr()
    fs => lfric_field%get_function_space()
    write(format_str, '("A, "", ("", I0, "") = "", ", I0, "F7.3")') fs%get_undf()
    call test_field%copy_from_lfric()
    test_data => test_field%get_data()
    write(output_unit, format_str) global_mpi%get_comm_rank(), test_data(:fs%get_undf())

  end subroutine print_test_field_32


  subroutine print_test_field_64( test_field )

    implicit none

    type(test_field_64_type), intent(inout) :: test_field

    class(field_parent_type),  pointer :: lfric_field
    type(function_space_type), pointer :: fs
    character(255)                     :: format_str
    real(real64),              pointer :: test_data(:)

    lfric_field => test_field%get_lfric_field_ptr()
    fs => lfric_field%get_function_space()
    write(format_str, '("(A, "" ("", I0, "") = "", ", I0, "F7.3)")') &
      fs%get_last_dof_owned()
    call test_field%copy_from_lfric()
    test_data => test_field%get_data()
    write(output_unit, format_str) trim(lfric_field%get_name()), &
                                   global_mpi%get_comm_rank(),   &
                                   test_data(:fs%get_last_dof_owned())

  end subroutine print_test_field_64

end module test_helpers_mod
