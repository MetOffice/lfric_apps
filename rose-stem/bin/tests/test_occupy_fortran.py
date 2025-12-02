from pathlib import Path
from textwrap import dedent

from ..modules.occupy_fortran import entry


def test_program(tmp_path: Path):
    """
    Ensures that globals may exist in programs. There can only ever be one
    program so globals are acceptable.
    """
    test_file = tmp_path / 'program.f90'
    test_file.write_text(
        dedent("""
            program test_program
              integer :: global_var
              real, save :: save_var
              integer :: implicit_save = 1
              real, pointer :: implicit_pointer => null()
              type :: test_type
                integer :: type_var
              end type test_type
              type(test_type) :: global_type
              type(test_type), save :: save_type
            end program test_program
            """)
    )

    dirty_list, clean_list, not_considered = entry([test_file])

    assert clean_list == [test_file]
    assert dirty_list == []
    assert not_considered == []


def test_parameters(tmp_path: Path):
    """
    Parameters are immutable so are intrinsically global. We don't mind those.
    """
    test_file = tmp_path / 'program.f90'
    test_file.write_text(
        dedent("""
        module parameter_mod
          integer, parameter :: int_param = 7
        contains
          subroutine param_sub
            real, parameter :: real_param = 12.3
          end subroutine param_sub
        end module parameter_mod
        """)
    )

    dirty_list, clean_list, not_considered = entry([test_file])

    assert clean_list == [test_file]
    assert dirty_list == []
    assert not_considered == []


def test_module(tmp_path: Path):
    """
    Ensures that globals are flagged in modules.

    todo: We have not currently come to a conclusion on how to handle pointer
          initialisation. So we ignore it for the moment.
    """
    test_file = tmp_path / 'module.f90'
    test_file.write_text(
        dedent("""
        module test_module
          integer, public :: global_var
          integer, public :: implicit_var = 2
          real, public, pointer :: implicit_pointer => null()
          type :: test_type
            integer :: type_var
            real, pointer :: type_pointer => null()
          end type test_type
          type(test_type), public :: global_type
        contains
          subroutine biscuits()
            integer, save :: save_local
            integer :: implicit_local = 3
            real, pointer :: implicit_pointer => null()
          end subroutine biscuits
        end module test_module
        """)
    )

    dirty_list, clean_list, not_considered = entry([test_file])

    assert clean_list == []
    assert not_considered == []
    assert len(dirty_list) == 1
    assert dirty_list[0].filename == test_file
    assert len(dirty_list[0].dirt) == 6

    assert dirty_list[0].dirt[0].line_number == 3
    assert dirty_list[0].dirt[0].fortran_type == 'integer'
    assert dirty_list[0].dirt[0].variable_name == 'global_var'

    assert dirty_list[0].dirt[1].line_number == 4
    assert dirty_list[0].dirt[1].fortran_type == 'integer'
    assert dirty_list[0].dirt[1].variable_name == 'implicit_var'

    # todo: This trips despite our indecision on pointers as it is
    #       plain and simple a global variable.
    #
    assert dirty_list[0].dirt[2].line_number == 5
    assert dirty_list[0].dirt[2].fortran_type == 'real'
    assert dirty_list[0].dirt[2].variable_name == 'implicit_pointer'

    assert dirty_list[0].dirt[3].line_number == 10
    assert dirty_list[0].dirt[3].fortran_type == 'test_type'
    assert dirty_list[0].dirt[3].variable_name == 'global_type'

    assert dirty_list[0].dirt[4].line_number == 13
    assert dirty_list[0].dirt[4].fortran_type == 'integer'
    assert dirty_list[0].dirt[4].variable_name == 'save_local'

    assert dirty_list[0].dirt[5].line_number == 14
    assert dirty_list[0].dirt[5].fortran_type == 'integer'
    assert dirty_list[0].dirt[5].variable_name == 'implicit_local'
