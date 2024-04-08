!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Provides a method for writing out a checksum
!>
!> This is a temporary solution based on the original jedi_lfric_tests implementation.
!>
module jedi_checksum_mod

  use checksum_alg_mod,                    only: checksum_alg
  use field_mod,                           only: field_type
  use jedi_lfric_tests_config_mod,         only: test_field
  use field_collection_mod,                only: field_collection_type

  implicit none

  private
  public output_checksum

contains

  !> @brief    Output the checksum.
  !>
  !> @param [in]    program_name The name of program.
  !> @param [inout] field_collection A collection that contains the field
  !>                                 to be output as a checksum.
  subroutine output_checksum( program_name, field_collection  )

    implicit none

    character(len=*), intent(in)               :: program_name
    type(field_collection_type), intent(inout) :: field_collection

    type(field_type), pointer :: working_field => null()

    !@todo: we should provide a means to checksum all fields in the field
    !       collection instead of just the test_field
    call field_collection%get_field( test_field, working_field )

    !---------------------------------------------------------------------------
    ! Model finalise
    !---------------------------------------------------------------------------
    ! Write checksums to file
    call checksum_alg( program_name, working_field, test_field )

  end subroutine output_checksum

end module jedi_checksum_mod
