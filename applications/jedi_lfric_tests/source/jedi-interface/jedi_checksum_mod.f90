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

  use sci_checksum_alg_mod,        only: checksum_alg
  use field_array_mod,             only: field_array_type
  use field_mod,                   only: field_type
  use driver_modeldb_mod,          only: modeldb_type
  use jedi_lfric_tests_config_mod, only: test_field
  use field_collection_mod,        only: field_collection_type

  implicit none

  private
  public output_checksum, output_linear_checksum

contains

  !> @brief    Output the checksum.
  !>
  !> @param [in]    program_name     The name of program.
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

  !> @brief    Output the checksum consistent with the linear model
  !>
  !> @param [in]    program_name The name of program.
  !> @param [inout] modeldb      The modelDB that contains the fields to be
  !>                             output as a checksum.
  subroutine output_linear_checksum( program_name, modeldb  )

    use constants_mod,          only: i_def
    use formulation_config_mod, only: moisture_formulation_dry
    use namelist_mod,           only: namelist_type

    implicit none

    character(len=*),      intent(in) :: program_name
    type(modeldb_type), intent(inout) :: modeldb

    ! Local
    type(field_collection_type), pointer :: moisture_fields
    type(field_collection_type), pointer :: prognostic_fields
    type(field_array_type ),     pointer :: mr_array
    type(field_type),            pointer :: mr(:)
    type(field_type),            pointer :: theta
    type(field_type),            pointer :: u
    type(field_type),            pointer :: rho
    type(namelist_type),         pointer :: formulation_nml
    integer(kind=i_def)                  :: moisture_formulation

    nullify(moisture_fields, prognostic_fields, mr_array)
    nullify(mr, theta, u, rho, formulation_nml)

    ! Get the fields to checksum
    prognostic_fields => modeldb%fields%get_field_collection("prognostic_fields")
    call prognostic_fields%get_field('theta', theta)
    call prognostic_fields%get_field('u', u)
    call prognostic_fields%get_field('rho', rho)

    moisture_fields => modeldb%fields%get_field_collection("moisture_fields")
    call moisture_fields%get_field("mr", mr_array)
    mr => mr_array%bundle

    ! Get configuration to inform if moisture is output
    formulation_nml => modeldb%configuration%get_namelist('formulation')
    call formulation_nml%get_value( 'moisture_formulation', moisture_formulation )

    ! Write checksums to file
    if (moisture_formulation /= moisture_formulation_dry) then
      call checksum_alg(program_name, rho, 'rho', theta, 'theta', u, 'u', &
                        field_bundle=mr, bundle_name='mr')
    else
      call checksum_alg(program_name, rho, 'rho', theta, 'theta', u, 'u')
    end if

  end subroutine output_linear_checksum

end module jedi_checksum_mod
