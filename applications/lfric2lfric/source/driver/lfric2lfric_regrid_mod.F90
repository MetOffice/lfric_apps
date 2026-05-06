!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Performs regridding of a field collection
!>
module lfric2lfric_regrid_mod

  use constants_mod,            only: str_def, i_def
  use driver_modeldb_mod,       only: modeldb_type
  use field_parent_mod,         only: field_parent_type
  use field_mod,                only: field_type
  use field_collection_mod,     only: field_collection_type
  use field_collection_iterator_mod, only: &
                                      field_collection_iterator_type
  use log_mod,                  only: log_event, &
                                      log_level_info, &
                                      log_scratch_space
  use model_clock_mod,          only: model_clock_type

  !------------------------------------
  ! lfric2lfric modules
  !------------------------------------
  use lfric2lfric_config_mod,         only: regrid_method_map,         &
                                            regrid_method_lfric2lfric, &
                                            regrid_method_oasis
  use lfric2lfric_map_regrid_mod,     only: lfric2lfric_map_regrid
  use lfric2lfric_no_regrid_mod,      only: lfric2lfric_no_regrid
  use lfric2lfric_oasis_regrid_mod,   only: lfric2lfric_oasis_regrid

  implicit none

  private
  public lfric2lfric_regrid

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief    Performs regridding of a field collection
  !> @details  Over one time step, all regridding is performed by algorithm
  !!           modules specific to the regrid method defined in the
  !!           configuration.
  !!           Fields to be regridded are extracted from the source field
  !!           collection. Corresponding source and target field pairs
  !!           are passed to the regridding algorithm, and written to fields
  !!           in the destination field collection.
  !> @param [in,out] modeldb                 The structure that holds model
  !!                                         state
  !> @param [in,out] oasis_clock             Clock for OASIS exchanges
  !> @param [in]     source_fields           Collection of fields to be regridded
  !> @param [in]     target_fields           Collection of regridded fields
  !> @param [in]     regrid_method           Method for regridding between the
  !>                                         source and destination meshes
  subroutine lfric2lfric_regrid( modeldb, oasis_clock,            &
                  source_fields, target_fields, regrid_method )

    implicit none

    type(modeldb_type),                   intent(inout) :: modeldb
    type(model_clock_type), allocatable,  intent(inout) :: oasis_clock
    type(field_collection_type), pointer, intent(in)    :: source_fields
    type(field_collection_type), pointer, intent(inout) :: target_fields
    integer(kind=i_def),                  intent(in)    :: regrid_method


    type(field_collection_iterator_type) :: iter

    class(field_parent_type), pointer :: field     => null()
    type(field_type),         pointer :: field_src => null()
    type(field_type),         pointer :: field_dst => null()

    character(len=str_def)   :: field_name

    ! Main loop over fields to be processed
    call iter%initialise(source_fields)
    do
      ! Locate the field to be processed in the field collections
      if ( .not.iter%has_next() ) exit
      field => iter%next()
      field_name = field%get_name()

      call source_fields%get_field(field_name, field_src)
      call target_fields%get_field(field_name, field_dst)

      write(log_scratch_space, '(A,A)') "Processing lfric field ", &
                                           trim(field_name)
      call log_event(log_scratch_space, log_level_info)

      ! Regrid source field depending on regrid method
      select case (regrid_method)
        case (regrid_method_map)
          call lfric2lfric_map_regrid(field_dst, field_src)

        case (regrid_method_lfric2lfric)
          write(log_scratch_space, '(A)') &
                              'Regrid method lfric2lfric not implemented yet'
          call log_event(log_scratch_space, log_level_info)

          call lfric2lfric_no_regrid(field_dst)

        case (regrid_method_oasis)
#ifdef MCT
          call lfric2lfric_oasis_regrid(modeldb, oasis_clock, &
                                        field_dst, field_src)
#endif
      end select
    end do

  end subroutine lfric2lfric_regrid

end module lfric2lfric_regrid_mod
