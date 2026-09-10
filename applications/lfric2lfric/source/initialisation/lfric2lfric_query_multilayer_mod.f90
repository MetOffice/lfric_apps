!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Query whether extra work is required for the multi-layer fields.
module lfric2lfric_query_multilayer_mod

  use constants_mod,            only: str_def, i_def, l_def
  use driver_modeldb_mod,       only: modeldb_type
  use field_mod,                only: field_type
  use field_collection_mod,     only: field_collection_type
  use field_parent_mod,         only: field_parent_type
  use field_collection_iterator_mod, only: &
                                      field_collection_iterator_type
  use function_space_mod,       only: function_space_type
  use log_mod,                  only: log_event,         &
                                      log_scratch_space, &
                                      log_level_info,    &
                                      log_level_error

  implicit none
  
  private
  public :: query_multilayer

contains

  !> @brief Set the tile_change logical to true or false
  !> @details tile_change is true if there is a change in the number of tiles
  !!         for any field between the source_fields and target_fields.
  subroutine query_multilayer( modeldb, source_fields, target_fields )

    implicit none

    type(modeldb_type),                   intent(inout) :: modeldb
    type(field_collection_type), pointer, intent(in)    :: source_fields
    type(field_collection_type), pointer, intent(in)    :: target_fields

    type(field_collection_iterator_type) :: iter

    class(field_parent_type),  pointer :: field
    type(field_type),          pointer :: field_src
    type(field_type),          pointer :: field_dst
    type(function_space_type), pointer :: fs_type

    character(len=str_def)   :: field_name
    integer(kind=i_def)      :: ndata_src, ndata_dst
    logical(kind=l_def)      :: tile_change

    tile_change = .false.
    
    ! Main loop over fields to be processed
    call iter%initialise(source_fields)
    do
      ! Locate the field to be processed in the field collections
      if ( .not.iter%has_next() ) exit
      field => iter%next()
      field_name = field%get_name()

      call source_fields%get_field(field_name, field_src)
      call target_fields%get_field(field_name, field_dst)
      
      ! Is this a multidata field and is ndata the same?
      fs_type => field_src%get_function_space()
      ndata_src = fs_type%get_ndata()
      fs_type => field_dst%get_function_space()
      ndata_dst = fs_type%get_ndata()
      if (ndata_dst /= ndata_src) then
        write(log_scratch_space, '(A)') "ndata not equal "
        call log_event(log_scratch_space, log_level_info)

        tile_change = .true.
      end if
    end do

    call modeldb%values%add_key_value('tile_change', tile_change)

  end subroutine query_multilayer

end module lfric2lfric_query_multilayer_mod
