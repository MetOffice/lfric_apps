!-------------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------

module conditional_collection_copy_mod

  !> @brief Make a deep copy of a subset of a field collection, based on a
  !!        list of fields provided by another collection
  !> @param[out] generic_fields_copied  New collection written with fields
  !!                                    saved at some point in time
  !> @param[in]  generic_fields_to_copy collection we want to save fields from
  !> @param[in]  field_list             list of fields to be saved
  subroutine conditional_collection_copy(generic_fields_copied,                &
                                         generic_fields_to_copy,               &
                                         field_list)

    use field_collection_mod,          only: field_collection_type
    use field_collection_iterator_mod, only: field_collection_iterator_type

    implicit none

    type(field_collection_type), intent(out) :: generic_fields_copied
    type(field_collection_type), intent(in)  :: generic_fields_to_copy
    type(field_collection_type), intent(in)  :: field_list

    ! Iterator for field collection
    type(field_collection_iterator_type) :: iterator

    ! One of the single fields out of the generic_fields_to_copy collection
    class(field_parent_type), pointer :: abstract_field_ptr
    type(field_type),         pointer :: single_generic_field

    ! The saved version of single_generic_field
    type(field_type) :: copied_generic_field

    logical(kind=l_def) :: l_copy

    nullify(abstract_field_ptr)
    nullify(single_generic_field)

    call generic_fields_copied%initialise(name='fields_copied')

    if ( generic_fields_to_copy%get_length() > 0 ) then

      call iterator%initialise(generic_fields_to_copy)

      do
        if ( .not.iterator%has_next() ) exit

        abstract_field_ptr => iterator%next()
        select type(abstract_field_ptr)
        type is (field_type)
          single_generic_field => abstract_field_ptr
        end select

        l_copy = field_list%field_exists(single_generic_field%get_name())

        if ( l_copy ) then

          ! We copy the field we want to save into a new field
          call single_generic_field%copy_field_properties(copied_generic_field)
          call invoke( setval_X(copied_generic_field, single_generic_field) )

          ! We add it to the field collection passed in
          call generic_fields_copied%add_field(copied_generic_field)

        end if ! ( l_copy )

      end do

    end if

  end subroutine conditional_collection_copy

end module conditional_collection_copy_mod