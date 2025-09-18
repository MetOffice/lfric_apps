!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the term_s
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Removes a field collection and its fields from the depository

module remove_field_collection_mod

  use constants_mod,        only : str_def
  use driver_modeldb_mod,   only : modeldb_type
  use field_collection_mod, only : field_collection_type
  use field_collection_iterator_mod, &
                            only : field_collection_iterator_type
  use field_mod,            only : field_type
  use field_parent_mod,     only : field_parent_type

  implicit none

  private
  public :: remove_field_collection

contains
  !> @brief Removes field collection and its fields from depository
  !> @param [in]  field_collection_name Name of field collection to be removed
  subroutine remove_field_collection( modeldb, field_collection_name )

    implicit none

    type( modeldb_type ), intent(inout), target  :: modeldb
    character(*),            intent(in)          :: field_collection_name

    type( field_collection_type ), pointer :: field_collection_ptr
    type( field_collection_type ), pointer :: depository
    type( field_collection_iterator_type ) :: iterator
    class( field_parent_type ), pointer    :: field_ptr
    character( str_def )                   :: name

    nullify( field_collection_ptr, depository )

    field_collection_ptr => &
       modeldb%fields%get_field_collection(trim(field_collection_name))
    depository => modeldb%fields%get_field_collection("depository")
    call iterator%initialise(field_collection_ptr)
    do
      if ( .not.iterator%has_next() ) exit
      field_ptr => iterator%next()

      select type(field_ptr)
      type is (field_type)
        name = trim(adjustl( field_ptr%get_name() ))
        call field_collection_ptr%remove_field(name)
        call depository%remove_field(name)
      end select
    end do
    field_ptr => null()

  end subroutine remove_field_collection

end module remove_field_collection_mod
