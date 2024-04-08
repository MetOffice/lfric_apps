!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Update the I/O field collection

!> @details Handles update of a field collection to include a set of fields
!>          defined by input meta data. This is for the purpose of performing
!>          file I/O via the XIOS read_state and write_state methods. The
!>          collection is updated so that fields are added if they are not
!>          present and removed if they are no longer needed.
module jedi_lfric_io_update_mod

  use constants_mod,                  only : l_def, i_def
  use log_mod,                        only : log_event, &
                                             LOG_LEVEL_INFO
  use mesh_mod,                       only : mesh_type
  use jedi_lfric_utils_mod,           only : add_real_field
  use jedi_lfric_field_meta_mod,      only : jedi_lfric_field_meta_type
  use field_collection_mod,           only : field_collection_type
  use field_collection_iterator_mod,  only : field_collection_iterator_type
  use field_parent_mod,               only : field_parent_type

  implicit none

  contains

  !> @brief        Sets up a field collection containing only fields defined in
  !>               the input meta-data
  !>
  !> @description  A field collection is updated to ensure that it contains
  !>               only the fields defined in the input-meta-data. The field
  !>               name and function_space to be used to create the fields is
  !>               defined in the meta-data. A logical is also stored in the
  !>               meta-data to define if the field is 2D. Both meshes are
  !>               passed in so that both 2D and 3D fields can be created.
  !>               Fields that are no-longer required are removed from the
  !>               collection as part of this routine.
  !>
  !> @param[inout] io_collection   the io_collection to update
  !> @param[in] mesh_3d            the <mesh_type> instance used for 3D fields
  !> @param[in] mesh_2d            the <mesh_type> instance used for 2D fields
  !> @param[in] field_meta_data    meta-data describing the fields to create
  subroutine update_io_field_collection( io_collection, mesh_3d, mesh_2d, &
                                         field_meta_data )

    implicit none

    type(field_collection_type), intent(inout)   :: io_collection
    type(mesh_type), pointer, intent(in)         :: mesh_3d
    type(mesh_type), pointer, intent(in)         :: mesh_2d
    type(jedi_lfric_field_meta_type), intent(in) :: field_meta_data

    ! Local
    class(field_parent_type), pointer     :: field_ptr => null()
    type(field_collection_iterator_type)  :: iter
    type(mesh_type), pointer              :: mesh
    logical(kind=l_def)                   :: remove_current_field
    logical(kind=l_def), allocatable      :: field_to_add(:)
    integer(kind=i_def)                   :: i

    call log_event( 'jedi_lfric: update_io_field_collection - start', LOG_LEVEL_INFO )

    ! Setup logical array to track which fields need to be added
    allocate( field_to_add(field_meta_data%get_n_variables()) )
    field_to_add = .true.

    ! If io_collection is not empty, loop over all elements to:
    ! 1. Remove fields that are not required,
    ! 2. Set field_to_add(i) = .false. if field is already
    !    present so it is not subsequently added
    if (io_collection%get_length() /= 0) then
      ! Initialise iterator (iter) on io_collection.
      call iter%initialise(io_collection)
      do
        if ( .not. iter%has_next() ) exit
        field_ptr => iter%next()

        ! Check to see if the field is required
        remove_current_field = .true.
        do i=1,field_meta_data%get_n_variables()
          if (field_ptr%get_name() == field_meta_data%get_variable_name(i)) then
            field_to_add(i) = .false.
            remove_current_field = .false.
            exit
          end if
        end do

        ! Remove fields that are not required
        if (remove_current_field) then
          call io_collection%remove_field(field_ptr%get_name())
        end if
      end do
    end if

    ! Add fields not found in the io_collection
    do i=1,field_meta_data%get_n_variables()
      if (field_to_add(i)) then
        if (field_meta_data%get_variable_is_2d(i)) then
          mesh => mesh_2d
        else
          mesh => mesh_3d
        end if
        call add_real_field( io_collection, mesh,                            &
                             field_meta_data%get_variable_function_space(i), &
                             field_meta_data%get_variable_name(i) )
      end if
    end do

    call log_event( 'jedi_lfric: update_io_field_collection - end', LOG_LEVEL_INFO )

  end subroutine update_io_field_collection

end module jedi_lfric_io_update_mod
