!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief  Handles initialisation of the IO file list for JEDI-LFRIC apps.
module jedi_lfric_init_files_mod

  use constants_mod,            only: i_def
  use lfric_xios_file_mod,      only: lfric_xios_file_type
  use linked_list_mod,          only: linked_list_type
  use jedi_lfric_file_meta_mod, only: jedi_lfric_file_meta_type

  implicit none

  private
  public :: jedi_lfric_init_files
contains

  !> @brief    Initialises IO file list for JEDI-LFRIC applications
  !>
  !> @details  A list of file objects is created based on file meta-data. The
  !>           file_list is a linked list of LFRic file objects that are used
  !>           in the XIOS file setup.
  !>
  !> @param[out] files_list  List to be populated with read/write files.
  !> @param[in]  file_meta   An object that includes the file meta data
  !>                         required to create the file_list.
  subroutine jedi_lfric_init_files(files_list, file_meta)

    implicit none

    type(linked_list_type),         intent(out) :: files_list
    type(jedi_lfric_file_meta_type), intent(in) :: file_meta(:)

    ! Local
    integer(i_def) :: i

    do i=1, size(file_meta)
      call files_list%insert_item(                   &
        lfric_xios_file_type(                        &
          file_name=file_meta(i)%file_name,          &
          xios_id=file_meta(i)%xios_id,              &
          io_mode=file_meta(i)%io_mode,              &
          freq=file_meta(i)%freq,                    &
          field_group_id=file_meta(i)%field_group_id &
        )                                            &
      )
    end do

  end subroutine jedi_lfric_init_files

end module jedi_lfric_init_files_mod
