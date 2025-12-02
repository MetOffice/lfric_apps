!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
module test_field_mod

  use, intrinsic :: iso_fortran_env, only : real32, real64

  use abstract_external_field_mod, only : abstract_external_field_type
  use constants_mod,               only : i_def
  use field_real32_mod,            only : field_real32_type, &
                                          field_real32_proxy_type
  use field_real64_mod,            only : field_real64_type, &
                                          field_real64_proxy_type
  use log_mod,                     only : log_event, log_level_error
  use field_parent_mod,            only : field_parent_type

  implicit none

  private

  type, public, extends(abstract_external_field_type) :: test_field_32_type
    private
    real(real32), pointer :: test_data(:) => null()
  contains
    private
    procedure, public :: copy_from_lfric => copy_32_from_lfric
    procedure, public :: copy_to_lfric => copy_32_to_lfric
    procedure, public :: get_data => get_data_32
    final destroy_32
  end type test_field_32_type

  interface test_field_32_type
    procedure test_field_32_constructor
  end interface test_field_32_type


  type, public, extends(abstract_external_field_type) :: test_field_64_type
    private
    real(real64), pointer :: test_data(:) => null()
  contains
    private
    procedure, public :: copy_from_lfric => copy_64_from_lfric
    procedure, public :: copy_to_lfric => copy_64_to_lfric
    procedure, public :: get_data => get_data_64
    final destroy_64
  end type test_field_64_type

  interface test_field_64_type
    procedure test_field_64_constructor
  end interface test_field_64_type

contains

  function test_field_32_constructor( lfric_field ) result(new_instance)

    implicit none

    class(field_real32_type), intent(in), pointer :: lfric_field
    type(test_field_32_type) :: new_instance

    type(field_real32_proxy_type)     :: proxy
    class(field_parent_type), pointer :: cast_field

    ! This little dance with pointers is needed to keep some compilers
    ! happy that you can, in fact, pass a child class to an argument
    ! expecting a parent class.
    !
    cast_field => lfric_field
    call new_instance%abstract_external_field_initialiser( cast_field )

    ! Mirror the data array
    !
    proxy = lfric_field%get_proxy()
    allocate( new_instance%test_data(size(proxy%data)) )

  end function test_field_32_constructor


  function test_field_64_constructor( lfric_field ) result(new_instance)

    implicit none

    class(field_real64_type), intent(in), pointer :: lfric_field
    type(test_field_64_type) :: new_instance

    type(field_real64_proxy_type)     :: proxy
    class(field_parent_type), pointer :: cast_field

    ! This little dance with pointers is needed to keep some compilers
    ! happy that you can, in fact, pass a child class to an argument
    ! expecting a parent class.
    !
    cast_field => lfric_field
    call new_instance%abstract_external_field_initialiser( cast_field )

    ! Mirror the data array
    !
    proxy = lfric_field%get_proxy()
    allocate( new_instance%test_data(size(proxy%data)) )

  end function test_field_64_constructor


  subroutine destroy_32( this )

    implicit none

    type(test_field_32_type), intent(inout) :: this

    if (associated(this%test_data)) then
      deallocate( this%test_data )
    end if

  end subroutine destroy_32


  subroutine destroy_64( this )

    implicit none

    type(test_field_64_type), intent(inout) :: this

    if (associated(this%test_data)) then
      deallocate( this%test_data )
    end if

  end subroutine destroy_64


  function get_data_32( this )

    implicit none

    class(test_field_32_type), intent(in) :: this
    real(real32), pointer :: get_data_32(:)

    get_data_32 => this%test_data

  end function get_data_32


  function get_data_64( this )

    implicit none

    class(test_field_64_type), intent(in) :: this
    real(real64), pointer :: get_data_64(:)

    get_data_64 => this%test_data

  end function get_data_64


  subroutine copy_32_from_lfric(self, return_code)

    implicit none

    class(test_field_32_type), intent(inout) :: self
    integer(i_def),  optional, intent(out)   :: return_code

    type(field_real32_type), pointer :: lfric_field
    type(field_real32_proxy_type)    :: lfric_proxy

    select type(lfric_field => self%get_lfric_field_ptr())
    class is (field_real32_type)
      lfric_proxy = lfric_field%get_proxy()
      self%test_data = lfric_proxy%data
    class default
      call log_event("Inconsistent test_field_32_type", log_level_error )
    end select

    if (present(return_code)) then
      return_code = 0
    end if

  end subroutine copy_32_from_lfric


  subroutine copy_64_from_lfric(self, return_code)

    implicit none

    class(test_field_64_type), intent(inout) :: self
    integer(i_def),  optional, intent(out)   :: return_code

    type(field_real64_type), pointer :: lfric_field
    type(field_real64_proxy_type)    :: lfric_proxy

    select type(lfric_field => self%get_lfric_field_ptr())
    class is (field_real64_type)
      lfric_proxy = lfric_field%get_proxy()
      self%test_data = lfric_proxy%data
    class default
      call log_event("Inconsistent test_field_64_type", log_level_error )
    end select

    if (present(return_code)) then
      return_code = 0
    end if

  end subroutine copy_64_from_lfric


  subroutine copy_32_to_lfric(self, return_code)

    implicit none

    class(test_field_32_type), intent(inout) :: self
    integer(i_def),  optional, intent(out)   :: return_code

    type(field_real32_type), pointer :: lfric_field
    type(field_real32_proxy_type)    :: lfric_proxy

    select type(lfric_field => self%get_lfric_field_ptr())
    class is (field_real32_type)
      lfric_proxy = lfric_field%get_proxy()
      lfric_proxy%data = self%test_data
    class default
      call log_event("Inconsistent test_field_32_type", log_level_error )
    end select

    if (present(return_code)) then
      return_code = 0
    end if

  end subroutine copy_32_to_lfric


  subroutine copy_64_to_lfric(self, return_code)

    implicit none

    class(test_field_64_type), intent(inout) :: self
    integer(i_def),  optional, intent(out)   :: return_code

    type(field_real64_type), pointer :: lfric_field
    type(field_real64_proxy_type)    :: lfric_proxy

    select type(lfric_field => self%get_lfric_field_ptr())
    class is (field_real64_type)
      lfric_proxy = lfric_field%get_proxy()
      lfric_proxy%data = self%test_data
    class default
      call log_event("Inconsistent test_field_32_type", log_level_error )
    end select

    if (present(return_code)) then
      return_code = 0
    end if

  end subroutine copy_64_to_lfric

end module test_field_mod
