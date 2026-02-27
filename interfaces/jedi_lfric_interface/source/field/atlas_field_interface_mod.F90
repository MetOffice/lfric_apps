!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing an Atlas-LFRic interface field class.
!>
!> @details The atlas field interface can hold a 64-bit real pointer to data
!> that is stored externally. The field can be copied to and from LFRic fields
!> stored in the lowest order W3 and Wtheta function-space including 2D fields.
!> Methods are provided to perform the copies that rely on access to the field
!> proxy. Only owned LFRic and Atlas points are updated.
!>
!> Adjoint of the data copy methods are also included. Noting that "=" is
!> ambiguous because it may represent a move or copy, it is assumed that a copy
!> is intended and add a rule: the adjoint of "falling out of scope" is
!> initialisation to zero. In the below adjoint code, where "=" is encountered
!> in the forward code:
!>
!> b=a
!>
!> we apply "+=" followed by initialisation to zero in the adjoint:
!>
!> a_hat=a_hat+b_hat
!> b_hat=0
!>
!> This results in non intuitive naming of functions. The adjoint of copying
!> "from" LFRic is an update "to" the LFRic field. The adjoint code below is
!> a standard line by line adjoint where the outlined rules are adopted.
!>
!> The convention noted above states that we need to include an additional
!> rule: The adjoint of "falling out of scope" is "initialisation to zero".
!> This means that in addition to the adjoint data copies, we also need to
!> include a method to zero the fields. To enable that, two further methods
!> (zero_lfric and zero_atlas) are included.
!>
module atlas_field_interface_mod

  use, intrinsic :: iso_fortran_env, only : real64

  use log_mod,                       only : log_event,          &
                                            log_scratch_space,  &
                                            LOG_LEVEL_ERROR
  use abstract_external_field_mod,   only : abstract_external_field_type
  use field_mod,                     only : field_type, field_proxy_type
  use field_parent_mod,              only : field_parent_type
  use fs_continuity_mod,             only : W3, Wtheta, name_from_functionspace
  use constants_mod,                 only : i_def, str_def, l_def

  implicit none

  private

  integer(i_def), parameter, public :: surface_level_present = 123
  integer(i_def), parameter, public :: surface_level_absent_copy_level_above = 456
  integer(i_def), parameter, public :: surface_level_absent_zero_level = 789

type, extends(abstract_external_field_type), public :: atlas_field_interface_type
  private

  !> The 64-bit floating point values of the field
  real( kind=real64 ), pointer :: atlas_data(:,:)
  !> Map that defines the order of the horizontal points in
  !> the atlas data to collumns in the LFRic field data
  integer(i_def),      pointer :: map_horizontal(:)
  !> The name of the atlas field
  character( len=str_def )     :: atlas_name
  !> Enumerator defining how the surface level is treated. There are three options
  !> present (no processing), absent (set to zero), absent (copy level above)
  integer(i_def)               :: surface_level_type
  !> Number of vertical points in the atlas data
  integer(i_def)               :: n_vertical
  !> Number of horizontal points in the atlas data
  integer(i_def)               :: n_horizontal
  !> Vertical start index for the LFRic data
  integer(i_def)               :: lfric_kstart
  !> Number of vertical points in the LFRic data
  integer(i_def)               :: n_vertical_lfric
  !> Vertical start index for the atlas data
  integer(i_def)               :: atlas_kstart
  !> Vertical end index for the atlas data
  integer(i_def)               :: atlas_kend
  !> Vertical fill direction for the atlas data
  integer(i_def)               :: atlas_kdirection

contains

  !> Field initialiser.
  procedure, public :: initialise => field_initialiser

  !> Copy atlas field from the LFRic field
  procedure, public :: copy_from_lfric

  !> Copy atlas field to the LFRic field
  procedure, public :: copy_to_lfric

  !> Adjoint of: (Copy atlas field from the LFRic field)
  procedure, public :: copy_from_lfric_ad

  !> Adjoint of: (Copy atlas field to the LFRic field)
  procedure, public :: copy_to_lfric_ad

  !> Zero the LFRic field
  procedure, public :: zero_lfric

  !> Zero the Atlas field
  procedure, public :: zero_atlas

  !> Get the LFRic field name
  procedure, public :: get_lfric_name

  !> Get the atlas field name
  procedure, public :: get_atlas_name

  !> Finalizer
  final             :: atlas_field_interface_destructor

end type atlas_field_interface_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!> @brief atlas_field_interface_type initialiser
!>
!> @param [in] atlas_data_ptr     Pointer to the atlas data
!> @param [in] map_horizontal_ptr Pointer to the horizontal map
!> @param [in] lfric_field_ptr    The LFRic field that atlas_field_interface
!>                                will copy to and from
!> @param [in] surface_level_type Enumerator defining how the lowest level
!>                                is treated.
!> @param [in] atlas_name         Optional name of the atlas field
!> @param [in] fill_direction_up  Optional logical set false if the atlas field
!>                                data is orientated from top-bottom
subroutine field_initialiser( self, atlas_data_ptr, map_horizontal_ptr, &
                              lfric_field_ptr, surface_level_type, &
                              atlas_name, fill_direction_up )

  implicit none

  class( atlas_field_interface_type ), intent(inout) :: self
  real( kind=real64 ),   pointer,  intent(in) :: atlas_data_ptr(:,:)
  integer(i_def),        pointer,  intent(in) :: map_horizontal_ptr(:)
  type(field_type),      pointer,  intent(in) :: lfric_field_ptr
  integer(i_def),        optional, intent(in) :: surface_level_type
  character( len=* ),    optional, intent(in) :: atlas_name
  logical( kind=l_def ), optional, intent(in) :: fill_direction_up

  ! locals
  logical( kind=l_def ), allocatable :: check_map(:)
  integer(i_def)                     :: ij
  integer(i_def)                     :: fs_enumerator
  integer(i_def)                     :: n_horizontal_lfric
  integer(i_def)                     :: n_vertical_lfric
  type( field_proxy_type )           :: field_proxy
  class(field_parent_type), pointer  :: cast_field
  logical( kind=l_def )              :: fill_direction_up_local

  ! Initialise the abstract parent
  !
  ! This little dance with pointers is needed to keep some compilers
  ! happy that you can, in fact, pass a child class to an argument
  ! expecting a parent class.
  !
  cast_field => lfric_field_ptr
  call self%abstract_external_field_initialiser( cast_field )

  ! Mandated inputs
  self%atlas_data => atlas_data_ptr
  self%map_horizontal => map_horizontal_ptr

  ! Optionals
  if ( present(atlas_name) ) then
    self%atlas_name = atlas_name
  else
    self%atlas_name = lfric_field_ptr%get_name()
  end if

  if ( present(surface_level_type) ) then
    ! Ensure the supplied option is valid
    if ( surface_level_type /= surface_level_present .AND. &
         surface_level_type /= surface_level_absent_copy_level_above .AND. &
         surface_level_type /= surface_level_absent_zero_level ) then
      write(log_scratch_space, '(A)') &
        "The supplied optional argument surface_level_type is not valid."
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    self%surface_level_type = surface_level_type
  else
    self%surface_level_type = surface_level_present
  end if

  if (present(fill_direction_up)) then
    fill_direction_up_local = fill_direction_up
  else
    fill_direction_up_local = .true.
  end if

  ! Setup data sizes and indices
  self%n_horizontal = size(atlas_data_ptr,dim=2)
  self%n_vertical = size(atlas_data_ptr,dim=1)

  ! Setup vertical indices for data mapping
  if ( self%surface_level_type == surface_level_absent_copy_level_above .OR. &
       self%surface_level_type == surface_level_absent_zero_level ) then
    self%lfric_kstart = 2
    self%n_vertical_lfric = self%n_vertical+1
  else
    self%lfric_kstart = 1
    self%n_vertical_lfric = self%n_vertical
  end if

  if ( fill_direction_up_local ) then
    self%atlas_kstart = 1
    self%atlas_kend = self%n_vertical
    self%atlas_kdirection = 1
  else
    self%atlas_kstart = self%n_vertical
    self%atlas_kend = 1
    self%atlas_kdirection = -1
  end if

  ! Check inputs are consistent

  ! 0. Check input field is supported and get associated sizes

  ! Get the number of points in the LFRic field
  fs_enumerator = lfric_field_ptr%which_function_space()
  field_proxy=lfric_field_ptr%get_proxy()
  select case ( fs_enumerator )
    case (Wtheta)
      n_vertical_lfric = field_proxy%vspace%get_nlayers() + 1
    case ( W3 )
      n_vertical_lfric = field_proxy%vspace%get_nlayers()
    case default
      write(log_scratch_space, '(3A)') "The ", &
        trim(name_from_functionspace(fs_enumerator)) , &
        " function space is not supported by the atlas_field_interface class."
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end select
  n_horizontal_lfric = field_proxy%vspace%get_last_dof_owned()/n_vertical_lfric

  ! 1. Vertical points
  if ( n_vertical_lfric /= self%n_vertical_lfric ) then
    write(log_scratch_space, '(A,I0,A,I0,A)') &
      "Field mismatch in atlas field interface constructor for the number of vertical points. The atlas field has ", &
      self%n_vertical_lfric, " points and the LFRic field has ", n_vertical_lfric, " points."
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

  ! 2. Horizotal points
  if ( n_horizontal_lfric /= self%n_horizontal ) then
    write(log_scratch_space, '(A,I0,A,I0,A,I0)') &
      "Field mismatch in atlas field interface constructor for the number of horizontal points. The atlas field has ", &
      self%n_horizontal, " points and the LFRic field has ", n_horizontal_lfric, " points."
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

  ! 3. Check the input horizontal map

  ! Check map size
  if ( self%n_horizontal /= size(self%map_horizontal) ) then
    write(log_scratch_space, '(A,I0,A,I0)') &
      "The data and map should be the same size: n_horizontal = ", &
      self%n_horizontal , " map_horizontal size = ", &
      size(self%map_horizontal)
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

  ! Check the input contains only unique points and that they are
  ! within the bounds of the field size
  allocate(check_map(self%n_horizontal))
  check_map = .false.
  do ij=1,self%n_horizontal
    ! Bounds check on input
    if ( self%map_horizontal(ij) > 0 .and. &
          self%map_horizontal(ij) <= self%n_horizontal ) then
      if ( check_map(self%map_horizontal(ij)) ) then
        write(log_scratch_space, '(A,I0,A)') &
          "The ", ij, " point in the map_horizontal has already been found."
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )

      else
        check_map(self%map_horizontal(ij)) = .true.
      end if
    else
      write(log_scratch_space, '(A,I0,A,I0,A,I0)') "The ", ij, &
        " point in the map_horizontal is out of bounds. Value provided is: ", &
        self%map_horizontal(ij), &
        ", which is outside the allowable bounds: 1-", self%n_horizontal
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
  end do

end subroutine field_initialiser

!> @brief Copy Atlas field from the LFRic field
!>
!> @param [out] return_code Optional code that is allways set to zero if
!>                          present
subroutine copy_from_lfric(self, return_code)

  implicit none

  class( atlas_field_interface_type ), intent(inout) :: self
  integer(i_def),              optional, intent(out) :: return_code

  type( field_proxy_type )    :: field_proxy
  integer(i_def)              :: lfric_kstart
  integer(i_def)              :: lfric_ij
  integer(i_def)              :: atlas_ij
  integer(i_def)              :: ij
  integer(i_def)              :: atlas_kstart
  integer(i_def)              :: atlas_kend
  integer(i_def)              :: atlas_kdirection
  integer(i_def)              :: n_vertical_lfric

  select type(field_ptr => self%get_lfric_field_ptr())
  class is (field_type)
    field_proxy = field_ptr%get_proxy()
  class default
    call log_event(                                                          &
      "Unexpected field type in atlas_field_interface_type%copy_from_lfric", &
      log_level_error                                                        &
    )
  end select  

  ! Get indices for atlas and LFRic data
  atlas_kstart = self % atlas_kstart
  atlas_kend = self % atlas_kend
  atlas_kdirection = self % atlas_kdirection

  lfric_kstart = self%lfric_kstart
  n_vertical_lfric = self%n_vertical_lfric

  ! Copy LFRic to Atlas
  do ij = 1,self%n_horizontal
    ! Get LFRic and Atlas index
    lfric_ij = (ij-1)*n_vertical_lfric
    atlas_ij = self%map_horizontal(ij)
    self%atlas_data(atlas_kstart:atlas_kend:atlas_kdirection,atlas_ij) &
      = field_proxy%data(lfric_ij+lfric_kstart:lfric_ij+n_vertical_lfric)
  end do

  if ( present(return_code) ) return_code = 0_i_def

end subroutine copy_from_lfric

!> @brief Copy Atlas field to the LFRic field
!>
!> @param [out] return_code Optional code that is allways set to zero if
!>                          present
subroutine copy_to_lfric( self, return_code )

  implicit none

  class( atlas_field_interface_type ), intent(inout) :: self
  integer(i_def),              intent(out), optional :: return_code

  type( field_proxy_type )    :: field_proxy
  integer(i_def)              :: lfric_kstart
  integer(i_def)              :: ij
  integer(i_def)              :: lfric_ij
  integer(i_def)              :: atlas_ij
  integer(i_def)              :: atlas_kstart
  integer(i_def)              :: atlas_kend
  integer(i_def)              :: atlas_kdirection
  integer(i_def)              :: n_vertical_lfric

  select type (field_ptr => self%get_lfric_field_ptr())
  class is (field_type)
    field_proxy = field_ptr%get_proxy()
  class default
    call log_event(                                                        &
      "Unexpected field type in atlas_field_interface_type%copy_to_lfric", &
      log_level_error                                                      &
    )
  end select  

  ! Get indices for atlas and LFRic data
  lfric_kstart = self%lfric_kstart
  n_vertical_lfric = self%n_vertical_lfric

  atlas_kstart = self % atlas_kstart
  atlas_kend = self % atlas_kend
  atlas_kdirection = self % atlas_kdirection

  ! Copy Atlas to lfric

  ! Fill data ommiting the lowest level if its not available
  do ij = 1,self%n_horizontal
    atlas_ij = self%map_horizontal(ij)
    lfric_ij = (ij-1)*n_vertical_lfric
    field_proxy%data(lfric_ij+lfric_kstart:lfric_ij+n_vertical_lfric) &
      = self%atlas_data(atlas_kstart:atlas_kend:atlas_kdirection,atlas_ij)
  end do

  ! Fill missing data if required
  if ( self%surface_level_type == surface_level_absent_zero_level ) then
    ! Set to zero
    do ij = 1,self%n_horizontal
      lfric_ij = (ij-1)*n_vertical_lfric
      field_proxy%data(lfric_ij+1) = 0.0_real64
    end do
  elseif ( self%surface_level_type == surface_level_absent_copy_level_above ) then
    ! Copy the level above
    do ij = 1,self%n_horizontal
      atlas_ij = self%map_horizontal(ij)
      lfric_ij = (ij-1)*n_vertical_lfric
      field_proxy%data(lfric_ij+1) &
          = self%atlas_data(atlas_kstart,atlas_ij)
    end do
  end if

  ! Set halo to dirty
  call field_proxy%set_dirty()

  if ( present(return_code) ) return_code = 0_i_def

end subroutine copy_to_lfric

!> @brief   Adjoint of copy_from_lfric
!>
!> @detail  This is a line by line adjoint of the copy_from_lfric method
!>          and as noted in the module detail results in data copies to the
!>          LFRic field and initialization of the Atlas field to zero.
subroutine copy_from_lfric_ad(self)

  implicit none

  class( atlas_field_interface_type ), intent(inout) :: self

  type( field_proxy_type )    :: field_proxy
  integer(i_def)              :: lfric_kstart
  integer(i_def)              :: ij
  integer(i_def)              :: lfric_ij
  integer(i_def)              :: atlas_ij
  integer(i_def)              :: atlas_kstart
  integer(i_def)              :: atlas_kend
  integer(i_def)              :: atlas_kdirection
  integer(i_def)              :: n_vertical_lfric

  select type (field_ptr => self%get_lfric_field_ptr())
  class is (field_type)
    field_proxy = field_ptr%get_proxy()
  class default
    call log_event(                                                             &
      "Unexpected field type in atlas_field_interface_type%copy_from_lfric_ad", &
      log_level_error                                                           &
    )
  end select  

  ! Get indices for atlas and LFRic data
  atlas_kstart = self % atlas_kstart
  atlas_kend = self % atlas_kend
  atlas_kdirection = self % atlas_kdirection

  lfric_kstart = self%lfric_kstart
  n_vertical_lfric = self%n_vertical_lfric

  ! Adjoint of copy LFRic to Atlas
  do ij = self%n_horizontal,1,-1
    atlas_ij = self%map_horizontal(ij)
    lfric_ij = (ij-1)*n_vertical_lfric
    field_proxy%data(lfric_ij+lfric_kstart:lfric_ij+n_vertical_lfric) &
        = field_proxy%data(lfric_ij+lfric_kstart:lfric_ij+n_vertical_lfric) &
        + self%atlas_data(atlas_kstart:atlas_kend:atlas_kdirection,atlas_ij)
  end do

  ! Initialise out of scope variable to zero
  self%atlas_data(:,:) = 0.0_real64

end subroutine copy_from_lfric_ad

!> @brief   Adjoint of copy_to_lfric
!>
!> @detail  This is a line by line adjoint of the copy_to_lfric method
!>          and as noted in the module detail results in data copies to the
!>          Atlas field and initialization of the LFRic field to zero.
!>
subroutine copy_to_lfric_ad( self )

  implicit none

  class( atlas_field_interface_type ), intent(inout) :: self

  type( field_proxy_type )    :: field_proxy
  integer(i_def)              :: lfric_kstart
  integer(i_def)              :: ij
  integer(i_def)              :: lfric_ij
  integer(i_def)              :: atlas_ij
  integer(i_def)              :: atlas_kstart
  integer(i_def)              :: atlas_kend
  integer(i_def)              :: atlas_kdirection
  integer(i_def)              :: n_vertical_lfric

  select type (field_ptr => self%get_lfric_field_ptr())
  class is (field_type)
    field_proxy = field_ptr%get_proxy()
  class default
    call log_event(                                                           &
      "Unexpected field type in atlas_field_interface_type%copy_to_lfric_ad", &
      log_level_error                                                         &
    )
  end select  

  ! Get indices for atlas and LFRic data
  lfric_kstart = self%lfric_kstart
  n_vertical_lfric = self%n_vertical_lfric

  atlas_kstart = self % atlas_kstart
  atlas_kend = self % atlas_kend
  atlas_kdirection = self % atlas_kdirection

  ! Adjoint of (copy Atlas to lfric)

  ! Adjoint of (fill missing data if required)

  ! Note: The adjoint of (set to zero) is omitted for:
  !       self%surface_level_type == surface_level_absent_zero_level. That is
  !       because the operation is redundant, in this case the adjoint is:
  !       self%atlas_data(:) = self%atlas_data(:) + 0.0.
  !
  if ( self%surface_level_type == surface_level_absent_copy_level_above ) then
    ! Adjoint of (copy the level above)
    do ij = self%n_horizontal,1,-1
      atlas_ij = self%map_horizontal(ij)
      lfric_ij = (ij-1)*n_vertical_lfric
      self%atlas_data(atlas_kstart,atlas_ij) &
                    = self%atlas_data(atlas_kstart,atlas_ij) &
                    + field_proxy%data(lfric_ij+1)
    end do
  end if

  ! Adjoint of (fill data ommiting the lowest level if its not available)
  do ij = self%n_horizontal,1,-1
    atlas_ij = self%map_horizontal(ij)
    lfric_ij = (ij-1)*n_vertical_lfric
    self%atlas_data(atlas_kstart:atlas_kend:atlas_kdirection,atlas_ij) &
        = self%atlas_data(atlas_kstart:atlas_kend:atlas_kdirection,atlas_ij) &
        + field_proxy%data(lfric_ij+lfric_kstart:lfric_ij+n_vertical_lfric)
  end do

  ! Initialise out of scope variable to zero
  field_proxy%data(1:self%n_vertical_lfric*self%n_horizontal) = 0.0_real64

  ! Set halo to dirty
  call field_proxy%set_dirty()

end subroutine copy_to_lfric_ad

!> @brief Zero the LFRic field
!>
subroutine zero_lfric( self )

  implicit none

  class( atlas_field_interface_type ), intent(inout) :: self

  type( field_proxy_type )    :: field_proxy

  select type (field_ptr => self%get_lfric_field_ptr())
  class is (field_type)
    field_proxy = field_ptr%get_proxy()

    ! Set all data to zero
    field_proxy%data(:) = 0.0_real64
  class default
    call log_event(                                                     &
      "Unexpected field type in atlas_field_interface_type%zero_lfric", &
      log_level_error                                                   &
    )
  end select  

end subroutine zero_lfric

!> @brief Zero the Atlas field
!>
subroutine zero_atlas( self )

  implicit none

  class( atlas_field_interface_type ), intent(inout) :: self

  ! Set owned data to zero
  self%atlas_data(:,:) = 0.0_real64

end subroutine zero_atlas

!> @brief Returns the lfric name of the field
!>
function get_lfric_name( self ) result( lfric_name )

  implicit none

  class( atlas_field_interface_type ), intent(in) :: self
  character( len=str_def )                        :: lfric_name

  select type(field_ptr => self%get_lfric_field_ptr())
  class is (field_type)
    lfric_name = field_ptr%get_name()
  class default
    call log_event(                                                         &
      "Unexpected field type in atlas_field_interface_type%get_lfric_name", &
      log_level_error                                                       &
    )
  end select

end function get_lfric_name

!> @brief Returns the atlas name of the field
!>
function get_atlas_name( self ) result( atlas_name )

  implicit none

  class( atlas_field_interface_type ), intent(in) :: self
  character( len=str_def )                        :: atlas_name

  atlas_name = self%atlas_name

end function get_atlas_name

!> @brief atlas_field_interface_type finalizer
!>
subroutine atlas_field_interface_destructor(self)!

  implicit none

  type(atlas_field_interface_type), intent(inout)    :: self

  ! Destroy the abstract parent
  call self%abstract_external_field_destructor()
  self%atlas_data => null()
  self%map_horizontal => null()

end subroutine atlas_field_interface_destructor

end module atlas_field_interface_mod
