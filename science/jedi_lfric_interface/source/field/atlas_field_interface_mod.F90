!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing an Atlas field class.
!>
!> @details The atlas field interface can hold a 64-bit real pointer to data
!> that is stored externaly. The field can be copied to and from LFRic fields
!> stored in the lowest order W3 and Wtheta function-space including 2D fields.
!> Methods are provided to perform the copies that rely on access to the field
!> proxy. Only owned LFRic points are updated.
!
module atlas_field_interface_mod

  use, intrinsic :: iso_fortran_env, only : real64

  use log_mod,                       only : log_event,          &
                                            log_scratch_space,  &
                                            LOG_LEVEL_ERROR
  use abstract_external_field_mod,   only : abstract_external_field_type
  use field_mod,                     only : field_type, field_proxy_type
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

!> @param [in] atlas_data_ptr pointer to the atlas data
!> @param [in] map_horizontal_ptr pointer to the horizontal map
!> @param [in] lfric_field_ptr the LFRic field that atlas_field_interface will
!>             copy to and from
!> @param [in] surface_level_type enumerator defining how the lowest level
!>             is treated.
!> @param [in] atlas_name Optional name of the atlas field
!> @param [in] fill_direction_up Optional logical set false if the atlas
!>               field data is orientated from top-bottom
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
  logical( kind=l_def )              :: fill_direction_up_local

  ! Initialise the abstract parent
  call self%abstract_external_field_initialiser( lfric_field_ptr )

  ! Mandated inputs
  self%atlas_data => atlas_data_ptr
  self%map_horizontal => map_horizontal_ptr

  ! Optionals
  if ( present(atlas_name) ) then
    self%atlas_name = atlas_name
  else
    self%atlas_name = lfric_field_ptr%get_name()
  endif

  if ( present(surface_level_type) ) then
    ! Ensure the supplied option is valid
    if ( surface_level_type /= surface_level_present .AND. &
         surface_level_type /= surface_level_absent_copy_level_above .AND. &
         surface_level_type /= surface_level_absent_zero_level ) then
      write(log_scratch_space, '(A)') &
        "The supplied optional argument surface_level_type is not valid."
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    endif
    self%surface_level_type = surface_level_type
  else
    self%surface_level_type = surface_level_present
  endif

  if (present(fill_direction_up)) then
    fill_direction_up_local = fill_direction_up
  else
    fill_direction_up_local = .true.
  endif

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
  endif

  if ( fill_direction_up_local ) then
    self%atlas_kstart = 1
    self%atlas_kend = self%n_vertical
    self%atlas_kdirection = 1
  else
    self%atlas_kstart = self%n_vertical
    self%atlas_kend = 1
    self%atlas_kdirection = -1
  endif

  ! Check inputs are consistent

  ! 0. Check input field is supported and get associated sizes

  ! get the number of points in the LFRic field
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
  endif

  ! 2. Horizotal points
  if ( n_horizontal_lfric /= self%n_horizontal ) then
    write(log_scratch_space, '(A,I0,A,I0,A,I0)') &
      "Field mismatch in atlas field interface constructor for the number of horizontal points. The atlas field has ", &
      self%n_horizontal, " points and the LFRic field has ", n_horizontal_lfric, " points."
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  endif

  ! 3. Check the input horizontal map

  ! Check map size
  if ( self%n_horizontal /= size(self%map_horizontal) ) then
    write(log_scratch_space, '(A,I0,A,I0)') &
      "The data and map should be the same size: n_horizontal = ", &
      self%n_horizontal , " map_horizontal size = ", &
      size(self%map_horizontal)
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  endif

  ! Check the input contains only unique points and that they are
  ! within the bounds of the field size
  allocate(check_map(self%n_horizontal))
  check_map = .false.
  do ij=1,self%n_horizontal
    ! bounds check on input
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

!> Copy Atlas field from the LFRic field
subroutine copy_from_lfric(self, return_code)

  implicit none

  class( atlas_field_interface_type ), intent(inout) :: self
  integer(i_def), intent(out), optional :: return_code

  type( field_type ), pointer :: lfric_field_ptr
  type( field_proxy_type )    :: field_proxy
  integer(i_def)              :: lfric_kstart
  integer(i_def)              :: lfric_ij
  integer(i_def)              :: atlas_ij
  integer(i_def)              :: atlas_kstart
  integer(i_def)              :: atlas_kend
  integer(i_def)              :: atlas_kdirection
  integer(i_def)              :: n_vertical_lfric

  lfric_field_ptr => self%get_lfric_field_ptr()
  field_proxy = lfric_field_ptr%get_proxy()

  ! get indices for atlas and LFRic data
  atlas_kstart = self % atlas_kstart
  atlas_kend = self % atlas_kend
  atlas_kdirection = self % atlas_kdirection

  lfric_kstart = self%lfric_kstart
  n_vertical_lfric = self%n_vertical_lfric

  ! copy LFRic to Atlas
  do lfric_ij = 1,self%n_horizontal
    atlas_ij = self%map_horizontal(lfric_ij)
    self%atlas_data(atlas_kstart:atlas_kend:atlas_kdirection,atlas_ij) &
      = field_proxy%data((lfric_ij-1)*n_vertical_lfric + lfric_kstart : lfric_ij*n_vertical_lfric)
  enddo

end subroutine copy_from_lfric

!> Copy Atlas field to the LFRic field
subroutine copy_to_lfric( self, return_code )

  implicit none

  class( atlas_field_interface_type ), intent(inout) :: self
  integer(i_def),              intent(out), optional :: return_code

  type( field_type ), pointer :: lfric_field_ptr
  type( field_proxy_type )    :: field_proxy
  integer(i_def)              :: lfric_kstart
  integer(i_def)              :: lfric_ij
  integer(i_def)              :: atlas_ij
  integer(i_def)              :: atlas_kstart
  integer(i_def)              :: atlas_kend
  integer(i_def)              :: atlas_kdirection
  integer(i_def)              :: n_vertical_lfric

  lfric_field_ptr => self%get_lfric_field_ptr()
  field_proxy = lfric_field_ptr%get_proxy()

  ! get indices for atlas and LFRic data
  lfric_kstart = self%lfric_kstart
  n_vertical_lfric = self%n_vertical_lfric

  atlas_kstart = self % atlas_kstart
  atlas_kend = self % atlas_kend
  atlas_kdirection = self % atlas_kdirection

  ! copy Atlas to lfric

  ! Fill data ommiting the lowest level if its not available
  do lfric_ij = 1,self%n_horizontal
    atlas_ij = self%map_horizontal(lfric_ij)
    field_proxy%data((lfric_ij-1)*n_vertical_lfric+lfric_kstart:lfric_ij*n_vertical_lfric) &
      = self%atlas_data(atlas_kstart:atlas_kend:atlas_kdirection,atlas_ij)
  enddo

  ! Fill missing data if required
  if ( self%surface_level_type == surface_level_absent_zero_level ) then
    ! Set to zero
    do lfric_ij = 1,self%n_horizontal
      field_proxy%data((lfric_ij-1)*n_vertical_lfric+1) = 0.0!!! set _rdef?
    end do
  elseif ( self%surface_level_type == surface_level_absent_copy_level_above ) then
    ! Copy the level above
    do lfric_ij = 1,self%n_horizontal
      field_proxy%data((lfric_ij-1)*n_vertical_lfric+1) = &
                      field_proxy%data((lfric_ij-1)*n_vertical_lfric+2)
    end do
  endif

  ! set halo to dirty
  call field_proxy%set_dirty()

end subroutine copy_to_lfric

!> Returns the lfric name of the field
function get_lfric_name( self ) result( lfric_name )

  implicit none

  class( atlas_field_interface_type ), intent(in) :: self
  type( field_type ), pointer                     :: lfric_field_ptr
  character( len=str_def )                        :: lfric_name

  lfric_field_ptr => self%get_lfric_field_ptr()
  lfric_name = lfric_field_ptr%get_name()

end function get_lfric_name

!> Returns the atlas name of the field
function get_atlas_name( self ) result( atlas_name )

  implicit none

  class( atlas_field_interface_type ), intent(in) :: self
  character( len=str_def )                        :: atlas_name

  atlas_name = self%atlas_name

end function get_atlas_name

!> Finalizer
subroutine atlas_field_interface_destructor(self)!

  implicit none

  type(atlas_field_interface_type), intent(inout)    :: self

  ! Destroy the abstract parent
  call self%abstract_external_field_destructor()
  self%atlas_data => null()
  self%map_horizontal => null()

end subroutine atlas_field_interface_destructor

end module atlas_field_interface_mod
