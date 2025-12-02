!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Infrastructure for selecting fields for use in create_physics_prognostics
!
module field_spec_mod

  use constants_mod, only : i_def, l_def, str_def, imdi
  use log_mod,       only : log_event, log_scratch_space, log_level_error, log_level_info
  use clock_mod,     only : clock_type

  implicit none

  !> @brief Dictionary of main field collections
  type :: main_coll_dict_type
    integer(i_def) :: derived
    integer(i_def) :: radiation
    integer(i_def) :: microphysics
    integer(i_def) :: electric
    integer(i_def) :: orography
    integer(i_def) :: turbulence
    integer(i_def) :: convection
    integer(i_def) :: cloud
    integer(i_def) :: surface
    integer(i_def) :: soil
    integer(i_def) :: snow
    integer(i_def) :: chemistry
    integer(i_def) :: aerosol
    integer(i_def) :: stph
    integer(i_def) :: lbc
    integer(i_def) :: none ! for non-physics fields, e.g., gungho
  contains
    procedure :: check => main_coll_check
  end type main_coll_dict_type

  !> @brief Map main collection enumerators to collections.
  type(main_coll_dict_type), parameter :: main_coll_dict &
    = main_coll_dict_type( &
    107,120,129,130,141,142,152,154,160,190,191,222,238,241,252,279)

  !> @brief Dictionary of advected field collections
  type :: adv_coll_dict_type
    integer(i_def) :: none       ! Not advected
    integer(i_def) :: all_adv    ! Adv_fields_all_outer
    integer(i_def) :: last_adv   ! Adv_fields_last_outer
    integer(i_def) :: all_con    ! Con_fields_all_outer
    integer(i_def) :: last_con   ! Con_fields_last_outer
  end type adv_coll_dict_type

  !> @brief Map advected field enumerators to collections.
  type(adv_coll_dict_type), parameter :: adv_coll_dict &
    = adv_coll_dict_type(387,391,395,399,412)

  !> @brief Dictionary of moisture field arrays
  type :: moist_arr_dict_type
    integer(i_def) :: none          ! field not part of a moisture array
    integer(i_def) :: mr            ! mr array
    integer(i_def) :: moist_dyn     ! moist_dyn array
    integer(i_def) :: moist_dyn_ref ! moist_dyn_ref array
  end type moist_arr_dict_type

   !> @brief Map moisture array enumerators to moisture array.
  type(moist_arr_dict_type), parameter :: moist_arr_dict &
    = moist_arr_dict_type(445,450,454,461)

  !> @brief Dictionary of time axes
    type :: time_axis_dict_type
    integer(i_def) :: none          ! field without time axis
    integer(i_def) :: lbc           ! lbc time axis
    integer(i_def) :: ls            ! ls time axis
  end type time_axis_dict_type

   !> @brief Map moisture array enumerators to moisture array.
  type(time_axis_dict_type), parameter :: time_axis_dict &
    = time_axis_dict_type(525,529,536)

  ! request function space discovery
  integer, parameter :: missing_fs = imdi

  !> @brief Metadata needed to construct a model field
  !> @details The members of the type are initialised to default values
  type :: field_spec_type
    character(str_def) :: name      = ''                  ! Field name
    integer(i_def)     :: main_coll = imdi                ! Enumerator of main field collection
    integer(i_def)     :: space     = missing_fs          ! Function space enumerator
    integer(i_def)     :: order_h   = 0                   ! Function space order in horizontal
    integer(i_def)     :: order_v   = 0                   ! Function space order in vertical
    integer(i_def)     :: adv_coll  = adv_coll_dict%none  ! Enumerator of advected field collection
    integer(i_def)     :: moist_arr = moist_arr_dict%none ! Enumerator of moisture array
    integer(i_def)     :: moist_idx = 0                   ! Index into moisture array
    integer(i_def)     :: time_axis = time_axis_dict%none ! Enumerator of time axis
    character(str_def) :: mult      = ''                  ! Name of multidata item, or blank string
    logical(l_def)     :: ckp       = .false.             ! Is it a checkpoint (prognostic) field?
    logical(l_def)     :: twod      = .false.             ! Is it two-dimensional?
    logical(l_def)     :: empty     = .false.             ! Is it empty (with an empty data array)?
    logical(l_def)     :: coarse    = .false.             ! Is it coarse?
    logical(l_def)     :: is_int    = .false.             ! Is it an integer field?
    logical(l_def)     :: legacy    = .false.             ! Is it a field with legacy checkpointing?
  end type field_spec_type

  private
  public :: field_spec_type, &
            main_coll_dict_type, main_coll_dict, &
            adv_coll_dict_type, adv_coll_dict, &
            moist_arr_dict, time_axis_dict, &
            processor_type, make_spec, if_advected, missing_fs, &
            space_has_xios_io

  !> @brief Base class for processor objects, operating on field specifiers
  type, abstract :: processor_type
    class(clock_type), pointer :: clock
  contains
    private
    ! accessors
    procedure, public :: get_clock => processor_get_clock
    procedure, public :: set_clock => processor_set_clock

    ! main interface
    procedure(apply_interface), public, deferred :: apply
  end type processor_type

  abstract interface
    !> @brief Apply a processor object to a field specifier.
    !> @param[in] self       Processor object
    !> @param[in] spec       Field specifier
    subroutine apply_interface(self, spec)
      import processor_type
      import field_spec_type
      class(processor_type), intent(in) :: self
      type(field_spec_type), intent(in) :: spec
    end subroutine apply_interface
  end interface

contains

  !> @brief Generic error handler for enumerator check
  !> @param[in] which              String identifying enumerator
  !> @param[in] enum               Faulty enumerator value
  subroutine enum_error(which, enum)
    implicit none
    character(*), intent(in) :: which
    integer(i_def), intent(in) :: enum
    write(log_scratch_space, &
      '("unexpected ", A, " enumerator: ", I3)') trim(which), enum
    call log_event(log_scratch_space, log_level_error)
  end subroutine enum_error

  !> @brief Check main collection enumerator
  !> @param[in] s         Main collection dictionary
  !> @param[in] enum      Enumerator value to be checked
  function main_coll_check(s, enum) result(ok)
    implicit none
    class(main_coll_dict_type), intent(in) :: s ! short for "self"
    integer(i_def), intent(in)             :: enum
    logical(l_def) :: ok
    ok = any ([s%derived, s%radiation, s%microphysics, s%electric, &
      s%orography, s%turbulence, s%convection, s%cloud, s%surface, &
      s%soil, s%snow, s%chemistry, s%aerosol, s%stph, s%lbc, s%none] == enum)
  end function main_coll_check

  !> @brief Getter for clock
  !> @param[in] self       Processor object
  !> @return               Clock returned
  function processor_get_clock(self) result(clock)
    implicit none
    class(processor_type), intent(in) :: self

    class(clock_type), pointer :: clock
    clock => self%clock
  end function processor_get_clock

  !> @brief Setter for clock
  !> @param[inout] self       Processor object
  !> @param[in]    clock      Model clock
  subroutine processor_set_clock(self, clock)
    implicit none
    class(processor_type), intent(inout) :: self
    class(clock_type), target, intent(in) :: clock

    self%clock => clock
  end subroutine processor_set_clock

  !> @brief Convenience function for creating field specifiers
  !> @param[in] name               Field name
  !> @param[in] main_coll          Enumerator of main fields collection
  !> @param[in, optional] space    Function space enumerator
  !> @param[in, optional] order_h  Function space order in horizontal direction
  !> @param[in, optional] order_v  Function space order in vertical direction
  !> @param[in, optional] adv_coll Enumerator of advected fields collection
  !> @param[in, optional] moist_arr  Moisture array enumerator
  !> @param[in, optional] moist_idx  Index into moisture array
  !> @param[in, optional] time_axis  Time axis enumerator
  !> @param[in, optional] mult     Name of multidata item, or blank string
  !> @param[in, optional] ckp      Is it a checkpoint (prognostic) field?
  !> @param[in, optional] twod     Is it two-dimensional?
  !> @param[in, optional] empty    Is it empty (with empty data array)?
  !> @param[in, optional] coarse   Is it on a coarse mesh?
  !> @param[in, optional] is_int   Is it an integer field?
  !> @param[in, optional] legacy   Is it a field with legacy checkpointing?
  !> @return                       Specifier returned
  function make_spec(name, main_coll, space, order_h, order_v, adv_coll, &
    moist_arr, moist_idx, time_axis, &
    mult, ckp, twod, empty, coarse, is_int, legacy) result(field_spec)
    implicit none
    character(*), intent(in) :: name
    integer(i_def), intent(in) :: main_coll
    integer(i_def), optional, intent(in) :: space
    integer(i_def), optional, intent(in) :: order_h
    integer(i_def), optional, intent(in) :: order_v
    integer(i_def), optional, intent(in) :: adv_coll
    integer(i_def), optional, intent(in) :: moist_arr
    integer(i_def), optional, intent(in) :: moist_idx
    integer(i_def), optional, intent(in) :: time_axis
    character(*), optional, intent(in) :: mult
    logical(l_def), optional, intent(in) :: ckp
    logical(l_def), optional, intent(in) :: twod
    logical(l_def), optional, intent(in) :: empty
    logical(l_def), optional, intent(in) :: coarse
    logical(l_def), optional, intent(in) :: is_int
    logical(l_def), optional, intent(in) :: legacy
    type(field_spec_type) :: field_spec

    field_spec%name = name
    field_spec%main_coll = main_coll
    if (present(space)) field_spec%space=space
    if (present(order_h)) field_spec%order_h=order_h
    if (present(order_v)) field_spec%order_v=order_v
    if (present(adv_coll)) field_spec%adv_coll=adv_coll
    if (present(moist_arr)) field_spec%moist_arr = moist_arr
    if (present(moist_idx)) field_spec%moist_idx = moist_idx
    if (present(time_axis)) field_spec%time_axis = time_axis
    if (present(mult)) field_spec%mult=mult
    if (present(ckp)) field_spec%ckp=ckp
    if (present(twod)) field_spec%twod=twod
    if (present(empty)) field_spec%empty=empty
    if (present(coarse)) field_spec%coarse=coarse
    if (present(is_int)) field_spec%is_int=is_int
    if (present(legacy)) field_spec%legacy=legacy

    if (.not. main_coll_dict%check(main_coll)) &
      call enum_error('main_coll', field_spec%main_coll)
  end function make_spec

  !> @brief If advected, return collection enumerator, otherwise NONE enumerator
  !> @details For use where the advected flag is known only at runtime.
  !> @param[in] advected   Is it for an advected field?
  !> @param[in] coll       Advected field collection enumerator to be used if needed
  !> @return               Collection returned
  function if_advected(advected, coll) result(adv_coll)
    use constants_mod, only : i_def
    implicit none
    logical(l_def), intent(in) :: advected
    integer(i_def), intent(in) :: coll

    integer(i_def) :: adv_coll
    adv_coll = coll
    if (.not. advected) adv_coll = adv_coll_dict%none
  end function if_advected

  !> @brief Return true if and only if a space is supported by XIOS
  !> @details In legacy mode, all the spaces are supported. In modern mode,
  !> W2 fields (like u) cannot be written directly. These will be checkpointed
  !> via the decomposition implemented by the routines split_complex_prognostics /
  !> combine_complex_prognostics in gungho_init_fields_mod.X90.
  !> @param[in] fs         Function space enumerator
  !> @param[in] legacy     Are we using legacy checkpoint domains (checkpoint_W2, etc.)?
  !> @return               True if and only if space is supported
  function space_has_xios_io(fs, legacy) result(flag)
    use fs_continuity_mod,              only : W1, W2
    implicit none

    integer(i_def), intent(in) :: fs ! function space enumerator
    logical(l_def), optional, intent(in) :: legacy ! using legacy io?

    logical(l_def) :: use_legacy
    logical(l_def) :: flag

    use_legacy = .false.
    if (present(legacy)) use_legacy = legacy

    select case (fs)
    case (W1)
      flag = .false. ! there is no legacy domain for W1
    case (W2)
      ! supported in legacy mode, otherwise not
      flag = use_legacy
    case default
      flag = .true.
  end select

  end function space_has_xios_io

end module field_spec_mod
