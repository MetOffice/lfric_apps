!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Auxiliary class for mapping field specifiers to
!> prognostic field collections and time axes
!>
module field_mapper_mod

  use constants_mod,                      only : i_def
  use log_mod,                            only : log_event,       &
                                                 log_level_error, &
                                                 log_scratch_space
  use field_mod,                          only : field_type
  use field_collection_mod,               only : field_collection_type
  use field_array_mod,                    only : field_array_type
  use field_spec_mod,                     only : main_coll_dict, &
                                                 adv_coll_dict, &
                                                 moist_arr_dict, &
                                                 time_axis_dict
  use fs_continuity_mod,                  only : W2, W3, Wtheta, W2H
  use boundaries_config_mod,              only : limited_area
  use gungho_time_axes_mod,               only : gungho_time_axes_type
  use lfric_xios_time_axis_mod,           only : time_axis_type

  implicit none

  private

  !> @brief Object type providing access to field collectios and
  !> time axes to a field maker
  type :: field_mapper_type

  private

    type(field_collection_type), pointer :: depository
    type(field_collection_type), pointer :: prognostic
    type(field_collection_type), pointer :: adv_all_outer
    type(field_collection_type), pointer :: adv_last_outer
    type(field_collection_type), pointer :: con_all_outer
    type(field_collection_type), pointer :: con_last_outer
    type(field_collection_type), pointer :: derived
    type(field_collection_type), pointer :: radiation
    type(field_collection_type), pointer :: microphysics
    type(field_collection_type), pointer :: electric
    type(field_collection_type), pointer :: orography
    type(field_collection_type), pointer :: turbulence
    type(field_collection_type), pointer :: convection
    type(field_collection_type), pointer :: cloud
    type(field_collection_type), pointer :: surface
    type(field_collection_type), pointer :: soil
    type(field_collection_type), pointer :: snow
    type(field_collection_type), pointer :: chemistry
    type(field_collection_type), pointer :: aerosol
    type(field_collection_type), pointer :: stph
    type(field_collection_type), pointer :: moisture
    type(field_collection_type), pointer :: lbc

    type(gungho_time_axes_type), pointer :: gungho_axes

  contains
    private

    ! accessors
    procedure, public :: get_depository
    procedure, public :: get_prognostic_fields
    procedure, public :: get_gungho_axes

    ! main interface
    procedure, public :: init
    procedure, public :: sanity_check
    procedure, public :: get_adv_coll_ptr
    procedure, public :: get_main_coll_ptr
    procedure, public :: get_moist_field_ptr
    procedure, public :: get_time_axis_ptr

    ! destructor - here to avoid gnu compiler bug
    final :: field_mapper_destructor
  end type field_mapper_type

  public field_mapper_type

contains

  !> @brief Accessor for depository collection
  !> @param[in] self       Field mapper object
  !> @return               Depository collection returned
  function get_depository(self) result(collection)
    implicit none
    class(field_mapper_type), intent(in) :: self

    class(field_collection_type), pointer :: collection
    collection => self%depository
  end function get_depository

  !> @brief Accessor for prognostic_fields collection
  !> @param[in] self       Field mapper object
  !> @return               Prognostic fields collection returned
  function get_prognostic_fields(self) result(collection)
    implicit none
    class(field_mapper_type), intent(in) :: self

    class(field_collection_type), pointer :: collection
    collection => self%prognostic
  end function get_prognostic_fields

  !> @brief Accessor for gungho axes registry
  !> @param[in] self       Field mapper object
  !> @return               Gungho axes object pointer returned
  function get_gungho_axes(self) result(axes)
    implicit none
    class(field_mapper_type), intent(in) :: self

    type(gungho_time_axes_type), pointer :: axes

    axes => self%gungho_axes
  end function get_gungho_axes

  !> @brief Initialise a field mapper object.
  !> @param[in,out] self Field mapper object
  !> @param[in,out] depository Main collection of all fields in memory
  !> @param[in,out] moisture_fields Collection of moisture field arrays
  !> @param[in,out] prognostic_fields The prognostic variables in the model
  !> @param[in,out] adv_tracer_all_outer Collection of fields that need to be advected every outer iteration
  !> @param[in,out] adv_tracer_last_outer Collection of fields that need to be advected at final outer iteration
  !> @param[in,out] con_tracer_all_outer Second collection of fields that need to be advected every outer iteration
  !> @param[in,out] con_tracer_last_outer Second collection of fields that need to be advected at final outer iteration
  !> @param[in,out] derived_fields Collection of FD fields derived from FE fields
  !> @param[in,out] radition_fields Collection of fields for radiation scheme
  !> @param[in,out] microphysics_fields Collection of fields for microphys scheme
  !> @param[in,out] electric_fields Collection of fields for electric scheme
  !> @param[in,out] orography_fields Collection of fields for orogr drag scheme
  !> @param[in,out] turbulence_fields Collection of fields for turbulence scheme
  !> @param[in,out] convection_fields Collection of fields for convection scheme
  !> @param[in,out] cloud_fields Collection of fields for cloud scheme
  !> @param[in,out] surface_fields Collection of fields for surface scheme
  !> @param[in,out] soil_fields Collection of fields for soil hydrology scheme
  !> @param[in,out] snow_fields Collection of fields for snow scheme
  !> @param[in,out] aerosol_fields Collection of fields for aerosol scheme
  !> @param[in,out] chemistry_fields Collection of fields for chemistry scheme
  !> @param[in,out] stph_fields Collection of fields for stph scheme
  !> @param[in,out] lbc_fields Collection of lbc fields
  !> @param[in,out] gungho_axes Registry of time axes
  subroutine init(self,    &
    depository_fields,     &
    moisture_fields,       &
    prognostic_fields,     &
    adv_tracer_all_outer,  &
    adv_tracer_last_outer, &
    con_tracer_all_outer,  &
    con_tracer_last_outer, &
    derived_fields,        &
    radiation_fields,      &
    microphysics_fields,   &
    electric_fields,       &
    orography_fields,      &
    turbulence_fields,     &
    convection_fields,     &
    cloud_fields,          &
    surface_fields,        &
    soil_fields,           &
    snow_fields,           &
    chemistry_fields,      &
    aerosol_fields,        &
    stph_fields,           &
    lbc_fields,            &
    gungho_axes)

    implicit none

    class(field_mapper_type), intent(inout) :: self
    type(field_collection_type), target, intent(inout) :: depository_fields
    type(field_collection_type), target, intent(inout) :: moisture_fields
    type(field_collection_type), target, intent(inout) :: prognostic_fields
    type(field_collection_type), target, intent(inout) :: adv_tracer_all_outer
    type(field_collection_type), target, intent(inout) :: adv_tracer_last_outer
    type(field_collection_type), target, intent(inout) :: con_tracer_all_outer
    type(field_collection_type), target, intent(inout) :: con_tracer_last_outer
    type(field_collection_type), target, intent(inout) :: derived_fields
    type(field_collection_type), target, intent(inout) :: radiation_fields
    type(field_collection_type), target, intent(inout) :: microphysics_fields
    type(field_collection_type), target, intent(inout) :: electric_fields
    type(field_collection_type), target, intent(inout) :: orography_fields
    type(field_collection_type), target, intent(inout) :: turbulence_fields
    type(field_collection_type), target, intent(inout) :: convection_fields
    type(field_collection_type), target, intent(inout) :: cloud_fields
    type(field_collection_type), target, intent(inout) :: surface_fields
    type(field_collection_type), target, intent(inout) :: soil_fields
    type(field_collection_type), target, intent(inout) :: snow_fields
    type(field_collection_type), target, intent(inout) :: chemistry_fields
    type(field_collection_type), target, intent(inout) :: aerosol_fields
    type(field_collection_type), target, intent(inout) :: stph_fields
    type(field_collection_type), target, intent(inout) :: lbc_fields

    type(gungho_time_axes_type), target, intent(inout) :: gungho_axes

    self%depository => depository_fields
    self%moisture => moisture_fields
    self%prognostic => prognostic_fields
    self%adv_all_outer => adv_tracer_all_outer
    self%adv_last_outer => adv_tracer_last_outer
    self%con_all_outer => con_tracer_all_outer
    self%con_last_outer => con_tracer_last_outer
    self%derived => derived_fields
    self%radiation => radiation_fields
    self%microphysics => microphysics_fields
    self%electric => electric_fields
    self%orography => orography_fields
    self%turbulence => turbulence_fields
    self%convection => convection_fields
    self%cloud => cloud_fields
    self%surface => surface_fields
    self%soil => soil_fields
    self%snow => snow_fields
    self%chemistry => chemistry_fields
    self%aerosol => aerosol_fields
    self%stph => stph_fields
    self%lbc => lbc_fields

    self%gungho_axes => gungho_axes

    ! Create collection of fields to be advected

#ifdef UM_PHYSICS
#endif
   call gungho_axes%initialise()

  end subroutine init

  !> @brief Post-initialisation sanity check
  !> @param[in] self Field mapper object
  subroutine sanity_check(self)
    implicit none
    class(field_mapper_type), intent(in) :: self

    type( field_type ), pointer :: theta => null()

    integer(i_def) :: theta_space

    call self%prognostic%get_field('theta', theta)
    theta_space = theta%which_function_space()

    if (theta_space /= Wtheta)then
      call log_event( 'Physics: requires theta variable to be in Wtheta',      &
                      log_level_error )
    end if

    if (theta%get_element_order_h() > 0 .or. &
        theta%get_element_order_v() > 0) then
      call log_event( 'Physics: requires lowest order elements',               &
                      log_level_error )
    end if
  end subroutine sanity_check

  !> @brief Map advection collection enumerator to collection pointer.
  !> @param[in] self     Field mapper object
  !> @param[in] adv_coll Advection collection enumerator
  !> @return             Collection returned
  function get_adv_coll_ptr(self, adv_coll) result(coll_ptr)
    implicit none
    class(field_mapper_type), intent(in) :: self
    integer(i_def), intent(in) :: adv_coll

    type(field_collection_type), pointer :: coll_ptr

    select case(adv_coll)
    case(adv_coll_dict%none)
      coll_ptr => null()
    case(adv_coll_dict%all_adv)
      coll_ptr => self%adv_all_outer
    case(adv_coll_dict%last_adv)
      coll_ptr => self%adv_last_outer
    case(adv_coll_dict%all_con)
      coll_ptr => self%con_all_outer
    case(adv_coll_dict%last_con)
      coll_ptr => self%con_last_outer
    case default
      coll_ptr => null()
      call log_event('unexpected advected collection enumerator', log_level_error)
    end select
  end function get_adv_coll_ptr

  !> @brief Map main collection enumerator to collection pointer.
  !> @param[in] self      Field mapper object
  !> @param[in] main_coll Main collection enumerator
  !> @return              Collection returned
  function get_main_coll_ptr(self, main_coll) result(coll_ptr)
    implicit none
    class(field_mapper_type), intent(in) :: self
    integer(i_def), intent(in) :: main_coll

    type(field_collection_type), pointer :: coll_ptr

    select case(main_coll)
    case (main_coll_dict%derived)
      coll_ptr => self%derived
    case (main_coll_dict%radiation)
      coll_ptr => self%radiation
    case (main_coll_dict%microphysics)
      coll_ptr => self%microphysics
    case (main_coll_dict%electric)
      coll_ptr => self%electric
    case (main_coll_dict%orography)
      coll_ptr => self%orography
    case (main_coll_dict%turbulence)
      coll_ptr => self%turbulence
    case (main_coll_dict%convection)
      coll_ptr => self%convection
    case (main_coll_dict%cloud)
      coll_ptr => self%cloud
    case (main_coll_dict%surface)
      coll_ptr => self%surface
    case (main_coll_dict%soil)
      coll_ptr => self%soil
    case (main_coll_dict%snow)
      coll_ptr => self%snow
    case (main_coll_dict%chemistry)
      coll_ptr => self%chemistry
    case (main_coll_dict%aerosol)
      coll_ptr => self%aerosol
    case (main_coll_dict%stph)
      coll_ptr => self%stph
    case (main_coll_dict%lbc)
      coll_ptr => self%lbc
    case (main_coll_dict%none)
      coll_ptr => null()
    case default
      coll_ptr => null()
      call log_event('unexpected main collection enumerator', log_level_error)
    end select
  end function get_main_coll_ptr

  !> @brief Accessor for moisture fields
  !> @param[in] self      Field mapper object
  !> @param[in] moist_arr Moisture array enumerator
  !> @param[in] moist_idx Moisture array index
  !> @return              field returned, null if not a moisture field
  function get_moist_field_ptr(self, moist_arr, moist_idx) result(field)
    implicit none
    class(field_mapper_type), intent(in) :: self

    integer(i_def), intent(in) :: moist_arr
    integer(i_def), intent(in) :: moist_idx

    type(field_array_type), pointer :: arr
    type(field_type), pointer :: field

    select case(moist_arr)
    case (moist_arr_dict%mr)
      call self%moisture%get_field("mr", arr)
    case (moist_arr_dict%moist_dyn)
      call self%moisture%get_field("moist_dyn", arr)
    case (moist_arr_dict%moist_dyn_ref)
      call self%moisture%get_field("moist_dyn_ref", arr)
    case (moist_arr_dict%none)
      arr => null()
    case default
      arr => null()
      call log_event('unexpected moisture array enumerator', log_level_error)
    end select

    if (associated(arr)) then
      field => arr%bundle(moist_idx)
    else
      field => null()
    end if
  end function get_moist_field_ptr

  !> @brief Accessor for time axes
  !> @param[in] self      Field mapper object
  !> @param[in] time_axis Time axis enumerator
  !> @return              Axis returned, null if field without time axis
  function get_time_axis_ptr(self, time_axis) result(axis)
    implicit none
    class(field_mapper_type), intent(in) :: self

    integer(i_def), intent(in) :: time_axis

    type(time_axis_type), pointer :: axis

    select case(time_axis)
    case (time_axis_dict%lbc)
      axis => self%gungho_axes%lbc_time_axis
    case (time_axis_dict%ls)
      axis => self%gungho_axes%lbc_time_axis
    case (time_axis_dict%none)
      axis => null()
    case default
      axis => null()
      call log_event('unexpected time axis enumerator', log_level_error)
    end select
  end function get_time_axis_ptr

  !> @brief Destructor for field mapper objects
  !> @param[inout] self       Field mapper object
  subroutine field_mapper_destructor(self)
    type(field_mapper_type), intent(inout) :: self
    ! empty
  end subroutine field_mapper_destructor

end module field_mapper_mod
