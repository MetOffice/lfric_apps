!-------------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Create fields for holding diagnostics
!> @details Initialise a diagnostic field

module initialise_diagnostics_mod

  use constants_mod,                   only: i_def, l_def, str_def
  use field_mod,                       only: field_type
  use mesh_mod,                        only: mesh_type
  use lfric_xios_diag_mod,             only: field_is_active
  use io_config_mod,                   only: diag_always_on_sampling, use_xios_io
  use field_from_metadata_mod,         only: init_field_from_metadata
  use field_parent_mod,                only: write_interface
  use lfric_xios_write_mod,            only: write_field_generic

  implicit none

  private

  ! field status indicators
  character(str_def), parameter :: activated                                  &
    = 'Activated'     ! needed as a dependency
  character(str_def), parameter :: deactivated                                &
    = 'Deactivated'   ! needed as a dependency
  character(str_def), parameter :: enabled                                    &
    = 'Enabled'       ! dynamically enabled, will be sampled
  character(str_def), parameter :: disabled                                   &
    = 'Disabled'      ! dynamically disabledm will not be sampled

  public :: init_diagnostic_field, diagnostic_to_be_sampled

contains

  !> @brief Return true if and only if XIOS will sample the field at the current timestep.
  !> @param[in]   unique_id   XIOS id of field
  !> @return                  Sampling on/off status of the field
  function diagnostic_to_be_sampled(unique_id) result(sampling_on)
    implicit none
    character(*), intent(in) :: unique_id
    logical(l_def) :: sampling_on
    if (diag_always_on_sampling .or. .not. use_xios_io) then
      sampling_on = .true. ! for testing
      ! This is used when xios is off to ensure existing behaviour of
      ! old nodal output is preserved
    else
      sampling_on = field_is_active(unique_id, .true.) ! derived from metadata
    end if
  end function diagnostic_to_be_sampled

  !> @brief Initialise a diagnostic field.
  !> @details If the field was requested, or if it is needed as a dependency,
  !> it will be created as an active field, otherwise with an empty data
  !> array (to save memory).
  !> Pass activate=.true. for fields needed as dependencies. If activate
  !> is not passed, the field's status will be derived from the XIOS metadata.
  !> Passing activate=.false. is equivalent to not passing activate.
  !> @post  The field name will be set equal to the XIOS id passed in.
  !> @param[out]          field            Field to initialise
  !> @param[in]           unique_id        XIOS id of field
  !> @param[in, optional] activate         Force-activate field
  !> @param[in, optional] force_mesh       Override derived mesh
  !> @param[in, optional] force_rad_levels Override derived radiation levels
  !> @param[in, optional] deactivate       Force-deactivate field
  !> @return                               Sampling status of the field
  function init_diagnostic_field(field, unique_id,                            &
    activate, force_mesh, force_rad_levels, deactivate) result (sampling_on)

    implicit none

    type(field_type),                   intent(out) :: field
    character(*),                       intent(in)  :: unique_id
    logical(l_def), optional,           intent(in)  :: activate
    logical(l_def), optional,           intent(in)  :: deactivate
    type(mesh_type), pointer, optional, intent(in)  :: force_mesh
    integer(kind=i_def), optional,      intent(in)  :: force_rad_levels

    logical(kind=l_def) :: is_activated
    logical(kind=l_def) :: is_deactivated
    logical(kind=l_def) :: sampling_on
    logical(kind=l_def) :: active
    character(str_def)  :: status

    procedure(write_interface), pointer :: write_behaviour => null()

    ! field sampling status
    sampling_on = diagnostic_to_be_sampled(unique_id)

    ! field activation status
    is_activated = .false.
    is_deactivated = .false.
    if (present(activate)) then
      if (activate) then
        is_activated = .true.
      end if
    end if
    if (present(deactivate)) then
      if (deactivate) then
        is_deactivated = .true.
      end if
    end if
    if (is_activated) then
      active = .true.
      status = activated
    else if (is_deactivated) then
      active = .false.
      status = deactivated
    else
      active = sampling_on
      if (active) then
        status = enabled
      else
        status = disabled
      end if
    end if

    call init_field_from_metadata( &
      field, unique_id, empty=(.not. active), &
      force_mesh=force_mesh, force_rad_levels=force_rad_levels, &
      diag_status=status)

    if (active) then
      write_behaviour => write_field_generic
      call field%set_write_behaviour(write_behaviour)
    end if
  end function init_diagnostic_field

end module initialise_diagnostics_mod
