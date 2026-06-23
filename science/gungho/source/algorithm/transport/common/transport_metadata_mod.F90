!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Contains metadata for transport.
module transport_metadata_mod

  use constants_mod,        only: i_def, l_def, r_def, str_def, r_tran
  use linked_list_data_mod, only: linked_list_data_type

  implicit none

  private

  ! Public types
  type, extends(linked_list_data_type), public :: transport_metadata_type

    private

    character(len=str_def) :: fname ! Name of the field (or field_group)
    integer(kind=i_def)    :: equation_form  ! Form of transport equation ( = advective, conservative, consistent)
    integer(kind=i_def)    :: splitting ! Horizontal/vertical splitting ( = none, strang_vhv/_hvh, hv, vh)
    integer(kind=i_def)    :: scheme    ! Transport scheme (= mol3d, ffsl3d, split)
    integer(kind=i_def)    :: horizontal_method ! Horizontal transport method (= mol, ffsl)
    integer(kind=i_def)    :: vertical_method ! Vertical transport method (= mol, sl/slice, ffsl)
    integer(kind=i_def)    :: horizontal_monotone      ! Horizontal monotone scheme
    integer(kind=i_def)    :: vertical_monotone        ! Vertical monotone scheme
    integer(kind=i_def)    :: vertical_monotone_order  ! Order of the vertical monotone scheme
    integer(kind=i_def)    :: special_edges_monotone   ! Horizontal monotone at special edges
    logical(kind=l_def)    :: enforce_min_value  ! enforce a min value (=T/F)
    real(kind=r_tran)      :: min_value          ! the min value to be enforced
    logical(kind=l_def)    :: log_space ! Do interpolation in log space
    logical(kind=l_def)    :: reversible ! Use a reversible transport scheme
    integer(kind=i_def)    :: ffsl_splitting ! Which FFSL splitting to use
    integer(kind=i_def)    :: ffsl_vertical_order ! Which FFSL order to use for vertical reconstructions

    ! Stored internal values, allowing some options to be temporarily changed
    integer(kind=i_def)    :: true_equation_form
    integer(kind=i_def)    :: true_horizontal_method
    integer(kind=i_def)    :: true_vertical_method
    integer(kind=i_def)    :: true_horizontal_monotone
    integer(kind=i_def)    :: true_vertical_monotone
    integer(kind=i_def)    :: true_splitting

    contains

    procedure, public :: get_name
    procedure, public :: get_equation_form
    procedure, public :: get_splitting
    procedure, public :: get_scheme
    procedure, public :: get_horizontal_method
    procedure, public :: get_vertical_method
    procedure, public :: get_monotone
    procedure, public :: get_enforce_min_value
    procedure, public :: get_min_value
    procedure, public :: get_log_space
    procedure, public :: get_horizontal_monotone
    procedure, public :: get_vertical_monotone
    procedure, public :: get_vertical_monotone_order
    procedure, public :: get_special_edges_monotone
    procedure, public :: get_monotone_order
    procedure, public :: get_reversible
    procedure, public :: get_ffsl_splitting
    procedure, public :: get_ffsl_vertical_order
    procedure, public :: update_metadata
    procedure, public :: reset_metadata

  end type transport_metadata_type

  !-----------------------------------------------------------------------------
  ! Constructors
  !-----------------------------------------------------------------------------
  !> Function to construct a transport_metadata object
  interface transport_metadata_type
    module procedure transport_metadata_constructor
  end interface

contains

    function transport_metadata_constructor(fname, equation_form, splitting, &
                                            scheme, horizontal_method,       &
                                            vertical_method,                 &
                                            horizontal_monotone,             &
                                            vertical_monotone,               &
                                            vertical_monotone_order,         &
                                            special_edges_monotone,          &
                                            enforce_min_value,               &
                                            min_value,                       &
                                            log_space,                       &
                                            reversible,                      &
                                            ffsl_splitting,                  &
                                            ffsl_vertical_order)             &
                                            result(self)

    implicit none

    type(transport_metadata_type) :: self

    character(len=str_def), intent(in) :: fname
    integer(kind=i_def),    intent(in) :: equation_form
    integer(kind=i_def),    intent(in) :: splitting
    integer(kind=i_def),    intent(in) :: scheme
    integer(kind=i_def),    intent(in) :: horizontal_method
    integer(kind=i_def),    intent(in) :: vertical_method
    integer(kind=i_def),    intent(in) :: horizontal_monotone
    integer(kind=i_def),    intent(in) :: vertical_monotone
    integer(kind=i_def),    intent(in) :: vertical_monotone_order
    integer(kind=i_def),    intent(in) :: special_edges_monotone
    logical(kind=l_def),    intent(in) :: enforce_min_value
    real(kind=r_tran),      intent(in) :: min_value
    logical(kind=l_def),    intent(in) :: log_space
    logical(kind=l_def),    intent(in) :: reversible
    integer(kind=i_def),    intent(in) :: ffsl_splitting
    integer(kind=i_def),    intent(in) :: ffsl_vertical_order

    self%fname                   = trim(fname)
    self%equation_form           = equation_form
    self%splitting               = splitting
    self%scheme                  = scheme
    self%horizontal_method       = horizontal_method
    self%vertical_method         = vertical_method
    self%horizontal_monotone     = horizontal_monotone
    self%vertical_monotone       = vertical_monotone
    self%vertical_monotone_order = vertical_monotone_order
    self%special_edges_monotone  = special_edges_monotone
    self%enforce_min_value       = enforce_min_value
    self%min_value               = min_value
    self%log_space               = log_space
    self%reversible              = reversible
    self%ffsl_splitting          = ffsl_splitting
    self%ffsl_vertical_order     = ffsl_vertical_order

    ! Set stored true values
    self%true_equation_form = equation_form
    self%true_horizontal_method = horizontal_method
    self%true_vertical_method = vertical_method
    self%true_horizontal_monotone = horizontal_monotone
    self%true_vertical_monotone = vertical_monotone
    self%true_splitting = splitting

  end function transport_metadata_constructor

  !> @brief Get the name of the metadata
  !> @param[in] self     The transport_metadata object
  !> @return             The Name of the metadata object
  function get_name(self) result(fname)

    implicit none

    class(transport_metadata_type), intent(in) :: self
    character(len=str_def)                     :: fname

    fname = self%fname

  end function get_name

  !> @brief Get the equation form
  !> @param[in] self     The transport_metadata object
  !> @return             The equation type (conservative, advective)
  function get_equation_form(self) result(equation_form)

    implicit none

    class(transport_metadata_type), intent(in) :: self
    integer(kind=i_def)                        :: equation_form

    equation_form = self%equation_form

  end function get_equation_form

  !> @brief Get the splitting type
  !> @param[in] self     The transport_metadata object
  !> @return             The splitting type
  function get_splitting(self) result(splitting)

    implicit none

    class(transport_metadata_type), intent(in) :: self
    integer(kind=i_def)                        :: splitting

    splitting = self%splitting

  end function get_splitting

  !> @brief Get the scheme type
  !> @param[in] self     The transport_metadata object
  !> @return             The scheme type
  function get_scheme(self) result(scheme)

    implicit none

    class(transport_metadata_type), intent(in) :: self
    integer(kind=i_def)                        :: scheme

    scheme = self%scheme

  end function get_scheme

  !> @brief Get the horizontal method
  !> @param[in] self     The transport_metadata object
  !> @return             The horizontal_method
  function get_horizontal_method(self) result(horizontal_method)

    implicit none

    class(transport_metadata_type), intent(in) :: self
    integer(kind=i_def)                        :: horizontal_method

    horizontal_method = self%horizontal_method

  end function get_horizontal_method

  !> @brief Get the vertical method
  !> @param[in] self     The transport_metadata object
  !> @return             The vertical method
  function get_vertical_method(self) result(vertical_method)

    implicit none

    class(transport_metadata_type), intent(in) :: self
    integer(kind=i_def)                        :: vertical_method

    vertical_method = self%vertical_method

  end function get_vertical_method

  !> @brief Get the monotone option
  !> @param[in] self       The transport_metadata object
  !> @param[in] direction  The direction of transport
  !> @return               The monotone option
  function get_monotone(self,direction) result(monotone)
    use transport_enumerated_types_mod, only: direction_h,  &
                                              direction_v,  &
                                              direction_3d
    implicit none

    class(transport_metadata_type), intent(in) :: self
    integer(kind=i_def),            intent(in) :: direction
    integer(kind=i_def)                        :: monotone

    select case (direction)
      case (direction_h, direction_3d)
         monotone = self%horizontal_monotone
      case (direction_v)
         monotone = self%vertical_monotone
    end select
  end function get_monotone

  !> @brief Get the monotone_order
  !> @param[in] self       The transport_metadata object
  !> @return               The monotone_order
  function get_monotone_order(self) result(monotone_order)
    implicit none
    class(transport_metadata_type), intent(in) :: self
    integer(kind=i_def)                        :: monotone_order

    monotone_order = self%vertical_monotone_order

  end function get_monotone_order

  !> @brief Get the horizontal_monotone option
  !> @param[in] self     The transport_metadata object
  !> @return             The horizontal_monotone option
  function get_horizontal_monotone(self) result(horizontal_monotone)

    implicit none

    class(transport_metadata_type), intent(in) :: self
    integer(kind=i_def)                        :: horizontal_monotone

    horizontal_monotone = self%horizontal_monotone

  end function get_horizontal_monotone

  !> @brief Get the vertical_monotone option
  !> @param[in] self     The transport_metadata object
  !> @return             The vertical_monotone option
  function get_vertical_monotone(self) result(vertical_monotone)

    implicit none

    class(transport_metadata_type), intent(in) :: self
    integer(kind=i_def)                        :: vertical_monotone

    vertical_monotone = self%vertical_monotone

  end function get_vertical_monotone

  !> @brief Get the vertical_monotone_order
  !> @param[in] self     The transport_metadata object
  !> @return             The vertical_monotone_order
  function get_vertical_monotone_order(self) result(vertical_monotone_order)

    implicit none

    class(transport_metadata_type), intent(in) :: self
    integer(kind=i_def)                        :: vertical_monotone_order

    vertical_monotone_order = self%vertical_monotone_order

  end function get_vertical_monotone_order

  !> @brief Get the special_edges_monotone option
  !> @param[in] self     The transport_metadata object
  !> @return             The special_edges_monotone option
  function get_special_edges_monotone(self) result(special_edges_monotone)

    implicit none

    class(transport_metadata_type), intent(in) :: self
    integer(kind=i_def)                        :: special_edges_monotone

    special_edges_monotone = self%special_edges_monotone

  end function get_special_edges_monotone

  !> @brief Get the enforce_min_value option
  !> @param[in] self     The transport_metadata object
  !> @return             The enforce_min_value switch
  function get_enforce_min_value(self) result(enforce_min_value)

    implicit none

    class(transport_metadata_type), intent(in) :: self
    logical(kind=l_def)                        :: enforce_min_value

    enforce_min_value = self%enforce_min_value

  end function get_enforce_min_value

  !> @brief Get the enforce_min_value option
  !> @param[in] self     The transport_metadata object
  !> @return             The min value enforced
  function get_min_value(self) result(min_value)

    implicit none

    class(transport_metadata_type), intent(in) :: self
    real(kind=r_tran)                          :: min_value

    min_value = self%min_value

  end function get_min_value

  !> @brief Get the log space option
  !> @param[in] self     The transport_metadata object
  !> @return             The log space switch
  function get_log_space(self) result(log_space)

    implicit none

    class(transport_metadata_type), intent(in) :: self
    logical(kind=l_def)                        :: log_space

    log_space = self%log_space

  end function get_log_space

  !> @brief Get the reversible option
  !> @param[in] self     The transport_metadata object
  !> @return             The reversible switch
  function get_reversible(self) result(reversible)

    implicit none

    class(transport_metadata_type), intent(in) :: self
    logical(kind=l_def)                        :: reversible

    reversible = self%reversible

  end function get_reversible

  !> @brief Get the splitting to use for consistent FFSL tracer transport
  !> @param[in] self     The transport_metadata object
  !> @return             The splitting
  function get_ffsl_splitting(self) result(ffsl_splitting)

    implicit none

    class(transport_metadata_type), intent(in) :: self
    integer(kind=i_def)                        :: ffsl_splitting

    ffsl_splitting = self%ffsl_splitting

  end function get_ffsl_splitting

  !> @brief Get the splitting to use for consistent FFSL tracer transport
  !> @param[in] self     The transport_metadata object
  !> @return             The splitting
  function get_ffsl_vertical_order(self) result(ffsl_vertical_order)

    implicit none

    class(transport_metadata_type), intent(in) :: self
    integer(kind=i_def)                        :: ffsl_vertical_order

    ffsl_vertical_order = self%ffsl_vertical_order

  end function get_ffsl_vertical_order

  !> @brief Update the metadata based on the outer loop
  !> @details Updates the metadata options for this outer loop, for instance
  !!          if reverting to advective form or not applying monotonicity
  !> @param[in,out] self     The transport_metadata object
  !> @param[in]     outer    The iteration of the semi-implicit outer loop
  !> @param[in]     adaptive_splitting  Whether to use adaptive splitting. If
  !!                         true, sets the splitting to the splitting argument.
  !> @param[in]     splitting The splitting to use if adaptive_splitting is true
  !> @param[in]     make_tracers_advective Whether to set the equation form to be
  !!                         advective when unable to compute departure points
  !!                         for consistent tracers
  subroutine update_metadata(self, outer, adaptive_splitting, splitting,       &
                             make_tracers_advective)

    use timestepping_config_mod,    only: outer_iterations
    use transport_config_mod,       only: si_outer_transport,                  &
                                          si_outer_transport_none,             &
                                          si_outer_transport_no_mono,          &
                                          si_outer_transport_advective,        &
                                          si_outer_transport_horizontal_sl,    &
                                          dry_field_name
    use transport_enumerated_types_mod,                                        &
                                    only: equation_form_advective,             &
                                          equation_form_consistent,            &
                                          split_method_ffsl,                   &
                                          split_method_sl,                     &
                                          monotone_qm_pos,                     &
                                          monotone_strict,                     &
                                          monotone_none

    implicit none

    class(transport_metadata_type), intent(inout) :: self
    integer(kind=i_def),            intent(in)    :: outer
    logical(kind=l_def),            intent(in)    :: adaptive_splitting
    integer(kind=i_def),            intent(in)    :: splitting
    logical(kind=l_def),            intent(in)    :: make_tracers_advective

    if (si_outer_transport /= si_outer_transport_none                          &
        .and. outer < outer_iterations                                         &
        .and. self%fname /= dry_field_name) then

      if (si_outer_transport == si_outer_transport_advective                   &
          .or. si_outer_transport == si_outer_transport_horizontal_sl) then
        ! Set equation form to be advective
        self%equation_form = equation_form_advective

        ! If scheme is FFSL, set scheme to be semi-Lagrangian
        if (self%horizontal_method == split_method_ffsl) then
          self%horizontal_method = split_method_sl
          ! Ensure correct monotone options
          if (self%horizontal_monotone == monotone_qm_pos ) then
            self%horizontal_monotone = monotone_strict
          end if
        end if

        if (si_outer_transport == si_outer_transport_advective                 &
            .and. self%vertical_method == split_method_ffsl) then
          self%vertical_method = split_method_sl
          ! Ensure correct monotone options
          if (self%vertical_monotone == monotone_qm_pos ) then
            self%vertical_monotone = monotone_strict
          end if
        end if

      else if (si_outer_transport == si_outer_transport_no_mono) then
        self%horizontal_monotone = monotone_none
        self%vertical_monotone = monotone_none
      end if

    end if

    if (adaptive_splitting) then
      self%splitting = splitting
    end if

    if ( make_tracers_advective .and.                                          &
         self%equation_form == equation_form_consistent ) then
      ! Set equation form to be advective
      self%equation_form = equation_form_advective
    end if

  end subroutine update_metadata

  !> @brief Reset the metadata to its original values
  !> @details Resets the transport metadata to its original values if it has
  !!          been temporarily changed, e.g. by not applying monotonicity or
  !!          using the advective form for a given outer loop
  !> @param[in,out] self     The transport_metadata object
  subroutine reset_metadata(self)

    implicit none

    class(transport_metadata_type), intent(inout) :: self

    self%equation_form = self%true_equation_form
    self%horizontal_method = self%true_horizontal_method
    self%vertical_method = self%true_vertical_method
    self%horizontal_monotone = self%true_horizontal_monotone
    self%vertical_monotone = self%true_vertical_monotone
    self%splitting = self%true_splitting

  end subroutine reset_metadata

end module transport_metadata_mod
