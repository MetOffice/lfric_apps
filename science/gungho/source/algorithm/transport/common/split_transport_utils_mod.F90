!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Contains utility routines relating to split transport schemes
module split_transport_utils_mod

  use constants_mod,                  only: i_def, l_def, r_tran, str_def
  use log_mod,                        only: log_event,                         &
                                            log_scratch_space,                 &
                                            LOG_LEVEL_ERROR
  use fs_continuity_mod,              only: W2, W2H, W2V
  use transport_enumerated_types_mod, only: direction_h,                       &
                                            direction_v,                       &
                                            direction_3d,                      &
                                            splitting_strang_hvh,              &
                                            splitting_strang_vhv,              &
                                            splitting_vhvhv_third,             &
                                            splitting_vhvhv_quarter,           &
                                            splitting_hv,                      &
                                            splitting_vh,                      &
                                            splitting_vhhv,                    &
                                            splitting_none

  implicit none

  private

  ! Index of the configuration relating to the designated *dry* reference field
  ! TODO: in future this could be stored in modeldb
  integer(kind=i_def), allocatable :: dry_config_index

  ! List of indices of the unique split step fractions
  ! TODO: in future this could be stored in modeldb
  integer(kind=i_def), allocatable :: fraction_idxs(:)

  ! Whether the shifted mesh is needed for transport of Wtheta variables
  ! TODO: in future this could be stored in modeldb
  logical(kind=l_def), allocatable :: need_shifted_mesh

  public :: get_dry_config
  public :: get_fraction_idx
  public :: get_fraction_from_idx
  public :: get_max_num_fractions
  public :: use_shifted_mesh
  public :: get_direction_w2_fs
  public :: get_splitting_factor
  public :: get_splitting_fraction
  public :: get_splitting_direction
  public :: get_splitting_name
  public :: get_num_split_steps
  public :: get_next_step_hori
  public :: get_first_hori_step

  public :: finalise_split_transport_utils

  private :: compute_fraction_idxs
  private :: splitting_error_message

contains

! ============================================================================ !
! CONFIGURATION UTILITIES
! ============================================================================ !
  !> @brief Finds the configuration index of the designated dry reference field
  !> @result dry_config_out  The configuration index of the dry reference field
  ! TODO: in future this should take modeldb as an argument
  function get_dry_config() result(dry_config_out)

    use transport_config_mod, only: field_names, profile_size, dry_field_name

    implicit none

    ! Internal arguments
    integer(kind=i_def) :: config, dry_config_out
    logical(kind=l_def) :: found_dry_config

    if (.not. allocated(dry_config_index)) then
      found_dry_config = .false.
      allocate(dry_config_index)
      ! Need to find dry config
      do config = 1, profile_size
        if ( trim(field_names(config)) == trim(dry_field_name) ) then
          found_dry_config = .true.
          dry_config_index = config
          exit
        end if
      end do

      if (.not. found_dry_config) then
        write(log_scratch_space, '(2A)') &
          'transport_utils: Cannot find config for ', trim(dry_field_name)
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
      end if
    end if

    ! Return index of dry config
    dry_config_out = dry_config_index

  end function get_dry_config

  !> @brief Finds the index for a particular splitting/step
  !> @param[in] splitting  Enumerator for particular splitting
  !> @param[in] step       Index of the splitting step
  !> @result    idx        The index of the splitting
  function get_fraction_idx(splitting, step) result(idx)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: splitting
    integer(kind=i_def), intent(in) :: step

    ! Internal variables
    integer(kind=i_def) :: i, idx, frac

    ! Set up list of splitting indices if this is not already done
    if (.not. allocated(fraction_idxs)) then
      call compute_fraction_idxs()
    end if

    ! Get splitting fraction
    frac = get_splitting_fraction(splitting, step)

    ! Find the index of the splitting in the list
    idx = 0
    do i = 1, size(fraction_idxs)
      if (fraction_idxs(i) == frac) then
        idx = i
        exit
      end if
    end do

    if (idx == 0) then
      call log_event(                                                          &
        'transport_utils: fraction_idx not found',                             &
        LOG_LEVEL_ERROR                                                        &
      )
    end if

  end function get_fraction_idx

  !> @brief Finds the split-step fraction given its unique index
  !> @param[in] idx        The index of the fraction in the list
  !> @result    fraction   The denominator for this step
  function get_fraction_from_idx(idx) result(frac)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: idx

    ! Internal variables
    integer(kind=i_def) :: frac

    ! Set up list of splitting indices if this is not already done
    if (.not. allocated(fraction_idxs)) then
      call compute_fraction_idxs()
    end if

    ! Get splitting fraction
    frac = fraction_idxs(idx)

  end function get_fraction_from_idx

  !> @brief Finds the maximum number of unique fractions for a split scheme
  !> @result max_num_fraction  The maximum number of unique split step fractions
  ! TODO: in future this should take modeldb as an argument
  function get_max_num_fractions() result(max_num_fractions)

    implicit none

    ! Internal arguments
    integer(kind=i_def) :: max_num_fractions

    ! Set up list of fraction indices if this is not already done
    if (.not. allocated(fraction_idxs)) then
      call compute_fraction_idxs()
    end if

    max_num_fractions = size(fraction_idxs)

  end function get_max_num_fractions

  !> @brief Finds whether this configuration needs a shifted mesh
  !> @result any_shifted  Whether shifted mesh is needed
  ! TODO: in future this should take modeldb as an argument
  function use_shifted_mesh() result(any_shifted)

    use extrusion_config_mod,           only: number_of_layers
    use formulation_config_mod,         only: moisture_formulation,            &
                                              moisture_formulation_dry
    use io_config_mod,                  only: write_conservation_diag
    use transport_config_mod,           only: profile_size,                    &
                                              field_names,                     &
                                              dry_field_name,                  &
                                              equation_form,                   &
                                              scheme,                          &
                                              horizontal_method,               &
                                              vertical_method
    use transport_enumerated_types_mod, only: equation_form_advective,         &
                                              equation_form_conservative,      &
                                              equation_form_consistent,        &
                                              split_method_ffsl,               &
                                              scheme_ffsl_3d,                  &
                                              scheme_split

    implicit none

    ! Internal arguments
    integer(kind=i_def) :: i
    logical(kind=l_def) :: any_shifted

    if (.not. allocated(need_shifted_mesh)) then

      allocate(need_shifted_mesh)
      need_shifted_mesh = .false.

      ! Need a shifted mesh if:
      ! (a) a variable uses the conservative or consistent transport equation
      ! (b) a Wtheta variable uses FFSL
      do i = 1, profile_size
        ! Check for a variable using conservative/consistent equation
        ! (but don't include "dry_field" which will never use shifted grid)
        if ( (equation_form(i) == equation_form_conservative .and. &
              field_names(i) /= dry_field_name) .or.               &
              equation_form(i) == equation_form_consistent ) then
            need_shifted_mesh = .true.
            exit
        end if

        ! Check if there is a transport scheme using FFSL
        select case (scheme(i))
          ! It could be either 3D FFSL or split scheme using FFSL
          case (scheme_ffsl_3d)
            if (equation_form(i) == equation_form_advective) then
              need_shifted_mesh = .true.
              exit
            end if

          case (scheme_split)
            if ( vertical_method(i) == split_method_ffsl &
                  .or. horizontal_method(i) == split_method_ffsl ) then
              if (equation_form(i) == equation_form_advective) then
                need_shifted_mesh = .true.
                exit
              end if
            end if

        end select
      end do

      if (write_conservation_diag .and. &
          moisture_formulation /= moisture_formulation_dry) then
        need_shifted_mesh = .true.
      end if

      ! TODO: total hack to enforce no shifted mesh in shallow water model
      if ((dry_field_name == 'q' .or. dry_field_name == 'geopot')              &
           .and. number_of_layers == 1) then
        need_shifted_mesh = .false.
      end if
    end if

    ! Return whether shifted mesh is needed
    any_shifted = need_shifted_mesh

  end function use_shifted_mesh

! ============================================================================ !
! UTILITIES THAT DESCRIBE SPLITTINGS
! ============================================================================ !

  !> @brief Returns the directional W2 function space for a given split step
  !> @param[in] splitting  Enumerator for particular splitting
  !> @param[in] step       Index of the splitting step
  !> @result    fs         Enumerator for the W2 function space
  function get_direction_w2_fs(splitting, step) result(fs)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: splitting
    integer(kind=i_def), intent(in) :: step

    ! Internal variables
    integer(kind=i_def) :: fs, direction

    direction = get_splitting_direction(splitting, step)

    select case (direction)
    case (direction_v)
      fs = W2V
    case (direction_h)
      fs = W2H
    case (direction_3d)
      fs = W2
    case default
      call log_event(                                                          &
        'transport_utils: splitting direction not implemented',                &
        LOG_LEVEL_ERROR                                                        &
      )
    end select

  end function get_direction_w2_fs


  !> @brief Returns the number of steps in a particular splitting
  !> @param[in] splitting        Enumerator for particular splitting
  !> @result    num_split_steps  Number of steps used in the splitting
  function get_num_split_steps(splitting) result(num_split_steps)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: splitting

    ! Internal variables
    integer(kind=i_def) :: num_split_steps

    ! Compute total number of split steps
    select case (splitting)
    case (splitting_strang_hvh, splitting_strang_vhv)
      num_split_steps = 3
    case (splitting_vhhv)
      num_split_steps = 4
    case (splitting_vhvhv_third, splitting_vhvhv_quarter)
      num_split_steps = 5
    case (splitting_vh, splitting_hv)
      num_split_steps = 2
    case (splitting_none)
      num_split_steps = 1
    case default
      call log_event(                                                          &
        'transport_utils, get_num_split_steps: splitting not implemented',     &
         LOG_LEVEL_ERROR                                                       &
      )
    end select

  end function get_num_split_steps


  !> @brief Returns the directional enumerator for a given split step
  !> @param[in] splitting  Enumerator for particular splitting
  !> @param[in] step       Index of the splitting step
  !> @result    direction  Enumerator for the split direction
  function get_splitting_direction(splitting, step) result(direction)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: splitting
    integer(kind=i_def), intent(in) :: step

    ! Internal variables
    integer(kind=i_def) :: direction

    ! Compute total number of split steps
    select case (splitting)
    case (splitting_strang_hvh)
      select case (step)
      case (1, 3)
        direction = direction_h
      case (2)
        direction = direction_v
      case default
        call splitting_error_message(splitting, step)
      end select
    case (splitting_strang_vhv)
      select case (step)
      case (1, 3)
        direction = direction_v
      case (2)
        direction = direction_h
      case default
        call splitting_error_message(splitting, step)
      end select
    case (splitting_vhhv)
      select case (step)
      case (1, 4)
        direction = direction_v
      case (2, 3)
        direction = direction_h
      case default
        call splitting_error_message(splitting, step)
      end select
    case (splitting_vhvhv_third, splitting_vhvhv_quarter)
      select case (step)
      case (1, 3, 5)
        direction = direction_v
      case (2, 4)
        direction = direction_h
      case default
        call splitting_error_message(splitting, step)
      end select
    case (splitting_hv)
      select case (step)
      case (1)
        direction = direction_h
      case (2)
        direction = direction_v
      case default
        call splitting_error_message(splitting, step)
      end select
    case (splitting_vh)
      select case (step)
      case (1)
        direction = direction_v
      case (2)
        direction = direction_h
      case default
        call splitting_error_message(splitting, step)
      end select
    case (splitting_none)
      direction = direction_3d
    case default
      call log_event(                                                          &
        'transport_utils, get_direction_w2_fs: splitting not implemented',     &
         LOG_LEVEL_ERROR                                                       &
      )
    end select

  end function get_splitting_direction

  !> @brief Returns the factor to multiply dt for this splitting step
  !> @param[in] splitting  Enumerator for particular splitting
  !> @param[in] step       Index of the splitting step
  !> @result    factor     The multiplying factor
  function get_splitting_factor(splitting, step) result(factor)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: splitting
    integer(kind=i_def), intent(in) :: step

    ! Internal variables
    real(kind=r_tran)   :: factor
    integer(kind=i_def) :: frac

    frac = get_splitting_fraction(splitting, step)
    factor = 1.0_r_tran / real(frac, kind=r_tran)

  end function get_splitting_factor

  !> @brief Returns the integer denominator for the length of this split step
  !> @param[in] splitting  Enumerator for particular splitting
  !> @param[in] step       Index of the splitting step
  !> @result    fraction   The denominator for this step
  function get_splitting_fraction(splitting, step) result(frac)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: splitting
    integer(kind=i_def), intent(in) :: step

    ! Internal variables
    integer(kind=i_def) :: frac

    select case (splitting)
    case (splitting_strang_hvh, splitting_strang_vhv)
      select case (step)
      case (1, 3)
        frac = 2
      case (2)
        frac = 1
      end select
    case (splitting_vhhv)
      frac = 2
    case (splitting_vhvhv_quarter)
      select case (step)
      case (2, 3, 4)
        frac = 2
      case (1, 5)
        frac = 4
      end select
    case (splitting_vhvhv_third)
      select case (step)
      case (2, 4)
        frac = 2
      case (1, 3, 5)
        frac = 3
      end select
    case (splitting_hv, splitting_vh, splitting_none)
      frac = 1
    case default
      call log_event(                                                          &
        'transport_utils, get_splitting_fraction: splitting not implemented',  &
         LOG_LEVEL_ERROR                                                       &
      )
    end select

  end function get_splitting_fraction

  !> @brief Returns whether the next step in a split scheme is horizontal
  !> @param[in] splitting       Enumerator for particular splitting
  !> @param[in] step            Index of the current step
  !> @result    next_is_hori    Whether the next step is horizontal
  function get_next_step_hori(splitting, step) result(next_is_hori)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: splitting
    integer(kind=i_def), intent(in) :: step

    ! Internal variables
    logical(kind=l_def) :: next_is_hori
    integer(kind=i_def) :: num_split_steps
    integer(kind=i_def) :: direction_next_step

    num_split_steps = get_num_split_steps(splitting)
    if (step == num_split_steps) then
      next_is_hori = .false.
    else
      direction_next_step = get_splitting_direction(splitting, step+1)
      next_is_hori = (                                                         &
          direction_next_step == direction_h                                   &
          .or. direction_next_step == direction_3d                             &
      )
    end if

  end function get_next_step_hori

  !> @brief Returns the index of the first horizontal step in the splitting
  !> @param[in] splitting        Enumerator for particular splitting
  !> @result    first_hori_step  The first horizontal step
  function get_first_hori_step(splitting) result(first_hori_step)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: splitting

    ! Internal variables
    integer(kind=i_def) :: first_hori_step
    integer(kind=i_def) :: step
    integer(kind=i_def) :: num_split_steps
    integer(kind=i_def) :: direction_next_step

    first_hori_step = 0

    num_split_steps = get_num_split_steps(splitting)

    do step = 1, num_split_steps
      direction_next_step = get_splitting_direction(splitting, step)
      if (direction_next_step == direction_h                                   &
          .or. direction_next_step == direction_3d) then
        first_hori_step = step
        EXIT
      end if
    end do

  end function get_first_hori_step

  !> @brief Private routine to set up the list of splitting indices. This lists
  !!        the unique splitting fractions for vertical/horizontal directions
  !> TODO: in future this should take modeldb as an argument
  subroutine compute_fraction_idxs()

    use transport_config_mod, only: profile_size, splitting

    implicit none

    ! Internal arguments
    integer(kind=i_def), allocatable :: tmp_fraction_idxs(:)
    integer(kind=i_def) :: i, j, k, frac, num_unique
    integer(kind=i_def) :: max_num_unique_steps
    integer(kind=i_def) :: num_steps
    logical(kind=l_def) :: found_fraction

    ! Determine maximum possible number of unique steps ------------------------
    ! Most cautious value is that each step for each variable is unique
    max_num_unique_steps = 0
    do i = 1, profile_size
      num_steps = get_num_split_steps(splitting(i))
      max_num_unique_steps = max_num_unique_steps + num_steps
    end do

    allocate(tmp_fraction_idxs(max_num_unique_steps))

    ! Find unique step fractions -----------------------------------------------
    num_unique = 0
    tmp_fraction_idxs(:) = 0

    ! Loop through all steps for all splittings
    do i = 1, profile_size
      do j = 1, get_num_split_steps(splitting(i))
        frac = get_splitting_fraction(splitting(i), j)

        ! Loop through list of unique fractions -- is this one already there?
        found_fraction = .false.
        do k = 1, num_unique
          if (frac == tmp_fraction_idxs(k)) then
            found_fraction = .true.
            exit
          end if
        end do

        ! If this fraction is not already in the list, add it
        if (.not. found_fraction) then
          num_unique = num_unique + 1
          tmp_fraction_idxs(num_unique) = frac
        end if
      end do
    end do

    ! Populate permanent list of vertical step fraction indices
    allocate(fraction_idxs(num_unique))
    fraction_idxs(:) = tmp_fraction_idxs(1:num_unique)

    deallocate(tmp_fraction_idxs)

  end subroutine compute_fraction_idxs

  !> @brief Returns the name of the splitting
  !> @param[in] splitting  Enumerator for particular splitting
  !> @result    name       The name of the splitting
  function get_splitting_name(splitting) result(name)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: splitting

    ! Internal variables
    character(len=str_def) :: name

    select case (splitting)
    case (splitting_strang_hvh)
      name = 'Strang HVH'
    case (splitting_strang_vhv)
      name = 'Strang VHV'
    case (splitting_vhhv)
      name = 'VHHV'
    case (splitting_vhvhv_third)
      name = 'VHVHV third'
    case (splitting_vhvhv_quarter)
      name = 'VHVHV quarter'
    case (splitting_hv)
      name = 'HV'
    case (splitting_vh)
      name = 'VH'
    case (splitting_none)
      name = 'None'
    case default
      call log_event(                                                          &
        'transport_utils, get_splitting_name: splitting not implemented',      &
         LOG_LEVEL_ERROR                                                       &
      )
      name = 'Unknown'
    end select

  end function get_splitting_name

  !> @brief Writes an error message when trying to access information about a
  !!        step for a splitting that does not have this many steps
  !> @param[in] splitting   Enumerator of splitting
  !> @param[in] step        Index of step
  subroutine splitting_error_message(splitting, step)

    implicit none

    integer(kind=i_def), intent(in) :: splitting
    integer(kind=i_def), intent(in) :: step

    character(len=str_def) :: splitting_name

    splitting_name = get_splitting_name(splitting)

    write(log_scratch_space, '(A,I6,A,A,A)')                                   &
      'Attempting to access information for step ', step, ' of splitting ',    &
      adjustl(trim(splitting_name)), ' but it does not have a step with this index'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)

  end subroutine splitting_error_message


  !> @brief Finalises this module, resetting the dry config variable
  subroutine finalise_split_transport_utils()

    implicit none

    if (allocated(dry_config_index)) deallocate(dry_config_index)
    if (allocated(need_shifted_mesh)) deallocate(need_shifted_mesh)
    if (allocated(fraction_idxs)) deallocate(fraction_idxs)

  end subroutine finalise_split_transport_utils

end module split_transport_utils_mod
