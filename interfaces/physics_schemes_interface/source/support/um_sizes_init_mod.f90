!----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief Controls the setting of UM high level variables which are either
!>         fixed in LFRic or derived from LFRic inputs.

module um_sizes_init_mod

  ! Other modules used
  use constants_mod,               only : i_um, r_um, rmdi, i_def, r_def

  implicit none

  private
  public :: um_sizes_init

contains

  !>@brief Initialise UM high levels variables which are either fixed in LFRic
  !>        or derived from LFRic inputs.
  !>@details Nothing in this file is ever likely to be promoted to the LFRic
  !>          namelist. Everything is either set from an LFRic variable
  !>          already in the namelist, or is "fixed" from the perspective
  !>          of LFRic (but cannot be made a parameter because it is required
  !>          to be variable in the UM and therefore declared as such in the
  !>          UM modules which contains it).
  !> @param[in] ncells  The number of cells in the horizontal domain that
  !>                    the UM code should loop over (i.e. not including halos)
  subroutine um_sizes_init(ncells)

    ! External subroutines
    use atm_fields_bounds_mod, only: atm_fields_bounds_init

    ! Sizes of fields
    use atm_step_local, only: rhc_row_length, rhc_rows
    use nlsizes_namelist_mod, only: row_length, rows
    use theta_field_sizes, only: t_i_length, t_j_length, &
                                 u_i_length, u_j_length, &
                                 v_i_length, v_j_length
    use tuning_segments_mod, only: bl_segment_size, precip_segment_size, &
                                   ussp_seg_size, gw_seg_size, &
                                   conv_gr_segment_size
    use physics_config_mod,  only : ls_ppn_segment, gw_segment, &
                                    bl_segment, ussp_segment, &
                                    configure_segments, conv_gr_segment
    use log_mod, only : log_event, log_scratch_space, LOG_LEVEL_ERROR

    implicit none

    integer(i_def),   intent(in)          :: ncells

    ! ----------------------------------------------------------------
    ! Model dimensions - contained in UM module nlsizes_namelist_mod
    ! ----------------------------------------------------------------
    ! Horizontal dimensions set to the value passed into this routine.
    ! This needs to match the number of cells passed to physics kernels.
    row_length = int( ncells, i_um )
    rows       = 1

    ! ----------------------------------------------------------------
    ! More model dimensions, this time from atm_step_local
    ! ----------------------------------------------------------------
    ! Dimensions of critical relative humidity array. Again, needs to
    ! match number of points passed to kernels.
    rhc_row_length = row_length
    rhc_rows       = rows

    ! ----------------------------------------------------------------
    ! Segment sizes for UM physics - contained in tuning_segments_mod
    ! ----------------------------------------------------------------
    ! These are set to 1 currently because only 1 grid-cell is passed to
    ! a kernel. However, multiple columns are passed to a kernel,
    ! these values will need to be set depending on how many columns
    ! a kernel is passed.

    !keep the first
    if (configure_segments .and. row_length .gt. 1) then
      select case (ls_ppn_segment)
        case (:-1)
          write(log_scratch_space,'(A)') &
                'Invalid value: specified large scale precipitation segment is -ve.'
          call log_event(log_scratch_space, LOG_LEVEL_ERROR)

        case (1:)
          ! Set the value from the namelist
          precip_segment_size = ls_ppn_segment

        case default
          ! Default behaviour is to set to row_length
          precip_segment_size = row_length

      end select
      select case (bl_segment)
        case (:-1)
          write(log_scratch_space,'(A)') &
                'Invalid value: specified boundary layer segment size is -ve.'
          call log_event(log_scratch_space, LOG_LEVEL_ERROR)

        case (1:)
          ! Set the value from the namelist
          bl_segment_size = bl_segment

        case default
          ! Default behaviour is to set to row_length
          bl_segment_size = row_length

      end select
      select case (gw_segment)
        case (:-1)
          write(log_scratch_space,'(A)') &
                'Invalid value: specified gravity wave segment is -ve.'
          call log_event(log_scratch_space, LOG_LEVEL_ERROR)

        case (1:)
           ! Set the value from the namelist
           gw_seg_size = gw_segment

        case default
          ! Default behaviour is to set to row_length
          gw_seg_size = row_length

      end select
      select case (ussp_segment)
        case (:-1)
        write(log_scratch_space,'(A)') &
              'Invalid value: specified ussp segment is -ve.'
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)

        case (1:)
          ! Set the value from the namelist
          ussp_seg_size = ussp_segment

        case default
          ! Default behaviour is to set to row_length
          ussp_seg_size = row_length

      end select
      select case (conv_gr_segment)
        case (:-1)
          write(log_scratch_space,'(A)') &
                'Invalid value: specified gr_convection segment is -ve.'
          call log_event(log_scratch_space, LOG_LEVEL_ERROR)

        case (1:)
          ! Set the value from the namelist
          conv_gr_segment_size = conv_gr_segment

        case default
          ! Default behaviour is to set to row_length
          conv_gr_segment_size = row_length

      end select
    else
      ! Default behaviour is to set to row_length
      precip_segment_size = row_length
      bl_segment_size     = row_length
      gw_seg_size         = row_length
      ussp_seg_size       = row_length
      conv_gr_segment_size = row_length
    end if

    ! Compute lengths in i and j direction. This is the earliest place that they
    ! are needed. They will be kept in the module from here onward.
    t_i_length = row_length
    t_j_length = rows
    u_i_length = row_length
    u_j_length = rows
    v_i_length = row_length
    v_j_length = rows

    ! Set the field bounds which are used by the UM code based on the
    ! information above.
    ! Hard-wired zeros are halo-sizes in UM code.
    ! Currently set to zero as we don't pass a stencil into the kernels
    ! but may change if we ever do.
    call atm_fields_bounds_init( 0_i_um, 0_i_um, 0_i_um, &
                                 0_i_um, row_length, rows, rows)

  end subroutine um_sizes_init

end module um_sizes_init_mod
