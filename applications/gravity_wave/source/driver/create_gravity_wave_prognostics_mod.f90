!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Create the prognostic fields and place them in the depository

!> @details Create the prognostic fields and place them both in the
!>          depository field collection and put pointers to them in the
!>          prognostic field collection

module create_gravity_wave_prognostics_mod

  use constants_mod,                  only : i_def
  use field_mod,                      only : field_type
  use field_parent_mod,               only : write_interface, &
                                             checkpoint_write_interface, &
                                             checkpoint_read_interface
  use finite_element_config_mod,      only : element_order_h, element_order_v
  use function_space_collection_mod,  only : function_space_collection
  use fs_continuity_mod,              only : W0, W2, W3, Wtheta
  use io_config_mod,                  only : write_diag, &
                                             use_xios_io, &
                                             checkpoint_write, &
                                             checkpoint_read
  use io_mod,                         only : checkpoint_read_netcdf, &
                                             checkpoint_write_netcdf
  use lfric_xios_read_mod,            only : checkpoint_read_xios
  use lfric_xios_write_mod,           only : write_field_generic, &
                                             checkpoint_write_xios

  use gravity_wave_constants_config_mod, &
                                      only : b_space, &
                                             b_space_w0, &
                                             b_space_w3, &
                                             b_space_wtheta
  use log_mod,                        only : log_event, &
                                             LOG_LEVEL_INFO, &
                                             LOG_LEVEL_ERROR
  use mesh_mod,                       only : mesh_type
  use pure_abstract_field_mod,        only : pure_abstract_field_type

  implicit none

  private
  public create_gravity_wave_prognostics

contains

  !> @brief Create the prognostic fields and place them in the depository
  !> @param [in]    mesh The primary mesh
  !> @param [inout] wind prognostic field
  !> @param [inout] pressure prognostic field
  !> @param [inout] buoyancy prognostic field
  subroutine create_gravity_wave_prognostics(mesh, wind, pressure, buoyancy)

    implicit none
    type( mesh_type ), intent(in), pointer :: mesh

    type( field_type), intent(inout) :: wind
    type( field_type), intent(inout) :: pressure
    type( field_type), intent(inout) :: buoyancy

    integer(i_def) :: buoyancy_space

    procedure(write_interface),       pointer :: tmp_write_ptr
    procedure(checkpoint_write_interface), pointer :: tmp_checkpoint_write_ptr
    procedure(checkpoint_read_interface),  pointer :: tmp_checkpoint_read_ptr

    ! Create prognostic fields
    select case(b_space)
      case(b_space_w0)
        buoyancy_space = W0
        call log_event( 'gravity_wave: Using W0 for buoyancy', LOG_LEVEL_INFO )
      case(b_space_w3)
        buoyancy_space = W3
        call log_event( 'gravity_wave: Using W3 for buoyancy', LOG_LEVEL_INFO )
      case(b_space_wtheta)
        buoyancy_space = Wtheta
        call log_event( 'gravity_wave: Using Wtheta for buoyancy', LOG_LEVEL_INFO )
      case default
        call log_event( 'gravity_wave: Invalid buoyancy space', LOG_LEVEL_ERROR )
    end select

    call wind%initialise( vector_space =                                       &
                       function_space_collection%get_fs(mesh, element_order_h, &
                                                        element_order_v, W2),  &
                       name="wind" )
    call buoyancy%initialise( vector_space =                                      &
               function_space_collection%get_fs(mesh, element_order_h,            &
                                                element_order_v, buoyancy_space), &
               name="buoyancy" )
    call pressure%initialise( vector_space =                                   &
                       function_space_collection%get_fs(mesh, element_order_h, &
                                                        element_order_v, W3),  &
                       name="pressure" )


    ! Set I/O behaviours for diagnostic output
    if (write_diag .and. use_xios_io) then
       ! Fields that are output on the XIOS face domain
       tmp_write_ptr => write_field_generic
       call wind%set_write_behaviour(tmp_write_ptr)
       call pressure%set_write_behaviour(tmp_write_ptr)
       call buoyancy%set_write_behaviour(tmp_write_ptr)
    end if

    ! Set I/O behaviours for checkpoint / restart
    if ( checkpoint_write .or. checkpoint_read) then
      if (use_xios_io) then
        ! Use XIOS for checkpoint / restart
        tmp_checkpoint_write_ptr => checkpoint_write_xios
        tmp_checkpoint_read_ptr => checkpoint_read_xios
        call log_event( 'GungHo: Using XIOS for checkpointing...', LOG_LEVEL_INFO )
      else
        ! Use old checkpoint and restart methods
        tmp_checkpoint_write_ptr => checkpoint_write_netcdf
        tmp_checkpoint_read_ptr => checkpoint_read_netcdf
        call log_event( 'GungHo: Using NetCDF for checkpointing...', LOG_LEVEL_INFO )
      end if

      call wind%set_checkpoint_write_behaviour(tmp_checkpoint_write_ptr)
      call pressure%set_checkpoint_write_behaviour(tmp_checkpoint_write_ptr)
      call buoyancy%set_checkpoint_write_behaviour(tmp_checkpoint_write_ptr)

      call wind%set_checkpoint_read_behaviour(tmp_checkpoint_read_ptr)
      call pressure%set_checkpoint_read_behaviour(tmp_checkpoint_read_ptr)
      call buoyancy%set_checkpoint_read_behaviour(tmp_checkpoint_read_ptr)

    end if

  end subroutine create_gravity_wave_prognostics

end module create_gravity_wave_prognostics_mod
