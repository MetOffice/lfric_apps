!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Creates fields that will hold data read from the tangent-linear start
!!        dump file.
module create_tl_prognostics_mod

  use field_mod,                      only : field_type
  use field_parent_mod,               only : read_interface, &
                                             write_interface
  use lfric_xios_read_mod,            only : read_field_generic
  use lfric_xios_write_mod,           only : write_field_generic
  use finite_element_config_mod,      only : element_order_h, element_order_v
  use function_space_collection_mod,  only : function_space_collection
  use field_collection_mod,           only : field_collection_type
  use fs_continuity_mod,              only : W3, Wtheta, W2H
  use log_mod,                        only : log_event,      &
                                             LOG_LEVEL_INFO, &
                                             LOG_LEVEL_ERROR
  use mesh_mod,                       only : mesh_type
  use initialization_config_mod,      only : read_w2h_wind

  implicit none

  private
  public :: create_tl_prognostics

contains

  !> @brief Initialise the fd_field collection.
  !> @details Add the required fields with appropriate read and write behaviour,
  !!        corresponding to the fields in the start dump.
  !> @param[in]     mesh                The mesh
  !> @param[in]     twod_mesh           The current 2d mesh
  !> @param[in,out] fd_field_collection The collection object to store fields in
  !> @param[in,out] depository          The depository field collection
  subroutine create_tl_prognostics( mesh, twod_mesh, fd_field_collection, &
                                    depository )

    implicit none

    type( mesh_type ), intent(in), pointer     :: mesh
    type( mesh_type ), intent(in), pointer     :: twod_mesh

    type(field_collection_type), intent(inout) :: fd_field_collection
    type(field_collection_type), intent(inout) :: depository

    procedure(read_interface),  pointer :: tmp_read_ptr  => null()
    procedure(write_interface), pointer :: tmp_write_ptr => null()

    ! FD field declarations
    type( field_type ) :: ew_wind_in_w3 ! U wind
    type( field_type ) :: ns_wind_in_w3 ! V wind
    type( field_type ) :: h_wind_in_w2h ! Horizontal wind (i.e. on W2h dofs)
    type( field_type ) :: dry_rho_in_w3 ! Dry rho
    type( field_type ) :: exner_in_w3
    ! Vertical theta levels
    type( field_type ) :: upward_wind_in_wtheta ! W wind
    type( field_type ) :: theta_in_wtheta ! Potential temp
    type( field_type ) :: mv_in_wtheta    ! Vapour mix ratio
    type( field_type ) :: mcl_in_wtheta   ! Cloud liquid mix ratio
    type( field_type ) :: mcf_in_wtheta   ! Clould ice mix ratio
    type( field_type ) :: mr_in_wtheta    ! Rain mix ratio

    call log_event( 'Creating Tangent linear initial perturbation prognostics...', &
                    LOG_LEVEL_INFO )

    if ( element_order_h > 0 .or. element_order_v > 0) then
      call log_event( 'Finite diff prognostics: requires lowest order elements', &
                       LOG_LEVEL_ERROR )
    end if

    ! Field collection created in linear_model.f90
    ! Create the field collection

    tmp_read_ptr => read_field_generic
    tmp_write_ptr => write_field_generic
    if ( read_w2h_wind ) then

       ! In this case we read in directly onto the W2H dofs
       call h_wind_in_w2h%initialise( vector_space =                   &
         function_space_collection%get_fs( mesh, element_order_h,      &
                                           element_order_v, W2H ), &
         name='h_u')
       call h_wind_in_w2h%set_read_behaviour(tmp_read_ptr)
       call h_wind_in_w2h%set_write_behaviour(tmp_write_ptr)
       call fd_field_collection%add_field(h_wind_in_w2h)

    else

      ! Create the fields, set the I/O behaviour and add to
      ! the field collection
      !========================================================================
      ! W3 fields - rho levels
      !========================================================================

      call ew_wind_in_w3%initialise( vector_space =                  &
        function_space_collection%get_fs( mesh, element_order_h,     &
                                          element_order_v, W3 ),     &
        name='ew_wind_in_w3')
      call ew_wind_in_w3%set_read_behaviour(tmp_read_ptr)
      call ew_wind_in_w3%set_write_behaviour(tmp_write_ptr)
      call fd_field_collection%add_field(ew_wind_in_w3)

      call ns_wind_in_w3%initialise( vector_space =                  &
        function_space_collection%get_fs( mesh, element_order_h,     &
                                          element_order_v, W3 ),     &
        name='ns_wind_in_w3')
      call ns_wind_in_w3%set_read_behaviour(tmp_read_ptr)
      call ns_wind_in_w3%set_write_behaviour(tmp_write_ptr)
      call fd_field_collection%add_field(ns_wind_in_w3)

    end if

    call dry_rho_in_w3%initialise( vector_space =                     &
         function_space_collection%get_fs( mesh, element_order_h,     &
                                           element_order_v, W3 ),     &
         name='rho')
    call dry_rho_in_w3%set_read_behaviour(tmp_read_ptr)
    call dry_rho_in_w3%set_write_behaviour(tmp_write_ptr)
    call fd_field_collection%add_field(dry_rho_in_w3)

    call exner_in_w3%initialise( vector_space =                       &
         function_space_collection%get_fs( mesh, element_order_h,     &
                                           element_order_v, W3 ),     &
         name='exner')
    call exner_in_w3%set_read_behaviour(tmp_read_ptr)
    call exner_in_w3%set_write_behaviour(tmp_write_ptr)
    call fd_field_collection%add_field(exner_in_w3)

    !========================================================================
    ! Wtheta fields - theta levels
    !========================================================================

    call upward_wind_in_wtheta%initialise( vector_space =                 &
         function_space_collection%get_fs( mesh, element_order_h,         &
                                           element_order_v, Wtheta ),     &
         name='v_u')
    call upward_wind_in_wtheta%set_read_behaviour(tmp_read_ptr)
    call upward_wind_in_wtheta%set_write_behaviour(tmp_write_ptr)
    call fd_field_collection%add_field(upward_wind_in_wtheta)

    call theta_in_wtheta%initialise( vector_space =                       &
         function_space_collection%get_fs( mesh, element_order_h,         &
                                           element_order_v, Wtheta ),     &
         name='theta')
    call theta_in_wtheta%set_read_behaviour(tmp_read_ptr)
    call theta_in_wtheta%set_write_behaviour(tmp_write_ptr)
    call fd_field_collection%add_field(theta_in_wtheta)

    call mv_in_wtheta%initialise( vector_space =                          &
         function_space_collection%get_fs( mesh, element_order_h,         &
                                           element_order_v, Wtheta ),     &
         name='m_v')
    call mv_in_wtheta%set_read_behaviour(tmp_read_ptr)
    call mv_in_wtheta%set_write_behaviour(tmp_write_ptr)
    call fd_field_collection%add_field(mv_in_wtheta)

    call mcl_in_wtheta%initialise( vector_space =                         &
         function_space_collection%get_fs( mesh, element_order_h,         &
                                           element_order_v, Wtheta ),     &
         name='m_cl')
    call mcl_in_wtheta%set_read_behaviour(tmp_read_ptr)
    call mcl_in_wtheta%set_write_behaviour(tmp_write_ptr)
    call fd_field_collection%add_field(mcl_in_wtheta)

    call mcf_in_wtheta%initialise( vector_space =                         &
         function_space_collection%get_fs( mesh, element_order_h,         &
                                           element_order_v, Wtheta ),     &
         name='m_s')
    call mcf_in_wtheta%set_read_behaviour(tmp_read_ptr)
    call mcf_in_wtheta%set_write_behaviour(tmp_write_ptr)
    call fd_field_collection%add_field(mcf_in_wtheta)

    call mr_in_wtheta%initialise( vector_space =                          &
         function_space_collection%get_fs( mesh, element_order_h,         &
                                           element_order_v, Wtheta ),     &
         name='m_r')
    call mr_in_wtheta%set_read_behaviour(tmp_read_ptr)
    call mr_in_wtheta%set_write_behaviour(tmp_write_ptr)
    call fd_field_collection%add_field(mr_in_wtheta)

  end subroutine create_tl_prognostics

end module create_tl_prognostics_mod
