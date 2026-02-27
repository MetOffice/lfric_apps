!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Initialise, define and finalise the linearisation state.

module linear_model_data_mod

  use constants_mod,                  only : i_def, r_def, l_def, str_def
  use pure_abstract_field_mod,        only : pure_abstract_field_type
  use field_array_mod,                only : field_array_type
  use field_mod,                      only : field_type
  use field_collection_mod,           only : field_collection_type
  use function_space_mod,             only : function_space_type
  use function_space_collection_mod,  only : function_space_collection
  use fs_continuity_mod,              only : W2, W3, WTheta, W2h
  use driver_modeldb_mod,             only : modeldb_type
  use gungho_time_axes_mod,           only : gungho_time_axes_type, &
                                             get_time_axes_from_collection
  use init_time_axis_mod,             only : setup_field
  use initialization_config_mod,      only : ls_option,           &
                                             ls_option_analytic,  &
                                             ls_option_file, &
                                             zero_w2v_wind
  use lfric_xios_time_axis_mod,       only : time_axis_type
  use lfric_xios_read_mod,            only : read_field_time_var
  use linear_data_algorithm_mod,      only : linear_copy_model_to_ls,  &
                                             linear_init_pert_random,  &
                                             init_ls_file_alg,         &
                                             linear_init_reference_ls, &
                                             linear_init_pert_analytical, &
                                             init_ls_file_alg,            &
                                             linear_init_pert_zero

  use linear_config_mod,              only : pert_option,          &
                                             pert_option_analytic, &
                                             pert_option_random,   &
                                             pert_option_file,     &
                                             pert_option_zero,     &
                                             ls_read_w2h
  use linked_list_mod,                only : linked_list_type
  use log_mod,                        only : log_event,         &
                                             log_scratch_space, &
                                             LOG_LEVEL_INFO,    &
                                             LOG_LEVEL_ERROR
  use mesh_mod,                       only : mesh_type
  use moist_dyn_mod,                  only : num_moist_factors
  use moist_dyn_factors_alg_mod,      only : moist_dyn_factors_alg
  use mr_indices_mod,                 only : nummr, &
                                             mr_names
  use linear_map_fd_alg_mod,          only : linear_map_fd_to_prognostics
  use set_any_dof_alg_mod,            only : set_any_dof_alg
  use reference_element_mod,          only : T
  use io_config_mod,                  only : checkpoint_read

  implicit none

  private
  public :: linear_create_ls,          &
            linear_create_ls_analytic, &
            linear_init_ls,            &
            linear_init_pert

contains

  !> @brief   Create the fields in the ls fields field collection.
  !> @details The fields in the ls fields field collection are setup
  !>          for use with a file update or analytical initialisation.
  !>
  !> @param[inout] modeldb   The working data set for a model run
  !> @param[in]    mesh      The current 3d mesh
  !>
  subroutine linear_create_ls( modeldb, mesh )

    implicit none

    type( modeldb_type ), target, intent(inout) :: modeldb

    type( mesh_type ), pointer, intent(in) :: mesh

    select case( ls_option )

      case( ls_option_analytic )

        call linear_create_ls_analytic( modeldb, mesh )

      case( ls_option_file )

        call linear_create_ls_file( modeldb, mesh )

      case default

        call log_event( "LS setup not available for requested ls_option ", &
                        LOG_LEVEL_ERROR)

    end select

  end subroutine linear_create_ls

  !> @brief   Create the fields in the ls fields field collection ready to be
  !>          initialised.
  !> @details The ls fields field collection is created ready to be updated
  !>          in the analytical test case. This method is also used by jedi-lfric
  !>
  !> @param[inout] modeldb   The working data set for a model run
  !> @param[in]    mesh      The current 3d mesh
  !>
  subroutine linear_create_ls_analytic( modeldb, mesh )

    implicit none

    type( modeldb_type ), target, intent(inout) :: modeldb

    type( mesh_type ), pointer, intent(in) :: mesh

    type( field_collection_type ), pointer :: depository
    type( field_collection_type ), pointer :: prognostics
    type( field_collection_type ), pointer :: ls_fields
    type( field_type ),            pointer :: ls_mr(:)
    type( field_type ),            pointer :: ls_moist_dyn(:)

    type(field_collection_type), pointer :: moisture_fields
    type(field_array_type),      pointer :: ls_mr_array
    type(field_array_type),      pointer :: ls_moist_dyn_array

    integer(i_def)            :: imr
    character(str_def)        :: name
    character(str_def)        :: moist_dyn_name
    logical(l_def), parameter :: checkpoint_restart_flag = .true.

    depository => modeldb%fields%get_field_collection("depository")
    prognostics => modeldb%fields%get_field_collection("prognostic_fields")

    moisture_fields => modeldb%fields%get_field_collection("moisture_fields")
    call moisture_fields%get_field("ls_mr",ls_mr_array)
    call moisture_fields%get_field("ls_moist_dyn", ls_moist_dyn_array)
    ls_mr => ls_mr_array%bundle
    ls_moist_dyn => ls_moist_dyn_array%bundle

    write(log_scratch_space,'(A,A)') "Create ls fields: "// &
                                     "Setting up ls field collection"
    call log_event(log_scratch_space, LOG_LEVEL_INFO)

    call modeldb%fields%add_empty_field_collection("ls_fields", table_len = 100)

    ls_fields => modeldb%fields%get_field_collection("ls_fields")


    call setup_field( ls_fields, depository, prognostics, "ls_rho", W3,       &
                      mesh, checkpoint_restart_flag )
    call setup_field( ls_fields, depository, prognostics, "ls_exner", W3,     &
                      mesh, checkpoint_restart_flag )
    call setup_field( ls_fields, depository, prognostics, "ls_theta", Wtheta, &
                      mesh, checkpoint_restart_flag )
    call setup_field( ls_fields, depository, prognostics, "ls_u", W2,         &
                       mesh, checkpoint_restart_flag )

    do imr = 1, nummr
      name = trim('ls_' // adjustl(mr_names(imr)) )
      call setup_field( ls_fields, depository, prognostics, name, Wtheta, &
                        mesh, checkpoint_restart_flag, mr=ls_mr, imr=imr )
    end do

    do imr = 1, num_moist_factors
      write(moist_dyn_name, "(A12, I1)") "ls_moist_dyn", imr
      name = trim(moist_dyn_name)
      call setup_field( ls_fields, depository, prognostics, name, Wtheta, &
                        mesh, checkpoint_restart_flag, mr=ls_moist_dyn,   &
                        imr=imr )
    end do

  end subroutine linear_create_ls_analytic

  !> @brief   Create the fields in the ls fields field collection to be setup
  !>          and updated by file.
  !> @details The ls fields field collection is created and is setup ready
  !>          for update via a file read with a time axis.
  !>
  !> @param[inout] modeldb   The working data set for a model run
  !> @param[in]    mesh      The current 3d mesh
  !>
  subroutine linear_create_ls_file( modeldb, mesh )

    implicit none

    type( modeldb_type ), target, intent(inout) :: modeldb

    type( mesh_type ), pointer, intent(in) :: mesh

    type( field_collection_type ), pointer :: depository
    type( field_collection_type ), pointer :: prognostics
    type( field_collection_type ), pointer :: ls_fields
    type( field_type ),            pointer :: ls_mr(:)
    type( field_type ),            pointer :: ls_moist_dyn(:)
    type( linked_list_type ),      pointer :: ls_times_list

    type(field_collection_type), pointer :: moisture_fields
    type(field_array_type),      pointer :: ls_mr_array
    type(field_array_type),      pointer :: ls_moist_dyn_array

    type(gungho_time_axes_type), pointer :: model_axes

    integer(i_def)     :: imr
    character(str_def) :: name
    character(str_def) :: moist_dyn_name

    type(time_axis_type), save  :: ls_time_axis
    logical(l_def),   parameter :: checkpoint_restart_flag = .false.
    logical(l_def),   parameter :: cyclic=.false.
    logical(l_def),   parameter :: interp_flag=.true.
    character(len=*), parameter :: axis_id="ls_axis"

    depository => modeldb%fields%get_field_collection("depository")
    prognostics => modeldb%fields%get_field_collection("prognostic_fields")

    ! Get model_axes out of modeldb
    model_axes => get_time_axes_from_collection(modeldb%values, "model_axes" )

    ls_times_list => model_axes%ls_times_list

    call modeldb%fields%add_empty_field_collection("ls_fields", table_len = 100)
    ls_fields => modeldb%fields%get_field_collection("ls_fields")

    moisture_fields => modeldb%fields%get_field_collection("moisture_fields")
    call moisture_fields%get_field("ls_mr",ls_mr_array)
    call moisture_fields%get_field("ls_moist_dyn", ls_moist_dyn_array)
    ls_mr => ls_mr_array%bundle
    ls_moist_dyn => ls_moist_dyn_array%bundle

    write(log_scratch_space,'(A,A)') "Create ls fields: "// &
                                     "Setting up ls field collection"
    call log_event(log_scratch_space, LOG_LEVEL_INFO)

    call ls_time_axis%initialise( "ls_time", file_id="ls", &
                                   yearly=cyclic,          &
                                   interp_flag = interp_flag )

    call setup_field( ls_fields, depository, prognostics, "ls_rho", W3,       &
                      mesh, checkpoint_restart_flag, time_axis=ls_time_axis )
    call setup_field( ls_fields, depository, prognostics, "ls_exner", W3,     &
                      mesh, checkpoint_restart_flag, time_axis=ls_time_axis )
    call setup_field( ls_fields, depository, prognostics, "ls_theta", Wtheta, &
                      mesh, checkpoint_restart_flag, time_axis=ls_time_axis )

    if ( ls_read_w2h ) then
      call setup_field( ls_fields, depository, prognostics, "ls_h_u", W2h,      &
                        mesh, checkpoint_restart_flag, time_axis=ls_time_axis )
      call setup_field( ls_fields, depository, prognostics, "ls_v_u", Wtheta,   &
                        mesh, checkpoint_restart_flag, time_axis=ls_time_axis )
    else
      call setup_field( ls_fields, depository, prognostics, "ls_u_in_w3", W3,      &
                        mesh, checkpoint_restart_flag, time_axis=ls_time_axis )
      call setup_field( ls_fields, depository, prognostics, "ls_v_in_w3", W3,      &
                        mesh, checkpoint_restart_flag, time_axis=ls_time_axis )
      call setup_field( ls_fields, depository, prognostics, "ls_w_in_wth", Wtheta, &
                        mesh, checkpoint_restart_flag, time_axis=ls_time_axis )
    end if

    call setup_field( ls_fields, depository, prognostics, "ls_u", W2,         &
                      mesh, checkpoint_restart_flag )

    do imr = 1, nummr-2

      name = trim('ls_' // adjustl(mr_names(imr)) )

      call setup_field( ls_fields, depository, prognostics, name, Wtheta, &
                        mesh, checkpoint_restart_flag, time_axis=ls_time_axis )
    enddo

    ! m_g and m_s
    do imr = nummr-1, nummr

      name = trim('ls_' // adjustl(mr_names(imr)) )

      call setup_field( ls_fields, depository, prognostics, name, Wtheta, &
                        mesh, checkpoint_restart_flag, mr=ls_mr, imr=imr )
    end do

    do imr = 1, num_moist_factors

      write(moist_dyn_name, "(A12, I1)") "ls_moist_dyn", imr
      name = trim(moist_dyn_name)

      call setup_field( ls_fields, depository, prognostics, name, Wtheta, &
                        mesh, checkpoint_restart_flag,                    &
                        mr=ls_moist_dyn, imr=imr )
    end do

    call ls_times_list%insert_item(ls_time_axis)

  end subroutine linear_create_ls_file

  !> @brief   Define the linearisation state values.
  !> @details At the present, these can only be defined from an analytical
  !!          solution.
  !> @param[in]    mesh        The current 3d mesh
  !> @param[in]    twod_mesh   The current 2d mesh
  !> @param[inout] modeldb     The working data set for a model run
  !>
  subroutine linear_init_ls( mesh, twod_mesh, modeldb )

    use gungho_step_mod,                only : gungho_step
    use sci_field_minmax_alg_mod,       only : log_field_minmax

    implicit none

    type( mesh_type ), pointer, intent(in) :: mesh
    type( mesh_type ), pointer, intent(in) :: twod_mesh

    type( modeldb_type ), target, intent(inout) :: modeldb

    integer(i_def)              :: i
    type( field_type ), pointer :: ls_field => null()
    type( field_collection_type ), pointer :: ls_fields
    integer(i_def), parameter   :: number_steps = 10

    type(field_collection_type), pointer :: moisture_fields => null()
    type(field_array_type), pointer      :: ls_mr_array => null()
    type(field_array_type), pointer      :: ls_moist_dyn_array => null()

    type(gungho_time_axes_type), pointer :: model_axes

    ! Get model_axes out of modeldb
    model_axes => get_time_axes_from_collection(modeldb%values, "model_axes" )

    ls_fields => modeldb%fields%get_field_collection("ls_fields")

    moisture_fields => modeldb%fields%get_field_collection("moisture_fields")
    call moisture_fields%get_field("ls_mr", ls_mr_array)
    call moisture_fields%get_field("ls_moist_dyn", ls_moist_dyn_array)

    select case( ls_option )

      case( ls_option_analytic )

        select case( pert_option )

          case( pert_option_random )

            ! Procedure to define the linearisation state from an analytical
            ! field
            ! 1. Define the analytical field in the gungho prognostics. (This
            !    is done in initialise_model_data in tl_test_driver).
            ! 2. Evolve the prognostic fields using gungho_step -
            !    this avoids the linearisation state being zero.
            ! 3. Copy the prognostic fields to the linearisation fields, and
            !    set the prognostic fields to zero.

            ! Evolve the prognostic fields.
            do i = 1, number_steps
              call gungho_step( mesh,      &
                                twod_mesh, &
                                modeldb,   &
                                modeldb%clock )
            end do

            ! Copy the prognostic fields to the LS and then zero the prognostics.
            call linear_copy_model_to_ls( modeldb )

        case( pert_option_analytic )

          call linear_init_reference_ls( modeldb )

        case( pert_option_file, pert_option_zero )
          call log_event("This pert_option not available with ls_option_analytic ", LOG_LEVEL_ERROR)

        case default

          call log_event("This pert_option not available", LOG_LEVEL_ERROR)

        end select

      case( ls_option_file )

        call init_ls_file_alg( model_axes%ls_times_list, &
                               modeldb%clock,                    &
                               ls_fields,                        &
                               ls_mr_array%bundle,               &
                               ls_moist_dyn_array%bundle )

      case default

        call log_event("This ls_option not available", LOG_LEVEL_ERROR)

    end select

    ! Print the min and max values of the linearisation fields.
    call ls_fields%get_field("ls_u", ls_field)
    call log_field_minmax( LOG_LEVEL_INFO, 'ls_u', ls_field )

    call ls_fields%get_field("ls_rho", ls_field)
    call log_field_minmax( LOG_LEVEL_INFO, 'ls_rho', ls_field )

    call ls_fields%get_field("ls_exner", ls_field)
    call log_field_minmax( LOG_LEVEL_INFO, 'ls_exner', ls_field )

    call ls_fields%get_field("ls_theta", ls_field)
    call log_field_minmax( LOG_LEVEL_INFO, 'ls_theta', ls_field )

    ls_field => ls_mr_array%bundle(1)
    call log_field_minmax( LOG_LEVEL_INFO, 'ls_mr', ls_field )

  end subroutine linear_init_ls

  !> @brief   Define the initial perturbation values.
  !> @details Define the initial perturbation - currently from random data
  !> @param[in]    mesh      The current 3d mesh
  !> @param[in]    twod_mesh The current 2d mesh
  !> @param[inout] modeldb   The working data set for a model run
  subroutine linear_init_pert( mesh, twod_mesh, modeldb )

    implicit none

    type( mesh_type ), pointer, intent(in) :: mesh
    type( mesh_type ), pointer, intent(in) :: twod_mesh

    type( modeldb_type ), target, intent(inout) :: modeldb
    type(field_collection_type), pointer :: prognostic_fields
    type( field_type ),          pointer :: u

    select case( pert_option )

      case( pert_option_random )

        call linear_init_pert_random( modeldb )

      case( pert_option_analytic )

        call linear_init_pert_analytical( mesh,      &
                                          twod_mesh, &
                                          modeldb )

      case( pert_option_file )

        call linear_map_fd_to_prognostics( modeldb )

      case( pert_option_zero )

        call linear_init_pert_zero( modeldb )

      case default

        call log_event("This pert_option not available", LOG_LEVEL_ERROR)

    end select

    if (.not. checkpoint_read .and. zero_w2v_wind) then
      prognostic_fields => modeldb%fields%get_field_collection(&
                                          "prognostic_fields")
      call prognostic_fields%get_field('u', u)
      call set_any_dof_alg(u, T, 0.0_r_def)
    end if

  end subroutine linear_init_pert

end module linear_model_data_mod