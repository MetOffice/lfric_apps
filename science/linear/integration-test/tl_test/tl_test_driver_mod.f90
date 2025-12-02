!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>@brief   Drives the execution of the tangent linear model tests.
!>@details The tests are initialised and finalised using a similar
!!         method to gungho, but with the addition of the linearisation state.
module tl_test_driver_mod

  use base_mesh_config_mod,       only : prime_mesh_name
  use calendar_mod,               only : calendar_type
  use constants_mod,              only : i_def, imdi, r_def
  use extrusion_mod,              only : TWOD
  use gungho_model_mod,           only : initialise_infrastructure, &
                                         initialise_model,          &
                                         finalise_infrastructure,   &
                                         finalise_model
  use gungho_init_fields_mod,     only : create_model_data,     &
                                         initialise_model_data, &
                                         finalise_model_data
  use driver_modeldb_mod,         only : modeldb_type
  use gungho_time_axes_mod,       only : gungho_time_axes_type
  use section_choice_config_mod,  only : stochastic_physics, &
                                         stochastic_physics_um
  use io_value_mod,               only : io_value_type
  use io_context_mod,             only : io_context_type
  use log_mod,                    only : log_event,         &
                                         LOG_LEVEL_ALWAYS
  use mesh_mod,                   only : mesh_type
  use mesh_collection_mod,        only : mesh_collection
  use model_clock_mod,            only : model_clock_type
  use linear_model_data_mod,      only : linear_create_ls,  &
                                         linear_init_ls
  use tl_test_kinetic_energy_gradient_mod, only : test_kinetic_energy_gradient
  use tl_test_advect_density_field_mod,    only : test_advect_density_field
  use tl_test_advect_theta_field_mod,      only : test_advect_theta_field
  use tl_test_vorticity_mod,               only : test_vorticity_advection
  use tl_test_project_eos_pressure_mod,    only : test_project_eos_pressure
  use tl_test_sample_eos_pressure_mod,     only : test_sample_eos_pressure
  use tl_test_hydrostatic_mod,             only : test_hydrostatic
  use tl_test_pressure_grad_bd_mod,        only : test_pressure_gradient_bd
  use tl_test_rk_alg_mod,                  only : test_rk_alg
  use tl_test_transport_control_mod,       only : test_transport_control
  use tl_test_rhs_sample_eos_mod,          only : test_rhs_sample_eos
  use tl_test_rhs_project_eos_mod,         only : test_rhs_project_eos
  use tl_test_rhs_alg_mod,                 only : test_rhs_alg
  use tl_test_semi_imp_alg_mod,            only : test_semi_imp_alg
  use tl_test_timesteps_alg_mod,           only : test_timesteps

  implicit none

  private
  public initialise,                  &
         finalise,                    &
         run_timesteps,               &
         run_kinetic_energy_gradient, &
         run_advect_density_field,    &
         run_advect_theta_field,      &
         run_vorticity_advection,     &
         run_project_eos_pressure,    &
         run_sample_eos_pressure,     &
         run_hydrostatic,             &
         run_pressure_gradient_bd,    &
         run_rk_alg,                  &
         run_rhs_alg,                 &
         run_rhs_project_eos,         &
         run_rhs_sample_eos,          &
         run_semi_imp_alg,            &
         run_transport_control

  type(mesh_type), pointer :: mesh              => null()
  type(mesh_type), pointer :: twod_mesh         => null()
  type(mesh_type), pointer :: aerosol_mesh      => null()
  type(mesh_type), pointer :: aerosol_twod_mesh => null()

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Sets up the required state in preparation for run.
  !>@param [in]     program_name An identifier given to the model being run
  !>@param [in,out] modeldb      The structure that holds model state
  !>
  subroutine initialise( program_name, modeldb, calendar )

    implicit none

    character(*),         intent(in)    :: program_name
    type(modeldb_type),   intent(inout) :: modeldb
    class(calendar_type), intent(in)    :: calendar

    type(gungho_time_axes_type)     :: model_axes
    type(io_value_type)             :: temp_corr_io_value
    type(io_value_type)             :: random_seed_io_value
    integer(i_def) :: random_seed_size
    real(r_def), allocatable :: real_array(:)

    call modeldb%values%initialise( 'values', 5 )

    ! Initialise infrastructure and setup constants
    !
    call initialise_infrastructure( program_name, modeldb )

    ! Add a place to store time axes in modeldb
    call modeldb%values%add_key_value('model_axes', model_axes)

    ! Get primary and 2D meshes for initialising model data
    mesh => mesh_collection%get_mesh(prime_mesh_name)
    twod_mesh => mesh_collection%get_mesh(mesh, TWOD)

    ! Assume aerosol mesh is the same as dynamics mesh
    aerosol_mesh => mesh
    aerosol_twod_mesh => twod_mesh

    ! gungho_init_field() expects these values to exist. The dependency of
    ! the linear application tests on this procedure will hopefully be resolved
    ! in the future, at which point this initialisation may be removed.
    !
    call temp_corr_io_value%init('temperature_correction_rate', [0.0_r_def])
    call modeldb%values%add_key_value( 'temperature_correction_io_value', &
                                       temp_corr_io_value )
    call modeldb%values%add_key_value( 'total_dry_mass', 0.0_r_def )
    call modeldb%values%add_key_value( 'total_energy', 0.0_r_def )
    call modeldb%values%add_key_value( 'total_energy_previous', 0.0_r_def )
    if ( stochastic_physics == stochastic_physics_um ) then
      ! Random seed for stochastic physics
      call random_seed(size = random_seed_size)
      allocate(real_array(random_seed_size))
      real_array(1:random_seed_size) = 0.0_r_def
      call random_seed_io_value%init("random_seed", real_array)
      call modeldb%values%add_key_value( 'random_seed_io_value', &
                                         random_seed_io_value )
      deallocate(real_array)
    end if

    ! Instantiate the fields stored in model_data
    call create_model_data( modeldb,      &
                            mesh,         &
                            twod_mesh,    &
                            aerosol_mesh, &
                            aerosol_twod_mesh )

    ! Instantiate the linearisation state
    call linear_create_ls( modeldb, mesh )

    ! Initialise the fields stored in the model_data prognostics. This needs
    ! to be done before initialise_model.
    call initialise_model_data( modeldb, mesh, twod_mesh )

    ! Model configuration initialisation
    call initialise_model( mesh,  &
                           modeldb )

    ! Initialise the linearisation state
    call linear_init_ls( mesh, twod_mesh, modeldb )

    ! Finalise model
    call finalise_model(modeldb)

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tests the tangent linear model for multiple timesteps
  !>@param [in,out] modeldb   The structure that holds model state
  subroutine run_timesteps(modeldb)

    implicit none

    type(modeldb_type), intent(inout) :: modeldb

    call test_timesteps( modeldb,   &
                         mesh,      &
                         twod_mesh, &
                         modeldb%clock )

  end subroutine run_timesteps

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tests the tangent linear model kinetic energy gradient kernel
  !>@param [in,out] modeldb   The structure that holds model state
  subroutine run_kinetic_energy_gradient(modeldb)

    implicit none

    type(modeldb_type), intent(inout) :: modeldb

    call test_kinetic_energy_gradient( modeldb, &
                                       mesh,    &
                                       twod_mesh )

  end subroutine run_kinetic_energy_gradient

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tests the tangent linear model for density advection kernel
  !>@param [in,out] modeldb   The structure that holds model state
  subroutine run_advect_density_field(modeldb)

    implicit none

    type(modeldb_type), intent(inout) :: modeldb

    call test_advect_density_field( modeldb, &
                                    mesh,    &
                                    twod_mesh )

  end subroutine run_advect_density_field

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tests the tangent linear model theta advection kernel
  !>@param [in,out] modeldb   The structure that holds model state
  subroutine run_advect_theta_field(modeldb)

    implicit none

    type(modeldb_type), intent(inout) :: modeldb

    call test_advect_theta_field( modeldb, &
                                  mesh,    &
                                  twod_mesh )

  end subroutine run_advect_theta_field

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tests the tangent linear model vorticity advection kernel
  !>@param [in,out] modeldb   The structure that holds model state
  subroutine run_vorticity_advection(modeldb)

    implicit none

    type(modeldb_type), intent(inout) :: modeldb

    call test_vorticity_advection( modeldb, &
                                   mesh,    &
                                   twod_mesh )

  end subroutine run_vorticity_advection

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tests the tangent linear model project pressure kernel
  !>@param [in,out] modeldb   The structure that holds model state
  subroutine run_project_eos_pressure(modeldb)

    implicit none

    type(modeldb_type), intent(inout) :: modeldb

    call test_project_eos_pressure( modeldb, &
                                    mesh,    &
                                    twod_mesh )

  end subroutine run_project_eos_pressure

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tests the tangent linear model sample pressure kernel
  !>@param [in,out] modeldb   The structure that holds model state
  subroutine run_sample_eos_pressure(modeldb)

    implicit none

    type(modeldb_type), intent(inout) :: modeldb

    call test_sample_eos_pressure( modeldb, &
                                   mesh,    &
                                   twod_mesh )

  end subroutine run_sample_eos_pressure

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tests the tangent linear model hydrostatic kernel
  !>@param [in,out] modeldb   The structure that holds model state
  subroutine run_hydrostatic(modeldb)

    implicit none

    type(modeldb_type), intent(inout) :: modeldb

    call test_hydrostatic( modeldb, &
                           mesh,    &
                           twod_mesh )

  end subroutine run_hydrostatic

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tests the tangent linear model pressure gradient bd kernel
  !>@param [in,out] modeldb   The structure that holds model state
  subroutine run_pressure_gradient_bd(modeldb)

    implicit none

    type(modeldb_type), intent(inout) :: modeldb

    call test_pressure_gradient_bd( modeldb, &
                                    mesh,    &
                                    twod_mesh )

  end subroutine run_pressure_gradient_bd

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tests the tangent linear model for runge kutta timestepping
  !>@param [in,out] modeldb   The structure that holds model state
  subroutine run_rk_alg(modeldb)

    implicit none

    type(modeldb_type), intent(inout) :: modeldb

    call test_rk_alg( modeldb, &
                      mesh)

  end subroutine run_rk_alg

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tests the tangent linear model transport control routine
  !>@param [in,out] modeldb   The structure that holds model state
  subroutine run_transport_control(modeldb)

    implicit none

    type(modeldb_type), intent(inout) :: modeldb

    call test_transport_control( modeldb, &
                                 mesh,    &
                                 twod_mesh )

  end subroutine run_transport_control

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tests the tangent linear model for semi-implicit timestepping
  !>@param [in,out] modeldb   The structure that holds model state
  subroutine run_semi_imp_alg(modeldb)

    implicit none

    type(modeldb_type), intent(inout) :: modeldb

    call test_semi_imp_alg( modeldb,   &
                            mesh,      &
                            twod_mesh )

  end subroutine run_semi_imp_alg

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tests the tangent linear model right-hand side
  !>@param [in,out] modeldb   The structure that holds model state
  subroutine run_rhs_alg(modeldb)

    implicit none

    type(modeldb_type), intent(inout) :: modeldb

    call test_rhs_alg( modeldb, &
                       mesh,    &
                       twod_mesh )

  end subroutine run_rhs_alg

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tests the tangent linear model for the project RHS EoS
  !>@param [in,out] modeldb   The structure that holds model state
  subroutine run_rhs_project_eos(modeldb)

    implicit none

    type(modeldb_type), intent(inout) :: modeldb

    call test_rhs_project_eos( modeldb, &
                               mesh,    &
                               twod_mesh )

  end subroutine run_rhs_project_eos

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tests the tangent linear model for the sample RHS EoS
  !>@param [in,out] modeldb   The structure that holds model state
  subroutine run_rhs_sample_eos(modeldb)

    implicit none

    type(modeldb_type), intent(inout) :: modeldb

    call test_rhs_sample_eos( modeldb, &
                              mesh,    &
                              twod_mesh )

  end subroutine run_rhs_sample_eos

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tidies up after a run.
  !>@param [in]     program_name An identifier given to the model being run
  !>@param [in,out] modeldb      The structure that holds model state
  subroutine finalise( program_name, modeldb )

    implicit none

    character(*),       intent(in)    :: program_name
    type(modeldb_type), intent(inout) :: modeldb

    call log_event( 'Finalising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    ! Destroy the fields stored in model_data
    call finalise_model_data( modeldb )

    ! Finalise infrastructure and constants
    call finalise_infrastructure(modeldb)

  end subroutine finalise

end module tl_test_driver_mod
