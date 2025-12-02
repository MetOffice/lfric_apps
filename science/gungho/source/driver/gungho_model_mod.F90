!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Handles initialisation and finalisation of infrastructure, constants
!>        and the gungho model simulations.
!>
module gungho_model_mod

  use add_mesh_map_mod,           only : assign_mesh_maps
  use sci_checksum_alg_mod,       only : checksum_alg
  use driver_fem_mod,             only : init_fem, final_fem, &
                                         init_function_space_chains
  use driver_io_mod,              only : init_io, final_io, &
                                         filelist_populator
  use driver_mesh_mod,            only : init_mesh
  use check_configuration_mod,    only : get_required_stencil_depth,    &
                                         check_any_shifted,             &
                                         check_any_wt_eqn_conservative, &
                                         check_configuration
  use energy_correction_config_mod, only : encorr_usage, encorr_usage_none
  use conservation_algorithm_mod, only : conservation_algorithm
  use constants_mod,              only : i_def, r_def, l_def, &
                                         PRECISION_REAL, r_second, str_def
  use convert_to_upper_mod,       only : convert_to_upper
  use create_gungho_prognostics_mod, only : process_gungho_prognostics
  use create_lbcs_mod,            only : process_lbc_fields
  use create_mesh_mod,            only : create_mesh
  use create_physics_prognostics_mod, only : &
                                            process_physics_prognostics
  use derived_config_mod,         only : set_derived_config, l_esm_couple, &
                                         l_couple_sea_ice, l_couple_ocean
  use extrusion_mod,              only : extrusion_type,              &
                                         uniform_extrusion_type,      &
                                         shifted_extrusion_type,      &
                                         double_level_extrusion_type, &
                                         TWOD, SHIFTED, DOUBLE_LEVEL
  use field_array_mod,            only : field_array_type
  use field_mod,                  only : field_type
  use field_spec_mod,             only : field_spec_type, processor_type
  use field_parent_mod,           only : write_interface
  use field_collection_mod,       only : field_collection_type
  use field_spec_mod,             only : field_spec_type, processor_type, &
                                         space_has_xios_io
  use lfric_xios_diag_mod,        only : set_variable
  use sci_geometric_constants_mod,       &
                                  only : get_chi_inventory, get_panel_id_inventory
  use gungho_extrusion_mod,       only : create_extrusion
  use driver_modeldb_mod,         only : modeldb_type
  use gungho_setup_io_mod,        only : init_gungho_files
  use init_altitude_mod,          only : init_altitude
  use inventory_by_mesh_mod,      only : inventory_by_mesh_type
  use lfric_string_mod,           only : split_string
  use lfric_xios_context_mod,     only : lfric_xios_context_type
  use lfric_xios_metafile_mod,    only : metafile_type, add_field, &
                                         CHECKPOINTING, RESTARTING
  use lfric_xios_write_mod,       only : create_checkpoint_list
  use linked_list_mod,            only : linked_list_type
  use log_mod,                    only : log_event,          &
                                         log_scratch_space,  &
                                         LOG_LEVEL_INFO,     &
                                         LOG_LEVEL_ERROR,    &
                                         LOG_LEVEL_TRACE,    &
                                         LOG_LEVEL_WARNING,  &
                                         LOG_LEVEL_ALWAYS
  use minmax_tseries_mod,         only : minmax_tseries,      &
                                         minmax_tseries_init, &
                                         minmax_tseries_final
  use mesh_mod,                   only : mesh_type
  use mesh_collection_mod,        only : mesh_collection
  use clock_mod,                  only : clock_type
  use model_clock_mod,            only : model_clock_type
  use moisture_conservation_alg_mod, &
                                  only : moisture_conservation_alg
  use namelist_collection_mod,    only : namelist_collection_type
  use namelist_mod,               only : namelist_type
  use mr_indices_mod,             only : nummr
  use no_timestep_alg_mod,        only : no_timestep_type
  use remove_duplicates_mod,      only : remove_duplicates
  use runtime_constants_mod,      only : create_runtime_constants, &
                                         final_runtime_constants
  use timestep_method_mod,        only : timestep_method_type, &
                                         get_timestep_method_from_collection
  use rk_alg_timestep_mod,        only : rk_timestep_type
  use semi_implicit_timestep_alg_mod, &
                                  only : semi_implicit_timestep_type
  use setup_orography_alg_mod,    only : setup_orography_alg
  use idealised_config_mod,       only : perturb_init, perturb_seed
  use initial_temperature_config_mod, only : perturb, perturb_random
#ifdef COUPLED
  use coupler_mod,                only : create_coupling_fields, &
                                         generate_coupling_field_collections
  use coupling_mod,               only : coupling_type,       &
                                         get_coupling_from_collection
#endif

#ifdef UM_PHYSICS
  use gas_calc_all_mod,            only : gas_calc_all
  use jules_control_init_mod,      only : jules_control_init
  use jules_physics_init_mod,      only : jules_physics_init
  use planet_constants_mod,        only : set_planet_constants
  use socrates_init_mod,           only : socrates_init
  use illuminate_alg_mod,          only : illuminate_alg
  use um_clock_init_mod,           only : um_clock_init
  use um_control_init_mod,         only : um_control_init
  use um_domain_init_mod,          only : um_domain_init
  use um_domain_init_mod,          only : um_domain_init
  use um_sizes_init_mod,           only : um_sizes_init
  use um_physics_init_mod,         only : um_physics_init
  use um_radaer_lut_init_mod,      only : um_radaer_lut_init
  use um_ukca_init_mod,            only : um_ukca_init
  use jules_timestep_alg_mod,      only : jules_timestep_type
#endif

  implicit none

  private

  logical(l_def) :: use_moisture

  !> @brief Processor class for persisting fields
  type, extends(processor_type) :: persistor_type
    type(metafile_type), allocatable :: ckp_out(:)
    type(metafile_type), allocatable :: ckp_inp(:)
  contains
    private

    ! main interface
    procedure, public :: init => persistor_init
    procedure, public :: apply => persistor_apply

    ! destructor - here to avoid gnu compiler bug
    final :: persistor_destructor
  end type persistor_type

  public initialise_infrastructure, &
         initialise_model,          &
         finalise_infrastructure,   &
         finalise_model
contains

  !> @brief  Initialise processor object for persisting LFRic fields
  !> @param[in]   self      Persistor object
  !> @param[in]   clock     Model clock
  subroutine persistor_init(self, clock)

    use initialization_config_mod, only: init_option, &
                                         init_option_checkpoint_dump
    use io_config_mod,             only: checkpoint_read, &
                                         checkpoint_write, &
                                         checkpoint_times, &
                                         end_of_run_checkpoint
    use files_config_mod,          only: checkpoint_stem_name

    implicit none

    class(persistor_type),       intent(inout) :: self
    class(clock_type), target,   intent(in)    :: clock

    character(str_def) :: ckp_out_id
    character(str_def), parameter :: ckp_inp_id = 'lfric_checkpoint_read'
    character(:), allocatable :: split_stem_name(:)
    integer(i_def) :: i
    real(r_second), allocatable :: all_checkpoint_times(:)

    call self%set_clock(clock)

    ! Add all the checkpoint write files to XIOS named for the checkpoint times.
    if(checkpoint_write)then
      select type(clock)
      type is (model_clock_type)
        call create_checkpoint_list( clock, checkpoint_times, &
                                     end_of_run_checkpoint, all_checkpoint_times )
        if(size(all_checkpoint_times) == 0) then
          call log_event("Checkpointing has been turned on but no " // &
                         "checkpoint times have been specified.", &
                         log_level_error)
        end if
        allocate(self%ckp_out(size(all_checkpoint_times)))
        split_stem_name = split_string( trim(checkpoint_stem_name), '/' )
        do i = 1, size(all_checkpoint_times)
          write(ckp_out_id,'(A,A,I10.10)') &
                trim(split_stem_name(size(split_stem_name))),"_", &
                clock%steps_from_seconds(all_checkpoint_times(i))
          call log_event("Persistor initialisation: Creating checkpoint file "// &
            trim(ckp_out_id), LOG_LEVEL_INFO)
          call self%ckp_out(i)%init(ckp_out_id)
        end do
      class default
        call log_event("Persistor initialisation: Unsupported clock type", &
          LOG_LEVEL_ERROR)
      end select
    end if

    if (checkpoint_read .or. init_option == init_option_checkpoint_dump) then
      allocate(self%ckp_inp(1))
      call self%ckp_inp(1)%init(ckp_inp_id)
    end if

  end subroutine persistor_init

  !> @brief     Persistor's apply method
  !> @details   Persist a field given by field specifier by adding it to the checkpoint files
  !> @param[in]   self      Persistor object
  !> @param[in]   spec      Field specifier
  subroutine persistor_apply(self, spec)

    use initialization_config_mod, only: init_option, &
                                         init_option_checkpoint_dump
    use io_config_mod,             only: checkpoint_read, &
                                         checkpoint_write, checkpoint_times
    use files_config_mod,          only: checkpoint_stem_name

    implicit none

    class(persistor_type), intent(in) :: self
    type(field_spec_type), intent(in) :: spec

    character(20), parameter :: operation = 'once'

    if (spec%ckp .and. space_has_xios_io(spec%space, spec%legacy)) then
      if (checkpoint_write) then
        call add_field(self%ckp_out, spec%name, mode=CHECKPOINTING, &
          operation=operation, id_as_name=.true., legacy=spec%legacy)
      end if
      if (checkpoint_read .or. init_option == init_option_checkpoint_dump) &
        call add_field(self%ckp_inp, spec%name, mode=RESTARTING, &
          operation=operation, id_as_name=.true., legacy=spec%legacy)
    end if

  end subroutine persistor_apply

  !> @brief Destructor of persistor object
  !> @param[inout] self  Persistor object
  subroutine persistor_destructor(self)
    type(persistor_type), intent(inout) :: self
    ! empty
  end subroutine persistor_destructor

  !> @brief Enable active state fields for checkpointing; sync xios axis dimensions;
  !>        set up scaled diagnostics fields
  !> @param[in] clock        The clock providing access to time information
  subroutine before_context_close(clock)

    use multidata_field_dimensions_mod, only: sync_multidata_field_dimensions
    use time_dimensions_mod,            only: sync_time_dimensions
    use boundaries_config_mod,          only: limited_area
    use formulation_config_mod,         only: use_physics
    use section_choice_config_mod,      only: stochastic_physics, &
                                              stochastic_physics_um
    use initialization_config_mod, only: init_option, &
                                         init_option_checkpoint_dump
    use io_config_mod,                  only: checkpoint_read, checkpoint_write

    implicit none
    class(clock_type), intent(in) :: clock

    type(persistor_type) :: persistor

    real(r_second) :: DT

    DT = clock%get_seconds_per_step()
    call set_variable("DT", DT, tolerant=.true.)

    call persistor%init(clock)
    call process_gungho_prognostics(persistor)
    ! Add the temperature_correction_rate to the appropriate files
    if(checkpoint_write) then
        if ( encorr_usage /= encorr_usage_none ) then
          call add_field( persistor%ckp_out, "temperature_correction_rate", mode=CHECKPOINTING, operation="once", &
                          id_as_name=.true.)
        end if
        if (stochastic_physics == stochastic_physics_um) then
          call add_field( persistor%ckp_out, "random_seed", mode=CHECKPOINTING, operation="once", &
                          id_as_name=.true.)
        end if
    end if
    if (checkpoint_read .or. init_option == init_option_checkpoint_dump) then
      if ( encorr_usage /= encorr_usage_none ) then
        call add_field( persistor%ckp_inp, "temperature_correction_rate", mode=RESTARTING, operation="once", &
                        id_as_name=.true.)
      end if
      if (stochastic_physics == stochastic_physics_um) then
        call add_field( persistor%ckp_inp, "random_seed", mode=RESTARTING, operation="once", &
                        id_as_name=.true.)
      end if
    end if


    if (limited_area) call process_lbc_fields(persistor)
    if (use_physics) then
      call process_physics_prognostics(persistor)
      call sync_multidata_field_dimensions()
      call sync_time_dimensions()
    end if
  end subroutine before_context_close

  !> @brief Initialisations that depend only on loaded namelists.
  !> @details To be called before the IO context opens so that the
  !> multidata dimensions are completely defined when the XIOS
  !> axis dimensions are being synched.
  !> @param[in] model_clock     The clock providing access to time information
  subroutine basic_initialisations(mesh,model_clock)

#ifdef UM_PHYSICS
    use formulation_config_mod,     only: use_physics
    use section_choice_config_mod,  only: radiation,          &
                                          radiation_socrates, &
                                          surface,            &
                                          surface_jules
#endif

    implicit none

    type( mesh_type ), intent(in), pointer :: mesh
    class(model_clock_type), intent(inout) :: model_clock

#ifdef UM_PHYSICS
    integer(i_def) :: ncells
    integer(i_def) :: ncells_ukca

    ! Set derived planet constants and presets
    call set_planet_constants()

    ! Initialise UM to run full field
    ncells_ukca = mesh%get_last_edge_cell()

    ! Initialise UM to run in columns
    ncells = 1_i_def

    if ( use_physics ) then

      if (radiation == radiation_socrates) then
        ! Initialisation for the Socrates radiation scheme
        call socrates_init()
      end if
      ! Initialisation of UM high-level variables
      call um_control_init()

      ! Initialisation of UM physics variables
      call um_sizes_init(ncells)

      ! Initialisation of UM clock
      call um_clock_init(model_clock)

      ! Initialisation of UM physics variables
      call um_physics_init()

      ! Read all the radaer lut namelist files
      call um_radaer_lut_init()

      ! Initialisation of Jules high-level variables
      call jules_control_init()

      if (surface == surface_jules) then
        ! Initialisation of Jules physics variables
        call jules_physics_init()
      end if

      ! Initialisation of UKCA physics variables
      call um_ukca_init(ncells_ukca, model_clock)

    end if
#endif
  end subroutine basic_initialisations

  !> @brief Initialises the infrastructure and sets up constants used by the
  !>        model.
  !> @param [in]               The name of the IO context of the model
  !> @param [in,out] modeldb   The full model database for the model run
  !>
  subroutine initialise_infrastructure( io_context_name, modeldb )

    use base_mesh_config_mod, only: GEOMETRY_PLANAR, &
                                    GEOMETRY_SPHERICAL

#ifdef UM_PHYSICS
    use formulation_config_mod,    only: use_physics
    use section_choice_config_mod, only: radiation, &
                                         radiation_socrates
#endif

    implicit none

    character(*),         intent(in)    :: io_context_name
    class(modeldb_type),  intent(inout) :: modeldb

    type(mesh_type), pointer :: mesh                => null()
    type(mesh_type), pointer :: shifted_mesh        => null()
    type(mesh_type), pointer :: double_level_mesh   => null()
    type(mesh_type), pointer :: twod_mesh           => null()
    type(mesh_type), pointer :: orography_twod_mesh => null()
    type(mesh_type), pointer :: orography_mesh      => null()

    procedure(filelist_populator), pointer :: files_init_ptr => null()

    type(field_type) :: surface_altitude

    class(extrusion_type),        allocatable :: extrusion
    type(uniform_extrusion_type), allocatable :: extrusion_2d

    type(shifted_extrusion_type),      allocatable :: extrusion_shifted
    type(double_level_extrusion_type), allocatable :: extrusion_double

    type(field_type), pointer :: chi(:) => null()

#ifdef COUPLED
    type(field_collection_type), pointer :: depository
    type(field_collection_type), pointer :: cpl_snd_2d
    type(field_collection_type), pointer :: cpl_rcv_2d
    type(field_collection_type), pointer :: cpl_snd_0d
    type(coupling_type), pointer         :: coupling_ptr
#endif
    type(inventory_by_mesh_type), pointer :: chi_inventory      => null()
    type(inventory_by_mesh_type), pointer :: panel_id_inventory => null()

    logical(l_def)                  :: mesh_already_exists
    integer(i_def)                  :: i, j, mesh_ctr
    character(str_def), allocatable :: base_mesh_names(:)
    character(str_def), allocatable :: meshes_to_shift(:)
    character(str_def), allocatable :: meshes_to_double(:)
    character(str_def), allocatable :: tmp_mesh_names(:)
    character(str_def), allocatable :: extra_io_mesh_names(:)
    character(str_def), allocatable :: chain_mesh_tags(:)
    character(str_def), allocatable :: multires_coupling_mesh_tags(:)
    character(str_def), allocatable :: twod_names(:)
    character(str_def), allocatable :: shifted_names(:)
    character(str_def), allocatable :: double_names(:)

    character(str_def), allocatable :: meshes_to_check(:)

#ifdef UM_PHYSICS
    type(field_collection_type),  pointer :: radiation_fields
#endif
    integer :: start_index, end_index

    character(str_def) :: prime_mesh_name
    character(str_def) :: orography_mesh_name

    logical(l_def) :: use_multires_coupling
    logical(l_def) :: l_multigrid
    logical(l_def) :: prepartitioned
    logical(l_def) :: apply_partition_check

    integer(i_def) :: geometry
    integer(i_def) :: extrusion_method
    integer(i_def) :: stencil_depth
    real(r_def)    :: domain_bottom
    real(r_def)    :: domain_height
    real(r_def)    :: scaled_radius
    integer(i_def) :: number_of_layers

    type(namelist_type), pointer :: base_mesh_nml
    type(namelist_type), pointer :: formulation_nml
    type(namelist_type), pointer :: extrusion_nml
    type(namelist_type), pointer :: planet_nml
    type(namelist_type), pointer :: multigrid_nml
    type(namelist_type), pointer :: multires_coupling_nml
#ifdef UM_PHYSICS
    real(r_def) :: dt
#endif


    integer(i_def), parameter :: one_layer = 1_i_def

    integer(i_def) :: random_seed_size, proc_rank, big_int, perturb_seed_in
    integer(i_def), allocatable :: ranseed(:)

    !=======================================================================
    ! 0.0 Extract configuration variables
    !=======================================================================

    call check_configuration(modeldb)

    base_mesh_nml   => modeldb%configuration%get_namelist('base_mesh')
    formulation_nml => modeldb%configuration%get_namelist('formulation')
    extrusion_nml   => modeldb%configuration%get_namelist('extrusion')
    planet_nml      => modeldb%configuration%get_namelist('planet')
    multires_coupling_nml => null()
    multigrid_nml         => null()

    call formulation_nml%get_value( 'l_multigrid', l_multigrid )
    call formulation_nml%get_value( 'use_multires_coupling', &
                                    use_multires_coupling )

    if ( use_multires_coupling ) then
      multires_coupling_nml => modeldb%configuration%get_namelist('multires_coupling')
      call multires_coupling_nml%get_value( 'multires_coupling_mesh_tags', &
                                            multires_coupling_mesh_tags )
      call multires_coupling_nml%get_value( 'orography_mesh_name', &
                                            orography_mesh_name )
    end if

    if ( l_multigrid ) then
      multigrid_nml => modeldb%configuration%get_namelist('multigrid')
      call multigrid_nml%get_value( 'chain_mesh_tags', chain_mesh_tags )
    end if

    call base_mesh_nml%get_value( 'prime_mesh_name', prime_mesh_name )
    call base_mesh_nml%get_value( 'geometry', geometry )
    call base_mesh_nml%get_value( 'prepartitioned', prepartitioned )
    call extrusion_nml%get_value( 'domain_height', domain_height )
    call extrusion_nml%get_value( 'method', extrusion_method )
    call extrusion_nml%get_value( 'number_of_layers', number_of_layers )
    call planet_nml%get_value( 'scaled_radius', scaled_radius )

    !-------------------------------------------------------------------------
    ! Initialise infrastructure
    !-------------------------------------------------------------------------
    write(log_scratch_space,'(A)')                        &
        'Application built with '//trim(PRECISION_REAL)// &
        '-bit real numbers'
    call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

    call set_derived_config( .true. )

    !=======================================================================
    ! 1.1 Determine the required meshes
    !=======================================================================

    ! 1.1a Meshes that require a prime/2d extrusion
    ! ---------------------------------------------------------
    if ( allocated(tmp_mesh_names) ) deallocate(tmp_mesh_names)

    allocate( tmp_mesh_names(100) )
    tmp_mesh_names(:) = ''
    tmp_mesh_names(1) = prime_mesh_name
    start_index = 2

    if ( l_multigrid ) then
      end_index = start_index + size(chain_mesh_tags) - 1
      tmp_mesh_names(start_index:end_index) = chain_mesh_tags
      start_index = end_index + 1
    end if

    if ( use_multires_coupling ) then
      end_index = start_index + size(multires_coupling_mesh_tags) - 1
      tmp_mesh_names(start_index:end_index) = multires_coupling_mesh_tags
      start_index = end_index + 1
    end if

    base_mesh_names = remove_duplicates( tmp_mesh_names(:)  )

    if ( l_multigrid .or. use_multires_coupling ) then
      meshes_to_check = remove_duplicates( tmp_mesh_names(2:) )
    end if
    deallocate(tmp_mesh_names)


    ! 1.1b Meshes the require a shifted extrusion
    ! ---------------------------------------------------------
    if ( check_any_shifted() ) then
      allocate(tmp_mesh_names(size(base_mesh_names)))

      mesh_ctr = 1
      tmp_mesh_names(1) = prime_mesh_name

      if ( use_multires_coupling ) then
        do i=1, size(multires_coupling_mesh_tags)
          ! Only add mesh if it has not already been added
          mesh_already_exists = .false.
          do j = 1, mesh_ctr
            if ( tmp_mesh_names(j) == multires_coupling_mesh_tags(i) ) then
              mesh_already_exists = .true.
              exit
            end if
          end do
          if ( .not. mesh_already_exists ) then
            mesh_ctr = mesh_ctr + 1
            tmp_mesh_names(mesh_ctr) = multires_coupling_mesh_tags(i)
          end if
        end do
      end if

      ! Transfer mesh names from temporary array to an array of appropriate size
      allocate(meshes_to_shift(mesh_ctr))
      do i=1, mesh_ctr
        meshes_to_shift(i) = tmp_mesh_names(i)
      end do
      deallocate(tmp_mesh_names)
    end if


    ! 1.1c Meshes that require a double-level extrusion
    ! ---------------------------------------------------------
    if ( check_any_shifted() ) then
      allocate(tmp_mesh_names(size(base_mesh_names)))

      mesh_ctr = 1
      tmp_mesh_names(1) = prime_mesh_name

      if ( use_multires_coupling ) then
        do i=1, size(multires_coupling_mesh_tags)
          ! Only add mesh if it has not already been added
          mesh_already_exists = .false.
          do j=1, mesh_ctr
            if ( tmp_mesh_names(j) == multires_coupling_mesh_tags(i) ) then
              mesh_already_exists = .true.
              exit
            end if
          end do
          if ( .not. mesh_already_exists ) then
            mesh_ctr = mesh_ctr + 1
            tmp_mesh_names(mesh_ctr) = multires_coupling_mesh_tags(i)
          end if
        end do
      end if

      ! Transfer mesh names from temporary array to an array of appropriate size
      if ( check_any_wt_eqn_conservative() ) then
        allocate(meshes_to_double(mesh_ctr))
        do i=1, mesh_ctr
          meshes_to_double(i) = tmp_mesh_names(i)
        end do
        deallocate(tmp_mesh_names)
      end if
    end if


    !=======================================================================
    ! 1.2 Generate required extrusions
    !=======================================================================

    ! 1.2a Extrusions for prime/2d meshes
    ! ---------------------------------------------------------
    select case (geometry)
    case (geometry_planar)
      domain_bottom = 0.0_r_def
    case (geometry_spherical)
      domain_bottom = scaled_radius
    case default
      call log_event("Invalid geometry for mesh initialisation", LOG_LEVEL_ERROR)
    end select

    allocate( extrusion, source=create_extrusion( extrusion_method, &
                                                  geometry,         &
                                                  number_of_layers, &
                                                  domain_height,    &
                                                  scaled_radius ) )

    extrusion_2d = uniform_extrusion_type( domain_bottom, &
                                           domain_bottom, &
                                           one_layer, TWOD )

    ! 1.2b Extrusions for shifted meshes
    ! ---------------------------------------------------------
    if ( allocated(meshes_to_shift) ) then
      if ( size(meshes_to_shift) > 0 ) then
        extrusion_shifted = shifted_extrusion_type(extrusion)
      end if
    end if

    ! 1.2c Extrusions for double-level meshes
    ! ---------------------------------------------------------
    if ( allocated(meshes_to_double) ) then
      if ( size(meshes_to_double) > 0 ) then
        extrusion_double = double_level_extrusion_type(extrusion)
      end if
    end if


    !=======================================================================
    ! 1.3 Initialise mesh objects and assign InterGrid maps
    !=======================================================================

    ! 1.3a Initialise prime/2d meshes
    ! ---------------------------------------------------------
    apply_partition_check = .false.
    if ( .not. prepartitioned .and. &
         ( l_multigrid .or. use_multires_coupling ) ) then
      apply_partition_check = .true.
    end if

    stencil_depth = get_required_stencil_depth()
    call init_mesh( modeldb%configuration,        &
                    modeldb%mpi%get_comm_rank(),  &
                    modeldb%mpi%get_comm_size(),  &
                    base_mesh_names,              &
                    extrusion,                    &
                    get_required_stencil_depth(), &
                    apply_partition_check )


    allocate( twod_names, source=base_mesh_names )
    do i=1, size(twod_names)
      twod_names(i) = trim(twod_names(i))//'_2d'
    end do
    call create_mesh( base_mesh_names, extrusion_2d, &
                      alt_name=twod_names )
    call assign_mesh_maps(twod_names)


    ! 1.3b Initialise shifted meshes
    ! ---------------------------------------------------------
    if (allocated(meshes_to_shift)) then
      if (size(meshes_to_shift) > 0) then

        allocate( shifted_names, source=meshes_to_shift )
        do i=1, size(shifted_names)
          shifted_names(i) = trim(shifted_names(i))//'_shifted'
        end do
        call create_mesh( meshes_to_shift,   &
                          extrusion_shifted, &
                          alt_name=shifted_names )
        call assign_mesh_maps(shifted_names)

      end if
    end if

    ! 1.3c Initialise double-level meshes
    ! ---------------------------------------------------------
    if (allocated(meshes_to_double)) then
      if (size(meshes_to_double) > 0) then

        allocate( double_names, source=meshes_to_double )
        do i=1, size(double_names)
          double_names(i) = trim(double_names(i))//'_double'
        end do
        call create_mesh( meshes_to_double, &
                          extrusion_double, &
                          alt_name=double_names )
        call assign_mesh_maps(double_names)

      end if
    end if


    !=======================================================================
    ! 2.0 Initialise FEM / Coordinates
    !=======================================================================
    chi_inventory => get_chi_inventory()
    panel_id_inventory => get_panel_id_inventory()

    call init_fem( mesh_collection, chi_inventory, panel_id_inventory )
    if ( l_multigrid ) then
      call init_function_space_chains( mesh_collection, chain_mesh_tags )
    end if


    mesh      => mesh_collection%get_mesh(prime_mesh_name)
    twod_mesh => mesh_collection%get_mesh(mesh, TWOD)
    call chi_inventory%get_field_array(mesh, chi)


    !=======================================================================
    ! 3.0 Initialise coupling
    !=======================================================================
#ifdef COUPLED
    if( l_esm_couple ) then
       call log_event("Initialising coupler", LOG_LEVEL_INFO)
       ! Add fields used in coupling
       depository => modeldb%fields%get_field_collection("depository")

       call create_coupling_fields( mesh,       &
                                    twod_mesh,  &
                                    modeldb )
       ! Define coupling interface

       ! Set up collections to hold 2d coupling fields
       call modeldb%fields%add_empty_field_collection("cpl_snd_2d" , &
                                                 table_len = 30)
       call modeldb%fields%add_empty_field_collection("cpl_rcv_2d" , &
                                                 table_len = 30)

       cpl_snd_2d => modeldb%fields%get_field_collection("cpl_snd_2d")
       cpl_rcv_2d => modeldb%fields%get_field_collection("cpl_rcv_2d")

       ! Set up collection to hold 0d (scalar) coupling fields
       call modeldb%fields%add_empty_field_collection("cpl_snd_0d" , &
                                                 table_len = 30)

       cpl_snd_0d => modeldb%fields%get_field_collection("cpl_snd_0d")

       call generate_coupling_field_collections(cpl_snd_2d, &
                                                cpl_rcv_2d, &
                                                cpl_snd_0d, &
                                                depository)

       ! Extract the coupling object from the modeldb key-value pair collection
       coupling_ptr => get_coupling_from_collection(modeldb%values, "coupling" )
       call coupling_ptr%define_partitions(modeldb%mpi,twod_mesh)
       call coupling_ptr%define_variables( cpl_snd_2d, &
                                           cpl_rcv_2d, &
                                           cpl_snd_0d )
       call coupling_ptr%end_definition( cpl_snd_2d, &
                                         cpl_rcv_2d, &
                                         cpl_snd_0d )

       ! Set flags for ocean and sea-ice coupling based on whether
       ! ice fraction and SST are coupled
       l_couple_sea_ice = cpl_rcv_2d%field_exists("lf_icefrc")
       l_couple_ocean = cpl_rcv_2d%field_exists("lf_ocn_sst")

       if (l_couple_sea_ice .and. .not. l_couple_ocean) then
          write(log_scratch_space, '(A)' ) "You are attempting to couple to seaice, without an ocean, this is not possible"
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
       endif
    endif
#else
  l_couple_sea_ice = .false.
  l_couple_ocean = .false.
#endif


    !=======================================================================
    ! 4.0 Initialise output
    !=======================================================================
    call basic_initialisations( mesh, modeldb%clock )

    call log_event("Initialising I/O context", LOG_LEVEL_INFO)

    files_init_ptr => init_gungho_files

    if ( use_multires_coupling ) then
      ! Compose list of meshes for I/O
      ! This is list of namelist mesh names minus the prime mesh name
      allocate(extra_io_mesh_names(size(multires_coupling_mesh_tags)-1))
      mesh_ctr = 1
      do i=1, size(multires_coupling_mesh_tags)
        if (multires_coupling_mesh_tags(i) /= prime_mesh_name) then
          extra_io_mesh_names(mesh_ctr) = multires_coupling_mesh_tags(i)
          mesh_ctr = mesh_ctr + 1
        end if
      end do

      call init_io( io_context_name, prime_mesh_name, modeldb, &
                    chi_inventory, panel_id_inventory,         &
                    populate_filelist=files_init_ptr,          &
                    alt_mesh_names=extra_io_mesh_names,        &
                    before_close=before_context_close )

    else
      call init_io( io_context_name, prime_mesh_name, modeldb, &
                    chi_inventory, panel_id_inventory,         &
                    populate_filelist=files_init_ptr,          &
                    before_close=before_context_close )
    end if


    !=======================================================================
    ! 5.0 Initialise orography
    !=======================================================================
    if ( use_multires_coupling ) then
      orography_mesh => mesh_collection%get_mesh(trim(orography_mesh_name))
      orography_twod_mesh => mesh_collection%get_mesh(orography_mesh, TWOD)
    else
      orography_mesh => mesh_collection%get_mesh(prime_mesh_name)
      orography_twod_mesh => mesh_collection%get_mesh(orography_mesh, TWOD)
    end if

    ! Set up surface altitude field - this will be used to generate orography
    ! for models with global land mass included (i.e GAL)
    call init_altitude( orography_twod_mesh, surface_altitude )

    call setup_orography_alg( base_mesh_names,                &
                              orography_mesh%get_mesh_name(), &
                              chi_inventory,                  &
                              panel_id_inventory,             &
                              surface_altitude )


    !=======================================================================
    ! 6.0 Initialise runtime constants
    !=======================================================================
    call create_runtime_constants()

#ifdef UM_PHYSICS
    if ( use_physics ) then
      ! Initialise time-varying trace gases
      call gas_calc_all()
      if (radiation == radiation_socrates) then
        ! Initialisation for the Socrates radiation scheme
        radiation_fields => modeldb%fields%get_field_collection("radiation_fields")
        dt = real(modeldb%clock%get_seconds_per_step(), r_def)
        call illuminate_alg( radiation_fields,                    &
                             modeldb%clock%get_step(),            &
                             dt)
      end if
      ! Initialisation of UM variables related to the mesh
      call um_domain_init(mesh)
    end if
#endif

    if ( perturb_init ) then
      if ( perturb_seed == 0 ) then
        call log_event("Gungho: Perturbation seed is zero", &
                       LOG_LEVEL_ERROR)
      end if
      perturb_seed_in = perturb_seed
    else
      perturb_seed_in = 0
    end if

    ! Set random seed for initial perturbation
    ! Note that the seed will be set again with another way
    ! for stochastic physics afterwards.
    if ( perturb_init .or. perturb == perturb_random ) then

      ! Get length of random seed array required by the compiler
      call random_seed(size = random_seed_size)
      allocate(ranseed(random_seed_size))

      ! Get processor rank for unique seed on each processor
      proc_rank = modeldb%mpi%get_comm_rank()

      ! Pick large number to multiply by proc_rank
      big_int = 1234567

      ! Allocate random seed across array with some added entropy from
      ! big_int and array index
      ranseed(:) = perturb_seed_in + proc_rank*big_int

      ! Set seed
      call random_seed(put = ranseed(1:random_seed_size))
    end if

    !=======================================================================
    ! Housekeeping
    !=======================================================================
    nullify(mesh, twod_mesh, shifted_mesh, double_level_mesh, chi, &
            chi_inventory, panel_id_inventory, files_init_ptr,     &
            orography_mesh, orography_twod_mesh)
    deallocate(base_mesh_names)
    if (allocated(meshes_to_shift))  deallocate(meshes_to_shift)
    if (allocated(meshes_to_double)) deallocate(meshes_to_double)

  end subroutine initialise_infrastructure

  !---------------------------------------------------------------------------
  !> @brief Initialises the gungho application
  !>
  !> @param[in] mesh  The primary mesh
  !> @param[in,out] modeldb The working data set for the model run
  !>
  subroutine initialise_model( mesh, modeldb )

    use timestepping_config_mod, only: method,                 &
                                       method_semi_implicit,   &
                                       method_rk,              &
                                       method_no_timestepping, &
                                       method_jules

    use io_config_mod,          only: write_conservation_diag, &
                                      write_minmax_tseries
    use formulation_config_mod, only: moisture_formulation, &
                                      moisture_formulation_dry
    implicit none

    type( mesh_type ),    intent(in),    pointer :: mesh
    type( modeldb_type ), intent(inout), target  :: modeldb

    type( field_collection_type ), pointer :: prognostic_fields => null()
    type( field_collection_type ), pointer :: moisture_fields => null()
    type( field_array_type ), pointer      :: mr_array
    type( field_type ),            pointer :: mr(:) => null()

    type( field_type), pointer :: theta => null()
    type( field_type), pointer :: u => null()
    type( field_type), pointer :: rho => null()
    type( field_type), pointer :: exner => null()

    class(timestep_method_type), pointer  :: timestep_method => null()

    use_moisture = ( moisture_formulation /= moisture_formulation_dry )

    ! Get pointers to field collections for use downstream
    prognostic_fields => modeldb%fields%get_field_collection( &
                                                         "prognostic_fields")

    moisture_fields => modeldb%fields%get_field_collection("moisture_fields")
    call moisture_fields%get_field("mr", mr_array)
    mr => mr_array%bundle

    ! Get pointers to fields in the prognostic/diagnostic field collections
    ! for use downstream
    call prognostic_fields%get_field('theta', theta)
    call prognostic_fields%get_field('u', u)
    call prognostic_fields%get_field('rho', rho)
    call prognostic_fields%get_field('exner', exner)

    if (write_minmax_tseries) then
      call minmax_tseries_init('u')
      call minmax_tseries(u, 'u', mesh)
    end if

    select case( method )
      case( method_semi_implicit )  ! Semi-Implicit
        ! Initialise the semi-implicit timestep method
        allocate( timestep_method, source=semi_implicit_timestep_type(modeldb) )
        ! Add to the model database
        call modeldb%values%add_key_value('timestep_method', &
                        timestep_method)
        ! Output initial conditions
        if ( write_conservation_diag ) then
         call conservation_algorithm( rho,              &
                                      u,                &
                                      theta,            &
                                      mr,               &
                                      exner )
         if ( use_moisture ) &
           call moisture_conservation_alg( rho,              &
                                           mr,               &
                                           'Before timestep' )
        end if
      case( method_rk )             ! RK
        ! Initialise the Runge-Kutta timestep method
        allocate( timestep_method, source=rk_timestep_type(modeldb) )
        ! Add to the model database
        call modeldb%values%add_key_value('timestep_method', &
                        timestep_method)

        ! Output initial conditions
        if ( write_conservation_diag ) then
          call conservation_algorithm( rho,              &
                                       u,                &
                                       theta,            &
                                       mr,               &
                                       exner )
         if ( use_moisture ) &
           call moisture_conservation_alg( rho,              &
                                           mr,               &
                                           'Before timestep' )
        end if
      case( method_no_timestepping )
        ! Initialise a null-timestep method
        allocate( timestep_method, source=no_timestep_type() )
        ! Add to the model database
        call modeldb%values%add_key_value('timestep_method', &
                        timestep_method)
        write( log_scratch_space, &
                    '(A, A)' ) 'CAUTION: Running with no timestepping. ' // &
                    ' Prognostic fields not evolved'
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
#ifdef UM_PHYSICS
      case( method_jules )  ! jules
        ! Initialise the jules timestep method
        allocate( timestep_method, source=jules_timestep_type(modeldb) )
        ! Add to the model database
        call modeldb%values%add_key_value('timestep_method', &
                      timestep_method)
#endif
      case default
        call log_event("Gungho: Incorrect time stepping option chosen, "// &
                        "stopping program! ",LOG_LEVEL_ERROR)
    end select

  end subroutine initialise_model

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Finalises infrastructure and constants used by the model.
  !>
  subroutine finalise_infrastructure(modeldb)

    implicit none

    class(modeldb_type), intent(inout) :: modeldb

    !-------------------------------------------------------------------------
    ! Finalise IO
    !-------------------------------------------------------------------------
    call final_io( modeldb )

    !-------------------------------------------------------------------------
    ! Finalise constants
    !-------------------------------------------------------------------------
    call final_runtime_constants()

    !-------------------------------------------------------------------------
    ! Finalise aspects of the grid
    !-------------------------------------------------------------------------
    call final_fem()

  end subroutine finalise_infrastructure

  !---------------------------------------------------------------------------
  !> @brief Finalise the gungho application
  !>
  !> @param[in,out] modeldb       The working data set for the model run
  !> @param[in]     program_name  An identifier given to the model run
  !>
  subroutine finalise_model( modeldb,       &
                             program_name )

    use io_config_mod, only: write_minmax_tseries

    implicit none

    type( modeldb_type ), target,     intent(inout) :: modeldb
    character(*),         optional,   intent(in)    :: program_name

    type( field_collection_type ), pointer :: diagnostic_fields => null()
    type( field_collection_type ), pointer :: moisture_fields => null()
    type( field_array_type ),      pointer :: mr_array
    type( field_type ),            pointer :: mr(:) => null()
    type( field_collection_type ), pointer :: fd_fields
    type( field_collection_type ), pointer :: prognostic_fields => null()

    type( field_type), pointer :: theta => null()
    type( field_type), pointer :: u => null()
    type( field_type), pointer :: rho => null()
    type( field_type), pointer :: exner => null()

    class(timestep_method_type), pointer :: timestep_method

    ! Pointer for setting I/O handlers on fields
    procedure(write_interface), pointer :: tmp_write_ptr => null()

    if ( present(program_name) ) then
      ! Get pointers to field collections for use downstream
      prognostic_fields => modeldb%fields%get_field_collection( &
                                                           "prognostic_fields")
      diagnostic_fields => modeldb%fields%get_field_collection( &
                                                           "diagnostic_fields")
      moisture_fields => modeldb%fields%get_field_collection("moisture_fields")
      call moisture_fields%get_field("mr", mr_array)
      mr => mr_array%bundle
      fd_fields => modeldb%fields%get_field_collection("fd_fields")

      ! Get pointers to fields in the prognostic/diagnostic field collections
      ! for use downstream
      call prognostic_fields%get_field('theta', theta)
      call prognostic_fields%get_field('u', u)
      call prognostic_fields%get_field('rho', rho)
      call prognostic_fields%get_field('exner', exner)

      ! Write checksums to file
      if (use_moisture) then
        call checksum_alg(program_name, rho, 'rho', theta, 'theta', u, 'u', &
                        field_bundle=mr, bundle_name='mr')
      else
        call checksum_alg(program_name, rho, 'rho', theta, 'theta', u, 'u')
      end if

      if (write_minmax_tseries) call minmax_tseries_final()
    end if

    ! Finalise the timestep method
    if ( modeldb%values%key_value_exists('timestep_method') ) then
      timestep_method => get_timestep_method_from_collection(     &
                    modeldb%values, "timestep_method")
      call timestep_method%finalise()
      call modeldb%values%remove_key_value('timestep_method')
    end if

  end subroutine finalise_model

end module gungho_model_mod
