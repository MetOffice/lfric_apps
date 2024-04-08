! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_lfric_driver_mod

USE constants_mod,              ONLY: i_def, r_def, l_def, r_second, str_def
USE log_mod,                    ONLY: log_event, log_scratch_space,            &
                                      LOG_LEVEL_INFO, LOG_LEVEL_ERROR,         &
                                      LOG_LEVEL_ALWAYS

! LFRic Modules
USE add_mesh_map_mod,           ONLY: assign_mesh_maps
USE create_mesh_mod,            ONLY: create_mesh
USE driver_collections_mod,     ONLY: init_collections, final_collections
USE driver_mesh_mod,            ONLY: init_mesh
USE driver_fem_mod,             ONLY: init_fem
USE driver_log_mod,             ONLY: init_logger, final_logger
USE derived_config_mod,         ONLY: set_derived_config
USE extrusion_mod,              ONLY: extrusion_type,         &
                                      uniform_extrusion_type, &
                                      TWOD
USE field_collection_mod,       ONLY: field_collection_type
USE field_mod,                  ONLY: field_type
USE geometric_constants_mod,    ONLY: get_chi_inventory, get_panel_id_inventory
USE gungho_extrusion_mod,       ONLY: create_extrusion
USE halo_comms_mod,             ONLY: initialise_halo_comms
USE inventory_by_mesh_mod,      ONLY: inventory_by_mesh_type
USE model_clock_mod,            ONLY: model_clock_type
USE io_context_mod,             ONLY: callback_clock_arg
USE lfric_xios_context_mod,     ONLY: lfric_xios_context_type, advance
USE lfric_xios_driver_mod,      ONLY: lfric_xios_initialise, &
                                      lfric_xios_finalise
USE lfricinp_setup_io_mod,      ONLY: io_config
USE linked_list_mod,            ONLY: linked_list_type
USE mesh_mod,                   ONLY: mesh_type
USE mesh_collection_mod,        ONLY: mesh_collection
USE namelist_collection_mod,    ONLY: namelist_collection_type
USE namelist_mod,               ONLY: namelist_type
USE lfricinp_runtime_constants_mod, ONLY: lfricinp_create_runtime_constants
USE step_calendar_mod,          ONLY: step_calendar_type

! Interface to mpi
USE mpi_mod,                    ONLY: global_mpi, create_comm, destroy_comm

! Configuration modules
USE base_mesh_config_mod,       ONLY: geometry_spherical, &
                                      geometry_planar

! lfricinp modules
USE lfricinp_um_parameters_mod, ONLY: fnamelen

IMPLICIT NONE

PRIVATE
PUBLIC :: lfricinp_initialise_lfric, lfricinp_finalise_lfric, lfric_fields,    &
          io_context

! Input namelist configuration
CHARACTER(LEN=fnamelen), PUBLIC :: lfric_nl_fname

CHARACTER(len=fnamelen) :: xios_id
! xios_ctx names needs to match iodef.xml file
CHARACTER(len=*), PARAMETER :: xios_ctx  = "gungho_atm"
CHARACTER(len=fnamelen) :: program_name

! MPI ranks
INTEGER(KIND=i_def), PUBLIC :: total_ranks
INTEGER(KIND=i_def), PUBLIC :: local_rank

INTEGER(KIND=i_def), PUBLIC :: comm = -999

TYPE(mesh_type), PUBLIC, pointer :: mesh      => null()
TYPE(mesh_type), PUBLIC, pointer :: twod_mesh => null()

! Container for all input fields
TYPE(field_collection_type) :: lfric_fields

TYPE(model_clock_type), PUBLIC, ALLOCATABLE :: model_clock
type(lfric_xios_context_type),  ALLOCATABLE :: io_context

CONTAINS

SUBROUTINE lfricinp_initialise_lfric(program_name_arg,                         &
                                     required_lfric_namelists,                 &
                                     start_date, time_origin,                  &
                                     first_step, last_step,                    &
                                     spinup_period, seconds_per_step)

! Description:
!  Initialises LFRic infrastructure, MPI, XIOS and halos.

IMPLICIT NONE

CHARACTER(LEN=*),    INTENT(IN) :: program_name_arg
CHARACTER(LEN=*),    INTENT(IN) :: required_lfric_namelists(:)
CHARACTER(LEN=*),    INTENT(IN) :: start_date, time_origin
INTEGER(KIND=i_def), INTENT(IN) :: first_step, last_step
REAL(r_second),      INTENT(IN) :: spinup_period
REAL(r_second),      INTENT(IN) :: seconds_per_step

TYPE(step_calendar_type), ALLOCATABLE :: model_calendar
TYPE(linked_list_type),   POINTER     :: file_list => null()

TYPE(field_type), POINTER :: chi(:) => null()
TYPE(field_type), POINTER :: panel_id => null()
TYPE(inventory_by_mesh_type), POINTER :: chi_inventory => null()
TYPE(inventory_by_mesh_type), POINTER :: panel_id_inventory => null()
PROCEDURE(callback_clock_arg), POINTER :: before_close => null()


TYPE(namelist_collection_type), SAVE :: configuration

TYPE(namelist_type), POINTER :: base_mesh_nml => null()
TYPE(namelist_type), POINTER :: planet_nml    => null()

CLASS(extrusion_type),        ALLOCATABLE :: extrusion
TYPE(uniform_extrusion_type), ALLOCATABLE :: extrusion_2d
CHARACTER(str_def),           ALLOCATABLE :: base_mesh_names(:)
CHARACTER(str_def),           ALLOCATABLE :: twod_names(:)

INTEGER(i_def), PARAMETER :: one_layer = 1_i_def
INTEGER(i_def) :: i

CHARACTER(str_def) :: prime_mesh_name

INTEGER(i_def) :: stencil_depth
INTEGER(i_def) :: geometry
REAL(r_def)    :: domain_bottom
REAL(r_def)    :: scaled_radius
LOGICAL(l_def) :: check_partitions

!=====================================================================

! Set module variables
program_name = program_name_arg
xios_id = TRIM(program_name) // "_client"

! Initialise MPI and create the default communicator: mpi_comm_world
CALL create_comm(comm)

! Initialise xios
CALL lfric_xios_initialise( program_name, comm, .false. )

! Save LFRic's part of the split communicator for later use, and
! set the total number of ranks and the local rank of the split
! communicator
CALL global_mpi%initialise(comm)
total_ranks = global_mpi%get_comm_size()
local_rank = global_mpi%get_comm_rank()

!Initialise halo functionality
CALL initialise_halo_comms( comm )

CALL configuration%initialise( program_name_arg, table_len=10 )
CALL load_configuration( lfric_nl_fname, required_lfric_namelists, &
                         configuration )

! Initialise logging system
CALL init_logger( comm, program_name )

CALL init_collections()

WRITE(log_scratch_space, '(2(A,I0))') 'total ranks = ', total_ranks, &
                         ', local_rank = ', local_rank
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)

! Sets variables used interally by the LFRic infrastructure.
CALL set_derived_config( .TRUE. )

CALL log_event('Initialising mesh', LOG_LEVEL_INFO)

! -------------------------------
! 0.0 Extract namelist variables
! -------------------------------
base_mesh_nml => configuration%get_namelist('base_mesh')
planet_nml    => configuration%get_namelist('planet')
CALL base_mesh_nml%get_value( 'prime_mesh_name', prime_mesh_name )
CALL base_mesh_nml%get_value( 'geometry', geometry )
CALL planet_nml%get_value( 'scaled_radius', scaled_radius )
base_mesh_nml => null()
planet_nml    => null()

!-------------------------------------------------------------------------
! 1.0 Create the meshes
!-------------------------------------------------------------------------
ALLOCATE(base_mesh_names(1))
base_mesh_names(1) = prime_mesh_name

!-------------------------------------------------------------------------
! 1.1 Create the required extrusions
!-------------------------------------------------------------------------
SELECT CASE ( geometry )
CASE ( GEOMETRY_PLANAR )
  domain_bottom = 0.0_r_def
CASE ( GEOMETRY_SPHERICAL )
  domain_bottom = scaled_radius
CASE default
  CALL log_event( "Invalid geometry for mesh initialisation", &
                  LOG_LEVEL_ERROR )
END SELECT

ALLOCATE( extrusion, source=create_extrusion() )
extrusion_2d = uniform_extrusion_type( domain_bottom, &
                                       domain_bottom, &
                                       one_layer, TWOD )

ALLOCATE( twod_names, source=base_mesh_names )
do i=1, size(twod_names)
  twod_names(i) = TRIM(twod_names(i))//'_2d'
end do

!-------------------------------------------------------------------------
! 1.2 Create the required meshes
!-------------------------------------------------------------------------
stencil_depth = 2_i_def
check_partitions = .false.
CALL init_mesh( configuration,              &
                local_rank, total_ranks,    &
                base_mesh_names, extrusion, &
                stencil_depth, check_partitions )

CALL create_mesh( base_mesh_names, extrusion_2d, &
                  alt_name=twod_names )
CALL assign_mesh_maps( twod_names )

!-------------------------------------------------------------------------
! 2.0 Create FEM specifics (function spaces and chi field)
!-------------------------------------------------------------------------
CALL log_event('Creating function spaces and chi', LOG_LEVEL_INFO)
chi_inventory => get_chi_inventory()
panel_id_inventory => get_panel_id_inventory()
CALL init_fem(mesh_collection, chi_inventory, panel_id_inventory)

! XIOS domain initialisation
mesh => mesh_collection%get_mesh(prime_mesh_name)
twod_mesh => mesh_collection%get_mesh(mesh, TWOD)
CALL chi_inventory%get_field_array(mesh, chi)
CALL panel_id_inventory%get_field(mesh, panel_id)
model_calendar = step_calendar_type(time_origin, start_date)
model_clock = model_clock_type( first_step, last_step, seconds_per_step, &
                                spinup_period )

allocate( io_context )
file_list => io_context%get_filelist()
CALL io_config%init_lfricinp_files(file_list)
call io_context%initialise( xios_ctx )
CALL io_context%initialise_xios_context( comm, chi, panel_id, &
                                         model_clock, model_calendar, before_close )
CALL advance(io_context, model_clock)

! Initialise runtime constants
CALL log_event('Initialising runtime constants', LOG_LEVEL_INFO)
CALL lfricinp_create_runtime_constants(mesh_collection,    &
                                       chi_inventory,      &
                                       panel_id_inventory, &
                                       seconds_per_step)

nullify(chi, panel_id, chi_inventory, panel_id_inventory)

END SUBROUTINE lfricinp_initialise_lfric

!------------------------------------------------------------------

SUBROUTINE load_configuration( lfric_nl, required_lfric_namelists, &
                               configuration )

! Description:
!  Reads lfric namelists and checks that all required namelists are present

USE configuration_mod, ONLY: read_configuration, ensure_configuration

IMPLICIT NONE

CHARACTER(*), INTENT(IN) :: lfric_nl

CHARACTER(*), INTENT(IN)  :: required_lfric_namelists(:)

TYPE(namelist_collection_type), INTENT(INOUT) :: configuration

LOGICAL              :: okay
LOGICAL, ALLOCATABLE :: success_map(:)
INTEGER              :: i

ALLOCATE(success_map(SIZE(required_lfric_namelists)))

CALL log_event('Loading '//TRIM(program_name)//' configuration ...',           &
               LOG_LEVEL_ALWAYS)

CALL read_configuration( lfric_nl, configuration )

okay = ensure_configuration(required_lfric_namelists, success_map)
IF (.NOT. okay) THEN
  WRITE(log_scratch_space, '(A)')                                              &
                         'The following required namelists were not loaded:'
  DO i = 1, SIZE(required_lfric_namelists)
    IF (.NOT. success_map(i))                                                  &
      log_scratch_space = TRIM(log_scratch_space) // ' '                       &
                          // required_lfric_namelists(i)
  END DO
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

DEALLOCATE(success_map)

END SUBROUTINE load_configuration

!-------------------------------------------------------------------------------

SUBROUTINE lfricinp_finalise_lfric()

! Description:
!  Call finalise routines for associated APIs and logging system

USE halo_comms_mod,            ONLY: finalise_halo_comms
USE log_mod,                   ONLY: log_event, LOG_LEVEL_INFO

! External libraries
USE xios,                      ONLY: xios_finalize


IMPLICIT NONE

CALL log_event( 'Calling lfric finalise routines', LOG_LEVEL_INFO )

! Finalise halos, XIOS, etc.
CALL finalise_halo_comms()
deallocate( io_context )
CALL lfric_xios_finalise()

CALL final_collections()

! Finalise the logging system. This has to be done before finallising MPI
! as logging is an MPI process.
!
CALL final_logger(program_name)

CALL global_mpi%finalise()
CALL destroy_comm()

END SUBROUTINE lfricinp_finalise_lfric

END MODULE lfricinp_lfric_driver_mod
