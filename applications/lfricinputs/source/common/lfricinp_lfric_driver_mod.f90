! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfricinp_lfric_driver_mod

use constants_mod,              only: i_def, r_def, l_def, r_second, str_def
use log_mod,                    only: log_event, log_scratch_space,            &
                                      LOG_LEVEL_INFO, LOG_LEVEL_ERROR,         &
                                      LOG_LEVEL_ALWAYS

! LFRic Modules
use add_mesh_map_mod,           only: assign_mesh_maps
use create_mesh_mod,            only: create_mesh
use driver_collections_mod,     only: init_collections, final_collections
use driver_mesh_mod,            only: init_mesh
use driver_fem_mod,             only: init_fem
use driver_log_mod,             only: init_logger, final_logger
use derived_config_mod,         only: set_derived_config
use event_mod,                  only: event_action
use event_actor_mod,            only: event_actor_type
use extrusion_mod,              only: extrusion_type,         &
                                      uniform_extrusion_type, &
                                      TWOD
use field_collection_mod,       only: field_collection_type
use field_mod,                  only: field_type
use sci_geometric_constants_mod,      &
                                only: get_chi_inventory, get_panel_id_inventory
use gungho_extrusion_mod,       only: create_extrusion
use halo_comms_mod,             only: initialise_halo_comms
use inventory_by_mesh_mod,      only: inventory_by_mesh_type
use model_clock_mod,            only: model_clock_type
use io_context_mod,             only: callback_clock_arg
use lfric_xios_context_mod,     only: lfric_xios_context_type
use lfric_xios_action_mod,      only: advance
use lfric_xios_driver_mod,      only: lfric_xios_initialise, &
                                      lfric_xios_finalise
use lfricinp_setup_io_mod,      only: io_config
use linked_list_mod,            only: linked_list_type
use mesh_mod,                   only: mesh_type
use mesh_collection_mod,        only: mesh_collection
use namelist_collection_mod,    only: namelist_collection_type
use namelist_mod,               only: namelist_type
use step_calendar_mod,          only: step_calendar_type

! Interface to mpi
use lfric_mpi_mod,              only: global_mpi, create_comm, destroy_comm, &
                                      lfric_comm_type

! Configuration modules
use base_mesh_config_mod,       only: geometry_spherical, &
                                      geometry_planar

! lfricinp modules
use lfricinp_um_parameters_mod, only: fnamelen

implicit none

private
public :: lfricinp_initialise_lfric, lfricinp_finalise_lfric, lfric_fields,    &
          io_context

! Input namelist configuration
character(len=fnamelen), public :: lfric_nl_fname

character(len=fnamelen) :: xios_id
! xios_ctx names needs to match iodef.xml file
character(len=*), parameter :: xios_ctx  = "gungho_atm"
character(len=fnamelen) :: program_name

! MPI ranks
integer(kind=i_def), public :: total_ranks
integer(kind=i_def), public :: local_rank

type(lfric_comm_type), public :: comm

type(mesh_type), public, pointer :: mesh      => null()
type(mesh_type), public, pointer :: twod_mesh => null()

! Container for all input fields
type(field_collection_type) :: lfric_fields

type(model_clock_type), public, allocatable :: model_clock
type(lfric_xios_context_type), target :: io_context

contains

subroutine lfricinp_initialise_lfric(program_name_arg,                         &
                                     required_lfric_namelists,                 &
                                     start_date, time_origin,                  &
                                     first_step, last_step,                    &
                                     spinup_period, seconds_per_step)

! Description:
!  Initialises LFRic infrastructure, MPI, XIOS and halos.

implicit none

character(len=*),    intent(in) :: program_name_arg
character(len=*),    intent(in) :: required_lfric_namelists(:)
character(len=*),    intent(in) :: start_date, time_origin
integer(kind=i_def), intent(in) :: first_step, last_step
real(r_second),      intent(in) :: spinup_period
real(r_second),      intent(in) :: seconds_per_step

type(step_calendar_type), allocatable :: model_calendar
type(linked_list_type),   pointer     :: file_list => null()

type(field_type), pointer :: chi(:) => null()
type(field_type), pointer :: panel_id => null()
type(inventory_by_mesh_type), pointer :: chi_inventory => null()
type(inventory_by_mesh_type), pointer :: panel_id_inventory => null()
procedure(callback_clock_arg), pointer :: before_close => null()
class(event_actor_type), pointer :: event_actor_ptr
procedure(event_action), pointer :: context_advance


type(namelist_collection_type), save :: configuration

type(namelist_type), pointer :: base_mesh_nml
type(namelist_type), pointer :: planet_nml
type(namelist_type), pointer :: extrusion_nml

class(extrusion_type),        allocatable :: extrusion
type(uniform_extrusion_type), allocatable :: extrusion_2d
character(str_def),           allocatable :: base_mesh_names(:)
character(str_def),           allocatable :: twod_names(:)

integer(i_def), parameter :: one_layer = 1_i_def
integer(i_def) :: i

character(str_def) :: prime_mesh_name

integer(i_def) :: stencil_depth
integer(i_def) :: geometry
real(r_def)    :: domain_bottom
real(r_def)    :: scaled_radius
logical(l_def) :: check_partitions
integer        :: extrusion_method
integer(i_def) :: number_of_layers
real(r_def)    :: domain_height

!=====================================================================

! Set module variables
program_name = program_name_arg
xios_id = trim(program_name) // "_client"

! Initialise MPI and create the default communicator: mpi_comm_world
call create_comm(comm)

! Initialise xios
call lfric_xios_initialise( program_name, comm, .false. )

! Save LFRic's part of the split communicator for later use, and
! set the total number of ranks and the local rank of the split
! communicator
call global_mpi%initialise(comm)
total_ranks = global_mpi%get_comm_size()
local_rank = global_mpi%get_comm_rank()

!Initialise halo functionality
call initialise_halo_comms( comm )

call configuration%initialise( program_name_arg, table_len=10 )
call load_configuration( lfric_nl_fname, required_lfric_namelists, &
                         configuration )

! Initialise logging system
call init_logger( comm, program_name )

call init_collections()

write(log_scratch_space, '(2(A,I0))') 'total ranks = ', total_ranks, &
                         ', local_rank = ', local_rank
call log_event(log_scratch_space, LOG_LEVEL_INFO)

! Sets variables used interally by the LFRic infrastructure.
call set_derived_config( .true. )

call log_event('Initialising mesh', LOG_LEVEL_INFO)

! -------------------------------
! 0.0 Extract namelist variables
! -------------------------------
base_mesh_nml => configuration%get_namelist('base_mesh')
planet_nml    => configuration%get_namelist('planet')
extrusion_nml => configuration%get_namelist('extrusion')

call base_mesh_nml%get_value( 'prime_mesh_name', prime_mesh_name )
call base_mesh_nml%get_value( 'geometry', geometry )
call planet_nml%get_value( 'scaled_radius', scaled_radius )
call extrusion_nml%get_value( 'method', extrusion_method )
call extrusion_nml%get_value( 'number_of_layers', number_of_layers )
call extrusion_nml%get_value( 'domain_height', domain_height )

!-------------------------------------------------------------------------
! 1.0 Create the meshes
!-------------------------------------------------------------------------
allocate(base_mesh_names(1))
base_mesh_names(1) = prime_mesh_name

!-------------------------------------------------------------------------
! 1.1 Create the required extrusions
!-------------------------------------------------------------------------
select case ( geometry )
case ( GEOMETRY_PLANAR )
  domain_bottom = 0.0_r_def
case ( GEOMETRY_SPHERICAL )
  domain_bottom = scaled_radius
case default
  call log_event( "Invalid geometry for mesh initialisation", &
                  LOG_LEVEL_ERROR )
end select

allocate( extrusion, source=create_extrusion( extrusion_method, &
                                              geometry,         &
                                              number_of_layers, &
                                              domain_height,    &
                                              scaled_radius ) )

extrusion_2d = uniform_extrusion_type( domain_bottom, &
                                       domain_bottom, &
                                       one_layer, TWOD )

allocate( twod_names, source=base_mesh_names )
do i=1, size(twod_names)
  twod_names(i) = trim(twod_names(i))//'_2d'
end do

!-------------------------------------------------------------------------
! 1.2 Create the required meshes
!-------------------------------------------------------------------------
stencil_depth = 2_i_def
check_partitions = .false.
call init_mesh( configuration,              &
                local_rank, total_ranks,    &
                base_mesh_names, extrusion, &
                stencil_depth, check_partitions )

call create_mesh( base_mesh_names, extrusion_2d, &
                  alt_name=twod_names )
call assign_mesh_maps( twod_names )

!-------------------------------------------------------------------------
! 2.0 Create FEM specifics (function spaces and chi field)
!-------------------------------------------------------------------------
call log_event('Creating function spaces and chi', LOG_LEVEL_INFO)
chi_inventory => get_chi_inventory()
panel_id_inventory => get_panel_id_inventory()
call init_fem(mesh_collection, chi_inventory, panel_id_inventory)

! XIOS domain initialisation
mesh => mesh_collection%get_mesh(prime_mesh_name)
twod_mesh => mesh_collection%get_mesh(mesh, TWOD)
call chi_inventory%get_field_array(mesh, chi)
call panel_id_inventory%get_field(mesh, panel_id)
model_calendar = step_calendar_type(time_origin, start_date)
model_clock = model_clock_type( first_step, last_step, seconds_per_step, &
                                max(spinup_period, 0.0_r_second) )

file_list => io_context%get_filelist()
call io_config%init_lfricinp_files(file_list)
call io_context%initialise( xios_ctx )
call io_context%initialise_xios_context( comm, chi, panel_id, &
                                         model_clock, model_calendar, before_close )
! Attach context advancement to the model's clock
context_advance => advance
event_actor_ptr => io_context
call model_clock%add_event( context_advance, event_actor_ptr)

call advance(io_context, model_clock)

nullify(chi, panel_id, chi_inventory, panel_id_inventory)

end subroutine lfricinp_initialise_lfric

!------------------------------------------------------------------

subroutine load_configuration( lfric_nl, required_lfric_namelists, &
                               configuration )

! Description:
!  Reads lfric namelists and checks that all required namelists are present

use configuration_mod, only: read_configuration, ensure_configuration

implicit none

character(*), intent(in) :: lfric_nl

character(*), intent(in)  :: required_lfric_namelists(:)

type(namelist_collection_type), intent(INOUT) :: configuration

logical              :: okay
logical, allocatable :: success_map(:)
integer              :: i

allocate(success_map(size(required_lfric_namelists)))

call log_event('Loading '//trim(program_name)//' configuration ...',           &
               LOG_LEVEL_ALWAYS)

call read_configuration( lfric_nl, configuration )

okay = ensure_configuration(required_lfric_namelists, success_map)
if (.not. okay) then
  write(log_scratch_space, '(A)')                                              &
                         'The following required namelists were not loaded:'
  do i = 1, size(required_lfric_namelists)
    if (.not. success_map(i))                                                  &
      log_scratch_space = trim(log_scratch_space) // ' '                       &
                          // required_lfric_namelists(i)
  end do
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end if

deallocate(success_map)

end subroutine load_configuration

!-------------------------------------------------------------------------------

subroutine lfricinp_finalise_lfric()

! Description:
!  Call finalise routines for associated APIs and logging system

use halo_comms_mod,            only: finalise_halo_comms
use log_mod,                   only: log_event, LOG_LEVEL_INFO

implicit none

call log_event( 'Calling lfric finalise routines', LOG_LEVEL_INFO )


! Finalise ...
! (note that the order of the calls of the following finalisers matters)
call final_collections()

call finalise_halo_comms()

call io_context%finalise_xios_context()
call lfric_xios_finalise()

call final_logger(program_name)

call global_mpi%finalise()
call destroy_comm()


end subroutine lfricinp_finalise_lfric

end module lfricinp_lfric_driver_mod
