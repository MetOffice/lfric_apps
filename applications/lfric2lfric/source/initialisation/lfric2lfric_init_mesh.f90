!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used
!-----------------------------------------------------------------------------
!>
!> @brief   Set up specified mesh(es) from global/local mesh input file(s).
!> @details This routine will create a mesh_object_type(s) from a
!>          specified mesh input file and extrusion.
!>
!>          The algorithm differs depending on whether the input files(s)
!>          are prepartitioned (local mesh files) or not (global mesh files).
!>
!>          The result will be:
!>            * A set of local mesh objects stored in the application local
!>              mesh collection object.
!>            * A set of mesh objects stored in the application mesh collection
!>              object.
!>
!>          Local mesh object names will use the same mesh name as given in the
!>          input file. Extruded meshes are allowed to have an alternative mesh
!>          name as several meshes could exist in memory based on the same local
!>          mesh object.
!>
module lfric2lfric_init_mesh_mod

  use add_mesh_map_mod,            only: assign_mesh_maps
  use constants_mod,               only: i_def, l_def, str_max_filename
  use check_global_mesh_mod,       only: check_global_mesh
  use check_local_mesh_mod,        only: check_local_mesh
  use create_mesh_mod,             only: create_mesh
  use extrusion_mod,               only: extrusion_type
  use load_global_mesh_mod,        only: load_global_mesh
  use load_local_mesh_mod,         only: load_local_mesh
  use load_local_mesh_maps_mod,    only: load_local_mesh_maps
  use log_mod,                     only: log_event,         &
                                         log_scratch_space, &
                                         log_level_info,    &
                                         log_level_error,   &
                                         log_level_debug
  use namelist_collection_mod,     only: namelist_collection_type
  use namelist_mod,                only: namelist_type
  use panel_decomposition_mod,     only: panel_decomposition_type
  use partition_mod,               only: partitioner_interface
  use runtime_partition_mod,       only: mesh_cubedsphere,       &
                                         mesh_planar,            &
                                         create_local_mesh_maps, &
                                         create_local_mesh
  use runtime_partition_lfric_mod, only: get_partition_parameters
  use global_mesh_collection_mod,  only: global_mesh_collection

  ! Lfric2lfric modules
  use lfric2lfric_config_mod,      only: regrid_method_map,                   &
                                         source_geometry_spherical,           &
                                         destination_geometry_spherical,      &
                                         source_topology_fully_periodic,      &
                                         destination_topology_fully_periodic

  implicit none

  private
  public :: init_mesh

contains
!=======================================


!===============================================================================
!> @brief  Generates mesh(es) from mesh input file(s) on a given extrusion.
!>
!> @param[in] configuration          Application configuration object.
!>                                   This configuration object should contain the
!>                                   following defined namelist objects:
!>                                      * partititioning
!> @param[in] local_rank             The MPI rank of this process.
!> @param[in] total_ranks            Total number of MPI ranks in this job.
!> @param[in] mesh_names             Mesh names to load from the mesh input file(s).
!> @param[in] extrusion              Extrusion object to be applied to meshes.
!> @param[in] stencil_depth          Required stencil depth for the application.
!> @param[in] regrid_method          Apply check for even partitions with the
!>                                   configured partition strategy if the
!>                                   regridding method is 'map'.
!>                                   (unpartitioned mesh input only)
!===============================================================================
subroutine init_mesh( configuration,           &
                      local_rank, total_ranks, &
                      mesh_names,              &
                      extrusion,               &
                      stencil_depth,           &
                      regrid_method )

  implicit none

  ! Arguments
  type(namelist_collection_type) :: configuration

  integer(kind=i_def),   intent(in) :: local_rank
  integer(kind=i_def),   intent(in) :: total_ranks
  character(len=*),      intent(in) :: mesh_names(2)
  class(extrusion_type), intent(in) :: extrusion
  integer(kind=i_def),   intent(in) :: stencil_depth
  integer(kind=i_def),   intent(in) :: regrid_method

  ! Parameters
  character(len=9), parameter :: routine_name = 'init_mesh'

  integer(kind=i_def), parameter :: dst = 1
  integer(kind=i_def), parameter :: src = 2

  ! Namelist variables
  type(namelist_type), pointer :: lfric2lfric_nml      => null()
  type(namelist_type), pointer :: src_partitioning_nml => null()
  type(namelist_type), pointer :: dst_partitioning_nml => null()

  ! partitioning namelist variables
  logical(l_def)                   :: generate_inner_halos(2)

  ! lfric2lfric namelist variables
  logical(kind=l_def)              :: prepartitioned

  character(len=str_max_filename)  :: meshfile_prefix(2)
  integer(kind=i_def)              :: geometry(2)
  integer(kind=i_def)              :: topology(2)
  integer(kind=i_def)              :: mesh_selection(2)

  ! Local variables
  character(len=str_max_filename)     :: mesh_file(2)

  procedure(partitioner_interface), pointer :: partitioner_src => null()
  procedure(partitioner_interface), pointer :: partitioner_dst => null()

  class(panel_decomposition_type), allocatable :: decomposition_src, &
                                                  decomposition_dst


  !============================================================================
  ! Extract and check configuration variables
  !============================================================================
  ! Read partitioning namelist for source and destination meshes
  src_partitioning_nml  => configuration%get_namelist('partitioning', &
                                                      'source')
  call src_partitioning_nml%get_value( 'generate_inner_halos', &
                                        generate_inner_halos(src) )

  dst_partitioning_nml  => configuration%get_namelist('partitioning', &
                                                      'destination')
  call dst_partitioning_nml%get_value( 'generate_inner_halos', &
                                        generate_inner_halos(dst) )

  ! Read lfric2lfric namelist
  lfric2lfric_nml    => configuration%get_namelist('lfric2lfric')

  call lfric2lfric_nml%get_value( 'prepartitioned_meshes', &
                                   prepartitioned )
  call lfric2lfric_nml%get_value( 'destination_meshfile_prefix', &
                                   meshfile_prefix(dst) )
  call lfric2lfric_nml%get_value( 'destination_geometry', &
                                   geometry(dst) )
  call lfric2lfric_nml%get_value( 'destination_topology', &
                                   topology(dst) )
  call lfric2lfric_nml%get_value( 'source_meshfile_prefix', &
                                   meshfile_prefix(src) )
  call lfric2lfric_nml%get_value( 'source_geometry', &
                                   geometry(src) )
  call lfric2lfric_nml%get_value( 'source_topology', &
                                   topology(src) )

  if ( regrid_method == regrid_method_map .and. &
     trim(meshfile_prefix(src)) /= trim(meshfile_prefix(dst)) ) then

    write( log_scratch_space, '(A)' )                                &
         'When using LFRic intermesh maps, source and destination '//&
         'meshes should be extracted from the same file.'
    call log_event(log_scratch_space, log_level_error)
  end if


  !===========================================================================
  ! Create local mesh objects:
  !  Two code pathes presented, either:
  !  1. The input files have been pre-partitioned.
  !     Meshes and are simply read from file and local mesh objects
  !     are populated.
  !  2. The input files have not been partitioned.
  !     Global meshes are loaded from file and partitioning is applied
  !     at runtime.  NOTE: This option is provided as legacy, and support
  !     is on a best endeavours basis.
  !===========================================================================
  if (prepartitioned) then

    !==========================================================================
    !  Read in local meshes / partition information / mesh maps
    !  direct from file.
    !==========================================================================
    !
    ! For this local rank, a mesh input file with a common base name
    ! of the following form should exist.
    !
    !   <input_basename>_<local_rank>_<total_ranks>.nc
    !
    ! Where 1 rank is assigned to each mesh partition.
    write(mesh_file(dst),'(A,2(I0,A))') &
        trim(meshfile_prefix(dst)) // '_', local_rank, '-', &
                                           total_ranks, '.nc'

    write(mesh_file(src),'(A,2(I0,A))') &
        trim(meshfile_prefix(src)) // '_', local_rank, '-',  &
                                           total_ranks, '.nc'

    ! Read in all local mesh data for this rank and
    ! initialise local mesh objects from them.
    !===========================================================
    ! Each partitioned mesh file will contain meshes of the
    ! same name as all other partitions.
    call log_event( 'Using pre-partitioned mesh file:', log_level_info )
    call log_event( '   '//trim(mesh_file(dst)), log_level_info )
    call log_event( "Loading local mesh(es)", log_level_info )

    if (mesh_file(dst) == mesh_file(src)) then
      call load_local_mesh( mesh_file(dst), mesh_names )
    else
      call load_local_mesh( mesh_file(dst), mesh_names(dst) )

      call log_event( 'Using pre-partitioned mesh file:', log_level_info )
      call log_event( '   '//trim(mesh_file(src)), log_level_info )
      call log_event( "Loading local mesh(es)", log_level_info )

      call load_local_mesh( mesh_file(src), mesh_names(src) )
    endif

    ! Apply configuration related checks to ensure that these
    ! meshes are suitable for the supplied application
    ! configuration.
    !===========================================================
    call check_local_mesh( configuration, &
                           stencil_depth, &
                           mesh_names )

    ! Load and assign mesh maps.
    !===========================================================
    ! Mesh map identifiers are determined by the source/target
    ! mesh IDs they relate to. As a result inter-grid mesh maps
    ! need to be loaded after the relevant local meshes have
    ! been loaded.
    if (regrid_method == regrid_method_map) then
      call load_local_mesh_maps( mesh_file(dst), mesh_names )
    end if
  else

    !==========================================================================
    ! Perform runtime partitioning of global meshes.
    !==========================================================================
    if ( geometry(src) == source_geometry_spherical .and. &
         topology(src) == source_topology_fully_periodic ) then
      mesh_selection(src) = mesh_cubedsphere
      call log_event( "Setting up cubed-sphere partition mesh(es)", &
                      log_level_debug )
    else
      mesh_selection(src) = mesh_planar
      call log_event( "Setting up planar partition mesh(es)", &
                      log_level_debug )
    end if

    if ( geometry(dst) == destination_geometry_spherical .and. &
         topology(dst) == destination_topology_fully_periodic ) then
      mesh_selection(dst) = mesh_cubedsphere
      call log_event( "Setting up cubed-sphere partition mesh(es)", &
                      log_level_debug )
    else
      mesh_selection(dst) = mesh_planar
      call log_event( "Setting up planar partition mesh(es)", &
                      log_level_debug )
    end if

    call log_event( "Setting up partition mesh(es)", log_level_info )
    write(mesh_file(src),'(A)') trim(meshfile_prefix(src)) // '.nc'
    write(mesh_file(dst),'(A)') trim(meshfile_prefix(dst)) // '.nc'

    ! Set constants that will control partitioning.
    !===========================================================
    call get_partition_parameters( src_partitioning_nml, &
                                   mesh_selection(src),  &
                                   total_ranks,          &
                                   decomposition_src,    &
                                   partitioner_src )

    call get_partition_parameters( dst_partitioning_nml, &
                                   mesh_selection(dst),  &
                                   total_ranks,          &
                                   decomposition_dst,    &
                                   partitioner_dst )


    ! Read in all global meshes from input file
    !===========================================================
    if (mesh_file(dst) == mesh_file(src)) then
      call load_global_mesh( mesh_file(dst), mesh_names )
    else
      call load_global_mesh( mesh_file(dst), mesh_names(dst) )
      call load_global_mesh( mesh_file(src), mesh_names(src) )
    endif

    if (regrid_method == regrid_method_map) then

      ! Apply configuration related checks to ensure that these
      ! meshes are suitable for the supplied application
      ! configuration.
      !===========================================================
      call check_global_mesh( configuration, mesh_names )

    end if

    ! Partition the global meshes
    !===========================================================
    call create_local_mesh( mesh_names(dst:dst),           &
                            local_rank, total_ranks,       &
                            decomposition_dst,             &
                            stencil_depth,                 &
                            generate_inner_halos(dst),    &
                            partitioner_dst )

    call create_local_mesh( mesh_names(src:src),           &
                            local_rank, total_ranks,       &
                            decomposition_src,             &
                            stencil_depth,                 &
                            generate_inner_halos(src),    &
                            partitioner_src )

    ! Read in the global intergrid mesh mappings,
    ! then create the associated local mesh maps
    !===========================================================
    if (regrid_method == regrid_method_map) then
      call create_local_mesh_maps( mesh_file(dst) )
    end if

  end if  ! prepartitioned

  !============================================================================
  ! Extrude the specified meshes from local mesh objects into
  ! mesh objects on the given extrusion.
  ! Alternative names are needed in case the source and destination
  ! mesh files use the same mesh name.
  !============================================================================
  call create_mesh( mesh_names, extrusion )

  !============================================================================
  ! Generate intergrid LiD-LiD maps and assign them to mesh objects.
  !============================================================================
  if (regrid_method == regrid_method_map) then
    call assign_mesh_maps(mesh_names)
  end if

end subroutine init_mesh

end module lfric2lfric_init_mesh_mod
