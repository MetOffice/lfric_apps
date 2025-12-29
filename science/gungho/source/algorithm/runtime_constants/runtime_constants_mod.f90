!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!> @brief A module that controls set-up of various run time constants.
!>
!> @details This module controls the set-up of various objects that are
!>          created at setup and are not changed thereafter but are needed
!>          throughout the algorithm layers.
module runtime_constants_mod

  use base_mesh_config_mod,    only: prime_mesh_name
  use constants_mod,           only: i_def, l_def
  use formulation_config_mod,  only: l_multigrid
  use log_mod,                 only: log_event, LOG_LEVEL_INFO, &
                                     LOG_LEVEL_ERROR
  use mesh_collection_mod,     only: mesh_collection
  use mesh_mod,                only: mesh_type
  use model_clock_mod,         only: model_clock_type
  use multigrid_config_mod,    only: chain_mesh_tags
  use timing_mod,              only: start_timing, stop_timing, tik, LPROF
  implicit none

  private

  ! Public functions to create and access the module contents
  public :: create_runtime_constants
  public :: final_runtime_constants

contains
  !> @brief Subroutine to create the runtime constants
  subroutine create_runtime_constants()

    use runge_kutta_init_mod,         only: runge_kutta_init
    use runtime_tools_mod,            only: init_hierarchical_mesh_id_list

    implicit none

    type(mesh_type),  pointer :: mesh

    ! Internal variables
    integer(kind=i_def)              :: i
    integer(kind=i_def), allocatable :: mg_mesh_ids(:)
    integer(kind=i_def)              :: num_mg_meshes
    integer(tik)                     :: id

    if ( LPROF ) call start_timing( id, 'runtime_constants_alg' )

    !==========================================================================!
    ! Turn all the meshes and coordinate fields into lists
    !==========================================================================!

    if ( l_multigrid ) then
      num_mg_meshes = SIZE(chain_mesh_tags)
      allocate(mg_mesh_ids(num_mg_meshes))
      do i = 1, num_mg_meshes
        mesh => mesh_collection%get_mesh(chain_mesh_tags(i))
        mg_mesh_ids(i) = mesh%get_id()
      end do
    end if

    !==========================================================================!
    ! Set up runtime_constants for each category
    !==========================================================================!

    if ( l_multigrid ) then
      ! mg_mesh_ids contains all mesh ids used in the multigrid chain
      ! including the primary mesh
      call init_hierarchical_mesh_id_list(mg_mesh_ids)
    else
      ! Just create a list with the primary mesh id in it
      mesh => mesh_collection%get_mesh(prime_mesh_name)
      call init_hierarchical_mesh_id_list( (/ mesh%get_id() /) )
    end if

    ! Set-up arrays for transport coefficients
    ! @TODO: can this be moved to somewhere more in line with the other
    ! code structure
    call runge_kutta_init()

    if ( LPROF ) call stop_timing( id, 'runtime_constants_alg' )

  end subroutine create_runtime_constants


  !> @brief Explicitly reclaim memory from module scope variables
  !
  subroutine final_runtime_constants()

    use dycore_constants_mod,        only: final_dycore_constants
    use sci_fem_constants_mod,       only: final_fem_constants
    use sci_geometric_constants_mod, only: final_geometric_constants
    use limited_area_constants_mod,  only: final_limited_area_constants
    use sci_mapping_constants_mod,   only: final_mapping_constants
    use physics_constants_mod,       only: final_physics_constants
    use solver_constants_mod,        only: final_solver_constants
    use transport_constants_mod,     only: final_transport_constants
    use runtime_tools_mod,           only: final_hierarchical_mesh_id_list

    implicit none

    call final_geometric_constants()
    call final_fem_constants()
    call final_physics_constants()
    call final_limited_area_constants()
    call final_mapping_constants()
    call final_dycore_constants()
    call final_transport_constants()
    call final_solver_constants()
    call final_hierarchical_mesh_id_list()


  end subroutine final_runtime_constants

end module runtime_constants_mod
