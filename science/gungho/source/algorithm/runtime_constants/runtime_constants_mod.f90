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
  use boundaries_config_mod,   only: limited_area
  use constants_mod,           only: i_def, r_def, str_def, l_def
  use extrusion_mod,           only: PRIME_EXTRUSION, SHIFTED, &
                                     TWOD, DOUBLE_LEVEL
  use field_mod,               only: field_type
  use formulation_config_mod,  only: l_multigrid, use_multires_coupling
  use inventory_by_mesh_mod,   only: inventory_by_mesh_type
  use io_config_mod,           only: subroutine_timers
  use log_mod,                 only: log_event, LOG_LEVEL_INFO, &
                                     LOG_LEVEL_ERROR
  use mesh_collection_mod,     only: mesh_collection_type
  use mesh_mod,                only: mesh_type
  use model_clock_mod,         only: model_clock_type
  use multigrid_config_mod,    only: chain_mesh_tags
  use multires_coupling_config_mod, &
                               only: multires_coupling_mesh_tags
  use runtime_tools_mod,       only: primary_mesh_label,      &
                                     shifted_mesh_label,      &
                                     double_level_mesh_label, &
                                     twod_mesh_label,         &
                                     multigrid_mesh_label,    &
                                     extra_mesh_label
  use timer_mod,               only: timer

  implicit none

  private

  ! Public functions to create and access the module contents
  public :: create_runtime_constants
  public :: final_runtime_constants

contains
  !>@brief Subroutine to create the runtime constants
  !> @param[in] mesh_collection      The collection of meshes to loop over
  !> @param[in] chi_inventory        Stores the coordinate fields by their mesh
  !> @param[in] panel_id_inventory   Stores the panel id fields by their mesh
  !> @param[in] model_clock          Model time.
  !> @param[in] create_rdef_div_flag Optional flag as to whether rdef div
  !!                                 operators need creating.
  subroutine create_runtime_constants(mesh_collection,       &
                                      chi_inventory,         &
                                      panel_id_inventory,    &
                                      model_clock,           &
                                      create_rdef_div_flag)

    ! Other runtime_constants modules
    use wt_advective_update_alg_mod,  only: wt_advective_update_set_num_meshes
    use fem_constants_mod,            only: create_fem_constants
    use geometric_constants_mod,      only: create_geometric_constants
    use intermesh_constants_mod,      only: create_intermesh_constants
    use limited_area_constants_mod,   only: create_limited_area_constants
    use physical_op_constants_mod,    only: create_physical_op_constants
    use reconstruct_w3_field_alg_mod, only: reconstruct_w3_field_alg_set_num_meshes
    use runge_kutta_init_mod,         only: runge_kutta_init
    use runtime_tools_mod,            only: init_mesh_id_list, &
                                            init_hierarchical_mesh_id_list

    implicit none

    type(mesh_collection_type),    intent(in) :: mesh_collection
    type(inventory_by_mesh_type),  intent(in) :: chi_inventory
    type(inventory_by_mesh_type),  intent(in) :: panel_id_inventory
    class(model_clock_type),       intent(in) :: model_clock
    logical(kind=l_def), optional, intent(in) :: create_rdef_div_flag

    type(mesh_type),  pointer :: mesh => null()
    type(field_type), pointer :: chi(:) => null()
    type(field_type), pointer :: panel_id => null()
    type(mesh_type),  pointer :: prime_extrusion_mesh => null()

    ! Internal variables
    character(str_def),  allocatable :: all_mesh_names(:)
    integer(kind=i_def)              :: num_meshes, i, j, k
    integer(kind=i_def), allocatable :: mesh_id_list(:)
    integer(kind=i_def), allocatable :: mg_mesh_ids(:)
    integer(kind=i_def), allocatable :: label_list(:)
    type(field_type),    allocatable :: chi_list(:,:)
    type(field_type),    allocatable :: panel_id_list(:)
    logical(kind=l_def)              :: is_multigrid_mesh
    logical(kind=l_def)              :: is_multires_mesh
    integer(kind=i_def)              :: num_mg_meshes
    logical(kind=l_def)              :: create_rdef_div_operators

    if ( subroutine_timers ) call timer('runtime_constants_alg')
    call log_event( "Gungho: creating runtime_constants", LOG_LEVEL_INFO )

    !==========================================================================!
    ! Turn all the meshes and coordinate fields into lists
    !==========================================================================!

    ! To loop through mesh collection, get all mesh names
    ! Then get mesh from collection using these names
    all_mesh_names = mesh_collection%get_mesh_names()
    num_meshes = SIZE(all_mesh_names)

    allocate(mesh_id_list(num_meshes))
    allocate(chi_list(3,num_meshes))
    allocate(panel_id_list(num_meshes))
    allocate(label_list(num_meshes))

    ! Populate these lists
    do i = 1, num_meshes
      mesh => mesh_collection%get_mesh(all_mesh_names(i))
      if (mesh%get_extrusion_id() == TWOD) then
        prime_extrusion_mesh => mesh_collection%get_mesh(mesh, PRIME_EXTRUSION)
        call chi_inventory%get_field_array(prime_extrusion_mesh, chi)
        call panel_id_inventory%get_field(prime_extrusion_mesh, panel_id)
      else
        call chi_inventory%get_field_array(mesh, chi)
        call panel_id_inventory%get_field(mesh, panel_id)
      end if

      ! Copy mesh id, chi field and panel_id into lists
      mesh_id_list(i) = mesh%get_id()
      call panel_id%copy_field_serial(panel_id_list(i))
      do j = 1, 3
        call chi(j)%copy_field_serial(chi_list(j,i))
      end do

      ! Label meshes based on their role
      select case (mesh%get_extrusion_id())
      case (SHIFTED)
        label_list(i) = shifted_mesh_label
      case (DOUBLE_LEVEL)
        label_list(i) = double_level_mesh_label
      case (TWOD)
        label_list(i) = twod_mesh_label
      case (PRIME_EXTRUSION)
        ! Now determine label based on namelist options
        if (all_mesh_names(i) == prime_mesh_name) then
          label_list(i) = primary_mesh_label
        else if (l_multigrid) then
          ! Search through chain mesh tags for matching mesh name
          is_multigrid_mesh = .false.
          is_multires_mesh = .false.
          do j = 1, SIZE(chain_mesh_tags)
            if (all_mesh_names(i) == chain_mesh_tags(j)) then
              ! Now check if it is used in multires-coupling: then it should
              ! be an "extra" mesh
              if (use_multires_coupling) then
                do k = 1, SIZE(multires_coupling_mesh_tags)
                  if (all_mesh_names(i) == multires_coupling_mesh_tags(k)) then
                    is_multires_mesh = .true.
                    exit
                  end if
                end do
              end if
              if (is_multires_mesh) then
                label_list(i) = extra_mesh_label
              else
                is_multigrid_mesh = .true.
                label_list(i) = multigrid_mesh_label
              end if
              exit
            end if
          end do
          ! If it was not a multigrid mesh, it must be an extra mesh
          if (.not. is_multigrid_mesh) label_list(i) = extra_mesh_label
        else
          label_list(i) = extra_mesh_label
        end if
      case default
        call log_event('Mesh extrusion not identified', LOG_LEVEL_ERROR)
      end select
    end do

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

    call init_mesh_id_list(mesh_id_list)

    if ( l_multigrid ) then
      ! mg_mesh_ids contains all mesh ids used in the multigrid chain
      ! including the primary mesh
      call init_hierarchical_mesh_id_list(mg_mesh_ids)
    else
      ! Just create a list with the primary mesh id in it
      mesh => mesh_collection%get_mesh(prime_mesh_name)
      call init_hierarchical_mesh_id_list( (/ mesh%get_id() /) )
    end if

    call create_geometric_constants(mesh_id_list,      &
                                    chi_list,          &
                                    panel_id_list,     &
                                    label_list         )

    ! Finite element constants should be created after geometric constants
    ! The chi fields set up in geometric constants are used here
    if (present(create_rdef_div_flag)) then
      create_rdef_div_operators = create_rdef_div_flag
    else
      create_rdef_div_operators = .false.
    end if

    call create_fem_constants(mesh_id_list,  &
                              chi_list,      &
                              panel_id_list, &
                              label_list,    &
                              model_clock,   &
                              create_rdef_div_operators)

    call create_physical_op_constants(mesh_id_list,  &
                                      chi_list,      &
                                      panel_id_list, &
                                      label_list,    &
                                      model_clock )

    if ( limited_area ) then
      call create_limited_area_constants(mesh_id_list, &
                                         chi_list,     &
                                         label_list )
    end if

    call create_intermesh_constants(mesh_collection, &
                                    use_multires_coupling)

    ! Set-up arrays for transport coefficients
    call runge_kutta_init()
    call reconstruct_w3_field_alg_set_num_meshes( num_meshes )
    call wt_advective_update_set_num_meshes( num_meshes )

    deallocate(mesh_id_list)
    deallocate(label_list)

    call log_event( "Gungho: created runtime_constants", LOG_LEVEL_INFO )
    if ( subroutine_timers ) call timer('runtime_constants_alg')

  end subroutine create_runtime_constants


  !> @brief Explicitly reclaim memory from module scope variables
  !
  subroutine final_runtime_constants()

    ! Other runtime_constants modules
    use fem_constants_mod,           only: final_fem_constants
    use geometric_constants_mod,     only: final_geometric_constants
    use intermesh_constants_mod,     only: final_intermesh_constants
    use limited_area_constants_mod,  only: final_limited_area_constants
    use physical_op_constants_mod,   only: final_physical_op_constants
    use runtime_tools_mod,           only: final_mesh_id_list, &
                                           final_hierarchical_mesh_id_list

    implicit none

    call final_geometric_constants()
    call final_fem_constants()
    call final_physical_op_constants()
    if ( limited_area ) call final_limited_area_constants()
    call final_intermesh_constants()
    call final_hierarchical_mesh_id_list()
    call final_mesh_id_list()


  end subroutine final_runtime_constants

end module runtime_constants_mod
