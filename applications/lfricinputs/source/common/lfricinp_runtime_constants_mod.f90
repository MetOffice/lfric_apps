! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_runtime_constants_mod
  use constants_mod,                     only: i_def, r_def, str_def
  use extrusion_mod,                     only: PRIME_EXTRUSION, TWOD
  use field_mod,                         only: field_type
  use inventory_by_mesh_mod,             only: inventory_by_mesh_type
  use runtime_tools_mod,                 only: primary_mesh_label,      &
                                               twod_mesh_label
  use log_mod,                           only: log_event, LOG_LEVEL_ERROR
  use mesh_mod,                          only: mesh_type
  use mesh_collection_mod,               only: mesh_collection_type

IMPLICIT NONE
PRIVATE
PUBLIC :: lfricinp_create_runtime_constants

CONTAINS

SUBROUTINE lfricinp_create_runtime_constants(mesh_collection,       &
                                             chi_inventory,         &
                                             panel_id_inventory,    &
                                             dt)
! Description:
  ! This routine contains a subset of the contents of create_runtime_constants
  ! as required by LFRic Inputs. Ideally we would use the centrally provided
  ! function, but this brings in infrastructure that is causing errors that need
  ! time to debug.

    use geometric_constants_mod,     only: create_geometric_constants
    use runtime_tools_mod,           only: init_mesh_id_list

    implicit none

    type(mesh_collection_type),   intent(in) :: mesh_collection
    type(inventory_by_mesh_type), intent(in) :: chi_inventory
    type(inventory_by_mesh_type), intent(in) :: panel_id_inventory
    real(r_def),                  intent(in) :: dt

    ! Internal variables
    character(str_def),  allocatable :: all_mesh_names(:)
    integer(kind=i_def)              :: num_meshes, i, j
    integer(kind=i_def), allocatable :: mesh_id_list(:)
    integer(kind=i_def), allocatable :: label_list(:)
    type(field_type),    allocatable :: chi_list(:,:)
    type(field_type),    allocatable :: panel_id_list(:)
    type(mesh_type),         pointer :: mesh => null()
    type(field_type),        pointer :: chi(:) => null()
    type(field_type),        pointer :: panel_id => null()
    type(mesh_type),         pointer :: prime_extrusion_mesh => null()


    !==========================================================================!
    ! Turn all the mesh IDs and coordinate fields into lists
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
      case (TWOD)
        label_list(i) = twod_mesh_label
      case (PRIME_EXTRUSION)
        label_list(i) = primary_mesh_label
      case default
        call log_event('Mesh extrusion not implemented in lfricinputs', LOG_LEVEL_ERROR)
      end select
    end do

    !==========================================================================!
    ! Set up runtime_constants for each category
    !==========================================================================!

    call init_mesh_id_list(mesh_id_list)

    call create_geometric_constants(mesh_id_list,      &
                                    chi_list,          &
                                    panel_id_list,     &
                                    label_list         )

    deallocate(mesh_id_list)
    deallocate(label_list)

END SUBROUTINE lfricinp_create_runtime_constants

END MODULE lfricinp_runtime_constants_mod
