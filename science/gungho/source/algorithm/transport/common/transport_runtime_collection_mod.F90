!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!>
!> @brief Holds and manages transport runtime objects created during a model run.
!>
!> @details A container which holds type definition of a collection of
!!          transport runtime objects. Provides functionality to return
!!          a runtime for a given underlying local mesh.
!!          TODO #3008: this should be a collection object based on a
!!          linked_list_type and not using array structures to manage the
!!          transport runtimes
module transport_runtime_collection_mod

  use constants_mod,              only : i_def, l_def, imdi
  use local_mesh_mod,             only : local_mesh_type
  use log_mod,                    only : log_event, log_scratch_space,    &
                                         LOG_LEVEL_ERROR
  use mesh_mod,                   only : mesh_type
  use transport_runtime_alg_mod,  only : transport_runtime_type
  use field_mod,                  only : field_type
  use model_clock_mod,            only : model_clock_type

  implicit none

  private

  integer(kind=i_def),          allocatable         :: local_mesh_id_list(:)
  type(transport_runtime_type), allocatable, target :: transport_runtime_list(:)

  public  :: init_transport_runtime_collection
  public  :: get_transport_runtime
  public  :: set_transport_runtime
  public  :: transport_runtime_collection_final
  private :: idx_from_local_mesh_id

contains

  !> @brief Initialises the arrays of transport runtimes
  !> TODO #3008 this should be replaced with the object constructor
  !> @param[in] local_mesh_ids  An array of local mesh IDs, each corresponding
  !!                            to a mesh with a transport_runtime object
  subroutine init_transport_runtime_collection(local_mesh_ids)

    implicit none

    integer(kind=i_def), intent(in) :: local_mesh_ids(:)
    integer(kind=i_def)             :: num_meshes

    num_meshes = size(local_mesh_ids)

    if (allocated(transport_runtime_list)) deallocate(transport_runtime_list)
    if (allocated(local_mesh_id_list)) deallocate(local_mesh_id_list)

    allocate(transport_runtime_list(num_meshes))
    allocate(local_mesh_id_list(num_meshes))

    local_mesh_id_list(:) = local_mesh_ids(:)

  end subroutine init_transport_runtime_collection

  !-----------------------------------------------------------------------------
  ! Get a runtime object
  !-----------------------------------------------------------------------------
  !> Function to get a transport_runtime instance
  !> @param[in] mesh  The mesh object for which to find the transport runtime
  function get_transport_runtime(mesh) result(transport_runtime)

    implicit none

    type(mesh_type),           intent(in) :: mesh

    type(transport_runtime_type), pointer :: transport_runtime
    type(local_mesh_type),        pointer :: local_mesh
    integer(kind=i_def)                   :: local_mesh_id, idx

    ! Obtain local mesh id
    local_mesh => mesh%get_local_mesh()
    local_mesh_id = local_mesh%get_id()
    idx = idx_from_local_mesh_id(local_mesh_id)

    transport_runtime => transport_runtime_list(idx)

    nullify(local_mesh)

  end function get_transport_runtime

  !-----------------------------------------------------------------------------
  ! Set a runtime object
  !-----------------------------------------------------------------------------
  !> Function to set an instance of a transport_runtime type
  !> @param[in] primal_wind_n_rdef    Wind field at start of current time step on the
  !!                                  prime extrusion
  !> @param[in] primal_wind_np1_rdef  Predictor of wind field at next time step on
  !!                                  the prime extrusion
  !> @param[in] model_clock           Time within the model
  !> @param[in] outer                 Outer (advection) iteration number
  !> @param[in] cheap_update          Logical flag for cheap transport update
  !> @param[in] shifted_mesh          (Optional) Shifted mesh
  subroutine set_transport_runtime(primal_wind_n_rdef,   &
                                   primal_wind_np1_rdef, &
                                   model_clock,          &
                                   outer,                &
                                   cheap_update,         &
                                   shifted_mesh)

    implicit none

    type(field_type),        target, intent(in) :: primal_wind_n_rdef
    type(field_type),        target, intent(in) :: primal_wind_np1_rdef
    class(model_clock_type), target, intent(in) :: model_clock
    integer(kind=i_def),             intent(in) :: outer
    logical(kind=l_def),             intent(in) :: cheap_update

    type(mesh_type),  pointer, intent(in), optional :: shifted_mesh

    integer(kind=i_def)            :: local_mesh_id, idx
    type(mesh_type),       pointer :: mesh => null()
    type(local_mesh_type), pointer :: local_mesh => null()

    mesh => primal_wind_n_rdef%get_mesh()
    local_mesh => mesh%get_local_mesh()
    local_mesh_id = local_mesh%get_id()
    idx = idx_from_local_mesh_id(local_mesh_id)

    if ( present( shifted_mesh ) ) then
      call transport_runtime_list(idx)%initialise(primal_wind_n_rdef,   &
                                                  primal_wind_np1_rdef, &
                                                  model_clock,          &
                                                  outer,                &
                                                  cheap_update,         &
                                                  shifted_mesh)
    else
      call transport_runtime_list(idx)%initialise(primal_wind_n_rdef,   &
                                                  primal_wind_np1_rdef, &
                                                  model_clock,          &
                                                  outer,                &
                                                  cheap_update )
    end if

  end subroutine set_transport_runtime


  !> Routine to finalise collection
  subroutine transport_runtime_collection_final()

    implicit none

    if (allocated(transport_runtime_list)) deallocate(transport_runtime_list)
    if (allocated(local_mesh_id_list)) deallocate(local_mesh_id_list)

  end subroutine transport_runtime_collection_final

  !> @brief Private function for getting array index corresponding to mesh_id
  !> TODO #3008: this should be removed
  !> @param[in] local_mesh_id   Identifier of local mesh to find
  function idx_from_local_mesh_id(local_mesh_id) result(idx)

    implicit none

    integer(kind=i_def), intent(in) :: local_mesh_id
    integer(kind=i_def)             :: idx, i

    idx = imdi
    do i = 1, size(local_mesh_id_list)
      if ( local_mesh_id == local_mesh_id_list(i) ) idx = i
    end do

    if ( idx == imdi ) then
      write(log_scratch_space, '(A,I4)') &
      'transport_runtime collection does not contain mesh: ', local_mesh_id
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if

  end function idx_from_local_mesh_id

end module transport_runtime_collection_mod
