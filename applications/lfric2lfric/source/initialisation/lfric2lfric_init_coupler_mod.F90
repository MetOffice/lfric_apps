!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Initialisation functionality for the coupled miniapp

!> @details Uses the coupler object to initialise the coupling
!>          partitions and variables

module lfric2lfric_init_coupler_mod

  use constants_mod,                  only : i_def, str_def, &
                                             l_def, radians_to_degrees
#ifdef MCT
  use coupling_mod,                   only : coupling_type,       &
                                             get_coupling_from_collection
#endif
  use driver_modeldb_mod,             only : modeldb_type
  use field_collection_mod,           only : field_collection_type
  use field_mod,                      only : field_type
  use function_space_collection_mod,  only : function_space_collection
  use fs_continuity_mod,              only : W3
  use mesh_mod,                       only : mesh_type
  use pure_abstract_field_mod,        only : pure_abstract_field_type
  use log_mod,                        only : log_event,         &
                                             log_scratch_space, &
                                             log_level_error

  implicit none

  contains

#ifdef MCT
  !> @details Initialises coupling partitions and variables for the
  !>          source fields
  !> @param[in]     mesh Representation of the mesh the code will run on
  !> @param[in]     coupler Name of the coupling component
  !> @param[in]     element_order_h The polynomial order of the set of
  !>                                compatible finite element in the
  !>                                horizontal
  !> @param[in]     element_order_v The polynomial order of the set of
  !>                                compatible finite element in the
  !>                                vertical
  !> @param[in,out] modeldb The structure that holds model state
  subroutine lfric2lfric_init_coupler_src(mesh, coupler, element_order_h, &
                                          element_order_v, modeldb)

    implicit none

    type(mesh_type), pointer, intent(in)     :: mesh
    character(*),             intent(in)     :: coupler
    integer(i_def),           intent(in)     :: element_order_h
    integer(i_def),           intent(in)     :: element_order_v
    type(modeldb_type),       intent(inout)  :: modeldb

    ! Field name
    character(str_def)                   :: name_cpl

    ! Coupling fields
    type(field_collection_type), pointer :: cpl_snd_2d => null()
    type(field_collection_type), pointer :: cpl_rcv_2d => null()
    type(field_collection_type), pointer :: cpl_snd_0d => null()

    type(field_type)                     :: field_cpl
    type(field_collection_type), pointer :: depository
    type(field_type), pointer            :: field
    class(pure_abstract_field_type), pointer :: field_ptr

    type(coupling_type), pointer         :: coupling_ptr


    ! Create prognostic fields
    ! Creates a field in the W3 function space (fully discontinuous field)
    name_cpl = "src_face"
    call field_cpl%initialise( vector_space = &
                    function_space_collection%get_fs(mesh,             &
                                element_order_h, element_order_v, W3), &
                                name=name_cpl)

    ! Add fields to modeldb
    depository => modeldb%fields%get_field_collection("depository")
    call depository%add_field(field_cpl)

    ! Set up collections to hold 2d coupling fields
    if (.not. modeldb%fields%field_collection_exists("cpl_snd_2d"))  &
      call modeldb%fields%add_empty_field_collection("cpl_snd_2d" ,  &
                                                table_len = 30)
    if (.not. modeldb%fields%field_collection_exists("dummy"))  &
      call modeldb%fields%add_empty_field_collection("dummy" ,  &
                                                table_len = 30)
    if (.not. modeldb%fields%field_collection_exists("cpl_snd_0d"))  &
      call modeldb%fields%add_empty_field_collection("cpl_snd_0d" ,  &
                                                table_len = 30)

    cpl_snd_2d => modeldb%fields%get_field_collection("cpl_snd_2d")
    cpl_rcv_2d => modeldb%fields%get_field_collection("dummy")
    cpl_snd_0d => modeldb%fields%get_field_collection("cpl_snd_0d")

    ! Set lfric fields as either the source entry
    call depository%get_field(trim(name_cpl), field)
    ! Add that field to the coupling send field collection for 2d fields
    field_ptr => field
    call cpl_snd_2d%add_reference_to_field(field_ptr)

    ! Extract the coupling object from the modeldb key-value pair collection
    coupling_ptr => get_coupling_from_collection(modeldb%values, trim(coupler))
    call coupling_ptr%define_partitions(modeldb%mpi,mesh)
    call coupling_ptr%define_variables( cpl_snd_2d, &
                                        cpl_rcv_2d, &
                                        cpl_snd_0d )

  end subroutine lfric2lfric_init_coupler_src


  !> @details Initialises coupling partitions and variables for the
  !>          destination fields
  !> @param[in]     mesh Representation of the mesh the code will run on
  !> @param[in]     coupler Name of the coupling component
  !> @param[in]     element_order_h The polynomial order of the set of
  !>                                compatible finite element in the
  !>                                horizontal
  !> @param[in]     element_order_v The polynomial order of the set of
  !>                                compatible finite element in the
  !>                                vertical
  !> @param[in,out] modeldb The structure that holds model state
  subroutine lfric2lfric_init_coupler_dst(mesh, coupler, element_order_h,    &
                                          element_order_v, modeldb)

    implicit none

    type(mesh_type), pointer, intent(in)     :: mesh
    character(*),             intent(in)     :: coupler
    integer(i_def),           intent(in)     :: element_order_h
    integer(i_def),           intent(in)     :: element_order_v
    type(modeldb_type),       intent(inout)  :: modeldb

    ! Field name
    character(str_def)                   :: name_cpl

    ! Coupling fields
    type(field_collection_type), pointer :: cpl_snd_2d => null()
    type(field_collection_type), pointer :: cpl_rcv_2d => null()
    type(field_collection_type), pointer :: cpl_snd_0d => null()

    type(field_type)                     :: field_cpl
    type(field_collection_type), pointer :: depository
    type(field_type), pointer            :: field
    class(pure_abstract_field_type), pointer :: field_ptr

    type(coupling_type), pointer         :: coupling_ptr


    ! Create prognostic fields
    ! Creates a field in the W3 function space (fully discontinuous field)
    name_cpl = "dst_face"
    call field_cpl%initialise( vector_space = &
                    function_space_collection%get_fs(mesh,             &
                                element_order_h, element_order_v, W3), &
                                name=name_cpl)

    ! Add fields to modeldb
    depository => modeldb%fields%get_field_collection("depository")
    call depository%add_field(field_cpl)

    ! Set up collections to hold 2d coupling fields
    if (.not. modeldb%fields%field_collection_exists("cpl_rcv_2d"))  &
      call modeldb%fields%add_empty_field_collection("cpl_rcv_2d" ,  &
                                                table_len = 30)
    if (.not. modeldb%fields%field_collection_exists("cpl_snd_0d"))  &
      call modeldb%fields%add_empty_field_collection("cpl_snd_0d" ,  &
                                                table_len = 30)
    if (.not. modeldb%fields%field_collection_exists("dummy"))       &
      call modeldb%fields%add_empty_field_collection("dummy" ,       &
                                                table_len = 30)

    cpl_snd_2d => modeldb%fields%get_field_collection("dummy")
    cpl_rcv_2d => modeldb%fields%get_field_collection("cpl_rcv_2d")
    cpl_snd_0d => modeldb%fields%get_field_collection("cpl_snd_0d")

    ! Set lfric fields as either the destination entry
    call depository%get_field(trim(name_cpl), field)
    ! Add that field to the coupling send field collection for 2d fields
    field_ptr => field
    call cpl_rcv_2d%add_reference_to_field(field_ptr)

    ! Extract the coupling object from the modeldb key-value pair collection
    coupling_ptr => get_coupling_from_collection(modeldb%values, trim(coupler))
    call coupling_ptr%define_partitions(modeldb%mpi,mesh)
    call coupling_ptr%define_variables( cpl_snd_2d, &
                                        cpl_rcv_2d, &
                                        cpl_snd_0d )

  end subroutine lfric2lfric_init_coupler_dst

  !> @details Ends the model coupling definition
  !> @param[in,out] modeldb The structure that holds model state
  !> @param[in]     coupler Name of the coupling component
  subroutine lfric2lfric_end_coupler_init(modeldb, coupler)

    implicit none

    type(modeldb_type),       intent(inout)  :: modeldb
    character(*),             intent(in)     :: coupler

    ! Coupling fields
    type(field_collection_type), pointer :: cpl_snd_2d => null()
    type(field_collection_type), pointer :: cpl_rcv_2d => null()
    type(field_collection_type), pointer :: cpl_snd_0d => null()

    type(coupling_type), pointer         :: coupling_ptr


    ! Get coupling field collections
    cpl_snd_2d => modeldb%fields%get_field_collection("cpl_snd_2d")
    cpl_rcv_2d => modeldb%fields%get_field_collection("cpl_rcv_2d")
    cpl_snd_0d => modeldb%fields%get_field_collection("cpl_snd_0d")

    ! Extract the coupling object from the modeldb key-value pair collection
    coupling_ptr => get_coupling_from_collection(modeldb%values, trim(coupler))
    call coupling_ptr%end_definition( cpl_snd_2d, &
                                      cpl_rcv_2d, &
                                      cpl_snd_0d )

  end subroutine lfric2lfric_end_coupler_init
#else
  subroutine lfric2lfric_init_coupler_src(mesh, coupler, element_order_h, &
                                          element_order_v, modeldb)

    implicit none

    type(mesh_type), pointer, intent(in)     :: mesh
    character(*),             intent(in)     :: coupler
    integer(i_def),           intent(in)     :: element_order_h
    integer(i_def),           intent(in)     :: element_order_v
    type(modeldb_type),       intent(inout)  :: modeldb

    write(log_scratch_space, '(A)' ) &
        "lfric2lfric_init_coupler_src: to use OASIS, " // &
        "cpp directive MCT must be set"
    call log_event( log_scratch_space, log_level_error )

  end subroutine lfric2lfric_init_coupler_src

  subroutine lfric2lfric_init_coupler_dst(mesh, coupler, element_order_h, &
                                          element_order_v, modeldb)

    implicit none

    type(mesh_type), pointer, intent(in)     :: mesh
    character(*),             intent(in)     :: coupler
    integer(i_def),           intent(in)     :: element_order_h
    integer(i_def),           intent(in)     :: element_order_v
    type(modeldb_type),       intent(inout)  :: modeldb

    write(log_scratch_space, '(A)' ) &
        "lfric2lfric_init_coupler_dst: to use OASIS, " // &
        "cpp directive MCT must be set"
    call log_event( log_scratch_space, log_level_error )

  end subroutine lfric2lfric_init_coupler_dst

  subroutine lfric2lfric_end_coupler_init(modeldb, coupler)

    implicit none

    type(modeldb_type),       intent(inout)  :: modeldb
    character(*),             intent(in)     :: coupler

    write(log_scratch_space, '(A)' ) &
        "lfric2lfric_end_coupler_init: to use OASIS, " // &
        "cpp directive MCT must be set"
    call log_event( log_scratch_space, log_level_error )

  end subroutine lfric2lfric_end_coupler_init
#endif

end module lfric2lfric_init_coupler_mod
