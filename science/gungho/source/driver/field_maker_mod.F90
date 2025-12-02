!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Infrastructure for creating fields for use in
!> creating prognostic fields
!
module field_maker_mod

  use constants_mod,                  only : i_def, l_def, str_def, r_def, &
                                             r_second
  use log_mod,                        only : log_event, log_scratch_space,     &
                                             log_level_error
  use field_mod,                      only : field_type
  use integer_field_mod,              only : integer_field_type
  use pure_abstract_field_mod,        only : pure_abstract_field_type
  use field_parent_mod,               only : write_interface, read_interface,  &
                                             checkpoint_write_interface,       &
                                             checkpoint_read_interface
  use field_collection_mod,           only : field_collection_type
  use field_mapper_mod,               only : field_mapper_type
  use field_from_metadata_mod,        only : init_field_from_metadata
  use function_space_mod,             only : function_space_type
  use function_space_collection_mod,  only : function_space_collection
  use mesh_mod,                       only : mesh_type
  use clock_mod,                      only : clock_type
  use model_clock_mod,                only : model_clock_type
  use field_spec_mod,                 only : main_coll_dict,                   &
                                             adv_coll_dict,                    &
                                             moist_arr_dict,                   &
                                             time_axis_dict,                   &
                                             field_spec_type,                  &
                                             processor_type,                   &
                                             missing_fs,                       &
                                             space_has_xios_io
  use lfric_xios_diag_mod,            only : field_is_valid
  use lfric_xios_time_axis_mod,       only : time_axis_type
  use lfric_xios_write_mod,           only : create_checkpoint_list
  use io_config_mod,                  only : use_xios_io, &
                                             checkpoint_write, &
                                             checkpoint_read, &
                                             checkpoint_times, &
                                             end_of_run_checkpoint
  use initialization_config_mod,      only : init_option,                      &
                                             init_option_checkpoint_dump
  use empty_data_mod,                 only : empty_real_data,                  &
                                             empty_integer_data
  use files_config_mod,               only : checkpoint_stem_name
  use lfric_string_mod,               only : split_string

#ifdef UM_PHYSICS
  use multidata_field_dimensions_mod, only : &
    get_ndata_val => get_multidata_field_dimension
#endif

  implicit none

  private
  public :: field_maker_type

  !> @brief Processor type for creating fields from specifiers
  type, extends(processor_type) :: field_maker_type
    type(mesh_type), pointer :: mesh
    type(mesh_type), pointer :: twod_mesh
    type(mesh_type), pointer :: coarse_mesh
    type(mesh_type), pointer :: coarse_twod_mesh
    type(field_mapper_type), pointer :: mapper
  contains
    private
    ! main interface
    procedure, public :: init => field_maker_init
    procedure, public :: apply => field_maker_apply

    ! destructor - here to avoid gnu compiler bug
    final :: field_maker_destructor
  end type field_maker_type

contains

!> Return true if and only if the given space supports XIOS IO
function has_xios_io(space, legacy) result(flag)
  implicit none

  type(function_space_type), pointer, intent(in) :: space
  logical(l_def), intent(in) :: legacy

  logical(l_def) :: flag

  if (associated(space)) then
    flag = space_has_xios_io(space%which(), legacy)
  else
    ! dynamic discovery requires XIOS
    flag = .true.
  end if

end function has_xios_io

  !> @brief Initialise field maker object
  !> @param[inout] self             Field maker object
  !> @param[in]    mesh             Mesh for spatial fields
  !> @param[in]    twod_mesh        Mesh for planar fields
  !> @param[in]    coarse_mesh      Coarse mesh for spatial fields
  !> @param[in]    coarse_twod_mesh Coarse mesh for planar fields
  !> @param[in]    mapper           Provides access to field collections
  !> @param[in]    clock            Model clock
  subroutine field_maker_init(self, mesh, twod_mesh, coarse_mesh, coarse_twod_mesh, &
                              mapper, clock)
    implicit none
    class(field_maker_type), intent(inout) :: self

    type(mesh_type), intent(in), pointer :: mesh
    type(mesh_type), intent(in), pointer :: twod_mesh
    type(mesh_type), intent(in), pointer :: coarse_mesh
    type(mesh_type), intent(in), pointer :: coarse_twod_mesh
    type(field_mapper_type), target, intent(in) :: mapper
    class(clock_type), intent(in) :: clock

    self%mesh => mesh
    self%twod_mesh => twod_mesh
    self%coarse_mesh => coarse_mesh
    self%coarse_twod_mesh => coarse_twod_mesh
    self%mapper => mapper
    call self%set_clock(clock)
  end subroutine field_maker_init

  !> @brief Destructor for field maker objects
  !> @param[inout] self       Field maker object
  subroutine field_maker_destructor(self)
    type(field_maker_type), intent(inout) :: self
    ! empty
  end subroutine field_maker_destructor

  !> @brief Apply a field maker object to a specifier, creating a field.
  !> @param[in] self       Field maker object
  !> @param[in] spec       Field specifier
  subroutine field_maker_apply(self, spec)
    implicit none
    class(field_maker_type), intent(in) :: self
    type(field_spec_type), intent(in) :: spec

    type(field_collection_type), pointer :: main_coll => null()
    type(field_collection_type), pointer :: adv_coll => null()
    type(field_collection_type), pointer :: depository => null()
    type(field_collection_type), pointer :: prognostic_fields => null()
    type(function_space_type), pointer :: space => null()
    type(field_type), pointer :: external_real_field => null()
    type(time_axis_type), pointer :: time_axis => null()
    type(integer_field_type), pointer :: external_int_field => null()
    logical(l_def) :: advected
    integer(i_def) :: ndata
    integer(i_def) :: i
    character(:), allocatable :: split_stem_name(:)
    character(str_def) :: field_prefix
    class(clock_type), pointer :: clock
    real(r_second), allocatable :: all_checkpoint_times(:)

#ifdef UM_PHYSICS
    ndata = get_ndata_val(spec%mult)
#else
    ndata = 1
#endif
    depository => self%mapper%get_depository()

    main_coll => self%mapper%get_main_coll_ptr(spec%main_coll)
    adv_coll => self%mapper%get_adv_coll_ptr(spec%adv_coll)

    advected = associated(adv_coll)
    if (.not. advected) adv_coll => depository ! arbirary, any collection will do

    if (spec%space == missing_fs) then
      space => null()                         ! to be inferred from metadata
      if (.not. use_xios_io) then
        call log_event('metadata-based space inference for field ' &
          // trim(spec%name) // ' requires use_xios_io to be true', log_level_error)
      end if
    else
      if (spec%coarse) then
        if (spec%twod) then
          space => function_space_collection%get_fs(self%coarse_twod_mesh, spec%order_h, &
                                                    spec%order_v, spec%space, ndata)
        else
          space => function_space_collection%get_fs(self%coarse_mesh, spec%order_h, &
                                                    spec%order_v, spec%space, ndata)
        end if
      else
        if (spec%twod) then
          space => function_space_collection%get_fs(self%twod_mesh, spec%order_h, &
                                                    spec%order_v, spec%space, ndata)
        else
          space => function_space_collection%get_fs(self%mesh, spec%order_h, &
                                                    spec%order_v, spec%space, ndata)
        end if
      end if
    end if

    ! Check whether checkpoint fields have been added to the XIOS context.
    if (use_xios_io .and. spec%ckp .and. space_has_xios_io(spec%space)) then
      if (checkpoint_write) then
        split_stem_name = split_string( &
          trim(checkpoint_stem_name), '/' )
        clock => self%get_clock()
        select type(clock)
        type is (model_clock_type)
          call create_checkpoint_list( clock, checkpoint_times, &
                                       end_of_run_checkpoint, all_checkpoint_times )
          do i = 1, size(all_checkpoint_times)
            write(field_prefix,'(A,A,I10.10,A)') &
                  trim(split_stem_name(size(split_stem_name))),"_", &
                  clock%steps_from_seconds(all_checkpoint_times(i)), "_"
            if (.not. field_is_valid(trim(field_prefix) // trim(spec%name))) then
              call log_event('Checkpoint field not enabled for ' &
                // trim(spec%name), log_level_error)
            end if
          end do
        class default
          call log_event('Unsupported clock type for checkpoint time handling', &
            log_level_error)
        end select
      end if
      if (checkpoint_read .or. init_option == init_option_checkpoint_dump) then
          if (.not. field_is_valid('restart_' // trim(spec%name))) then
            call log_event('restart field not enabled for ' &
              // trim(spec%name), log_level_error)
          end if
      end if
    end if

    prognostic_fields => self%mapper%get_prognostic_fields()
    if (spec%is_int) then
      if (spec%moist_arr /= moist_arr_dict%none) &
        call log_event('integer field ' // trim(spec%name) &
          // ' cannot be a moisture field', &
          log_level_error)
      external_int_field => null()
      if (spec%time_axis /= time_axis_dict%none) then
        call log_event('integer field ' // trim(spec%name) &
          // ' cannot have a time axis', &
          log_level_error)
      end if
      time_axis => null()
      call add_integer_field( &
        main_coll, &
        depository, &
        prognostic_fields, &
        adv_coll, &
        spec%name, &
        space, &
        spec%order_h, &
        spec%order_v, &
        spec%empty, &
        external_int_field, &
        time_axis, &
        spec%legacy, &
        spec%ckp, &
        advected)
    else
      if (spec%moist_arr /= moist_arr_dict%none .and. spec%moist_idx == 0) then
        call log_event('moisture field ' // trim(spec%name) &
          // ' needs a non-zero array index', &
          log_level_error)
      end if
      external_real_field => &
        self%mapper%get_moist_field_ptr(spec%moist_arr, spec%moist_idx)
      time_axis => &
        self%mapper%get_time_axis_ptr(spec%time_axis)
      call add_real_field( &
        main_coll, &
        depository, &
        prognostic_fields, &
        adv_coll, &
        spec%name, &
        space, &
        spec%order_h, &
        spec%order_v, &
        spec%empty, &
        external_real_field, &
        time_axis, &
        spec%legacy, &
        spec%ckp, &
        advected)
    end if
    nullify(main_coll, adv_coll, depository, prognostic_fields, &
      space, external_real_field, external_int_field, time_axis)
  end subroutine field_maker_apply

  !>@brief Add field to field collection and set its write,
  !>       checkpoint-restart and advection behaviour
  !> @param[in,out] field_collection  Field collection that 'name' will be added to
  !> @param[in,out] depository        Collection of all fields
  !> @param[in,out] prognostic_fields Collection of checkpointed fields
  !> @param[in,out] advected_fields   Collection of fields to be advected
  !> @param[in]     name              Name of field to be added to collection
  !> @param[in]     vector_space      Function space of field to set behaviour for
  !> @param[in]     order_h           Horizontal function space order (for space discovery)
  !> @param[in]     order_v           Vertical function space order (for space discovery)
  !> @param[in]     empty             Flag whether this field is empty
  !> @param[in]     external_field    Pointer to external field or null
  !> @param[in]     time_axis         Pointer to time axis or null
  !> @param[in]     legacy            Flag whether this field uses legacy IO
  !> @param[in]     checkpoint_flag   Optional flag to allow checkpoint-
  !>                                   restart behaviour of field to be set
  !> @param[in]     advection_flag    Optional flag whether this field is to be advected
   subroutine add_real_field(field_collection, &
                              depository, prognostic_fields, advected_fields, &
                              name, vector_space, order_h, order_v, empty, external_field, &
                              time_axis, legacy, checkpoint_flag, advection_flag)

    use io_config_mod,           only : use_xios_io, &
                                        write_diag, checkpoint_write, &
                                        checkpoint_read
    use initialization_config_mod, only: init_option,               &
                                         init_option_checkpoint_dump
    use lfric_xios_read_mod,     only : read_field_generic
    use lfric_xios_write_mod,    only : write_field_generic
    use io_mod,                  only : checkpoint_write_netcdf, &
                                        checkpoint_read_netcdf

    implicit none

    character(*), intent(in)                       :: name
    type(field_collection_type), pointer, &
                                 intent(inout)     :: field_collection
    type(field_collection_type), intent(inout)     :: depository
    type(field_collection_type), intent(inout)     :: prognostic_fields
    type(field_collection_type), intent(inout)     :: advected_fields
    type(function_space_type), pointer, intent(in) :: vector_space
    integer(i_def), intent(in)                     :: order_h
    integer(i_def), intent(in)                     :: order_v
    logical(l_def), intent(in)                     :: empty
    type(field_type), pointer, intent(in)          :: external_field
    type(time_axis_type), pointer, intent(in)      :: time_axis
    logical(l_def), intent(in)                     :: legacy
    logical(l_def), optional, intent(in)           :: checkpoint_flag
    logical(l_def), optional, intent(in)           :: advection_flag
    !Local variables
    type(field_type), target                       :: new_field
    type(field_type), pointer                      :: new_field_ptr => null()
    type(field_type), pointer                      :: field_ptr => null()
    class(pure_abstract_field_type), pointer       :: tmp_ptr => null()
    type(function_space_type), pointer             :: window_size_space
    logical(l_def)                                 :: checkpointed
    logical(l_def)                                 :: advected

    ! pointers for xios write interface
    procedure(write_interface), pointer :: write_behaviour => null()
    procedure(read_interface),  pointer :: read_behaviour => null()
    procedure(checkpoint_write_interface), pointer :: checkpoint_write_behaviour => null()
    procedure(checkpoint_read_interface), pointer  :: checkpoint_read_behaviour => null()

    ! Create the new field
    if (associated(external_field)) then
      new_field_ptr => external_field
    else
      new_field_ptr => new_field
    end if

    ! deal with time axis
    if (associated(time_axis)) then
      ! pre-initialise field with window size, add to time axis;
      ! cf. init_time_axis_mod / setup_field
      if (associated(vector_space)) then
        ! re-create function space, overriding ndata with window size
        window_size_space => &
          function_space_collection%get_fs(     &
            vector_space%get_mesh(),            &
            vector_space%get_element_order_h(), &
            vector_space%get_element_order_v(), &
            vector_space%which(),               &
            time_axis%get_window_size())
        if (empty) then
          call new_field_ptr%initialise( window_size_space, name=trim(name), &
            override_data = empty_real_data )
        else
          call new_field_ptr%initialise( window_size_space, name=trim(name) )
        end if
      else
        ! discover space, overriding ndata with window size
        call init_field_from_metadata( &
           new_field_ptr, trim(name), force_order_h=order_h, force_order_v=order_v,&
            force_ndata = time_axis%get_window_size(), empty=empty )
      end if
      call time_axis%add_field(new_field_ptr)
    end if

    ! regular field initialisation
    if (associated(vector_space)) then
      if (empty) then
        call new_field_ptr%initialise( vector_space, name=trim(name), &
          override_data = empty_real_data )
      else
        call new_field_ptr%initialise( vector_space, name=trim(name) )
      end if
    else
      call init_field_from_metadata( &
         new_field_ptr, trim(name), force_order_h=order_h, &
         force_order_v=order_v, empty=empty )
    end if

    ! Set advection flag
    if (present(advection_flag)) then
      advected = advection_flag
    else
      advected = .false.
    end if

    ! Set checkpoint flag
    if (present(checkpoint_flag)) then
      checkpointed = checkpoint_flag
    else
      checkpointed = .false.
    end if

    ! Set read and write behaviour
    if (use_xios_io) then
        write_behaviour => write_field_generic
        read_behaviour  => read_field_generic
        if (has_xios_io(vector_space, legacy) &
            .and. (write_diag .or. (checkpoint_write .and. checkpointed))) &
          call new_field_ptr%set_write_behaviour(write_behaviour)
        if (has_xios_io(vector_space, legacy) &
            .and. (checkpoint_read .or. init_option == init_option_checkpoint_dump) &
            .and. checkpointed) &
          call new_field_ptr%set_read_behaviour(read_behaviour)
    else
        checkpoint_write_behaviour => checkpoint_write_netcdf
        checkpoint_read_behaviour  => checkpoint_read_netcdf
        call new_field_ptr%set_checkpoint_write_behaviour(checkpoint_write_behaviour)
        call new_field_ptr%set_checkpoint_read_behaviour(checkpoint_read_behaviour)
    endif

    ! Add the field to the depository, unless it is external
    if (associated(external_field)) then
      tmp_ptr => external_field
    else
      call depository%add_field(new_field)
      call depository%get_field(name, field_ptr)
      tmp_ptr => field_ptr
    end if
    ! Put a pointer to the field in the required collection
    if (associated(field_collection)) &
      call field_collection%add_reference_to_field( tmp_ptr )
    ! If checkpointing the field, put a pointer to it in the prognostics collection
    if ( checkpointed ) then
      call prognostic_fields%add_reference_to_field( tmp_ptr )
    endif
    ! If advecting the field, put a pointer to it in the advected collection
    if ( advected ) then
      call advected_fields%add_reference_to_field( tmp_ptr )
    endif

  end subroutine add_real_field

  !>@brief Add integer field to field collection and set its write,
  !>       checkpoint-restart and advection behaviour
  !> @param[in,out] field_collection  Field collection that 'name' will be added to
  !> @param[in,out] depository        Collection of all fields
  !> @param[in,out] prognostic_fields Collection of checkpointed fields
  !> @param[in,out] advected_fields   Collection of fields to be advected
  !> @param[in]     name              Name of field to be added to collection
  !> @param[in]     vector_space      Function space of field to set behaviour for
  !> @param[in]     order_h           Horizontal function space order (for space discovery)
  !> @param[in]     order_v           Vertitcal function space order (for space discovery)
  !> @param[in]     empty             Flag whether this field is empty
  !> @param[in]     external_field    Pointer to external field or null
  !> @param[in]     time_axis         Pointer to time axis or null
  !> @param[in]     legacy            Flag whether this field uses legacy checkpointing
  !> @param[in]     checkpoint_flag   Optional flag to allow checkpoint-
  !>                                   restart behaviour of field to be set
  !> @param[in]     advection_flag    Optional flag whether this field is to be advected
  subroutine add_integer_field(field_collection, &
                              depository, prognostic_fields, advected_fields, &
                              name, vector_space, order_h, order_v, empty, &
                              external_field, time_axis, legacy, &
                              checkpoint_flag, advection_flag)

    use io_config_mod,           only : use_xios_io, &
                                        write_diag, checkpoint_write, &
                                        checkpoint_read
    use initialization_config_mod, only: init_option,               &
                                         init_option_checkpoint_dump
    use lfric_xios_read_mod,     only : read_field_generic
    use lfric_xios_write_mod,    only : write_field_generic
    use io_mod,                  only : checkpoint_write_netcdf, &
                                        checkpoint_read_netcdf

    implicit none

    character(*), intent(in)                       :: name
    type(field_collection_type), pointer, &
                                 intent(inout)     :: field_collection
    type(field_collection_type), intent(inout)     :: depository
    type(field_collection_type), intent(inout)     :: prognostic_fields
    type(field_collection_type), intent(inout)     :: advected_fields
    type(function_space_type), pointer, intent(in) :: vector_space
    integer(i_def), intent(in)                     :: order_h
    integer(i_def), intent(in)                     :: order_v
    logical(l_def), intent(in)                     :: empty
    type(integer_field_type), pointer, intent(in)  :: external_field
    type(time_axis_type), pointer, intent(in)      :: time_axis
    logical(l_def), intent(in)                     :: legacy
    logical(l_def), optional, intent(in)           :: checkpoint_flag
    logical(l_def), optional, intent(in)           :: advection_flag
    !Local variables
    type(integer_field_type), target               :: new_field
    type(integer_field_type), pointer              :: new_field_ptr => null()
    type(integer_field_type), pointer              :: field_ptr => null()
    class(pure_abstract_field_type), pointer       :: tmp_ptr => null()
    logical(l_def)                                 :: checkpointed
    logical(l_def)                                 :: advected

    ! pointers for xios write interface
    procedure(write_interface), pointer :: write_behaviour => null()
    procedure(read_interface),  pointer :: read_behaviour => null()
    procedure(checkpoint_write_interface), pointer :: checkpoint_write_behaviour => null()
    procedure(checkpoint_read_interface), pointer  :: checkpoint_read_behaviour => null()

    ! Create the new field
    if (associated(external_field)) then
      new_field_ptr => external_field
    else
      new_field_ptr => new_field
    end if
    if (associated(time_axis)) then
      call log_event('impossible error: integer field ' // trim(name) &
        // ' with non-null time axis', &
      log_level_error)
    end if
    if (associated(vector_space)) then
      if (empty) then
        call new_field_ptr%initialise( vector_space, name=trim(name), &
          override_data = empty_integer_data )
      else
        call new_field_ptr%initialise( vector_space, name=trim(name) )
      end if
    else
      call init_field_from_metadata( &
        new_field_ptr, trim(name), force_order_h=order_h, &
        force_order_v=order_v, empty=empty )
    end if

    ! Set advection flag
    if (present(advection_flag)) then
      advected = advection_flag
    else
      advected = .false.
    end if

    ! Set checkpoint flag
    if (present(checkpoint_flag)) then
      checkpointed = checkpoint_flag
    else
      checkpointed = .false.
    end if

    ! Set read and write behaviour
    if (use_xios_io) then
      write_behaviour => write_field_generic
      read_behaviour  => read_field_generic
      if (has_xios_io(vector_space, legacy) &
          .and. (write_diag .or. (checkpoint_write .and. checkpointed))) &
        call new_field_ptr%set_write_behaviour(write_behaviour)
      if (has_xios_io(vector_space, legacy) &
           .and. (checkpoint_read .or. init_option == init_option_checkpoint_dump) &
           .and. checkpointed) &
        call new_field_ptr%set_read_behaviour(read_behaviour)
    else
      checkpoint_write_behaviour => checkpoint_write_netcdf
      checkpoint_read_behaviour  => checkpoint_read_netcdf
      call new_field_ptr%set_checkpoint_write_behaviour(checkpoint_write_behaviour)
      call new_field_ptr%set_checkpoint_read_behaviour(checkpoint_read_behaviour)
    endif

    ! Add the field to the depository, unless it is external
    if (associated(external_field)) then
      tmp_ptr => external_field
    else
      call depository%add_field(new_field)
      call depository%get_field(name, field_ptr)
      tmp_ptr => field_ptr
    end if

    ! Put a pointer to the field in the required collection
    if (associated(field_collection)) &
      call field_collection%add_reference_to_field( tmp_ptr )
    ! If checkpointing the field, put a pointer to it in the prognostics collection
    if ( checkpointed ) then
      call prognostic_fields%add_reference_to_field( tmp_ptr )
    endif
    ! If advecting the field, put a pointer to it in the advected collection
    if ( advected ) then
      call advected_fields%add_reference_to_field( tmp_ptr )
    endif

  end subroutine add_integer_field

  end module field_maker_mod
