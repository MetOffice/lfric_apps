!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Initialisation functionality for the oasis coupling

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

  public :: lfric2lfric_init_coupler_src, &
            lfric2lfric_init_coupler_dst, &
            lfric2lfric_end_coupler_init, &
            read_nfields,                 &
            parse_namcouple

  contains

#ifdef MCT
  !> @details Initialises coupling partitions and variables for the
  !>          source fields
  !> @param[in]     name_cpl_list   List of field names to be coupled
  !> @param[in]     coupler         Name of the coupling component
  !> @param[in]     element_order_h The polynomial order of the set of
  !>                                compatible finite element in the
  !>                                horizontal
  !> @param[in]     element_order_v The polynomial order of the set of
  !>                                compatible finite element in the
  !>                                vertical
  !> @param[in]     mesh Representation of the mesh the code will run on
  !> @param[in,out] modeldb The structure that holds model state
  subroutine lfric2lfric_init_coupler_src(name_cpl_list, coupler, element_order_h, &
                                          element_order_v, mesh, modeldb)

    implicit none

    character(len=str_def), allocatable, intent(in)    :: name_cpl_list(:)
    character(*),                        intent(in)    :: coupler
    integer(i_def),                      intent(in)    :: element_order_h
    integer(i_def),                      intent(in)    :: element_order_v
    type(mesh_type), pointer,            intent(in)    :: mesh
    type(modeldb_type),                  intent(inout) :: modeldb
    
    ! Coupling fields
    type(field_collection_type), pointer :: cpl_snd_2d => null()
    type(field_collection_type), pointer :: cpl_rcv_2d => null()
    type(field_collection_type), pointer :: cpl_snd_0d => null()

    type(field_type)                     :: field_cpl
    type(field_collection_type), pointer :: depository
    type(field_type), pointer            :: field
    class(pure_abstract_field_type), pointer :: field_ptr

    type(coupling_type), pointer         :: coupling_ptr
    character(str_def)                   :: name_cpl
    integer :: i

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

    depository => modeldb%fields%get_field_collection("depository")

    ! Loop over names in the coupling list (as defined by the namcouple)
    do i =1, size(name_cpl_list)

      name_cpl = name_cpl_list(i)
       
      ! Create prognostic fields
      ! Creates a field in the W3 function space (fully discontinuous field)
      call field_cpl%initialise( vector_space = &
                    function_space_collection%get_fs(mesh,             &
                                element_order_h, element_order_v, W3), &
                                name=name_cpl)

      ! Add fields to modeldb 
      call depository%add_field(field_cpl)

      ! Set lfric fields as either the source entry
      call depository%get_field(TRIM(name_cpl), field)

      ! Add that field to the coupling send field collection for 2d fields
      field_ptr => field
      call cpl_snd_2d%add_reference_to_field(field_ptr)

    end do

    ! Extract the coupling object from the modeldb key-value pair collection
    coupling_ptr => get_coupling_from_collection(modeldb%values, trim(coupler))
    call coupling_ptr%define_partitions(modeldb%mpi,mesh)
    call coupling_ptr%define_variables( cpl_snd_2d, &
                                        cpl_rcv_2d, &
                                        cpl_snd_0d )
  end subroutine lfric2lfric_init_coupler_src

  !> @details Initialises coupling partitions and variables for the
  !>          destination fields
  !> @param[in]     name_cpl_list   List of field names to be coupled
  !> @param[in]     coupler         Name of the coupling component
  !> @param[in]     element_order_h The polynomial order of the set of
  !>                                compatible finite element in the
  !>                                horizontal
  !> @param[in]     element_order_v The polynomial order of the set of
  !>                                compatible finite element in the
  !>                                vertical
  !> @param[in]     mesh Representation of the mesh the code will run on
  !> @param[in,out] modeldb The structure that holds model state 
  subroutine lfric2lfric_init_coupler_dst(name_cpl_list, coupler, element_order_h, &
                                          element_order_v, mesh, modeldb)

    implicit none

    character(len=str_def), allocatable, intent(in)    :: name_cpl_list(:)
    character(*),                        intent(in)    :: coupler
    integer(i_def),                      intent(in)    :: element_order_h
    integer(i_def),                      intent(in)    :: element_order_v
    type(mesh_type), pointer,            intent(in)    :: mesh
    type(modeldb_type),                  intent(inout) :: modeldb

    ! Coupling fields
    type(field_collection_type), pointer :: cpl_snd_2d => null()
    type(field_collection_type), pointer :: cpl_rcv_2d => null()
    type(field_collection_type), pointer :: cpl_snd_0d => null()

    type(field_type)                     :: field_cpl
    type(field_collection_type), pointer :: depository
    type(field_type), pointer            :: field
    class(pure_abstract_field_type), pointer :: field_ptr

    type(coupling_type), pointer         :: coupling_ptr
    character(str_def)                   :: name_cpl
    integer :: i

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

    depository => modeldb%fields%get_field_collection("depository")

    do i =1, size(name_cpl_list)

      name_cpl = name_cpl_list(i)
    
      ! Create prognostic fields
      ! Creates a field in the W3 function space (fully discontinuous field)

      call field_cpl%initialise( vector_space = &
                     function_space_collection%get_fs(mesh,            &
                                element_order_h, element_order_v, W3), &
                                name=name_cpl)

      ! Add fields to modeldb
      call depository%add_field(field_cpl)

      ! Set lfric fields as either the destination entry
      call depository%get_field(trim(name_cpl), field)
      ! Add that field to the coupling send field collection for 2d fields
      field_ptr => field
      call cpl_rcv_2d%add_reference_to_field(field_ptr)

    end do

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

  !> @brief   Read the total number of fields from the namcouple
  !> @details Read the line below $NFIELDS. This is the number of fields
  !!          and is equivalent to the number of entries.
  !> @param[inout] total_fields total number of fields to be coupled
  subroutine read_nfields(total_fields)
    
    implicit none

    integer :: unit, ios, total_fields
    character(len=str_def) :: line

    ! Open the namcouple file (that already exists)
    open(newunit=unit, file='namcouple', status='old', iostat=ios)

    ! Error message if unable to open
    if (ios /= 0) then
      write( log_scratch_space, '(A,I3)' ) &
        'lfric2lfric_parse_namcouple: Could not open namcouple file, iostat=', ios
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    ! Find the number of entries from the namcouple
    ! This is the line below $NFIELDS
    do
      read(unit,'(A)',iostat=ios) line

      ! Error message if end of file found before $NFIELDS
      if (ios /= 0) then
        write( log_scratch_space, '(A)' ) &
        'lfric2lfric_parse_namcouple: $NFIELDS not found in namcouple file'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      ! Read nfields
      if (trim(adjustl(line)) == '$NFIELDS') then
        do 
          read(unit, *, iostat=ios) line
          if (line(1:1) /= '#') then
            backspace(unit)
            read(unit, *, iostat=ios) total_fields
            exit
          end if
        end do

        ! Error message if unable to read
        if (ios /= 0) then
          write( log_scratch_space, '(A)' ) &
          'lfric2lfric_parse_namcouple: not able to read $NFIELDS from namcouple'
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
       
        exit
      endif

    enddo

    close(unit)

  end subroutine read_nfields

  !> @brief   Populate src_name_cpl_list and dst_name_cpl_list
  !> @details Create the lists of fields to be initialised by the oasis
  !!          coupler, by reading from the namcouple file.
  !> @param[inout] src_name_cpl_list List of source fields to be initialised
  !> @param[inout] dst_name_cpl_list List of destination fields to be
  !!               initialised
  subroutine parse_namcouple(src_name_cpl_list, dst_name_cpl_list)

    implicit none

    character(len=str_def), allocatable, intent(inout) :: src_name_cpl_list(:)
    character(len=str_def), allocatable, intent(inout) :: dst_name_cpl_list(:)

    type :: coupling_field
      ! Coupling field matches the format of the entry in the namcouple
      ! e.g. {src} {dst} 1 {freq} 1 lfric2lfric_restart_oa.nc EXPORTED
      character(len=str_def) :: src_name
      character(len=str_def) :: dst_name
      integer                :: int1
      integer                :: freq
      integer                :: int2
      character(str_def)     :: file
      character(str_def)     :: exported
    end type coupling_field

    type(coupling_field), allocatable :: fields(:)
    integer :: unit, ios, total_fields, nfields
    character(len=str_def) :: line

    ! Get the total number of fields to couple
    call read_nfields(total_fields)
    
    ! Allocate the array sizes to be total_fields
    allocate(fields(total_fields))
    allocate(src_name_cpl_list(total_fields))
    allocate(dst_name_cpl_list(total_fields))

    ! Open the namcouple file (that already exists)
    open(newunit=unit, file='namcouple', status='old', iostat=ios)

    ! Error message if unable to open
    if (ios /= 0) then
      write( log_scratch_space, '(A,I3)' ) &
        'lfric2lfric_parse_namcouple: Could not open namcouple file, iostat=', ios
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    
    ! Read the entries in the namcouple (denoted by 'EXPORTED') and
    ! populate src_name_cpl_list and dst_name_cpl_list
    nfields = 0
    do
      read(unit,'(A)',iostat=ios) line
      if (ios /= 0) exit

      if (index(line,'EXPORTED') > 0) then
        nfields = nfields + 1

        ! Error message if the number of entries read is greater than the total_fields
        if (nfields > total_fields) then
           write( log_scratch_space, '(A,I3,A,I3)' ) &
           'lfric2lfric_parse_namcouple: actual number of entries read from namcouple', &
           nfields, 'greater than the expected number of entries ', &
           total_fields
           call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if

        ! Read the entry into the coupling_field type
        read(line, *, iostat=ios) fields(nfields)%src_name, &
                                  fields(nfields)%dst_name, &
                                  fields(nfields)%int1,     &
                                  fields(nfields)%freq,     &
                                  fields(nfields)%int2,     &
                                  fields(nfields)%file,     &
                                  fields(nfields)%exported
        ! Error message if unable to read entry
        if (ios /= 0) then
          write( log_scratch_space, '(A, I3)' ) &
          'lfric2lfric_parse_namcouple: Not able to read namcouple entry number', nfields
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if

        ! populate the lists with the src_name and dst_name
        src_name_cpl_list(nfields) = trim(fields(nfields)%src_name)
        dst_name_cpl_list(nfields) = trim(fields(nfields)%dst_name)
      endif
    enddo

    close(unit)

    ! Error message if the number of entries read is different to total_fields
    if (nfields /= total_fields) then
      write( log_scratch_space, '(A,I3,A,I3)' ) &
           'lfric2lfric_parse_namcouple: actual number of entries read from namcouple', &
           nfields, 'not equal to the expected number of entries ', &
           total_fields
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
  
  end subroutine parse_namcouple

#else
  subroutine lfric2lfric_init_coupler_src(name_cpl_list, coupler, element_order_h, &
                                          element_order_v, mesh, modeldb)

    implicit none

    character(*),             intent(in)     :: coupler
    integer(i_def),           intent(in)     :: element_order_h
    integer(i_def),           intent(in)     :: element_order_v
    type(mesh_type), pointer, intent(in)     :: mesh
    type(modeldb_type),       intent(inout)  :: modeldb
    character(len=str_def), allocatable, intent(inout) :: name_cpl_list(:)

    write(log_scratch_space, '(A)' ) &
        "lfric2lfric_init_coupler_src: to use OASIS, " // &
        "cpp directive MCT must be set"
    call log_event( log_scratch_space, log_level_error )

  end subroutine lfric2lfric_init_coupler_src

  subroutine lfric2lfric_init_coupler_dst(name_cpl_list, coupler, element_order_h, &
                                          element_order_v, mesh, modeldb)

    implicit none

    character(*),             intent(in)     :: coupler
    integer(i_def),           intent(in)     :: element_order_h
    integer(i_def),           intent(in)     :: element_order_v
    type(modeldb_type),       intent(inout)  :: modeldb
    type(mesh_type), pointer, intent(in)     :: mesh
    character(len=str_def), allocatable, intent(inout) :: name_cpl_list(:)
    
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

  subroutine parse_namcouple(src_name_cpl_list, dst_name_cpl_list)
    
    implicit none

    character(len=str_def), allocatable, intent(inout) :: src_name_cpl_list(:)
    character(len=str_def), allocatable, intent(inout) :: dst_name_cpl_list(:)

    write(log_scratch_space, '(A)' ) &
        "lfric2lfric_parse_namcouple: to use OASIS, " // &
        "cpp directive MCT must be set"
    call log_event( log_scratch_space, log_level_error )
    
  end subroutine parse_namcouple

#endif

end module lfric2lfric_init_coupler_mod
