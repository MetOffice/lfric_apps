!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> &brief coupling related routines for use in coupled configuration

module coupler_mod

  use accumulate_send_fields_2d_mod,  only: accumulate_send_fields_2d
  use coupling_mod,                   only: coupling_type,       &
                                            get_coupling_fields, &
                                            get_coupling_from_collection
  use driver_water_constants_mod,     only: T_freeze_h2o_sea
  use jules_sea_seaice_config_mod,    only: therm_cond_sice => kappai,         &
                                            therm_cond_sice_snow => kappai_snow
  use field_mod,                      only: field_type
  use field_parent_mod,               only: field_parent_type
  use pure_abstract_field_mod,        only: pure_abstract_field_type
  use function_space_mod,             only: function_space_type
  use fs_continuity_mod,              only: W3, Wtheta
  use function_space_collection_mod,  only: function_space_collection
  use field_collection_iterator_mod,  only: field_collection_iterator_type
  use field_collection_mod,           only: field_collection_type
  use driver_modeldb_mod,             only: modeldb_type
  use constants_mod,                  only: i_def, r_def, str_def, l_def, &
                                            i_halo_index, imdi, rmdi
  use timestepping_config_mod,        only: dt
  use lfric_mpi_mod,                  only: global_mpi
  use log_mod,                        only: log_event,       &
                                            LOG_LEVEL_INFO,  &
                                            LOG_LEVEL_DEBUG, &
                                            LOG_LEVEL_ERROR, &
                                            log_scratch_space
  use mesh_mod,                       only: mesh_type
  use field_parent_mod,               only: write_interface, read_interface,  &
                                            checkpoint_write_interface,       &
                                            checkpoint_read_interface
  use process_send_fields_0d_mod,     only: process_send_fields_0d
  use process_send_fields_2d_mod,     only: process_send_fields_2d,           &
                                            cpl_reset_field
  use coupler_exchange_2d_mod,        only: coupler_exchange_2d_type
  use coupler_exchange_0d_mod,        only: coupler_send_0d
  use coupler_update_prognostics_mod, only: coupler_update_prognostics
  use coupling_mod,                   only: cpl_cat
  use process_recv_fields_2d_mod,     only: process_recv_fields_2d
  use esm_couple_config_mod,          only: l_esm_couple_test

#if defined(UM_PHYSICS)
  use jules_control_init_mod,         only: n_sea_tile, first_sea_tile
  ! Note: n_sea_ice_tile has to be retrieved from jules_sea_seaice_config_mod
  !       and not jules_control_init_mod as the coupler is initialised before
  !       jules
  use jules_sea_seaice_config_mod,    only: n_sea_ice_tile => nice
#endif

  implicit none

#if !defined(UM_PHYSICS)
  !
  ! Dummy variables required when NOT running with UM_PHYSICS
  !
  integer(i_def),parameter              :: n_sea_tile = imdi
  integer(i_def),parameter              :: first_sea_tile = imdi
  integer(i_def),parameter              :: n_sea_ice_tile = imdi
#endif

  private

  ! Prefix of fields that will be coupled to/from LFRic
  character(len=3), parameter  :: cpl_prefix = "lf_"
  ! Name of the first level of a multi-category field
  character(len=2), parameter  :: cpl_fixed_catno = "01"

  !routines
  public zero_coupling_send_fields, create_coupling_fields, &
         cpl_snd, cpl_fld_update, cpl_rcv,                  &
         generate_coupling_field_collections, set_cpl_name

  contains

  !>@brief Resets the "send" coupling fields to 0 - both those in the send
  !>       collection and their counterparts in the depository
  !>@param[inout] modeldb The working data set for a model run
  !>
  subroutine zero_coupling_send_fields(modeldb)
    implicit none
    type( modeldb_type ), intent(inout) :: modeldb

    type( field_collection_type ), pointer       :: cpl_rcv_2d
    type( field_collection_type ), pointer       :: depository

    ! Field in the form that comes out of a field collection iterator
    class( field_parent_type ), pointer          :: field_iter
    ! Pointer to an actual field
    type( field_type ), pointer                  :: field_ptr
    ! Iterator over a field collection
    type( field_collection_iterator_type)        :: iter
    ! Number of accumulation steps
    real(r_def), pointer                         :: acc_step

    depository => modeldb%fields%get_field_collection("depository")
    cpl_rcv_2d => modeldb%fields%get_field_collection("cpl_rcv_2d")

    ! Reset the accumulation step counter
    call modeldb%values%get_value( 'accumulation_steps', acc_step )
    acc_step = 0.0_r_def

    call iter%initialise(cpl_rcv_2d)
    do
      if (.not.iter%has_next())exit
      field_iter => iter%next()
      select type(field_iter)
      type is (field_type)
        field_ptr => field_iter
        write(log_scratch_space,'(2A)') &
              "zero_coupling_send_fields: set initial value for ", &
              trim(adjustl(field_ptr%get_name()))
        call log_event(log_scratch_space,LOG_LEVEL_DEBUG)
        call cpl_reset_field(field_ptr, depository)
      class default
        write(log_scratch_space, '(2A)')"Error: zero_coupling_send_fields: ", &
                                trim(field_ptr%get_name())//" is NOT field_type"
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end select
    end do

  end subroutine zero_coupling_send_fields

  !>@brief Adds fields used in coupling to depository and prognosic_fields
  !>       collections
  !>@details These fields are the raw fields that haven't been processed
  !>         for coupling, so are the accumulations etc. that are passed
  !>         from timestep to timestep.
  !>@param [in]    mesh      The full 3d mesh
  !>@param [in]    twod_mesh Mesh on which coupling fields are defined
  !>@param [inout] modeldb   The working data set for a model run
  !>
  subroutine create_coupling_fields( mesh,       &
                                     twod_mesh,  &
                                     modeldb )
    implicit none

    type( mesh_type ),             intent(in), pointer :: mesh
    type( mesh_type ),             intent(in), pointer :: twod_mesh
    type(modeldb_type),            intent(inout)       :: modeldb

    type( field_collection_type ),             pointer :: depository
    type( field_collection_type ),             pointer :: prognostic_fields
    !
    !vactor space for coupling field
    type(function_space_type), pointer :: vector_space
    type(function_space_type), pointer :: sice_space
    type(function_space_type), pointer :: threed_space
    type(function_space_type), pointer :: wtheta_space
    real(r_def) :: acc_step
    !checkpoint flag for coupling field
    logical(l_def)                     :: checkpoint_restart_flag

    write(log_scratch_space, * ) "create_coupling_fields: add coupling fields to repository"
    call log_event( log_scratch_space, LOG_LEVEL_DEBUG )

    depository => modeldb%fields%get_field_collection("depository")
    prognostic_fields => modeldb%fields%get_field_collection("prognostic_fields")

   vector_space=> function_space_collection%get_fs( twod_mesh, 0, 0, W3 )
   sice_space  => function_space_collection%get_fs( twod_mesh, 0, 0, W3,       &
                                                               n_sea_ice_tile )
    threed_space => function_space_collection%get_fs(mesh, 0, 0, W3 )
    wtheta_space => function_space_collection%get_fs(mesh, 0, 0, Wtheta)

    ! Coupling uses some sea-ice fractions from the previous step. These,
    ! therefore, need to be stored in the depository for later use.
    checkpoint_restart_flag = .false.
    ! Sea ice fractions used in the previous timestep (not updated by coupling)
    ! (stored as a fraction of the marine portion of the grid box)
    call add_cpl_field(depository, prognostic_fields, &
         'sea_ice_frac_previous', sice_space, checkpoint_restart_flag)
    call add_cpl_field(depository, prognostic_fields, &
         'sea_frac_previous', vector_space, checkpoint_restart_flag)
    ! Raw sea-ice fractions that are unaltered by post processing.
    ! These are used in scaling data when going from atmosphere to ocean.
    ! r_sea_ice_frac_raw is its reciprocal.
    call add_cpl_field(depository, prognostic_fields, &
         'r_sea_ice_frac_raw', sice_space, checkpoint_restart_flag)

    ! Raw ocean fractions that are unaltered by post processing.
    ! These are used in scaling data when going from atmosphere to ocean.
    ! r_ocean_fraction is its reciprocal.
    call add_cpl_field(depository, prognostic_fields, &
         'r_ocean_fraction', vector_space, checkpoint_restart_flag)

    ! The following fields contain accumulations, but will be written to
    ! the restart file as means (ready for coupling). Create a variable to
    ! hold the number of timesteps that make up the accumulations, so it can
    ! be used to convert accumulations to means
    acc_step = 0.0_r_def
    call modeldb%values%add_key_value('accumulation_steps', acc_step)

    !coupling fields
    !sending-depository

    ! these need to be in restart file as they are the accumulations used to
    ! generate the fields passed to the ocean or river model
    checkpoint_restart_flag = .true.
    call add_cpl_field(depository, prognostic_fields, &
         'lf_taux', vector_space, checkpoint_restart_flag)

    call add_cpl_field(depository, prognostic_fields, &
         'lf_tauy', vector_space, checkpoint_restart_flag)

    call add_cpl_field(depository, prognostic_fields, &
         'lf_w10', vector_space, checkpoint_restart_flag)

    call add_cpl_field(depository, prognostic_fields, &
         'lf_solar', vector_space, checkpoint_restart_flag)

    call add_cpl_field(depository, prognostic_fields, &
         'lf_heatflux', vector_space, checkpoint_restart_flag)

    call add_cpl_field(depository, prognostic_fields, &
         'lf_train', vector_space, checkpoint_restart_flag)

    call add_cpl_field(depository, prognostic_fields, &
         'lf_tsnow', vector_space, checkpoint_restart_flag)

    call add_cpl_field(depository, prognostic_fields, &
         'lf_rsurf', vector_space, checkpoint_restart_flag)

    call add_cpl_field(depository, prognostic_fields, &
         'lf_rsub', vector_space, checkpoint_restart_flag)

    ! The following fields are taken care of elsewhere (theoretically)
    ! but we might need duplicates for coupling restarts.
    call add_cpl_field(depository, prognostic_fields, &
         'lf_evap', vector_space, checkpoint_restart_flag)

    call add_cpl_field(depository, prognostic_fields, &
         'lf_topmelt', sice_space, checkpoint_restart_flag)

    call add_cpl_field(depository, prognostic_fields, &
         'lf_iceheatflux', sice_space, checkpoint_restart_flag)

    call add_cpl_field(depository, prognostic_fields, &
         'lf_sublimation', sice_space, checkpoint_restart_flag)

    call add_cpl_field(depository, prognostic_fields, &
         'lf_iceskint', sice_space, checkpoint_restart_flag)

    call add_cpl_field(depository, prognostic_fields, &
         'lf_pensolar', sice_space, checkpoint_restart_flag)

    ! The following fields don't need to be in checkpoint files as they are
    ! calculated instantaneously from snow depth just before coupling
    call add_cpl_field(depository, prognostic_fields, &
         'lf_antarctic', vector_space, .false.)

    call add_cpl_field(depository, prognostic_fields, &
         'lf_greenland', vector_space, .false.)

    !receiving - depository
    vector_space => function_space_collection%get_fs( twod_mesh, 0, 0, W3, ndata=1 )


    ! These do not need to be in the restart file because they come FROM
    ! the ocean/seaice model!
    checkpoint_restart_flag = .false.

    call add_cpl_field(depository, prognostic_fields, &
         'lf_ocn_sst', vector_space, checkpoint_restart_flag)

    call add_cpl_field(depository, prognostic_fields, &
         'lf_icefrc', sice_space, checkpoint_restart_flag)

    call add_cpl_field(depository, prognostic_fields, &
         'lf_icetck', sice_space, checkpoint_restart_flag)

    call add_cpl_field(depository, prognostic_fields, &
         'lf_icelayert', sice_space, checkpoint_restart_flag)

    call add_cpl_field(depository, prognostic_fields, &
         'lf_conductivity', sice_space, checkpoint_restart_flag)

    call add_cpl_field(depository, prognostic_fields, &
         'lf_snow_depth', sice_space, checkpoint_restart_flag)

    call add_cpl_field(depository, prognostic_fields, &
         'lf_pond_frac', sice_space, checkpoint_restart_flag)

    call add_cpl_field(depository, prognostic_fields, &
         'lf_pond_depth', sice_space, checkpoint_restart_flag)

    call add_cpl_field(depository, prognostic_fields, &
         'lf_sunocean', vector_space, checkpoint_restart_flag)

    call add_cpl_field(depository, prognostic_fields, &
         'lf_svnocean', vector_space, checkpoint_restart_flag)

  end subroutine create_coupling_fields


  !>@brief Top level routine for setting coupling fields and sending data
  !>@param [in,out] modeldb The structure that holds model state
  !>
  subroutine cpl_snd( modeldb )

    implicit none

    type( modeldb_type ), intent(inout)          :: modeldb

    !local variables
    ! Field collection containing all fields
    type( field_collection_type ), pointer       :: depository
    ! Field collection containing 2d fields to be sent to coupler
    type( field_collection_type ), pointer       :: cpl_snd_2d
    ! Field collection containing 0d fields to be senr to coupler
    type( field_collection_type ), pointer       :: cpl_snd_0d
    ! Pointer to the coupling object
    type(coupling_type), pointer                 :: coupling_ptr
    ! Pointer to the index used for ordering coupled fields
    integer(i_def), pointer                      :: local_index(:)
    ! Flag that indicates whether pre-dump preparation is required
    logical(l_def)                               :: ldump_prep
    !pointer to a field (parent)
    class( field_parent_type ), pointer          :: field
    !pointer to a field
    type( field_type ), pointer                  :: field_ptr
    !field pointer to field to be written
    type( field_type ), pointer                  :: dep_fld
    !iterator
    type( field_collection_iterator_type)        :: iter
    !pointer to sea ice fractions
    type( field_type ),         pointer          :: ice_frac_ptr
    ! External field used for sending data to the coupler
    type(coupler_exchange_2d_type)               :: coupler_exchange_2d
    !name of the field
    character(str_def)                           :: sname
    ! Coupled field variable id
    integer(i_def)                               :: var_id
    !Ice sheet mass scalar
    real( r_def )                                :: ice_mass
    ! Returned error code
    integer(i_def)                               :: ierror
    ! Number of accumulation steps
    real(r_def), pointer                         :: acc_step

    depository => modeldb%fields%get_field_collection("depository")
    cpl_snd_2d => modeldb%fields%get_field_collection("cpl_snd_2d")
    cpl_snd_0d => modeldb%fields%get_field_collection("cpl_snd_0d")

    ! Extract the coupling object from the modeldb key-value pair collection
    coupling_ptr => get_coupling_from_collection(modeldb%values, "coupling" )
    local_index => coupling_ptr%get_local_index()

    ldump_prep = .false.

    ! increment accumulation step
    call modeldb%values%get_value( 'accumulation_steps', acc_step )
    acc_step = acc_step + 1.0_r_def

    ! Send fields using 2d coupling
    call iter%initialise(cpl_snd_2d)
    do
      if ( .not. iter%has_next() ) exit
      field => iter%next()
      sname = trim(adjustl(field%get_name()))
      select type(field)
      type is (field_type)
        field_ptr => field
        call accumulate_send_fields_2d(field_ptr, depository, acc_step, &
                                       modeldb%clock, ldump_prep)
        ! Create a coupling external field
        call coupler_exchange_2d%initialise(field_ptr, local_index)
        call coupler_exchange_2d%set_time(modeldb%clock)
        if( coupler_exchange_2d%is_coupling_time() ) then
          ! Process the accumulations to make the fields that will be coupled
          call process_send_fields_2d(field_ptr, &
                                      modeldb)
          ! Call through to coupler_send_2d in coupler_exchange_2d_mod
          call coupler_exchange_2d%copy_from_lfric(ierror)
        else
          ! No coupling at this time-step
          ierror=1
          write(log_scratch_space, '(3A)' ) "cpl_snd: field ", &
                         trim(sname), " NOT exchanged on this timestep"
          call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
        end if
        call coupler_exchange_2d%clear()

        ! If coupling was successful then reset field to 0 ready to start
        ! accumulation for the next exchange
        if (ierror == 0) then
          call cpl_reset_field(field, depository)
        endif

        ! Write out the depository version of the field
        sname = trim(adjustl(field%get_name()))
        call depository%get_field(trim(sname), dep_fld)
        call dep_fld%write_field(trim(sname))

      class default
        write(log_scratch_space, '(2A)' ) "PROBLEM cpl_snd: field ", &
              trim(field%get_name())//" is NOT field_type"
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end select

    end do

    ! Reset the accumulation step count (there is only one step count for all
    ! fields, so this just checks that the last field was successfully coupled)
    if (ierror == 0) then
      acc_step = 0.0_r_def
    end if

    ! Send fields using 0d coupling (scalars)
    call iter%initialise(cpl_snd_0d)
    do
      if ( .not. iter%has_next() ) exit
      field => iter%next()
      sname = field%get_name()
      select type(field)
      type is (field_type)
        field_ptr => field
        var_id = field%get_cpl_id(1)
        call process_send_fields_0d(field_ptr, modeldb, ice_mass)
        call coupler_send_0d(ice_mass, sname, var_id, modeldb%clock )
      class default
        write(log_scratch_space, '(2A)' ) "PROBLEM cpl_snd: field ", &
              trim(field%get_name())//" is NOT field_type"
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end select

    end do

  end subroutine cpl_snd


  !>@brief Top level routine for updating coupling fields
  !>@param [in,out] modeldb The structure that holds model state
  !>
  subroutine cpl_fld_update( modeldb )
    implicit none
    type( modeldb_type ), intent(inout)          :: modeldb

    !local variables
    ! Field collection containing all fields
    type( field_collection_type ), pointer       :: depository
    ! Field collection containing 2d fields to be sent to coupler
    type( field_collection_type ), pointer       :: cpl_snd_2d
    ! Flag that indicates whether pre-dump preparation is required
    logical(l_def)                               :: ldump_prep
    !pointer to a field (parent)
    class( field_parent_type ), pointer          :: field
    !pointer to a field
    type( field_type ), pointer                  :: field_ptr
    !field pointer to field to be written
    type( field_type ), pointer                  :: dep_fld
    !iterator
    type( field_collection_iterator_type)        :: iter
    !name of the field
    character(str_def)                           :: name
    ! Number of accumulation steps
    real(r_def), pointer                         :: acc_step

    depository => modeldb%fields%get_field_collection("depository")
    cpl_snd_2d => modeldb%fields%get_field_collection("cpl_snd_2d")

    ldump_prep = .true.

    call modeldb%values%get_value( 'accumulation_steps', acc_step )
    acc_step = acc_step + 1.0_r_def

    ! We need to loop over each output field and ensure it gets updated
    call iter%initialise(cpl_snd_2d)
    do
      if(.not.iter%has_next())exit
      field => iter%next()
      select type(field)
      type is (field_type)
        field_ptr => field
        call accumulate_send_fields_2d(field_ptr, depository, acc_step, &
                                       modeldb%clock, ldump_prep)
        ! Write out the depository version of the field
        name = trim(adjustl(field%get_name()))
        call depository%get_field(trim(name), dep_fld)
        call dep_fld%write_field(trim(name))
      class default
        write(log_scratch_space, '(2A)' ) "PROBLEM cpl_fld_update: field ", &
                         trim(field%get_name())//" is NOT field_type"
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end select
    end do

    ! All fields have been updated and checkpointed,
    ! so reset the accumulation step count
    acc_step = 0.0_r_def

  end subroutine cpl_fld_update


  !>@brief Top level routine for receiving data
  !>@param [in,out] modeldb The structure that holds model state
  !>
  subroutine cpl_rcv( modeldb )
    implicit none
    type( modeldb_type ), intent(inout)          :: modeldb

    !local variables
    ! Field collection containing all fields
    type( field_collection_type ), pointer       :: depository
     ! Field collection containing 2d fields to be received from the coupler
    type( field_collection_type ), pointer       :: cpl_rcv_2d
    ! Pointer to the coupling object
    type(coupling_type), pointer                 :: coupling_ptr
    ! Pointer to the index used for ordering coupled fields
    integer(i_def), pointer                      :: local_index(:)
    !pointer to a field (parent)
    class( field_parent_type ), pointer          :: field
    !pointer to a field
    type( field_type ), pointer                  :: field_ptr
    ! External field used for receiving data from Oasis
    type(coupler_exchange_2d_type)               :: coupler_exchange_2d
    !iterator
    type( field_collection_iterator_type)        :: iter
    !flag for processing data that has just been exchanged
    ! (set to 1 once data has been successfully passed through the coupler)
    integer(i_def)                               :: ierror

    depository => modeldb%fields%get_field_collection("depository")
    cpl_rcv_2d => modeldb%fields%get_field_collection("cpl_rcv_2d")
    ! Extract the coupling object from the modeldb key-value pair collection
    coupling_ptr => get_coupling_from_collection(modeldb%values, "coupling" )
    local_index => coupling_ptr%get_local_index()

    ! Set defaults
    ierror = 0

    call iter%initialise(cpl_rcv_2d)
    do
      if (.not.iter%has_next())exit
      field => iter%next()
      select type(field)
      type is (field_type)
        field_ptr => field
        ! Create a coupling external field and call copy_to_lfric
        ! to receive the coupling field from Oasis
        call coupler_exchange_2d%initialise(field_ptr, local_index)
        call coupler_exchange_2d%set_time(modeldb%clock)
        ! Call through to coupler_receive_2d in coupler_exchange_2d_mod
        call coupler_exchange_2d%copy_to_lfric(ierror)
        call coupler_exchange_2d%clear()
      class default
        write(log_scratch_space, '(2A)' ) "PROBLEM cpl_rcv: field ", &
                        trim(field%get_name())//" is NOT field_type"
           call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end select
    end do

    if (ierror == 0 .and. l_esm_couple_test) then
      write(log_scratch_space, '(2A)' ) "Skipping updating of prognostics ",&
                            "from coupler (due to l_esm_couple_test=.true.)"
      call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
      ierror = 1
    end if

    if(ierror == 0) then

      ! If exchange is successful then process the data that has
      ! come through the coupler
      call process_recv_fields_2d(cpl_rcv_2d, depository, &
                                  n_sea_ice_tile, T_freeze_h2o_sea,  &
                                  therm_cond_sice, therm_cond_sice_snow )

      ! Update the prognostics
      call iter%initialise(cpl_rcv_2d)
      do
        if(.not.iter%has_next())exit
        field => iter%next()
        call cpl_rcv_2d%get_field( &
                               trim(field%get_name()), field_ptr)
        call field_ptr%write_field(trim(field%get_name()))
        call coupler_update_prognostics(field_ptr, depository)
      end do

    end if

  end subroutine cpl_rcv


  !>@brief Adds field used in coupling code to depository and prognostic_fields
  !>collection
  !>@param [in,out] depository         Field collection - all fields
  !>@param [in,out] prognostic_fields  Prognostic_fields collection
  !>@param [in]     name               Name of the fields to be added
  !>@param [in]     vector_space       Function space of field
  !>@param [in]     checkpoint_flag    Flag to allow checkpoint and
  !>                                   restart behaviour of field to be set
  subroutine add_cpl_field(depository,        &
                           prognostic_fields, &
                           name,              &
                           vector_space,      &
                           checkpoint_flag)
    use io_config_mod,           only : use_xios_io, &
                                        write_diag, checkpoint_write, &
                                        checkpoint_read
    use lfric_xios_read_mod,     only : read_field_generic
    use lfric_xios_write_mod,    only : write_field_generic
    use io_mod,                  only : checkpoint_write_netcdf, &
                                       checkpoint_read_netcdf

    implicit none

    character(*), intent(in)                       :: name
    type(field_collection_type), intent(inout)     :: depository
    type(field_collection_type), intent(inout)     :: prognostic_fields
    type(function_space_type), pointer, intent(in) :: vector_space
    logical(l_def), optional, intent(in)           :: checkpoint_flag

    !Local variables
    !field to initialize
    type(field_type)                               :: new_field
    !pointer to a field
    type(field_type), pointer                      :: field_ptr
    class(pure_abstract_field_type), pointer       :: abstract_field
    !flag for field checkpoint
    logical(l_def)                                 :: checkpointed

    ! pointers for xios write interface
    procedure(write_interface), pointer :: write_behaviour
    procedure(read_interface),  pointer :: read_behaviour
    procedure(checkpoint_write_interface), pointer :: checkpoint_write_behaviour
    procedure(checkpoint_read_interface), pointer  :: checkpoint_read_behaviour

    call new_field%initialise( vector_space, name=trim(name) )

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
      if (write_diag .or. checkpoint_write) &
        call new_field%set_write_behaviour(write_behaviour)
      if (checkpoint_read .and. checkpointed) &
        call new_field%set_read_behaviour(read_behaviour)
    else
      checkpoint_write_behaviour => checkpoint_write_netcdf
      checkpoint_read_behaviour  => checkpoint_read_netcdf
      call new_field%set_checkpoint_write_behaviour(checkpoint_write_behaviour)
      call new_field%set_checkpoint_read_behaviour(checkpoint_read_behaviour)
    endif

    ! Add the field to the depository
    call depository%add_field(new_field)
    call depository%get_field(name, field_ptr)
    ! If checkpointing the field, put a pointer to it
    ! in the prognostics collection
    if ( checkpointed ) then
      abstract_field => field_ptr
      call prognostic_fields%add_reference_to_field( abstract_field )
    endif

  end subroutine add_cpl_field

  !>@brief Populate the coupling field collections from information
  !>       from the coupling configuration
  !>@param [inout] cpl_snd_2d Collection of 2d coupling fields for sending
  !>@param [inout] cpl_rcv_2d Collection of 2d coupling fields for receiving
  !>@param [inout] cpl_snd_0d Collection of scalars (0d data) for sending
  !>@param [in]    depository Field collection of all fields
  !
  subroutine generate_coupling_field_collections(cpl_snd_2d, &
                                                 cpl_rcv_2d, &
                                                 cpl_snd_0d, &
                                                 depository)
    implicit none
    type( field_collection_type ), intent(inout) :: cpl_snd_2d
    type( field_collection_type ), intent(inout) :: cpl_rcv_2d
    type( field_collection_type ), intent(inout) :: cpl_snd_0d
    type( field_collection_type ), intent(in)    :: depository

    ! Names of coupling send fields listed in the coupling configuration
    character(str_def), allocatable              :: send_field_names(:)
    ! Names of coupling receive fields listed in the coupling configuration
    character(str_def), allocatable              :: recv_field_names(:)
    ! Loop variable over the configuration names
    integer(i_def)                               :: nv
    ! Field name
    character(str_def)                           :: name
    ! Index of the coupling prefix  in the name
    integer(i_def)                               :: index_prefix
    ! Index of the multi-category suffix in the name
    integer(i_def)                               :: index_cat
    ! Index of the category "01" in the suffix of the name
    integer(i_def)                               :: index_01
    ! Pointer to a field extracted from the depository
    type(field_type), pointer                    :: field1
    ! Temporary field used
    type(field_type)                             :: field2
    ! Pointer to a field to add a reference to a collection
    class(pure_abstract_field_type), pointer     :: abstract_field

    ! Get the contents of the coupling configuration
    ! Note: this allocates the two arrays - so they will need deallocating later
    call get_coupling_fields(send_field_names, recv_field_names)

    ! Loop over all entries in the namcouple file, looking for lfric fields
    ! as either the source or destination entry
    do nv = 1, size(send_field_names)

      name = trim(adjustl(send_field_names(nv)))
      ! Check the source entry in the namcouple file to see if it is an
      ! lfric field (i.e. starts with the substring: cpl_prefix)
      index_prefix = index(name, cpl_prefix)
      if (index_prefix > 0) then
        ! Check if the field is a multi-category field (i.e. contains cpl_cat)
        index_cat = index(name, cpl_cat)
        ! Multiple category entry (so will need a multidata lfric field)
        if (index_cat > 0) then
          ! Check if name ends in cpl_cat substring followed by 2 characters
          if (len(trim(name)) - index_cat+1 - len(cpl_cat) .ne. 2 ) then
            write(log_scratch_space,'(3A)')             &
               "generate_coupling_field_collections : ", &
               "incorrect send variable name in coupling config: ", name
            call log_event(log_scratch_space, LOG_LEVEL_ERROR)
          end if
          ! Create a multidata field if there are multiple categories
          ! (only create field once - i.e. when name contains cpl_fixed_catno)
          index_01 = index(name, cpl_fixed_catno)
          if (index_01 > 0) then
            ! Extract the field from the depository that has the same name as
            ! the entry in the namcouple (minus the cpl_cat suffix)
            call depository%get_field( name(1:index_cat-1), field1)
            ! Make a copy, so we can apply transient calculations that
            ! we need for coupling, but don't want in the prognostic field
            call field1%copy_field_properties(field2, name(1:index_cat-1))
            ! Add that field to the coupling-send field collection
            call cpl_snd_2d%add_field(field2)
          endif
        ! Create a single field
        else
          ! Deal with fields that will be coupled through a scalar (0d)
          if( (name == 'lf_greenland') .or. &
              (name == 'lf_antarctic') ) then
            ! Extract the field from the depository that has the same name as
            ! the entry in the namcouple
            call depository%get_field(trim(name), field1)
            ! Make a copy, so we can apply transient calculations that
            ! we need for coupling, but don't want in the prognostic field
            call field1%copy_field_properties(field2, trim(name))
            ! Add that field to the coupling send field collection for 0d fields
            call cpl_snd_0d%add_field(field2)
          else
            ! Deal with fields that will be coupled through 2d coupling

            ! Extract the field from the depository that has the same name as
            ! the entry in the namcouple
            call depository%get_field( trim(name),  field1)
            ! Make a copy, so we can apply transient calculations that
            ! we need for coupling, but don't want in the prognostic field
            call field1%copy_field_properties(field2, trim(name))
            ! Add that field to the coupling send field collection for 2d fields
            call cpl_snd_2d%add_field(field2)
          end if
        endif
      endif

      name = trim(adjustl(recv_field_names(nv)))
      ! Check the destination entry in the namcouple file to see if it is an
      ! lfric field (i.e. starts with the substring: cpl_prefix)
      index_prefix = index(name, cpl_prefix)
      if (index_prefix > 0) then
        index_cat = index(name, cpl_cat)
        ! if name ends in cpl_cat substring, check it is the correct length
        if (index_cat > 0 .and. &
           (index_cat - 1 + len(cpl_cat) + len(cpl_fixed_catno) .ne. &
                                                   len(trim(name)))) then
          write(log_scratch_space,'(3A)')             &
               "generate_coupling_field_collections : ", &
               "incorrect receive variable name in coupling config: ", name
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        endif
        if (index_cat > 0) then  ! multiple category
          index_01 = index(name, cpl_fixed_catno)
          if (index_01 > 0) then
            call depository%get_field(name(1:index_cat-1), field1)
            abstract_field => field1
            call cpl_rcv_2d%add_reference_to_field(abstract_field)
          endif
        else
          call depository%get_field(trim(name), field1)
          abstract_field => field1
          call cpl_rcv_2d%add_reference_to_field(abstract_field)
        endif
      endif

    enddo

    if(allocated(send_field_names)) deallocate(send_field_names)
    if(allocated(recv_field_names)) deallocate(recv_field_names)

  end subroutine generate_coupling_field_collections

  !>@brief Puts the coupling component name into modeldb
  !>@param[inout] modeldb The working data set for a model run
  !>@param[in]    name The name of this coupling component
  !
  subroutine set_cpl_name(modeldb, name)
    type( modeldb_type ), intent(inout) :: modeldb
    character(*),         intent(in)    :: name
    call modeldb%values%add_key_value('cpl_name', name)
  end subroutine set_cpl_name

end module coupler_mod
