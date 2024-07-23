!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> &brief coupling related routines for use in coupled configuration

module coupler_mod
#ifdef MCT
  use mod_oasis,                      only: oasis_init_comp,        &
                                            oasis_get_localcomm, oasis_abort,  &
                                            oasis_terminate, oasis_enddef,     &
                                            oasis_def_var, oasis_def_partition,&
                                            oasis_out, prism_ok, nnamcpl,      &
                                            namsrcfld, namdstfld, oasis_in,    &
                                            oasis_get_ncpl, oasis_get_freqs,   &
                                            prism_real
#endif
  use accumulate_send_fields_2d_mod,  only: accumulate_send_fields_2d
  use driver_water_constants_mod,     only: T_freeze_h2o_sea
  use surface_config_mod,             only: therm_cond_sice, &
                                            therm_cond_sice_snow
  use field_mod,                      only: field_type
  use field_parent_mod,               only: field_parent_type
  use pure_abstract_field_mod,        only: pure_abstract_field_type
  use function_space_mod,             only: function_space_type
  use finite_element_config_mod,      only: element_order
  use fs_continuity_mod,              only: W3, Wtheta
  use psykal_lite_mod,                only: invoke_nodal_coordinates_kernel
  use function_space_collection_mod,  only: function_space_collection
  use field_collection_iterator_mod,  only: field_collection_iterator_type
  use field_collection_mod,           only: field_collection_type
  use sort_mod,                       only: bubble_sort
  use constants_mod,                  only: i_def, r_def, str_def, l_def, &
                                            i_halo_index, imdi, rmdi
  use timestepping_config_mod,        only: dt
  use log_mod,                        only: log_event,       &
                                            LOG_LEVEL_INFO,  &
                                            LOG_LEVEL_DEBUG, &
                                            LOG_LEVEL_ERROR, &
                                            log_scratch_space
  use mesh_mod,                       only: mesh_type
  use model_clock_mod,                only: model_clock_type
  use mpi_mod,                        only: global_mpi
  use field_parent_mod,               only: write_interface, read_interface,  &
                                            checkpoint_write_interface,       &
                                            checkpoint_read_interface
  use process_send_fields_0d_mod,     only: process_send_fields_0d
  use process_send_fields_2d_mod,     only: process_send_fields_2d,           &
                                            cpl_reset_field,                  &
                                            initialise_extra_coupling_fields, &
                                            acc_step, ldump_prep,             &
                                            r_sea_ice_frac_raw
  use coupler_exchange_2d_mod,        only: coupler_exchange_2d_type
  use coupler_exchange_0d_mod,        only: coupler_send_0d
  use coupler_update_prognostics_mod, only: coupler_update_prognostics,       &
                                            initialise_snow_mass
  use process_recv_fields_2d_mod,     only: process_recv_fields_2d
  use derived_config_mod,             only: l_esm_couple
  use esm_couple_config_mod,          only: l_esm_couple_test

#if defined(UM_PHYSICS)
  use jules_control_init_mod,         only: n_sea_tile, first_sea_tile
  ! Note: n_sea_ice_tile has to be retrieved from surface_config_mod and not
  !       jules_control_init_mod as the coupler is initialised before jules
  use surface_config_mod,             only: n_sea_ice_tile
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

  !maximum number of components lfric can send the same data
  integer(i_def), parameter             :: nmax = 8

  !Index to sort data for sending
  integer(i_def), allocatable           :: slocal_index(:)

  !name of component in OASIS
  character(len=80),     parameter      :: cpl_name = 'lfric'
#ifdef MCT
  !OASIS component id
  integer(i_def)                        :: il_comp_id
  !OASIS partition id for icesheets
  integer(i_def)                        :: il_part_id
  !keeps info about level
  character(len=2)                      :: cpl_catno
#endif
  !prefix for lfric fields in namcouple
  character(len=3), parameter           :: cpl_prefix = "lf_"
  !prefix for field category (level)
  character(len=4), parameter           :: cpl_cat = "_cat"
  !name of the first level for multi data level field
  character(len=2), parameter           :: cpl_fixed_catno = "01"
  !format for writing the category number
  character(len=6), parameter           :: cpl_fmt = "(i2.2)"

  !routines
  public cpl_init_fields, cpl_fields, cpl_snd, cpl_fld_update, cpl_rcv
  public cpl_initialize, cpl_define, cpl_finalize

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Science routines for coupling
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !>@brief Initializes coupling fields (for sending) to 0
  !>
  !> @param [in,out] cpl_rcv_2d field collection containing fields for sending
  !> @param[in] depository field collection - all fields
  !
  subroutine cpl_init_fields(cpl_rcv_2d, depository)
   implicit none
   type( field_collection_type ), intent(in) :: cpl_rcv_2d
   type( field_collection_type ), intent(in) :: depository
   !local variables
   !iterator over fields in cpl_rcv_2d collection
   class( field_parent_type ), pointer          :: cfield_iter   => null()
   !poiter to a coupling field
   type( field_type ),         pointer          :: cfield        => null()
   !iterator
   type( field_collection_iterator_type)        :: iter

   ! Initilaise accumulation step counter
   acc_step = 0.0_r_def
   ldump_prep = .false.

   call iter%initialise(cpl_rcv_2d)
   do
     if (.not.iter%has_next())exit
     cfield_iter => iter%next()
     select type(cfield_iter)
       type is (field_type)
          call cpl_rcv_2d%get_field(trim(cfield_iter%get_name()), cfield)
          write(log_scratch_space,'(2A)') &
                "cpl_init_fields: set initial value for ", &
                trim(adjustl(cfield%get_name()))
          call log_event(log_scratch_space,LOG_LEVEL_DEBUG)
          call cpl_reset_field(cfield, depository)
          cfield   => null()
       class default
         write(log_scratch_space, '(2A)') "Problem cpl_init_fields: field ", &
                               trim(cfield%get_name())//" is NOT field_type"
         call log_event( log_scratch_space, LOG_LEVEL_ERROR )
       end select
   end do

   nullify(cfield_iter)

  end subroutine cpl_init_fields

  !>@brief Adds fields used in coupling to depository and prognosic_fields
  !>       collections
  !> @details These fields are the raw fields that haven't been processed
  !>          for coupling, so are the accumulations etc. that are passed
  !>          from timestep to timestep.
  !> @param [in]     twod_mesh    mesh on which coupling fields are defined (W3)
  !> @param [in,out] depository   field collection - all fields
  !> @param [in,out] prognostic_fields field collection - prognostic fields
  !
  subroutine cpl_fields( mesh, twod_mesh, depository, prognostic_fields )
   implicit none

   type( mesh_type ),             intent(in), pointer :: mesh
   type( mesh_type ),             intent(in), pointer :: twod_mesh
   type( field_collection_type ), intent(inout)       :: depository
   type( field_collection_type ), intent(inout)       :: prognostic_fields
   !
   !vactor space for coupling field
   type(function_space_type), pointer :: vector_space => null()
   type(function_space_type), pointer :: sice_space => null()
   type(function_space_type), pointer :: threed_space => null()
   type(function_space_type), pointer :: wtheta_space => null()
   !checkpoint flag for coupling field
   logical(l_def)                     :: checkpoint_restart_flag

   write(log_scratch_space, * ) "cpl_fields: add coupling fields to repository"
   call log_event( log_scratch_space, LOG_LEVEL_DEBUG )

   vector_space=> function_space_collection%get_fs( twod_mesh, 0, W3 )
   sice_space  => function_space_collection%get_fs( twod_mesh, 0, W3,          &
                                                               n_sea_ice_tile )
   threed_space => function_space_collection%get_fs(mesh, 0, W3 )
   wtheta_space => function_space_collection%get_fs(mesh, 0, Wtheta)
   !coupling fields
   !sending-depository

   ! these need to be in restart file as they are the accumulations used to
   ! generate the fields passed to the ocean or river model
   checkpoint_restart_flag = .true.
   call add_cpl_field(depository, prognostic_fields, &
        'lf_taux',   vector_space, checkpoint_restart_flag)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_tauy',   vector_space, checkpoint_restart_flag)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_w10',   vector_space, checkpoint_restart_flag)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_solar',   vector_space, checkpoint_restart_flag)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_heatflux',   vector_space, checkpoint_restart_flag)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_train',   vector_space, checkpoint_restart_flag)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_tsnow',   vector_space, checkpoint_restart_flag)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_rsurf',   vector_space, checkpoint_restart_flag)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_rsub',   vector_space, checkpoint_restart_flag)

   ! The following fields are taken care of elsewhere (theoretically)
   ! but we might need duplicates for coupling restarts.
   call add_cpl_field(depository, prognostic_fields, &
        'lf_evap',   vector_space, checkpoint_restart_flag)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_topmelt',   sice_space, checkpoint_restart_flag)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_iceheatflux',sice_space, checkpoint_restart_flag)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_sublimation',sice_space, checkpoint_restart_flag)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_iceskint',sice_space, checkpoint_restart_flag)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_pensolar',sice_space, checkpoint_restart_flag)

   ! The following fields don't need to be in checkpoint files as they are
   ! calculated instantaneously from snow depth just before coupling
   call add_cpl_field(depository, prognostic_fields, &
        'lf_antarctic', vector_space, .false.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_greenland', vector_space, .false.)

   !receiving - depository
   vector_space => function_space_collection%get_fs( twod_mesh, 0, W3, ndata=1 )


   ! These do not need to be in the restart file because they come FROM
   ! the ocean/seaice model!
   checkpoint_restart_flag = .false.

   call add_cpl_field(depository, prognostic_fields, &
        'lf_ocn_sst', vector_space, checkpoint_restart_flag)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_icefrc',   sice_space, checkpoint_restart_flag)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_icetck',   sice_space, checkpoint_restart_flag)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_icelayert',sice_space, checkpoint_restart_flag)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_conductivity',sice_space, checkpoint_restart_flag)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_snow_depth',sice_space, checkpoint_restart_flag)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_pond_frac',sice_space, checkpoint_restart_flag)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_pond_depth',sice_space, checkpoint_restart_flag)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_sunocean', vector_space, checkpoint_restart_flag)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_svnocean', vector_space, checkpoint_restart_flag)

   ! Special 3D fields needed when converting incoming ocean u/v
   ! from W3 (cell centres) to W2 (cell faces)
   call add_cpl_field(depository, prognostic_fields, &
        'sea_u_3d', threed_space, checkpoint_restart_flag)

   call add_cpl_field(depository, prognostic_fields, &
        'sea_v_3d', threed_space, checkpoint_restart_flag)

   call add_cpl_field(depository, prognostic_fields, &
        'sea_w_3d', wtheta_space, checkpoint_restart_flag)

  end subroutine cpl_fields

  !>@brief Top level routine for setting coupling fields and sending data
  !>
  !> @param [in,out] cpl_snd_2d field collection with fields sent to another
  !>                            component as 2d data
  !> @param [in,out] cpl_snd_2d field collection with fields sent to another
  !>                            component as 0d data (scalars)
  !> @param [in]    depository  field collection - all fields
  !> @param [in]    model_clock Time within the model.
  !
  subroutine cpl_snd(cpl_snd_2d, cpl_snd_0d, depository, model_clock)

    implicit none

    type( field_collection_type ), intent(in)    :: cpl_snd_2d
    type( field_collection_type ), intent(in)    :: cpl_snd_0d
    type( field_collection_type ), intent(in)    :: depository
    class(model_clock_type),       intent(in)    :: model_clock

    !local variables
    !pointer to a field (parent)
    class( field_parent_type ), pointer          :: field   => null()
    !pointer to a field
    type( field_type ), pointer                  :: field_ptr   => null()
    !field pointer to field to be written
    type( field_type ), pointer                  :: dep_fld
    !iterator
    type( field_collection_iterator_type)        :: iter
    !pointer to sea ice fractions
    type( field_type ),         pointer          :: ice_frac_ptr   => null()
    ! External field used for sending data to Oasis
    type(coupler_exchange_2d_type)               :: coupler_exchange_2d
    !name of the field
    character(str_def)                           :: sname
    ! Oasis variable id
    integer(i_def)                               :: var_id
    !Ice sheet mass scalar
    real( r_def )                                :: ice_mass
    ! Returned OASIS error code
    integer(i_def)                               :: ierror

    ldump_prep = .false.

    ! increment accumulation step
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
          call accumulate_send_fields_2d(field_ptr, depository, model_clock)
          ! Create a coupling external field
          call coupler_exchange_2d%initialise(field_ptr, slocal_index)
          call coupler_exchange_2d%set_time(model_clock)
          if( coupler_exchange_2d%is_coupling_time() ) then
            ! Process the accumulations to make the fields that will be coupled
            call process_send_fields_2d(field_ptr, depository, model_clock)
            ! Call through to coupler_send_2d in coupler_exchange_2d_mod
            call coupler_exchange_2d%copy_from_lfric(ierror)
          else
            ! No coupling at this time-step
            ierror=1
            write(log_scratch_space, '(3A)' ) "cpl_snd: field ", &
                           trim(sname), " NOT exchanged on this timestep"
            call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
          end if

          ! If coupling was successful then reset field to 0 ready to start
          ! accumulation for the next exchange
          if (ierror == 0) then
            call cpl_reset_field(field, depository)
            acc_step = 0.0_r_def
          endif

          ! Write out the depository version of the field
          call depository%get_field(trim(sname), dep_fld)
          call dep_fld%write_field(trim(sname))

        class default
          write(log_scratch_space, '(2A)' ) "PROBLEM cpl_snd: field ", &
                trim(field%get_name())//" is NOT field_type"
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end select

    end do

    ! Send fields using 0d coupling (scalars)
    call iter%initialise(cpl_snd_0d)
    do
      if ( .not. iter%has_next() ) exit
      field => iter%next()
      select type(field)
        type is (field_type)
          field_ptr => field
          sname = field%get_name()
          var_id = field%get_cpl_id(1)
          call process_send_fields_0d(field_ptr, depository, model_clock, &
                                      ice_mass)
          call coupler_send_0d(ice_mass, sname, var_id, model_clock )
        class default
          write(log_scratch_space, '(2A)' ) "PROBLEM cpl_snd: field ", &
                trim(field%get_name())//" is NOT field_type"
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end select

    end do
    ice_frac_ptr   => null()
    nullify(field)

  end subroutine cpl_snd


  !>@brief Top level routine for updating coupling fields
  !>
  !> @param [in,out] cpl_snd_2d field collection with fields sent to another
  !>                         component
  !> @param [in]    depository field collection - all fields
  !> @param [in]    model_clock Time within the model.
  !
  subroutine cpl_fld_update(cpl_snd_2d, depository, model_clock)
   implicit none
   type( field_collection_type ), intent(in)    :: cpl_snd_2d
   type( field_collection_type ), intent(in)    :: depository
   class(model_clock_type),       intent(in)    :: model_clock
   !local variables
   !pointer to a field (parent)
   class( field_parent_type ), pointer          :: field   => null()
   !pointer to a field
   type( field_type ), pointer                  :: field_ptr   => null()
   !field pointer to field to be written
   type( field_type ), pointer                  :: dep_fld
   !iterator
   type( field_collection_iterator_type)        :: iter
   !name of the field
   character(str_def)                           :: sname

   ldump_prep = .true.
   acc_step = acc_step + 1.0_r_def

   ! We need to loop over each output field and ensure it gets updated
   call iter%initialise(cpl_snd_2d)
   do
      if(.not.iter%has_next())exit
      field => iter%next()
      select type(field)
        type is (field_type)
          field_ptr => field
          call process_send_fields_2d(field_ptr, depository, model_clock)

          ! Write out the depository version of the field
          sname = trim(adjustl(field%get_name()))
          call depository%get_field(trim(sname), dep_fld)
          call dep_fld%write_field(trim(sname))
        class default
          write(log_scratch_space, '(2A)' ) "PROBLEM cpl_fld_update: field ", &
                         trim(field%get_name())//" is NOT field_type"
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end select
   end do

   acc_step = 0.0_r_def

   nullify(field)

  end subroutine cpl_fld_update


  !>@brief Top level routine for receiving data
  !>
  !> @param [in,out] cpl_rcv_2d field collection with names of the fields received
  !>                          from another component
  !> @param [in] depository   field collection - all fields
  !> @param [in] model_clock  Time within the model.
  !
  subroutine cpl_rcv(cpl_rcv_2d, depository, model_clock)
   implicit none
   type( field_collection_type ), intent(in)    :: cpl_rcv_2d
   type( field_collection_type ), intent(in)    :: depository
   class(model_clock_type),       intent(in)    :: model_clock
   !local variables
   !pointer to a field (parent)
   class( field_parent_type ), pointer          :: field   => null()
   !pointer to a field
   type( field_type ), pointer                  :: field_ptr => null()
   ! External field used for receiving data from Oasis
   type(coupler_exchange_2d_type)               :: coupler_exchange_2d
   !iterator
   type( field_collection_iterator_type)        :: iter
   !flag for processing data that has just been exchanged
   ! (set to 1 once data has been successfully passed through the coupler)
   integer(i_def)                               :: ierror

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
          call coupler_exchange_2d%initialise(field_ptr, slocal_index)
          call coupler_exchange_2d%set_time(model_clock)
          ! Call through to coupler_receive_2d in coupler_exchange_2d_mod
          call coupler_exchange_2d%copy_to_lfric(ierror)
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
      ierror = 0
   end if

   if(ierror == 0) then

      ! If exchange is successful then process the data that has
      ! come through the coupler
      call process_recv_fields_2d(cpl_rcv_2d, depository,             &
                                  n_sea_ice_tile, T_freeze_h2o_sea, &
                                  therm_cond_sice, therm_cond_sice_snow )

      ! Update the prognostics
      call iter%initialise(cpl_rcv_2d)
      do
         if(.not.iter%has_next())exit
         field => iter%next()
         call cpl_rcv_2d%get_field(trim(field%get_name()), field_ptr)
         call field_ptr%write_field(trim(field%get_name()))
         call coupler_update_prognostics(field_ptr, depository)
      end do

   end if

   nullify(field)

  end subroutine cpl_rcv


  !> @brief Adds field used in coupling code to depository and prognostic_fields
  !> collection
  !>
  !> @param [in,out] depository   field collection - all fields
  !> @param [in,out] prognostic_fields  prognostic_fields collection
  !> @param [in]    name name of the fields to be added
  !> @param [in]    vector_space Function space of field to set behaviour for
  !> @param [in]    checkpoint_flag Flag to allow checkpoint and
  !>                                restart behaviour of field to be set
  !
  subroutine add_cpl_field(depository, prognostic_fields, &
                               name, vector_space, &
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
   type(field_type), pointer                      :: field_ptr => null()
   class(pure_abstract_field_type), pointer       :: tmp_ptr => null()
   !flag for field checkpoint
   logical(l_def)                                 :: checkpointed

   ! pointers for xios write interface
   procedure(write_interface), pointer :: write_behaviour => null()
   procedure(read_interface),  pointer :: read_behaviour => null()
   procedure(checkpoint_write_interface), pointer ::                           &
                                           checkpoint_write_behaviour => null()
   procedure(checkpoint_read_interface), pointer  ::                           &
                                            checkpoint_read_behaviour => null()

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
     tmp_ptr => field_ptr
     call prognostic_fields%add_reference_to_field( tmp_ptr )
   endif

  end subroutine add_cpl_field

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Technical routines for coupling
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !>@brief Initializes OASIS coupler
  !>
  !> @param [out] comm_out Communicator returned from OASIS to run the model in
  !> @param [in]  comm_in  Input communicator that OASIS can split
  !
  subroutine cpl_initialize(comm_out, comm_in)
   implicit none
   integer(i_def),   intent(out) :: comm_out
   integer(i_def),   intent(in)  :: comm_in
#ifdef MCT
   integer(i_def)                :: ierror ! error return by OASIS
   call oasis_init_comp (il_comp_id, cpl_name, ierror, commworld=comm_in)

   if (ierror .NE. prism_ok) then
     call oasis_abort(il_comp_id, trim(cpl_name), 'cpl_initialize')
   endif

   call oasis_get_localcomm ( comm_out, ierror)

   if (ierror .NE. prism_ok) then
      call oasis_abort(il_comp_id, trim(cpl_name), 'cpl_initialize')
   endif
#else
   comm_out = -1
   write(log_scratch_space, * ) &
        "cpl_initialize: to use OASIS cpp directive MCT must be set"
   call log_event( log_scratch_space, LOG_LEVEL_ERROR )
#endif
  end subroutine cpl_initialize


  !>@brief Defines grid for coupling and initializes lfric component in OASIS
  !>
  !> @param [in]     twod_mesh  2D mesh on which fields are defined (W3)
  !> @param [in]     chi        Input coordinate field
  !> @param [in]     depository model depository with all fields
  !> @param [in,out] cpl_snd_2d field collection with fields to send
  !> @param [in,out] cpl_rcv_2d field collection with fields to receive
  !> @param [in,out] cpl_snd_0d field collection with fields to send as scalars
  !
  subroutine cpl_define( twod_mesh, chi, depository, &
                         cpl_snd_2d, cpl_rcv_2d, cpl_snd_0d )
   implicit none

   type( mesh_type ),  intent(in), pointer     :: twod_mesh
   type( field_type ), intent(in)              :: chi(:)
   type( field_collection_type ), intent(in)   :: depository
   type( field_collection_type ), intent(inout):: cpl_snd_2d
   type( field_collection_type ), intent(inout):: cpl_rcv_2d
   type( field_collection_type ), intent(inout):: cpl_snd_0d

#ifdef MCT
   !index for different do loops
   integer(i_def)                              :: i
   !number of levels in the mesh
   integer(i_def)                              :: num_levels
   !coordinates
   type( field_type )                          :: coord_output(3)
   !function space for coupling field
   type(function_space_type), pointer          :: fld_cpld_fs
   type(function_space_type), pointer          :: sice_space => null()
   !global index for the mesh
   integer(i_halo_index), pointer              :: global_index(:)
   !global index for the first mesh level
   integer(i_def), allocatable                 :: sglobal_index(:)
   !partition for OASIS
   integer(i_def), allocatable                 :: ig_paral(:)
   !partition for OASIS for icesheets
   integer(i_def)                              :: ig_paral_isheets(3)
   !rank/bundles of coupling fields
   integer(i_def)                              :: il_var_nodims(2)
   !dimension of coupled fields
   integer(i_def)                              :: var_shape(2)
   !dimension of icesheet mass fields
   integer(i_def)                              :: imass_shape(1)
   !error return by oasis routine
   integer(i_def)                              :: ierror
   !loop index
   integer(i_def)                              :: nv
   !index of cpl_prefix in the string
   integer(i_def)                              :: index_prefix
   !index of cpl_prefix in the string (receive)
   integer(i_def)                              :: indr
   !index of category (cpl_cat) in the string
   integer(i_def)                              :: index_cat
   !index of cpl_catno in the string
   integer(i_def)                              :: index_01
   !number of data levls
   integer(i_def)                              :: ndata
   !pointer to a field from depository - moving to field1, so field will need to be removed
   type(field_type), pointer                   :: field, field1
   ! temporary field
   type(field_type)                            :: field2
   !pointer to a field to add a reference to a collection
   class(pure_abstract_field_type), pointer    :: tmp_ptr => null()
   !pointer to a field type
   class( field_parent_type ), pointer         :: field_itr => null()
   !pointer to a field
   type( field_type ), pointer                 :: field_ptr => null()
   !ID for transient fields (receive)
   integer(i_def)                              :: oasis_rvar_id
   !name for transient fields (receive)
   character(str_def)                          :: rvar_name
   !name with level information for transient field (receive)
   character(str_def)                          :: rvar_name_lev
   !name of the source entry in the namcouple file
   character(str_def)                          :: source_name
   !name of the destination entry in the namcouple file
   character(str_def)                          :: dest_name
   !ID for transient fields (send)
   integer(i_def)                              :: oasis_svar_id
   !name for transient fields (send)
   character(str_def)                          :: svar_name
   !name with level information for transient field (send)
   character(str_def)                      :: svar_name_lev
   !length of the string used to determine if variable has multiple levels
   integer(i_def)                              :: islgth
   !iterator
   type( field_collection_iterator_type)       :: iter
   !rank number of current PE
   integer(i_def) :: local_rank
   ! Number of coupling components the data will be sent to
   integer(i_def)                              :: ncpl
   ! Coupling frequency of each model
   integer(i_def)                              :: cpl_freqs(nmax)
   !length of the snd_field/rcv_field
   integer(i_def)                              :: icpl_size

   nullify( fld_cpld_fs, global_index )

   num_levels = twod_mesh%get_nlayers()

   if (num_levels > 1) then
      write(log_scratch_space,'(2A)') "cpl_define: only 2D mesh can be used", &
         " to define grid for OASIS"
         call log_event( log_scratch_space, LOG_LEVEL_ERROR )
   endif

   fld_cpld_fs => function_space_collection%get_fs( twod_mesh,     &
                                                    element_order, &
                                                    W3 )

   !fields holding output coordinates
   do i = 1,3
     call coord_output(i)%initialise( vector_space = fld_cpld_fs)
   end do

   !Convert field to physical nodal output & sample chi on nodal points
   call invoke_nodal_coordinates_kernel(coord_output, chi)

   icpl_size = fld_cpld_fs%get_last_dof_owned()

   allocate(sglobal_index(icpl_size))
   allocate(slocal_index(icpl_size))

   global_index => fld_cpld_fs%get_global_dof_id()

   if (maxval(global_index) > int(huge(i_def), i_halo_index)) then
      write(log_scratch_space,'(3A)') "cpl_define: global index", &
         " outside default intager range"
         call log_event( log_scratch_space, LOG_LEVEL_ERROR )
   endif

   sglobal_index(1:icpl_size) =                                   &
           int(global_index(1:int(icpl_size, i_halo_index)), i_def)

   do i = 1, icpl_size
     slocal_index(i) = i
   enddo

   !sort global index to improve OASIS performance
   call bubble_sort(icpl_size, sglobal_index, slocal_index)

   !oasis partition
   il_var_nodims(1) = 1 ! rank of coupling field
   il_var_nodims(2) = 1 ! number of bundles in coupling field (always 1)

   allocate(ig_paral(2+icpl_size))
   ig_paral(1) = 4
   ig_paral(2) = icpl_size

   do i = 1, icpl_size
     ig_paral(i + 2) = sglobal_index(i) + 1
   enddo

   var_shape(1) = 1
   var_shape(2) = icpl_size
   imass_shape(1) = 1

   call oasis_def_partition (il_comp_id, ig_paral, ierror)

   !Set up a special partition for 0D icesheet coupling
   local_rank  = global_mpi%get_comm_rank()
   ig_paral_isheets(1)=0
   ig_paral_isheets(2)=0
   if (local_rank == 0 ) then
     ig_paral_isheets(3)=1
   else
     ig_paral_isheets(3)=0
   endif
   call oasis_def_partition (il_part_id, ig_paral_isheets, ierror)

   ! Loop over all entries in the namcouple file, looking for lfric fields
   ! as either the source or destination entry
   do nv = 1, nnamcpl

      source_name = trim(adjustl(namsrcfld(nv)))
      ! Check the source entry in the namcouple file to see if it is an
      ! lfric field (i.e. starts with the substring: cpl_prefix)
      index_prefix = index(source_name, cpl_prefix)
      if (index_prefix > 0) then
        ! Check if the field is a multi-category field (i.e. contains cpl_cat)
        index_cat = index(source_name, cpl_cat)
        ! Multiple category entry (so will need a multidata lfric field)
        if (index_cat > 0) then
          ! Check if name ends in cpl_cat substring followed by 2 characters
          if (len(trim(source_name)) - index_cat+1 - len(cpl_cat) .ne. 2 ) then
             write(log_scratch_space,'(2A)') &
               "cpl_define : incorrect variable name in namcouple: ", &
               source_name
             call log_event(log_scratch_space, LOG_LEVEL_ERROR)
          end if
          ! Create a multidata field if there are multiple categories
          ! (only create field once - i.e. when name contains cpl_fixed_catno)
          index_01 = index(source_name, cpl_fixed_catno)
          if (index_01 > 0) then
            ! Extract the field from the depository that has the same name as
            ! the entry in the namcouple (minus the cpl_cat suffix)
            call depository%get_field( source_name(1:index_cat-1), field1)
            ! Make a copy, so we can apply transient calculations that
            ! we need for coupling, but don't want in the prognostic field
            call field1%copy_field_properties(field2, source_name(1:index_cat-1))
            ! Add that field to the coupling-send field collection
            call cpl_snd_2d%add_field(field2)
          endif
        ! Create a single field
        else
          if( (source_name == 'lf_greenland') .or. &
              (source_name == 'lf_antarctic') ) then
            ! Extract the field from the depository that has the same name as
            ! the entry in the namcouple
            call depository%get_field(trim(source_name), field1)
            ! Make a copy, so we can apply transient calculations that
            ! we need for coupling, but don't want in the prognostic field
            call field1%copy_field_properties(field2, trim(source_name))
            ! Add that field to the coupling send field collection for 0d fields
            call cpl_snd_0d%add_field(field2)
          else
            ! Extract the field from the depository that has the same name as
            ! the entry in the namcouple
            call depository%get_field( trim(source_name),  field1)
            ! Make a copy, so we can apply transient calculations that
            ! we need for coupling, but don't want in the prognostic field
            call field1%copy_field_properties(field2, trim(source_name))
            ! Add that field to the coupling send field collection for 2d fields
            call cpl_snd_2d%add_field(field2)
          end if
        endif
      endif

      dest_name = trim(adjustl(namdstfld(nv)))
      ! Check the destination entry in the namcouple file to see if it is an
      ! lfric field (i.e. starts with the substring: cpl_prefix)
      indr = index(dest_name, cpl_prefix)
      if (indr > 0) then
        index_cat = index(dest_name, cpl_cat)
        islgth = len(trim(dest_name))
        if (index_cat > 0 .and. &
           (index_cat - 1 + len(cpl_cat) + len(cpl_fixed_catno) .ne. islgth)) then
           !has _cat in name, but no number after it
           write(log_scratch_space,'(3A)') " cpl_define :", &
                 " incorrect variable name in namcouple (ends with _cat): ", &
                                                           trim(namdstfld(nv))
           call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        endif
        if (index_cat > 0) then  ! multiple category
           index_01 = index(dest_name, cpl_fixed_catno)
           if (index_01 > 0) then
              call depository%get_field(dest_name(1:index_cat-1), field)
              tmp_ptr => field
              call cpl_rcv_2d%add_reference_to_field(tmp_ptr)
           endif
        else
              call depository%get_field(trim(dest_name), field)
              tmp_ptr => field
              call cpl_rcv_2d%add_reference_to_field(tmp_ptr)
        endif
      endif

   enddo

   call iter%initialise(cpl_rcv_2d)
   do
      if (.not.iter%has_next())exit
      field_itr => iter%next()
      rvar_name = trim(adjustl(field_itr%get_name()))
      fld_cpld_fs => field_itr%get_function_space()
      ndata = fld_cpld_fs%get_ndata()
      if (ndata > 1) then
        do i = 1, ndata
          write(cpl_catno, cpl_fmt) i
          rvar_name_lev = trim(rvar_name)//cpl_cat//cpl_catno
          call oasis_def_var( oasis_rvar_id, trim(rvar_name_lev), il_comp_id, &
                il_var_nodims, oasis_in, var_shape, prism_real, ierror)
          call field_itr%set_cpl_id(oasis_rvar_id, i)

          write(log_scratch_space, '(A)' ) &
                    "cpl_define: field "//trim(rvar_name_lev)//" receive"
          call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
        enddo
      else
         call oasis_def_var( oasis_rvar_id, trim(rvar_name), il_comp_id, &
               il_var_nodims, oasis_in, var_shape, prism_real, ierror)
         call field_itr%set_cpl_id(oasis_rvar_id, 1)

         write(log_scratch_space, '(A)' ) &
                      "cpl_define: field "//trim(rvar_name)//" receive"
         call log_event( log_scratch_space, LOG_LEVEL_DEBUG )

      endif
      field_itr   => null()
   end do

   call iter%initialise(cpl_snd_2d)
   do
      if (.not.iter%has_next())exit
      field_itr => iter%next()
      svar_name     = trim(adjustl(field_itr%get_name()))
      fld_cpld_fs => field_itr%get_function_space()
      ndata = fld_cpld_fs%get_ndata()
      if (ndata > 1) then
         do i = 1, ndata
           write(cpl_catno, cpl_fmt) i
           svar_name_lev = trim(svar_name)//cpl_cat//cpl_catno
           call oasis_def_var( oasis_svar_id, trim(svar_name_lev), il_comp_id, &
                 il_var_nodims, oasis_out, var_shape, prism_real, ierror)
           call field_itr%set_cpl_id(oasis_svar_id, i)

           write(log_scratch_space, '(A)' ) &
                       "cpl_define: field "//trim(svar_name_lev)//" send"
           call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
        enddo
      else
         call oasis_def_var( oasis_svar_id, trim(svar_name), il_comp_id, &
               il_var_nodims, oasis_out, var_shape, prism_real, ierror)
         call field_itr%set_cpl_id(oasis_svar_id, 1)

         write(log_scratch_space, '(A)' ) &
                          "cpl_define: field "//trim(svar_name)//" send"
         call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
      endif
      field_itr   => null()
      field_ptr   => null()
   end do

   call iter%initialise(cpl_snd_0d)
   do
      if (.not.iter%has_next())exit
      field_itr => iter%next()
      svar_name     = trim(adjustl(field_itr%get_name()))
      call oasis_def_var( oasis_svar_id, trim(svar_name), il_part_id, &
            il_var_nodims, oasis_out, imass_shape, prism_real, ierror)
      call field_itr%set_cpl_id(oasis_svar_id, 1)

      write(log_scratch_space, '(A)' ) &
                       "cpl_define: field "//trim(svar_name)//" send"
      call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
   end do

   call oasis_enddef (ierror)

   ! Check each coupling field as the coupling frequency must be the same for
   ! all components the field is sent to
   call iter%initialise(cpl_snd_2d)
   do
     if (.not.iter%has_next())exit
     field_itr => iter%next()
     svar_name     = trim(adjustl(field_itr%get_name()))
     fld_cpld_fs => field_itr%get_function_space()
     ndata = fld_cpld_fs%get_ndata()
     if (ndata > 1) then
       do i = 1, ndata
         write(cpl_catno, cpl_fmt) i
         svar_name_lev = trim(svar_name)//cpl_cat//cpl_catno
         oasis_svar_id = field_itr%get_cpl_id(i)
         call oasis_get_ncpl(oasis_svar_id, ncpl, ierror)
         call oasis_get_freqs(oasis_svar_id, oasis_out, ncpl, &
                              cpl_freqs(1:ncpl), ierror)
         if (maxval(cpl_freqs(1:ncpl)) /= minval(cpl_freqs(1:ncpl))) then
           write(log_scratch_space, '(3A)' ) "ERROR: coupling field ", &
                  trim(svar_name_lev),                                 &
                  " has different coupling frequencies for different components"
           call log_event( log_scratch_space, LOG_LEVEL_ERROR )
         endif
       enddo
     else
       oasis_svar_id = field_itr%get_cpl_id(1)
       call oasis_get_ncpl(oasis_svar_id, ncpl, ierror)
       call oasis_get_freqs(oasis_svar_id, oasis_out, ncpl, &
                            cpl_freqs(1:ncpl), ierror)
       if (maxval(cpl_freqs(1:ncpl)) /= minval(cpl_freqs(1:ncpl))) then
         write(log_scratch_space, '(3A)' ) "ERROR: coupling field ", &
                trim(svar_name),                                     &
                " has different coupling frequencies for different components"
         call log_event( log_scratch_space, LOG_LEVEL_ERROR )
       endif
     endif
   end do
   call iter%initialise(cpl_snd_0d)
   do
     if (.not.iter%has_next())exit
     field_itr => iter%next()
     svar_name     = trim(adjustl(field_itr%get_name()))
     oasis_svar_id = field_itr%get_cpl_id(1)

     call oasis_get_ncpl(oasis_svar_id, ncpl, ierror)
     call oasis_get_freqs(oasis_svar_id, oasis_out, ncpl, &
                          cpl_freqs(1:ncpl), ierror)
     if (maxval(cpl_freqs(1:ncpl)) /= minval(cpl_freqs(1:ncpl))) then
       write(log_scratch_space, '(3A)' ) "ERROR: coupling field ", &
              trim(svar_name),                                     &
              " has different coupling frequencies for different components"
       call log_event( log_scratch_space, LOG_LEVEL_ERROR )
     endif
   end do

   sice_space  => function_space_collection%get_fs(twod_mesh, 0, W3,           &
                                                          ndata=n_sea_ice_tile)

   ! Initialize extra coupling variables
   call initialise_extra_coupling_fields( fld_cpld_fs, sice_space )
   call initialise_snow_mass( sice_space )
   call r_sea_ice_frac_raw%initialise( vector_space = sice_space, &
                                       name = "r_sea_ice_frac_raw" )

   nullify(field)
   nullify(global_index)
   deallocate(sglobal_index)
   deallocate(ig_paral)
#else
   write(log_scratch_space, * ) &
                      "cpl_define: to use OASIS cpp directive MCT must be set"
   call log_event( log_scratch_space, LOG_LEVEL_ERROR )
#endif

  end subroutine cpl_define


  !>@brief Finalizes coupler
  !
  subroutine cpl_finalize()
   implicit none
   integer(i_def) :: ierror           ! error flag from OASIS
   ! finalize OASIS only if coupled configuration
   if ( l_esm_couple ) then
#ifdef MCT
      ierror = prism_ok
      call oasis_terminate(ierror)
      if (ierror .NE. prism_ok) then
          write(log_scratch_space,'(A, I4)') "lfric: oasis_terminate error: ", &
                                                                        ierror
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
          call oasis_abort(il_comp_id, 'finalise','abort1')
      else
          write(log_scratch_space,'(A)') " lfric : cpl_finalize OK"
          call log_event( log_scratch_space, LOG_LEVEL_INFO )
      endif
#else
   ierror = 1
   write(log_scratch_space, * ) &
         "cpl_finalize: to use OASIS cpp directive MCT must be set"
   call log_event( log_scratch_space, LOG_LEVEL_ERROR )
#endif
   endif

  end subroutine cpl_finalize

end module coupler_mod
