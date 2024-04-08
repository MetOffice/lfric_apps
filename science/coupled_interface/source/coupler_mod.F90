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
                                            prism_real
#endif
  use driver_water_constants_mod,     only: T_freeze_h2o_sea
  use surface_config_mod,             only: therm_cond_sice, &
                                            therm_cond_sice_snow
  use field_mod,                      only: field_type, field_proxy_type
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
  use constants_mod,                  only: i_def, r_def, i_halo_index, l_def, &
                                            imdi, rmdi
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
  use coupler_diagnostics_mod,        only: cpl_diagnostics, cpl_reset_field, &
                                            initialise_extra_coupling_fields, &
                                            acc_step, ldump_prep
  use coupler_external_field_mod,     only: coupler_external_field_type, &
                                            initialise_send_fields,      &
                                            set_snow_mass_fields
  use coupler_update_prognostics_mod, only: coupler_update_prognostics,       &
                                            initialise_snow_mass
  use process_o2a_algorithm_mod,      only: process_o2a_algorithm
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

  !length of the snd_field/rcv_field
  integer(i_def)                        :: icpl_size

  !Max length of coupling field names.
  !UM uses 20, but NEMO names can be
  !much longer and OASIS caters for a
  !max length of 80 so we use that
  integer(i_def),        parameter      :: slength = 80

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
  character(len=2)                      :: cpl_lev
#endif
  !prefix for lfric fields in namcouple
  character(len=3), parameter           :: cpl_prefix = "lf_"
  !prefix for field category (level)
  character(len=4), parameter           :: cpl_cat = "_cat"
  !name of the first level for multi data level field
  character(len=2), parameter           :: cpl_flev = "01"
  !this is len of cpl_flev
  character(len=6), parameter           :: cpl_fmt = "(i2.2)"

  !routines
  public cpl_finalize, cpl_initialize, cpl_define, cpl_init_fields, &
         cpl_snd, cpl_rcv, cpl_fld_update
  public cpl_fields

  contains

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

  !>@brief Initializes coupling fields (for sending) to 0
  !>
  !> @param [in,out] dcpl_rcv field collection containing fields for sending
  !
  subroutine cpl_init_fields(dcpl_rcv)
   implicit none
   type( field_collection_type ), intent(in) :: dcpl_rcv
   !local variables
   !iterator over fields in dcpl_rcv collection
   class( field_parent_type ), pointer          :: cfield_iter   => null()
   !poiter to a coupling field
   type( field_type ),         pointer          :: cfield        => null()
   !iterator
   type( field_collection_iterator_type)        :: iter

   ! Initilaise accumulation step counter
   acc_step = 0.0
   ldump_prep = .false.

   call iter%initialise(dcpl_rcv)
   do
     if (.not.iter%has_next())exit
     cfield_iter => iter%next()
     select type(cfield_iter)
       type is (field_type)
          call dcpl_rcv%get_field(trim(cfield_iter%get_name()), cfield)
          write(log_scratch_space,'(2A)') &
                "cpl_init_fields: set initial value for ", &
                trim(adjustl(cfield%get_name()))
          call log_event(log_scratch_space,LOG_LEVEL_DEBUG)
          call cpl_reset_field(cfield)
          cfield   => null()
       class default
         write(log_scratch_space, '(2A)') "Problem cpl_init_fields: field ", &
                               trim(cfield%get_name())//" is NOT field_type"
         call log_event( log_scratch_space, LOG_LEVEL_ERROR )
       end select
   end do

   nullify(cfield_iter)

  end subroutine cpl_init_fields

  !>@brief Defines grid for coupling and initializes lfric component in OASIS
  !>
  !> @param [in]    twod_mesh  2D mesh on which fields are defined (W3)
  !> @param [in]    chi        Input coordinate field
  !> @param [in]    depository model depository with all fields
  !> @param [in,out] dcpl_snd   field collection with fields to receive
  !> @param [in,out] dcpl_rcv   field collection with fields to send
  !
  subroutine cpl_define( twod_mesh, chi, depository, dcpl_snd, dcpl_rcv )
   implicit none

   type( mesh_type ),  intent(in), pointer     :: twod_mesh
   type( field_type ), intent(in)              :: chi(:)
   type( field_collection_type ), intent(inout):: dcpl_snd
   type( field_collection_type ), intent(inout):: dcpl_rcv
   type( field_collection_type ), intent(in)   :: depository

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
   !field proxy
   type( field_proxy_type )                    :: field_proxy
   !proxies for coordinates
   type( field_proxy_type ), target            :: proxy_coord_output(3)
   !loop index
   integer(i_def)                              :: nv
   !index of cpl_prefix in the string (send)
   integer(i_def)                              :: inds
   !index of cpl_prefix in the string (receive)
   integer(i_def)                              :: indr
   !index of category (cpl_cat) in the string
   integer(i_def)                              :: indc
   !index of cpl_lev in the string
   integer(i_def)                              :: ind01
   !number of data levls
   integer(i_def)                              :: ndata
   !pointer to a field from depository
   type(field_type), pointer                   :: field => null()
   !pointer to a field to add a reference to a collection
   class(pure_abstract_field_type), pointer    :: tmp_ptr => null()
   !pointer to a field type
   class( field_parent_type ), pointer         :: field_itr => null()
   !pointer to a field
   type( field_type ), pointer                 :: field_ptr => null()
   !ID for transient fields (receive)
   integer(i_def)                              :: oasis_rvar_id
   !name for transient fields (receive)
   character(len=slength)                      :: rvar_name
   !name with level information for transient field (receive)
   character(len=slength)                      :: rvar_name_lev
   !ID for transient fields (send)
   integer(i_def)                              :: oasis_svar_id
   !name for transient fields (send)
   character(len=slength)                      :: svar_name
   !name with level information for transient field (send)
   character(len=slength)                      :: svar_name_lev
   !length of the string used to determine if variable has multiple levels
   integer(i_def)                              :: islgth
   !iterator
   type( field_collection_iterator_type)       :: iter
   !rank number of current PE
   integer(i_def) :: local_rank

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
     !Get proxies for coordinates so we can access them
     proxy_coord_output(i) = coord_output(i)%get_proxy()
   end do

   !Convert field to physical nodal output & sample chi on nodal points
   call invoke_nodal_coordinates_kernel(coord_output, chi)

   icpl_size = proxy_coord_output(1)%vspace%get_last_dof_owned()

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

   !add fields to cpl_snd and cpl_rcv collection
   do nv = 1, nnamcpl
      inds = index(namsrcfld(nv), cpl_prefix)
      if (inds > 0) then
        indc = index(namsrcfld(nv), cpl_cat)
        islgth = len(trim(namsrcfld(nv)))
        if (indc > 0 .and. &
                 (indc - 1 + len(cpl_cat) + len(cpl_flev) .ne. islgth)) then
           !has _cat in name, but no number after it
           write(log_scratch_space,'(3A)') &
            " cpl_define :",               &
            " incorrect variable name in namcouple (ends with _cat): ", &
                                                      trim(namsrcfld(nv))
           call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        endif
        if (indc > 0) then  ! multiple category
           ind01 = index(namsrcfld(nv), cpl_flev)
           if (ind01 > 0) then
              rvar_name = trim(namsrcfld(nv))
              call depository%get_field(rvar_name(1:indc-1), field)
              tmp_ptr => field
              call dcpl_snd%add_reference_to_field(tmp_ptr)
           endif
        else    ! single category
           call depository%get_field(trim(namsrcfld(nv)), field)
           tmp_ptr => field
           call dcpl_snd%add_reference_to_field(tmp_ptr)
        endif
      endif

      indr = index(namdstfld(nv), cpl_prefix)
      if (indr > 0) then
        indc = index(namdstfld(nv), cpl_cat)
        islgth = len(trim(namdstfld(nv)))
        if (indc > 0 .and. &
           (indc - 1 + len(cpl_cat) + len(cpl_flev) .ne. islgth)) then
           !has _cat in name, but no number after it
           write(log_scratch_space,'(3A)') " cpl_define :", &
                 " incorrect variable name in namcouple (ends with _cat): ", &
                                                           trim(namdstfld(nv))
           call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        endif
        if (indc > 0) then  ! multiple category
           ind01 = index(namdstfld(nv), cpl_flev)
           if (ind01 > 0) then
              svar_name = trim(namdstfld(nv))
              call depository%get_field(svar_name(1:indc-1), field)
              tmp_ptr => field
              call dcpl_rcv%add_reference_to_field(tmp_ptr)
           endif
        else
              call depository%get_field(trim(namdstfld(nv)), field)
              tmp_ptr => field
              call dcpl_rcv%add_reference_to_field(tmp_ptr)
        endif
      endif
   enddo

   call iter%initialise(dcpl_rcv)
   do
      if (.not.iter%has_next())exit
      field_itr => iter%next()
      rvar_name     = trim(adjustl(field_itr%get_name()))
      call dcpl_rcv%get_field(rvar_name, field_ptr)
      field_proxy = field_ptr%get_proxy()
      ndata = field_proxy%vspace%get_ndata()
      if (ndata > 1) then
        do i = 1, ndata
            write(cpl_lev, cpl_fmt) i
            rvar_name_lev = trim(rvar_name)//cpl_cat//cpl_lev
            call oasis_def_var( oasis_rvar_id, rvar_name_lev, il_comp_id, &
                  il_var_nodims, oasis_in, var_shape, prism_real, ierror)
            call field_ptr%set_cpl_id(oasis_rvar_id, i)

            write(log_scratch_space, '(A)' ) &
                      "cpl_define: field "//trim(rvar_name_lev)//" receive"
            call log_event( log_scratch_space, LOG_LEVEL_DEBUG )

         enddo
      else
         call oasis_def_var( oasis_rvar_id, rvar_name, il_comp_id, &
               il_var_nodims, oasis_in, var_shape, prism_real, ierror)
         call dcpl_rcv%get_field(rvar_name, field_ptr)
         call field_ptr%set_cpl_id(oasis_rvar_id, 1)

         write(log_scratch_space, '(A)' ) &
                      "cpl_define: field "//trim(rvar_name)//" receive"
         call log_event( log_scratch_space, LOG_LEVEL_DEBUG )

      endif
      field_itr   => null()
      field_ptr   => null()
   end do

   call iter%initialise(dcpl_snd)
   do
      if (.not.iter%has_next())exit
      field_itr => iter%next()
      svar_name     = trim(adjustl(field_itr%get_name()))
      call dcpl_snd%get_field(svar_name, field_ptr)
      field_proxy = field_ptr%get_proxy()
      ndata = field_proxy%vspace%get_ndata()
      if (ndata > 1) then
         do i = 1, ndata
            write(cpl_lev, cpl_fmt) i
            svar_name_lev = trim(svar_name)//cpl_cat//cpl_lev
            call oasis_def_var( oasis_svar_id, svar_name_lev, il_comp_id, &
                  il_var_nodims, oasis_out, var_shape, prism_real, ierror)
            call field_ptr%set_cpl_id(oasis_svar_id, i)

            write(log_scratch_space, '(A)' ) &
                        "cpl_define: field "//trim(svar_name_lev)//" send"
            call log_event( log_scratch_space, LOG_LEVEL_DEBUG )

         enddo
      else if ((svar_name == 'lf_greenland') .OR. &
               (svar_name == 'lf_antarctic')) then
         call oasis_def_var( oasis_svar_id, svar_name, il_part_id, &
               il_var_nodims, oasis_out, imass_shape, prism_real, ierror)
         call field_ptr%set_cpl_id(oasis_svar_id, 1)

         write(log_scratch_space, '(A)' ) &
                          "cpl_define: field "//trim(svar_name)//" send"
         call log_event( log_scratch_space, LOG_LEVEL_DEBUG )

      else
         call oasis_def_var( oasis_svar_id, svar_name, il_comp_id, &
               il_var_nodims, oasis_out, var_shape, prism_real, ierror)
         call field_ptr%set_cpl_id(oasis_svar_id, 1)

         write(log_scratch_space, '(A)' ) &
                          "cpl_define: field "//trim(svar_name)//" send"
         call log_event( log_scratch_space, LOG_LEVEL_DEBUG )

      endif
      field_itr   => null()
      field_ptr   => null()
   end do

   call oasis_enddef (ierror)

   sice_space  => function_space_collection%get_fs(twod_mesh, 0, W3,           &
                                                          ndata=n_sea_ice_tile)

   ! Initialize extra coupling variables
   call initialise_extra_coupling_fields( fld_cpld_fs, sice_space )
   call initialise_send_fields( fld_cpld_fs, sice_space )
   call initialise_snow_mass( sice_space )

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

  !>@brief Adds coupling fields used in the model to depository and
  !> prognosic_fields collections
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

   ! these need to be in restart file to pass to the ocean or river model!
   checkpoint_restart_flag = .true.
   call add_cpl_field(depository, prognostic_fields, &
        'lf_taux',   vector_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_tauy',   vector_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_solar',   vector_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_heatflux',   vector_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_train',   vector_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_tsnow',   vector_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_w10',   vector_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_rsurf',   vector_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_rsub',   vector_space, checkpoint_restart_flag, twod=.true.)

   ! The following fields are taken care of elsewhere (theoretically)
   ! but we might need duplicates for coupling restarts.
   call add_cpl_field(depository, prognostic_fields, &
        'lf_evap',   vector_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_topmelt',   sice_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_iceheatflux',sice_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_pensolar',sice_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_sublimation',sice_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_iceskint',sice_space, checkpoint_restart_flag, twod=.true.)

   ! The following fields don't need to be in checkpoint files as they are
   ! calculated instantaneously from snow depth just before coupling
   call add_cpl_field(depository, prognostic_fields, &
        'lf_antarctic', vector_space, .false., twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_greenland', vector_space, .false., twod=.true.)

   !receiving - depository
   vector_space => function_space_collection%get_fs( twod_mesh, 0, W3, ndata=1 )


   ! These do not need to be in the restart file because they come FROM
   ! the ocean/seaice model!
   checkpoint_restart_flag = .false.

   call add_cpl_field(depository, prognostic_fields, &
        'lf_ocn_sst', vector_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_icefrc',   sice_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_icetck',   sice_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_icelayert',sice_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_conductivity',sice_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_snow_depth',sice_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_pond_frac',sice_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_pond_depth',sice_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_sunocean', vector_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_svnocean', vector_space, checkpoint_restart_flag, twod=.true.)

   ! Special 3D fields needed when converting incoming ocean u/v
   ! from W3 (cell centres) to W2 (cell faces)
   call add_cpl_field(depository, prognostic_fields, &
        'sea_u_3d', threed_space, checkpoint_restart_flag)

   call add_cpl_field(depository, prognostic_fields, &
        'sea_v_3d', threed_space, checkpoint_restart_flag)

   call add_cpl_field(depository, prognostic_fields, &
        'sea_w_3d', wtheta_space, checkpoint_restart_flag)

  end subroutine cpl_fields

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

  !>@brief Top level routine for setting coupling fields and sending data
  !>
  !> @param [in,out] dcpl_snd field collection with fields sent to another
  !>                         component
  !> @param [in]    depository field collection - all fields
  !> @param [in]    model_clock Time within the model.
  !
  subroutine cpl_snd(dcpl_snd, depository, model_clock)

    implicit none

    type( field_collection_type ), intent(in)    :: dcpl_snd
    type( field_collection_type ), intent(in)    :: depository
    class(model_clock_type),       intent(in)    :: model_clock

    !local variables
    !pointer to a field (parent)
    class( field_parent_type ), pointer          :: field   => null()
    !pointer to a field
    type( field_type ), pointer                  :: field_ptr   => null()
    !iterator
    type( field_collection_iterator_type)        :: iter
    !pointer to sea ice fractions
    type( field_type ),         pointer          :: ice_frac_ptr   => null()
    ! External field used for sending data to Oasis
    type(coupler_external_field_type)            :: coupler_external_field
    !name of the field
    character(len=slength)                       :: sname

    ldump_prep = .false.

    ! increment accumulation step
    acc_step = acc_step + 1.0

    call set_snow_mass_fields(depository)

    call iter%initialise(dcpl_snd)
    do
      if ( .not. iter%has_next() ) exit
      field => iter%next()

      select type(field)
        type is (field_type)
          field_ptr => field
          call cpl_diagnostics(field_ptr, depository, model_clock)
          ! Create a coupling external field and call copy_from_lfric
          ! to send the coupling field to Oasis
          call coupler_external_field%initialise_cpl_external_field(field_ptr, &
                              nmax, icpl_size, slength, slocal_index)
          call coupler_external_field%set_coupling_time(model_clock)
          ! Call through to cpl_field_send in coupler_external_field_mod
          call coupler_external_field%copy_from_lfric()
          sname = trim(adjustl(field%get_name()))
          if((sname /= 'lf_greenland') .AND. (sname /= 'lf_antarctic')) &
             call field_ptr%write_field(trim(field%get_name()))
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
  !> @param [in,out] dcpl_snd field collection with fields sent to another
  !>                         component
  !> @param [in]    depository field collection - all fields
  !> @param [in]    model_clock Time within the model.
  !
  subroutine cpl_fld_update(dcpl_snd, depository, model_clock)
   implicit none
   type( field_collection_type ), intent(in)    :: dcpl_snd
   type( field_collection_type ), intent(in)    :: depository
   class(model_clock_type),       intent(in)    :: model_clock
   !local variables
   !pointer to a field (parent)
   class( field_parent_type ), pointer          :: field   => null()
   !pointer to a field
   type( field_type ), pointer                  :: field_ptr   => null()
   !iterator
   type( field_collection_iterator_type)        :: iter
   !pointer to sea ice fractions
   type( field_type ),         pointer          :: ice_frac_ptr        => null()
   !name of the field
   character(len=slength)          :: sname

   ldump_prep = .true.
   acc_step = acc_step + 1.0

   ! We need to loop over each output field and ensure it gets updated
   call iter%initialise(dcpl_snd)
   do
      if(.not.iter%has_next())exit
      field => iter%next()
      select type(field)
        type is (field_type)
          field_ptr => field
          call cpl_diagnostics(field_ptr, depository, model_clock)

          sname = trim(adjustl(field%get_name()))
          if((sname /= 'lf_greenland') .AND. (sname /= 'lf_antarctic')) &
             call field_ptr%write_field(trim(field%get_name()))
        class default
          write(log_scratch_space, '(2A)' ) "PROBLEM cpl_fld_update: field ", &
                         trim(field%get_name())//" is NOT field_type"
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end select
   end do

   acc_step = 0.0

   ice_frac_ptr   => null()
   nullify(field)

  end subroutine cpl_fld_update

  !>@brief Top level routine for receiving data
  !>
  !> @param [in,out] dcpl_rcv field collection with names of the fields received
  !>                          from another component
  !> @param [in] depository   field collection - all fields
  !> @param [in] model_clock  Time within the model.
  !
  subroutine cpl_rcv(dcpl_rcv, depository, model_clock)
   implicit none
   type( field_collection_type ), intent(in)    :: dcpl_rcv
   type( field_collection_type ), intent(in)    :: depository
   class(model_clock_type),       intent(in)    :: model_clock
   !local variables
   !pointer to a field (parent)
   class( field_parent_type ), pointer          :: field   => null()
   !pointer to a field
   type( field_type ), pointer                  :: field_ptr => null()
   ! External field used for receiving data from Oasis
   type(coupler_external_field_type)            :: coupler_external_field
   !iterator
   type( field_collection_iterator_type)        :: iter
   !flag for processing data that has just been exchanged
   ! (set to 1 once data has been successfully passed through the coupler)
   integer(i_def)                               :: exchange_flag

   ! Set defaults
   exchange_flag = 0

   call iter%initialise(dcpl_rcv)
   do
      if (.not.iter%has_next())exit
      field => iter%next()
      select type(field)
        type is (field_type)
          field_ptr => field
          ! Create a coupling external field and call copy_to_lfric
          ! to receive the coupling field from Oasis
          call coupler_external_field%initialise_cpl_external_field(field_ptr, &
                              nmax, icpl_size, slength, slocal_index)
          call coupler_external_field%set_coupling_time(model_clock)
          ! Call through to cpl_field_receive in coupler_external_field_mod
          call coupler_external_field%copy_to_lfric(exchange_flag)
        class default
          write(log_scratch_space, '(2A)' ) "PROBLEM cpl_rcv: field ", &
                        trim(field%get_name())//" is NOT field_type"
             call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end select
   end do

   if (exchange_flag == 1 .and. l_esm_couple_test) then
      write(log_scratch_space, '(2A)' ) "Skipping updating of prognostics ",&
                            "from coupler (due to l_esm_couple_test=.true.)"
      call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
      exchange_flag = 0
   end if

   if(exchange_flag == 1) then

      ! If exchange is successful then process the data that has
      ! come through the coupler
      call process_o2a_algorithm(dcpl_rcv, depository,             &
                                 n_sea_ice_tile, T_freeze_h2o_sea, &
                                 therm_cond_sice, therm_cond_sice_snow )

      ! Update the prognostics
      call iter%initialise(dcpl_rcv)
      do
         if(.not.iter%has_next())exit
         field => iter%next()
         call dcpl_rcv%get_field(trim(field%get_name()), field_ptr)
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
  !> @param [in]    twod            Optional flag to determine if this is a
  !>                                2D field
  !
  subroutine add_cpl_field(depository, prognostic_fields, &
                               name, vector_space, &
                               checkpoint_flag, twod)
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
   logical(l_def), optional, intent(in)           :: twod

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
end module coupler_mod
