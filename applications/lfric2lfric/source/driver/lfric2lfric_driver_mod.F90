!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Drives the execution of the lfric2lfric miniapp.
!>
module lfric2lfric_driver_mod

  use constants_mod,          only: str_def, i_def, l_def, r_second
  use driver_fem_mod,         only: final_fem
  use driver_io_mod,          only: final_io
  use driver_modeldb_mod,     only: modeldb_type
  use field_collection_iterator_mod, only: &
                                    field_collection_iterator_type
  use field_collection_mod,   only: field_collection_type
  use field_mod,              only: field_type
  use field_parent_mod,       only: field_parent_type
  use function_space_mod,     only: function_space_type
  use lfric_xios_context_mod, only: lfric_xios_context_type
  use lfric_xios_read_mod,    only: read_checkpoint
  use lfric_xios_write_mod,   only: write_checkpoint
  use log_mod,                only: log_event, &
                                    log_level_info, &
                                    log_scratch_space
  use mesh_collection_mod,    only: mesh_collection
  use mesh_mod,               only: mesh_type
  use namelist_mod,           only: namelist_type
  use sci_checksum_alg_mod,   only: checksum_alg

  !------------------------------------
  ! lfric2lfric modules
  !------------------------------------
  use lfric2lfric_config_mod,         only: regrid_method_map,         &
                                            regrid_method_lfric2lfric, &
                                            regrid_method_oasis
  use lfric2lfric_infrastructure_mod, only: initialise_infrastructure, &
                                            context_dst, context_src,  &
                                            source_collection_name,    &
                                            target_collection_name
  use lfric2lfric_map_regrid_mod,     only: lfric2lfric_map_regrid
  use lfric2lfric_no_regrid_mod,      only: lfric2lfric_no_regrid
  use lfric2lfric_oasis_regrid_mod,   only: lfric2lfric_oasis_regrid

  implicit none

  private
  public initialise, run, finalise

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief          Sets up required state in preparation for run.
  !> @details        Calls the `initialise_infrastructure` subroutine that
  !!                 checks the configuration namelist, initialises meshes,
  !!                 extrusions, XIOS contexts and files, field collections
  !!                 and fields.
  !> @param [in,out] modeldb                 The structure holding model state
  subroutine initialise( modeldb )

    implicit none

    type(modeldb_type), intent(inout) :: modeldb

    call initialise_infrastructure( modeldb )

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief    Performs regridding of a field collection
  !> @details  Populates fields in the source field collection with data
  !!           read from the source XIOS context file, regrids the field
  !!           collection to a destination mesh, switches to the destination
  !!           XIOS context and writes to a checkpoint file.
  !!           TODO: #262 Porting algorithms and kernels
  !!           Over one time step, all regridding is performed by algorithm
  !!           modules specific to the regrid method defined in the
  !!           configuration.
  !!           Fields to be regridded are extracted from the source field
  !!           collection and located in the source dump file using
  !!           `read_checkpoint`. Corresponding source and target field pairs
  !!           are passed to the regridding algorithm, written to fields
  !!           in the destination field collection and then written to
  !!           an outfile.
  !> @param [in,out] modeldb                 The structure that holds model
  !!                                         state
  subroutine run( modeldb )

    implicit none

    type(modeldb_type), intent(inout) :: modeldb

    ! LFRic-XIOS constants
    integer(kind=i_def), parameter :: start_timestep = 1_i_def

    ! Namelist variables
    character(len=str_def) :: start_dump_filename
    character(len=str_def) :: checkpoint_stem_name
    integer(kind=i_def)    :: regrid_method

    ! Local parameters
    type(namelist_type), pointer :: files_nml
    type(namelist_type), pointer :: lfric2lfric_nml

    type(field_collection_type), pointer :: source_fields
    type(field_collection_type), pointer :: target_fields

    type(field_collection_iterator_type) :: iter

    class(field_parent_type), pointer :: field     => null()
    type(field_type),         pointer :: field_src => null()
    type(field_type),         pointer :: field_dst => null()

    type(lfric_xios_context_type), pointer :: io_context

    character(len=str_def)   :: field_name

    real(r_second) :: checkpoint_times(1)

    ! Namelist pointers
    files_nml       => modeldb%configuration%get_namelist('files')
    lfric2lfric_nml => modeldb%configuration%get_namelist('lfric2lfric')

    ! Extract configuration variables
    call files_nml%get_value( 'start_dump_filename', start_dump_filename )
    call files_nml%get_value( 'checkpoint_stem_name', checkpoint_stem_name )
    call lfric2lfric_nml%get_value( 'regrid_method', regrid_method )

    ! Point to source and target field collections
    source_fields => modeldb%fields%get_field_collection(source_collection_name)
    target_fields => modeldb%fields%get_field_collection(target_collection_name)

    call read_checkpoint(source_fields,      &
                         start_timestep,     &
                         start_dump_filename )

    ! Main loop over fields to be processed
    call iter%initialise(source_fields)
    do
      ! Locate the field to be processed in the field collections
      if ( .not.iter%has_next() ) exit
      field => iter%next()
      field_name = field%get_name()

      call source_fields%get_field(field_name, field_src)
      call target_fields%get_field(field_name, field_dst)

      write(log_scratch_space, '(A,A)') "Processing lfric field ", &
                                           trim(field_name)
      call log_event(log_scratch_space, log_level_info)

      ! Regrid source field depending on regrid method
      select case (regrid_method)
        case (regrid_method_map)
          call lfric2lfric_map_regrid(field_dst, field_src)

        case (regrid_method_lfric2lfric)
          write(log_scratch_space, '(A)') &
                              'Regrid method lfric2lfric not implemented yet'
          call log_event(log_scratch_space, log_level_info)

          call lfric2lfric_no_regrid(field_dst)

        case (regrid_method_oasis)
#ifdef MCT
          call lfric2lfric_oasis_regrid(modeldb, &
                                        field_dst, field_src)
#endif
      end select

      ! Free memory of the processed field
      call field_src%field_final()
    end do

    ! Write output
    call modeldb%io_contexts%get_io_context(context_dst, io_context)
    call io_context%set_current()

    checkpoint_times(1) = modeldb%clock%seconds_from_steps(modeldb%clock%get_step())
    call write_checkpoint( target_fields, modeldb%values, modeldb%clock, &
                           checkpoint_stem_name, checkpoint_times )

    ! Write checksum
    call checksum_alg("lfric2lfric", field_collection=target_fields)

  end subroutine run

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief   Tidies up after a run.
  !>
  !> @param [in]     program_name An identifier given to the model being run
  !> @param [in,out] modeldb      The structure that holds model state
  subroutine finalise( program_name, modeldb )

    implicit none

    character(len=*),   intent(in)     :: program_name
    type(modeldb_type), intent(inout)  :: modeldb


    !--------------------------------------------------------------------------
    ! Model finalise
    !--------------------------------------------------------------------------

    call log_event( program_name//': Miniapp completed', log_level_info )

    !-------------------------------------------------------------------------
    ! Driver layer finalise
    !-------------------------------------------------------------------------

    ! Finalise IO
    call final_io( modeldb )

    call final_fem()

  end subroutine finalise

end module lfric2lfric_driver_mod
