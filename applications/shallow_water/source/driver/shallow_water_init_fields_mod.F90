!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Container for shallow water model run working data set including
!!        methods to initialise, copy and finalise the data set.
!!
!> @details This module provides a type to hold all the model fields and methods
!!          to initialise (create and read), copy and finalise (write and
!!          destroy) the data contained within the type.
!!
module shallow_water_init_fields_mod

  use constants_mod,                        only: r_second
  use driver_modeldb_mod,                   only: modeldb_type
  use field_mod,                            only: field_type
  use field_collection_mod,                 only: field_collection_type
  use files_config_mod,                     only: checkpoint_stem_name
  use log_mod,                              only: log_event,      &
                                                  LOG_LEVEL_INFO
  use io_config_mod,                        only: checkpoint_read,  &
                                                  checkpoint_write
  use lfric_xios_read_mod,                  only: read_checkpoint
  use lfric_xios_write_mod,                 only: write_checkpoint, &
                                                  write_state
  use mesh_mod,                             only: mesh_type
  use create_shallow_water_prognostics_mod, only: &
                                                 create_shallow_water_prognostics
  use swe_init_fields_alg_mod,              only: swe_init_fields_alg, &
                                                  swe_init_surface_alg

  implicit none

  private

  public :: create_model_data,     &
            finalise_model_data,   &
            initialise_model_data, &
            output_model_data

contains

  !=============================================================================
  !> @brief Create the fields contained in model_data.
  !> @param[in,out] modeldb    The working data set for a model run
  !> @param[in]     mesh       Mesh to initialise variables on
  subroutine create_model_data( modeldb, &
                                mesh )

    implicit none

    type(modeldb_type),    intent(inout) :: modeldb
    type(mesh_type), pointer, intent(in) :: mesh

    type( field_collection_type ), pointer :: depository => null()
    type( field_collection_type ), pointer :: prognostic_fields => null()

    !-------------------------------------------------------------------------
    ! Instantiate the fields
    !-------------------------------------------------------------------------

    depository => modeldb%fields%get_field_collection("depository")
    prognostic_fields => modeldb%fields%get_field_collection("prognostics")

    ! Create prognostics
    call create_shallow_water_prognostics( mesh,                         &
                                           depository,                   &
                                           prognostic_fields )

  end subroutine create_model_data

  !=======================================================================
  !> @brief Initialises the working data set depending on namelist config.
  !> @param[in,out] modeldb    The working data set for a model run
  !> @param[in]     mesh       Mesh to initialise variables on
  subroutine initialise_model_data( modeldb, &
                                    mesh )

    implicit none

    type(modeldb_type),    intent(inout) :: modeldb
    type(mesh_type), pointer, intent(in) :: mesh

    type( field_collection_type ), pointer :: depository => null()
    type( field_collection_type ), pointer :: prognostic_fields => null()

    type(field_type), pointer :: s_geopot

    depository => modeldb%fields%get_field_collection("depository")
    prognostic_fields => modeldb%fields%get_field_collection("prognostics")

    call prognostic_fields%get_field('s_geopot', s_geopot)
    ! Initialise the surface geopotential
    call swe_init_surface_alg(s_geopot)

    ! Initialise prognostic fields
    if (checkpoint_read) then               ! Recorded check point to start from
      call read_checkpoint(depository, &
                       modeldb%clock%get_first_step() - 1, checkpoint_stem_name)
    else                                      ! No check point to start from
      call swe_init_fields_alg(mesh,                &
                               s_geopot,            &
                               prognostic_fields)
    end if

  end subroutine initialise_model_data

  !=============================================================================
  !> @brief Writes out a checkpoint and dump file dependent on namelist
  !!        options.
  !> @param[in,out] modeldb The working data set for the model run
  subroutine output_model_data( modeldb )

    implicit none

    type(modeldb_type),    intent(inout) :: modeldb
    real(r_second) :: checkpoint_times(1)

    type( field_collection_type ), pointer :: prognostic_fields => null()

    ! Get pointers to field collections for use downstream
    prognostic_fields => modeldb%fields%get_field_collection("prognostics")

    !=================== Write fields to checkpoint files ====================!
    if ( checkpoint_write ) then
      checkpoint_times(1) = modeldb%clock%seconds_from_steps(modeldb%clock%get_step())
      call write_checkpoint( prognostic_fields, modeldb%values, modeldb%clock, &
                             checkpoint_stem_name, checkpoint_times )
    end if

  end subroutine output_model_data

  !=============================================================================
  !> @brief Routine to destroy all the field collections in the working data set
  !> @param[in,out] modeldb The working data set for a model run
  subroutine finalise_model_data( modeldb )

    implicit none

    type(modeldb_type),    intent(inout) :: modeldb

    type( field_collection_type ), pointer :: depository => null()
    type( field_collection_type ), pointer :: prognostic_fields => null()

    depository => modeldb%fields%get_field_collection("depository")
    prognostic_fields => modeldb%fields%get_field_collection("prognostics")

    ! Clear all the fields in each field collection
    call depository%clear()
    call prognostic_fields%clear()

    call log_event( 'finalise_model_data: all fields have been cleared', &
                     LOG_LEVEL_INFO )

  end subroutine finalise_model_data

end module shallow_water_init_fields_mod
