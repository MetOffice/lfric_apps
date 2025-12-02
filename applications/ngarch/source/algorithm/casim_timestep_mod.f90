!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Timestep module that only runs the CASIM microphysics scheme
!> @details This module is included as a template for how to create timesteps
!>          for the ngarch project. This timestep can be selected in the
!>          configuration as ngarch>method="casim". To set up new timesteps
!>          to be selectable by configuration they need to be set up in
!>          ngarch/source/driver/override_timestep_mod.f90 and in
!>          ngarch/rose-meta/lfric-ngarch/HEAD/rose-meta.conf.
module casim_timestep_mod

  use constants_mod,                 only: i_def, r_def
  use field_mod,                     only: field_type
  use field_array_mod,               only: field_array_type
  use sci_field_bundle_builtins_mod, only: clone_bundle, &
                                           set_bundle_scalar
  use field_collection_mod,          only: field_collection_type
  use driver_modeldb_mod,            only: modeldb_type
  use log_mod,                       only: log_event,         &
                                           log_scratch_space, &
                                           LOG_LEVEL_INFO,    &
                                           LOG_LEVEL_ERROR,   &
                                           LOG_LEVEL_TRACE
  use mr_indices_mod,                only: nummr
  use timestep_method_mod,           only: timestep_method_type
  use mesh_mod,                      only: mesh_type

  use casim_alg_mod, only: casim_alg

  implicit none

  type, extends(timestep_method_type), public :: casim_timestep_type
    private

  contains
    private

    procedure, public  :: step     => casim_alg_step
    procedure, public  :: finalise => casim_alg_final

  end type casim_timestep_type

contains

  !> @brief Executes a timestep in which only the CASIM scheme is run.
  !> @param [in]  modeldb  The structure that holds model state
  subroutine casim_alg_step( self, modeldb )
    implicit none

    class( casim_timestep_type ), intent(inout) :: self
    type( modeldb_type ), intent(in), target    :: modeldb

    type( field_collection_type ), pointer :: collection, turbulence_fields
    type( field_type ),            pointer :: theta, rho, dtheta_mphys
    type( field_type )                     :: dmr_mphys(nummr), dcfl, dcff, dbcf
    type( field_array_type ),      pointer :: mr
    type( mesh_type ),             pointer :: mesh

    type(field_collection_type),   pointer :: microphysics_fields
    type(field_collection_type),   pointer :: derived_fields
    type(field_collection_type),   pointer :: cloud_fields
    type(field_collection_type),   pointer :: aerosol_fields

    collection => modeldb%fields%get_field_collection( "moisture_fields" )
    microphysics_fields => modeldb%fields%get_field_collection("microphysics_fields")
    derived_fields => modeldb%fields%get_field_collection("derived_fields")
    cloud_fields => modeldb%fields%get_field_collection("cloud_fields")
    aerosol_fields => modeldb%fields%get_field_collection("aerosol_fields")

    call collection%get_field( "mr", mr )
    call clone_bundle( mr%bundle, dmr_mphys, nummr )
    call set_bundle_scalar( 0.0_r_def, dmr_mphys, nummr )

    call microphysics_fields%get_field( 'dtheta_mphys', dtheta_mphys )

    collection => modeldb%fields%get_field_collection( "depository" )
    call collection%get_field( "theta", theta )
    call collection%get_field( "rho", rho )
    mesh => theta%get_mesh()

    turbulence_fields => modeldb%fields%get_field_collection( "turbulence_fields" )

    ! Call the algorithm
    call log_event( "Running CASIM", LOG_LEVEL_INFO )
    call casim_alg( mr%bundle, theta, rho,                  &
                    derived_fields,                         &
                    microphysics_fields,                    &
                    cloud_fields,                           &
                    aerosol_fields,                         &
                    turbulence_fields, mesh,                &
                    dmr_mphys, dtheta_mphys, dcfl, dcff, dbcf )
    call log_event( "CASIM completed", LOG_LEVEL_INFO )

  end subroutine casim_alg_step

  subroutine casim_alg_final( self )
    implicit none
    class( casim_timestep_type ), intent(inout) :: self
  end subroutine casim_alg_final

end module casim_timestep_mod