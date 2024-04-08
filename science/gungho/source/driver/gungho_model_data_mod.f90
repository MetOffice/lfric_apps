!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Container for Gung Ho model run working data set.
!>
module gungho_model_data_mod

  use field_mod,            only : field_type
  use field_collection_mod, only : field_collection_type

  implicit none

  private

  !> Holds the working data set for a model run and other working state.
  !>
  type, public :: model_data_type

    private

    !> @name Fields needed to time-step the model.
    !> @{

    !> Tracers that should be advected
    type( field_collection_type ), public   :: adv_tracer_last_outer
    type( field_collection_type ), public   :: adv_tracer_all_outer
    !> Second group of tracers that should be advected
    type( field_collection_type ), public   :: con_tracer_last_outer
    type( field_collection_type ), public   :: con_tracer_all_outer
    !> FD fields derived from FE fields for use in physics time-stepping schemes
    type( field_collection_type ), public   :: derived_fields
    !> LBC fields - lateral boundary conditions to run a limited area model
    type( field_collection_type ), public   :: lbc_fields
    !> Fields owned by the radiation scheme
    type( field_collection_type ), public   :: radiation_fields
    !> Fields owned by the microphysics scheme
    type( field_collection_type ), public   :: microphysics_fields
    !> Fields owned by the electric scheme
    type( field_collection_type ), public   :: electric_fields
    !> Fields owned by the orographic drag schemes
    type( field_collection_type ), public   :: orography_fields
    !> Fields owned by the turbulence scheme
    type( field_collection_type ), public   :: turbulence_fields
    !> Fields owned by the convection schemes
    type( field_collection_type ), public   :: convection_fields
    !> Fields owned by the cloud schemes
    type( field_collection_type ), public   :: cloud_fields
    !> Fields owned by the surface exchange scheme
    type( field_collection_type ), public   :: surface_fields
    !> Fields owned by the soil hydrology scheme
    type( field_collection_type ), public   :: soil_fields
    !> Fields owned by the snow scheme
    type( field_collection_type ), public   :: snow_fields
    !> Fields owned by the chemistry schemes
    type( field_collection_type ), public   :: chemistry_fields
    !> Fields owned by the aerosol schemes
    type( field_collection_type ), public   :: aerosol_fields
    !> Fields owned by the stochastic physics schemes
    type( field_collection_type ), public   :: stph_fields
    !> Fields that contain coupling data
    type( field_collection_type ), public :: cpl_snd
    type( field_collection_type ), public :: cpl_rcv
    !> @}

    !> FD fields used to read initial conditions from LFRic-Input files
    type( field_collection_type ), public   :: fd_fields

    !> Fields used to store data read in from ancillary files
    type( field_collection_type ), public   :: ancil_fields

    !> Fields for the tangent linear linearisation state
    type( field_collection_type ), public   :: ls_fields

    contains

  end type model_data_type

end module gungho_model_data_mod
