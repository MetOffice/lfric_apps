!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Container for everything that makes up the model state
!!
!> @details This module holds all the information required to run an instance
!>          of the model (such as an ensemble member or adjoint etc.). That
!>          includes all the scientific and technical state. Nothing
!>          should be held in global space, unless it is truely global
!!
module gungho_modeldb_mod

  use driver_modeldb_mod,    only: driver_modeldb_type => modeldb_type
  use gungho_model_data_mod, only: model_data_type

  implicit none

  private

  !> Holds the technical and scientific model state for a model run
  !>
  !> @todo We are in the middle of migrating all the functionality this class
  !>       provides into its parent. Once done this class will be removed.
  !>
  type, extends(driver_modeldb_type) :: modeldb_type

    private

    !> Stores all the fields used by the model
    type( model_data_type ), public :: model_data

    contains

  end type modeldb_type

  public modeldb_type

contains

end module gungho_modeldb_mod
