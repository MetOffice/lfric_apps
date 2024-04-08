!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Holds the tangent linear version of the transport_runtime object
module tl_transport_runtime_collection_mod

  use transport_runtime_alg_mod, only: transport_runtime_type

  implicit none

  private

  ! TODO #3008: remove this single runtime and replace with runtime_collection
  type(transport_runtime_type), public :: tl_transport_runtime

end module tl_transport_runtime_collection_mod
