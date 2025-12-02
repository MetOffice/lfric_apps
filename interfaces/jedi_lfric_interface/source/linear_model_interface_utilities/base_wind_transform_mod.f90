!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Abstract base class for JEDI <-> LFRic wind transformations.
!> @details We have two options for variable transformation between the JEDI analysis wind variables and the LFRic
!>          prognostic wind variables. Two derived types of this abstract type encapsulate this subtlety.

module base_wind_transform_mod

  use field_collection_mod, only: field_collection_type

  implicit none

  private

  type, public, abstract :: base_wind_transform_type

    private

    contains

    procedure(scalar_to_vector_interface),       public, deferred :: initialise

    procedure(process_interface),                public, deferred :: process
    procedure(initialise_for_adjoint_interface), public, deferred :: initialise_for_adjoint
    procedure(adj_process_interface),            public, deferred :: adj_process

    procedure(scalar_to_vector_interface),       public, deferred :: scalar_to_vector
    procedure(adj_scalar_to_vector_interface),   public, deferred :: adj_scalar_to_vector
    procedure(vector_to_scalar_interface),       public, deferred :: vector_to_scalar
    procedure(adj_vector_to_scalar_interface),   public, deferred :: adj_vector_to_scalar

  end type base_wind_transform_type

  abstract interface

    subroutine initialise_interface( self, fields )
      import base_wind_transform_type, field_collection_type
      class(base_wind_transform_type), intent(inout) :: self
      type(field_collection_type),     intent(in)    :: fields
    end subroutine initialise_interface

    subroutine process_interface( self, fields )
      import base_wind_transform_type, field_collection_type
      implicit none
      class(base_wind_transform_type), intent(inout) :: self
      type(field_collection_type),     intent(in)    :: fields
    end subroutine process_interface

    subroutine initialise_for_adjoint_interface(self)
      import base_wind_transform_type
      class(base_wind_transform_type), intent(inout) :: self
    end subroutine initialise_for_adjoint_interface

    subroutine adj_process_interface( self, fields )
      import base_wind_transform_type, field_collection_type
      class(base_wind_transform_type), intent(inout) :: self
      type(field_collection_type),     intent(in)    :: fields
    end subroutine adj_process_interface

    subroutine scalar_to_vector_interface( self, fields )
      import base_wind_transform_type, field_collection_type
      class(base_wind_transform_type), intent(inout) :: self
      type(field_collection_type),     intent(in)    :: fields
    end subroutine scalar_to_vector_interface

    subroutine adj_scalar_to_vector_interface( self, fields )
      import base_wind_transform_type, field_collection_type
      class(base_wind_transform_type), intent(inout) :: self
      type(field_collection_type),     intent(in)    :: fields
    end subroutine adj_scalar_to_vector_interface

    subroutine vector_to_scalar_interface( self, fields )
      import base_wind_transform_type, field_collection_type
      class(base_wind_transform_type), intent(inout) :: self
      type(field_collection_type),     intent(in)    :: fields
    end subroutine vector_to_scalar_interface

    subroutine adj_vector_to_scalar_interface( self, fields )
      import base_wind_transform_type, field_collection_type
      class(base_wind_transform_type), intent(inout) :: self
      type(field_collection_type),     intent(in)    :: fields
    end subroutine adj_vector_to_scalar_interface

  end interface

end module base_wind_transform_mod
