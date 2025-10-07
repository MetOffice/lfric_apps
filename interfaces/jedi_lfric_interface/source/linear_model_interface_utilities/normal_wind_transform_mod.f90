!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Basic variable transformation between the JEDI analysis
!>        wind variables and the LFRic prognostic wind variables.

module normal_wind_transform_mod

  use base_wind_transform_mod, only: base_wind_transform_type
  use field_collection_mod,    only: field_collection_type
  use field_mod,               only: field_type

  implicit none

  private

  type, public, extends(base_wind_transform_type) :: normal_wind_transform_type

    private

    contains

    procedure, public :: initialise

    procedure, public :: process
    procedure, public :: initialise_for_adjoint
    procedure, public :: adj_process

    procedure, public :: scalar_to_vector
    procedure, public :: adj_scalar_to_vector
    procedure, public :: vector_to_scalar
    procedure, public :: adj_vector_to_scalar

  end type normal_wind_transform_type

  contains

  !> @brief Do-nothing initialiser.
  !> @param[in] fields Field collection containing wind fields.
  subroutine initialise( self, fields )
    implicit none
    class(normal_wind_transform_type), intent(inout) :: self
    type(field_collection_type),       intent(in)    :: fields
  end subroutine initialise

  !> @brief Do-nothing process routine.
  !> @param[in] fields Field collection containing wind fields.
  subroutine process( self, fields )
    implicit none
    class(normal_wind_transform_type), intent(inout) :: self
    type(field_collection_type),       intent(in)    :: fields
  end subroutine process

  !> @brief Do-nothing adjoint initialiser.
  subroutine initialise_for_adjoint(self)
    implicit none
    class(normal_wind_transform_type), intent(inout) :: self
  end subroutine initialise_for_adjoint

  !> @brief Do-nothing adjoint process routine.
  !> @param[in] fields Field collection containing wind fields.
  subroutine adj_process( self, fields )
    implicit none
    class(normal_wind_transform_type), intent(inout) :: self
    type(field_collection_type),       intent(in)    :: fields
  end subroutine adj_process

  !> @brief Transform JEDI analysis wind variables to LFRic prognostic wind variables.
  !> @param[in] fields Field collection containing wind fields.
  subroutine scalar_to_vector( self, fields )

    use interpolation_alg_mod, only: interp_w3wth_to_w2_alg

    implicit none

    class(normal_wind_transform_type), intent(inout) :: self
    type(field_collection_type),       intent(in)    :: fields

    type(field_type), pointer :: u_in_w3
    type(field_type), pointer :: v_in_w3
    type(field_type), pointer :: w_in_wth
    type(field_type), pointer :: u_in_w2

    call fields%get_field( 'u_in_w3', u_in_w3 )
    call fields%get_field( 'v_in_w3', v_in_w3 )
    call fields%get_field( 'w_in_wth', w_in_wth )
    call fields%get_field( 'u', u_in_w2 )

    call interp_w3wth_to_w2_alg( u_in_w2, u_in_w3, v_in_w3, w_in_wth )

  end subroutine scalar_to_vector

  !> @brief (Adjoint of) transform JEDI analysis wind variables to LFRic prognostic wind variables.
  !> @param[in] fields Field collection containing wind fields.
  subroutine adj_scalar_to_vector( self, fields )

    use adj_interpolation_alg_mod, only: adj_interp_w3wth_to_w2_alg

    implicit none

    class(normal_wind_transform_type), intent(inout) :: self
    type(field_collection_type),       intent(in)    :: fields

    type(field_type), pointer :: u_in_w3
    type(field_type), pointer :: v_in_w3
    type(field_type), pointer :: w_in_wth
    type(field_type), pointer :: u_in_w2

    call fields%get_field( 'u_in_w3', u_in_w3 )
    call fields%get_field( 'v_in_w3', v_in_w3 )
    call fields%get_field( 'w_in_wth', w_in_wth )
    call fields%get_field( 'u', u_in_w2 )

    call adj_interp_w3wth_to_w2_alg( u_in_w2, u_in_w3, v_in_w3, w_in_wth )

  end subroutine adj_scalar_to_vector

  !> @brief Transform LFRic prognostic wind variables to JEDI analysis wind variables.
  !> @param[in] fields Field collection containing wind fields.
  subroutine vector_to_scalar( self, fields )

    use interpolation_alg_mod, only: interp_w2_to_w3wth_alg

    implicit none

    class(normal_wind_transform_type), intent(inout) :: self
    type(field_collection_type),       intent(in)    :: fields

    type(field_type), pointer :: u_in_w3
    type(field_type), pointer :: v_in_w3
    type(field_type), pointer :: w_in_wth
    type(field_type), pointer :: u_in_w2

    call fields%get_field( 'u_in_w3', u_in_w3 )
    call fields%get_field( 'v_in_w3', v_in_w3 )
    call fields%get_field( 'w_in_wth', w_in_wth )
    call fields%get_field( 'u', u_in_w2 )

    call interp_w2_to_w3wth_alg( u_in_w2, u_in_w3, v_in_w3, w_in_wth )

  end subroutine vector_to_scalar

  !> @brief (Adjoint of) transform LFRic prognostic wind variables to JEDI analysis wind variables.
  !> @param[in] fields Field collection containing wind fields.
  subroutine adj_vector_to_scalar( self, fields )

    use adj_interpolation_alg_mod, only: adj_interp_w2_to_w3wth_alg

    implicit none

    class(normal_wind_transform_type), intent(inout) :: self
    type(field_collection_type),       intent(in)    :: fields

    type(field_type), pointer :: u_in_w3
    type(field_type), pointer :: v_in_w3
    type(field_type), pointer :: w_in_wth
    type(field_type), pointer :: u_in_w2

    call fields%get_field( 'u_in_w3', u_in_w3 )
    call fields%get_field( 'v_in_w3', v_in_w3 )
    call fields%get_field( 'w_in_wth', w_in_wth )
    call fields%get_field( 'u', u_in_w2 )

    call adj_interp_w2_to_w3wth_alg( u_in_w2, u_in_w3, v_in_w3, w_in_wth )

  end subroutine adj_vector_to_scalar

end module normal_wind_transform_mod
