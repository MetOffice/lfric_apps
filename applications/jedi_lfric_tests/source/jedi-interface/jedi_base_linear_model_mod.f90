!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing a base class than handles the linear model (tlm)
!>
!> @details This module includes a base class that defines the forecast methods
!>          for the Tangent Linear (TL) and Adjoint (AD) models. These methods
!>          rely on deferred interface methods: init, step and final. This
!>          enables the possibility to instantiate multiple extended linear
!>          model instances. An additional deferred method "set_trajectory" is
!>          included. This is used via the "traj" post processor to enable
!>          creation of the linear state trajectory used by the linear model.
!>
module jedi_base_linear_model_mod

  use jedi_geometry_mod,       only : jedi_geometry_type
  use jedi_increment_mod,      only : jedi_increment_type
  use jedi_lfric_datetime_mod, only : jedi_datetime_type
  use jedi_lfric_duration_mod, only : jedi_duration_type
  use jedi_state_mod,          only : jedi_state_type

  implicit none

  private

type, public, abstract :: jedi_base_linear_model_type
  private
contains
  private

  procedure(set_trajectory_interface), public, deferred :: set_trajectory

  procedure(model_initTL_interface),   public, deferred :: model_initTL
  procedure(model_stepTL_interface),   public, deferred :: model_stepTL
  procedure(model_finalTL_interface),  public, deferred :: model_finalTL

  procedure(model_initAD_interface),   public, deferred :: model_initAD
  procedure(model_stepAD_interface),   public, deferred :: model_stepAD
  procedure(model_finalAD_interface),  public, deferred :: model_finalAD

  ! Forecast methods
  procedure, public :: forecastTL
  procedure, public :: forecastAD

end type jedi_base_linear_model_type

abstract interface
  subroutine set_trajectory_interface( self, jedi_state )
    import jedi_base_linear_model_type, jedi_state_type
    implicit none
    class( jedi_base_linear_model_type ), intent(inout) :: self
    type( jedi_state_type ),              intent(inout) :: jedi_state
  end subroutine set_trajectory_interface

  subroutine model_initTL_interface(self, increment)
    import jedi_base_linear_model_type, jedi_increment_type
    class( jedi_base_linear_model_type ), intent(inout) :: self
    type( jedi_increment_type ),          intent(inout) :: increment
  end subroutine model_initTL_interface

  subroutine model_stepTL_interface(self, increment)
    import jedi_base_linear_model_type, jedi_increment_type
    class( jedi_base_linear_model_type ), target, intent(inout) :: self
    type( jedi_increment_type ),                  intent(inout) :: increment
  end subroutine model_stepTL_interface

  subroutine model_finalTL_interface(self, increment)
    import jedi_base_linear_model_type, jedi_increment_type
    class( jedi_base_linear_model_type ), intent(inout) :: self
    type( jedi_increment_type ),          intent(inout) :: increment
  end subroutine model_finalTL_interface

  subroutine model_initAD_interface(self, increment)
    import jedi_base_linear_model_type, jedi_increment_type
    class( jedi_base_linear_model_type ), intent(inout) :: self
    type( jedi_increment_type ),          intent(inout) :: increment
  end subroutine model_initAD_interface

  subroutine model_stepAD_interface(self, increment)
    import jedi_base_linear_model_type, jedi_increment_type
    class( jedi_base_linear_model_type ), target, intent(inout) :: self
    type( jedi_increment_type ),                  intent(inout) :: increment
  end subroutine model_stepAD_interface

  subroutine model_finalAD_interface(self, increment)
    import jedi_base_linear_model_type, jedi_increment_type
    class( jedi_base_linear_model_type ), intent(inout) :: self
    type( jedi_increment_type ),          intent(inout) :: increment
  end subroutine model_finalAD_interface
end interface

!------------------------------------------------------------------------------
! Contained functions/subroutines
!------------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
! OOPS defined forecastTL and forecastAD methods
!------------------------------------------------------------------------------

!> @brief    Run a Tangent Linear forecast using the model init, step and final
!>
!> @param [inout] increment       The Increment object to propagate
!> @param [inout] forecast_length The duration of the forecastTL
subroutine forecastTL( self, increment, forecast_length )

  implicit none

  class( jedi_base_linear_model_type ), intent(inout) :: self
  type( jedi_increment_type ),          intent(inout) :: increment
  type( jedi_duration_type ),              intent(in) :: forecast_length

  ! Local
  type( jedi_datetime_type ) :: end_time

  ! End time
  end_time = increment%valid_time() + forecast_length

  ! Initialize the model
  call self%model_initTL( increment )
  ! Initialize the post processor and call first process
  ! call post_processor%pp_init( increment )
  ! call post_processor%process( increment )

  ! Loop until end_time
  do while ( end_time > increment%valid_time() )
    call self%model_stepTL( increment )
    ! call post_processor%process( increment )
  end do

  ! Finalize model and post processor
  ! call post_processor%pp_final( increment )
  call self%model_finalTL( increment )

end subroutine forecastTL

!> @brief    Run a Adjoint of the Tangent Linear forecast using the model init,
!>           step and final
!>
!> @param [inout] increment       The Increment object to propagate
!> @param [in]    forecast_length The duration of the forecastAD
subroutine forecastAD( self, increment, forecast_length )

  implicit none

  class( jedi_base_linear_model_type ), intent(inout) :: self
  type( jedi_increment_type ),          intent(inout) :: increment
  type( jedi_duration_type ),              intent(in) :: forecast_length

  ! Local
  type( jedi_datetime_type ) :: begin_time

  ! Begin time
  begin_time = increment%valid_time() + forecast_length*(-1)

  ! Initialize the model
  call self%model_initAD( increment )
  ! Initialize the post processor and call first process
  ! call post_processor%pp_init( increment )
  ! call post_processor%process( increment )

  ! Loop until begin_time
  do while ( begin_time < increment%valid_time() )
    call self%model_stepAD( increment )
    ! call post_processor%process( increment )
  end do

  ! Finalize model and post processor
  ! call post_processor%pp_final( increment )
  call self%model_finalAD( increment )

end subroutine forecastAD

end module jedi_base_linear_model_mod
