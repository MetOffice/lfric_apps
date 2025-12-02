!----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief Controls the setting of UM clock

module um_clock_init_mod

  ! LFRic modules
  use clock_mod,     only : clock_type
  use constants_mod, only : r_um

  implicit none

  private
  public :: um_clock_init

contains

  !> @brief Initialise UM clock
  !> @param[in] clock  The LFRic model clock object
  subroutine um_clock_init( clock )

    ! UM modules containing things that need setting
    use timestep_mod, only: timestep, recip_timestep

    implicit none

    class(clock_type), intent(in) :: clock

    ! Timestep used in UM code - contained in UM timestep_mod.
    ! Set from LFRic input timestep.
    timestep = real(clock%get_seconds_per_step(), r_um)
    recip_timestep = 1.0_r_um / timestep

  end subroutine um_clock_init

end module um_clock_init_mod
