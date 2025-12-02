!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Provides functionality for analytic orography related classes.
!>
!> @details Contains abstract orography type definition and abstract interface
!>          which set up templates to calculate analytic orography profiles.
!>          This functionality is later used to assign analytic orography
!>          profiles for Schar and Witch-of-Agnesi mountains in spherical and
!>          Cartesian coordinates.
!-------------------------------------------------------------------------------
module analytic_orography_mod

  use constants_mod, only : r_def

  implicit none

  private

  !-----------------------------------------------------------------------------
  !> @brief Abstract type definition for analytic orography profiles.
  !-----------------------------------------------------------------------------
  type, public, abstract :: analytic_orography_type

    private

  contains

    procedure(analytic_orography_function), deferred, pass(self) :: analytic_orography

  end type analytic_orography_type

  !-----------------------------------------------------------------------------
  !> @brief Abstract interface definition for function which calculates surface
  !>        heights in a selected analytic orography profile.
  !-----------------------------------------------------------------------------
  abstract interface

    function analytic_orography_function(self, chi_1, chi_2) result(chi_surf)

      import :: analytic_orography_type
      import :: r_def

      implicit none

      class(analytic_orography_type), intent(in) :: self
      real(kind=r_def),               intent(in) :: chi_1, chi_2
      real(kind=r_def)                           :: chi_surf

    end function analytic_orography_function

  end interface

  !-----------------------------------------------------------------------------
  !> @brief Stores functionality for a selected analytic orography profile.
  !-----------------------------------------------------------------------------
  class(analytic_orography_type), allocatable, public :: orography_profile

end module analytic_orography_mod


