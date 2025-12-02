!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Collection of routines that are needed for physics.

!> @detail Collection of routines that are needed for physics. These may be
!>         replaced/overloaded by specific schemes as they are brought in from
!>         the UM

module physics_common_mod

  use constants_mod,                 only: r_def
  implicit none
  private

  public qsaturation

contains

  ! Function to return the saturation mr over water
  ! Based on Tetens' formula
  ! QS=3.8/(P*EXP(-17.2693882*(T-273.15)/(T-35.86))-6.109)
  function qsaturation (T, p)
    implicit none
    real(kind=r_def), intent(in) :: T     ! Temperature in Kelvin
    real(kind=r_def), intent(in) :: p     ! Pressure in mb

    real(kind=r_def)             :: Qsaturation

    real(kind=r_def),  parameter :: tk0c = 273.15_r_def      ! Temperature of freezing in Kelvin
    real(kind=r_def),  parameter :: qsa1 = 3.8_r_def         ! Top constant in qsat equation
    real(kind=r_def),  parameter :: qsa2 = -17.2693882_r_def ! Constant in qsat equation in Kelvin
    real(kind=r_def),  parameter :: qsa3 = 35.86_r_def      ! Constant in qsat equation
    real(kind=r_def),  parameter :: qsa4 = 6.109_r_def       ! Constant in qsat equation in mbar

    if (T > qsa3 .and. p * exp (qsa2 * (t - tk0c) / (T - qsa3)) > qsa4) then
      qsaturation=qsa1/(p*exp(qsa2*(t-tk0c)/(T-qsa3))-qsa4)
    else
      qsaturation=999.0_r_def
    end if
  end function qsaturation

end module physics_common_mod
