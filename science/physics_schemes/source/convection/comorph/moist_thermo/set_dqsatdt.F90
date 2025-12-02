! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module set_dqsatdt_mod

implicit none

contains

! Estimates the gradient of the saturation water vapour mixing
! ratio curve as a function of temperature T,
! by assuming d/dT of the saturation vapour pressure follows
! the Claussius-Clapeyron equation

!----------------------------------------------------------------
! Routine for liquid at all temperatures
!----------------------------------------------------------------
subroutine set_dqsatdt_liq( n_points, temperature, qsat,                       &
                            dqsatdt )

use comorph_constants_mod, only: R_dry, R_vap, real_cvprec, one
use lat_heat_mod, only: set_l_con

implicit none

! Number of points
integer, intent(in) :: n_points

! Temperature
real(kind=real_cvprec), intent(in) :: temperature(n_points)

! Already calculated saturation vapour mixing ratio
! w.r.t. liquid water
real(kind=real_cvprec), intent(in) :: qsat(n_points)

! Partial d/dT of qsat
real(kind=real_cvprec), intent(out) :: dqsatdt(n_points)

! Latent heat of condensation
real(kind=real_cvprec) :: L_con(n_points)

! Loop counter
integer :: ic

! Compute latent heat of condensation
call set_l_con( n_points, temperature, L_con )

! Compute dqsat/dT:
do ic = 1, n_points
  dqsatdt(ic) = qsat(ic) * ( one + (R_vap/R_dry) * qsat(ic) )                  &
                         * L_con(ic)                                           &
              / ( R_vap * temperature(ic) * temperature(ic) )
end do

return
end subroutine set_dqsatdt_liq


!----------------------------------------------------------------
! Routine for ice at all temperatures
!----------------------------------------------------------------
subroutine set_dqsatdt_ice( n_points, temperature, qsat,                       &
                            dqsatdt )

use comorph_constants_mod, only: R_dry, R_vap, real_cvprec, one
use lat_heat_mod, only: set_l_sub

implicit none

! Number of points
integer, intent(in) :: n_points

! Temperature
real(kind=real_cvprec), intent(in) :: temperature(n_points)

! Already calculated saturation vapour mixing ratio
! w.r.t. ice
real(kind=real_cvprec), intent(in) :: qsat(n_points)

! Partial d/dT of qsat
real(kind=real_cvprec), intent(out) :: dqsatdt(n_points)

! Latent heat of sublimation
real(kind=real_cvprec) :: L_sub(n_points)

! Loop counter
integer :: ic

! Compute latent heat of condensation and multiply by it
call set_l_sub( n_points, temperature, L_sub )

! Compute dqsat/dT:
do ic = 1, n_points
  dqsatdt(ic) = qsat(ic) * ( one + (R_vap/R_dry) * qsat(ic) )                  &
                         * L_sub(ic)                                           &
              / ( R_vap * temperature(ic) * temperature(ic) )
end do

return
end subroutine set_dqsatdt_ice


end module set_dqsatdt_mod
