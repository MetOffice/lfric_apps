! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module sat_adjust_mod

implicit none

contains

! Subroutine to adjust the temperature, vapour and liquid cloud
! so-as to obtain liquid saturation.
! Simple method based on linearising q_sat about the initial T
subroutine sat_adjust( n_points, l_update_q, pressure,                         &
                       temperature, q_vap, q_cl, cp_tot )

use comorph_constants_mod, only: real_cvprec, one
use lat_heat_mod, only: lat_heat_incr, set_l_con,                              &
                        i_phase_change_con
use set_qsat_mod, only: set_qsat_liq
use set_dqsatdt_mod, only: set_dqsatdt_liq

implicit none

! Number of points
integer, intent(in) :: n_points

! Flag for whether to update the vapour and liquid-cloud
! mixing-ratios (this routine can be called to just
! compute the temperature after adjustment to saturation,
! without changing q_vap or q_cl)
logical, intent(in) :: l_update_q

! Pressure
real(kind=real_cvprec), intent(in) :: pressure(n_points)

! Temperature, water-vapour and liquid-cloud mixing ratios
! to be updated
real(kind=real_cvprec), intent(in out) :: temperature(n_points)
real(kind=real_cvprec), intent(in out) :: q_vap(n_points)
real(kind=real_cvprec), intent(in out) :: q_cl(n_points)

! Parcel total heat capacity (also gets updated)
real(kind=real_cvprec), intent(in out) ::cp_tot(n_points)

! Saturation mixing-ratio at the input temperature
real(kind=real_cvprec) :: qsat(n_points)

! Rate of change of qsat with temperature
real(kind=real_cvprec) :: dqsatdt(n_points)

! Amount of condensation required (negative for evaporation)
real(kind=real_cvprec) :: dq(n_points)

! Latent heat of condensation (including temperature dependence)
real(kind=real_cvprec) :: l_con(n_points)

! Loop counter
integer :: ic


! Calculate saturation mixing ratio at the initial temperature
call set_qsat_liq( n_points, temperature, pressure, qsat )

! Calculate rate of change of qsat with T for linearisation
call set_dqsatdt_liq( n_points, temperature, qsat, dqsatdt )

! Set latent heat of condensation as a function of temperature
call set_l_con( n_points, temperature, l_con )

! Estimate amount of condensation dq required to
! reach saturation:
! We have:
! Latent heating:
!     dT = (Lc/cp) dq                        (5)
! Condensation to adjust q_vap to qsat:
!     dq = q_vap - q_sat( T + dT )
! Linearising q_sat:
!     dq = q_vap - q_sat(T) - dqsat/dT dT
! Substituting (5) for dT and rearranging:
!     dq( 1 + dqsat/dT (Lc/cp) ) = q_vap - q_sat(T)
! Therefore:
!     dq = ( q_vap - q_sat(T) )
!        / ( 1 + dqsat/dT (Lc/cp) )
do ic = 1, n_points
  dq(ic) = ( q_vap(ic) - qsat(ic) )                                            &
         / ( one + dqsatdt(ic) * (l_con(ic)/cp_tot(ic)) )
end do

! If doing evaporation to reach saturation, don't allow
! removing more liquid cloud than is present
do ic = 1, n_points
  dq(ic) = max( dq(ic), -q_cl(ic) )
end do

! Do accurate updating of temperature and cp_tot
call lat_heat_incr( n_points, n_points, i_phase_change_con,                    &
                    cp_tot, temperature, dq=dq )

if ( l_update_q ) then
  ! Update vapour and liquid cloud
  do ic = 1, n_points
    q_vap(ic) = q_vap(ic) - dq(ic)
    q_cl(ic)  = q_cl(ic)  + dq(ic)
  end do
end if


return
end subroutine sat_adjust

end module sat_adjust_mod
