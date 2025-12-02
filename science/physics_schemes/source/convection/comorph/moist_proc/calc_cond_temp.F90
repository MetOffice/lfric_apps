! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_cond_temp_mod

implicit none

contains

! Subroutine to calculate the temperature of hydrometeors.
! Method is to back-track out the coefficients for the
! hydrometeor temperature formula from those for the vapour
! exchange formula.

subroutine calc_cond_temp( n_points, nc, index_ic, l_full_do, l_ice,           &
                           qsat_liq_ref, qsat_ice_ref, dqsatdt,                &
                           kq_cond, imp_temp, imp_q_vap,                       &
                           coefs_cond, cond_temp )

use comorph_constants_mod, only: real_cvprec, zero
use phase_change_coefs_mod, only: n_coefs, i_q, i_t, i_0

implicit none

! Number of points
integer, intent(in) :: n_points

! Points where calculation needs to be done
integer, intent(in) :: nc
integer, intent(in) :: index_ic(nc)
! Flag for whether to do full-field do-loops instead of indirect indexing
logical, intent(in) :: l_full_do

! Flag for whether this is an ice species
logical, intent(in) :: l_ice

! Saturation water vapour mixing ratio qsat and dqsat/dT at the
! reference temperature, w.r.t. liquid water and ice.
real(kind=real_cvprec), intent(in) :: qsat_liq_ref(n_points)
real(kind=real_cvprec), intent(in) :: qsat_ice_ref(n_points)
real(kind=real_cvprec), intent(in) :: dqsatdT(n_points)

! Vapour exchange coefficient
real(kind=real_cvprec), intent(in) :: kq_cond(n_points)

! Estimated temperature and water vapour mixing ratio after
! implicit solution of phase-changes (minus reference values)
real(kind=real_cvprec), intent(in) :: imp_temp(n_points)
real(kind=real_cvprec), intent(in) :: imp_q_vap(n_points)

! Coefficients for condensation / evaporation
real(kind=real_cvprec), intent(in) :: coefs_cond                               &
                                      ( n_points, n_coefs )

! Diagnosed temperature of the hydrometeors
! (minus the reference temperature)
real(kind=real_cvprec), intent(in out) :: cond_temp(n_points)

! Loop counters
integer :: ic, ic2


! If calculations required at the majority of points
if ( l_full_do ) then
  ! Full-field calculation

  ! Note: both the coefs and kq_cond include factors of delta_t,
  ! and these cancel out below:

  do ic = 1, n_points
    cond_temp(ic) = (                                                          &
        ( kq_cond(ic) - coefs_cond(ic,i_q) ) * imp_q_vap(ic)                   &
      - coefs_cond(ic,i_t) * imp_temp(ic)                                      &
      - coefs_cond(ic,i_0)                                                     &
                    ) / dqsatdT(ic)
  end do

  ! If this is an ice species, add on extra term
  if ( l_ice ) then
    do ic = 1, n_points
      cond_temp(ic) = cond_temp(ic)                                            &
        - kq_cond(ic) * ( qsat_ice_ref(ic) - qsat_liq_ref(ic) ) / dqsatdT(ic)
    end do
  end if

  ! Final factor of 1/kq_cond (except where no condensate present; kq_cond=0)
  do ic = 1, n_points
    if ( kq_cond(ic) > zero ) then
      cond_temp(ic) = cond_temp(ic) / kq_cond(ic)
    else
      cond_temp(ic) = zero
    end if
  end do

  ! Calculations only needed at a minority of points
else
  ! Indirect indexing version of the above calculation

  do ic2 = 1, nc
    ic = index_ic(ic2)
    cond_temp(ic) = (                                                          &
        ( kq_cond(ic) - coefs_cond(ic,i_q) ) * imp_q_vap(ic)                   &
      - coefs_cond(ic,i_t) * imp_temp(ic)                                      &
      - coefs_cond(ic,i_0)                                                     &
                    ) / dqsatdT(ic)
  end do

  if ( l_ice ) then
    do ic2 = 1, nc
      ic = index_ic(ic2)
      cond_temp(ic) = cond_temp(ic)                                            &
        - kq_cond(ic) * ( qsat_ice_ref(ic) - qsat_liq_ref(ic) ) / dqsatdT(ic)
    end do
  end if

  do ic2 = 1, nc
    ic = index_ic(ic2)
    cond_temp(ic) = cond_temp(ic) / kq_cond(ic)
  end do

end if


return
end subroutine calc_cond_temp

end module calc_cond_temp_mod
