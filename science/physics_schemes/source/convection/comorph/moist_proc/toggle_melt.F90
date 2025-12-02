! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module toggle_melt_mod

implicit none

contains

! Subroutine to alter the coefficients for the implicit solution
! of T,q so-as to turn melting on or off
subroutine toggle_melt( n_points, nc, index_ic, l_switch_on,                   &
                        L_sub, L_fus, cp_tot,                                  &
                        coefs_melt, coefs_cond, coefs_cond_m,                  &
                        coefs_temp, coefs_q_vap, l_melt,                       &
                        coefs_cond_mlt, i_ice, l_trunc_mlt )

use comorph_constants_mod, only: real_cvprec, n_cond_species,                  &
                                 n_cond_species_liq
use phase_change_coefs_mod, only: n_coefs, i_q, i_t, i_0

implicit none

! Number of points
integer, intent(in) :: n_points

! List of points where melting wants switching on or off
integer, intent(in) :: nc
integer, intent(in) :: index_ic(nc)

! Flag for whether to turn melting on or off
logical, intent(in) :: l_switch_on

! Latent heats of sublimation and fusion
real(kind=real_cvprec), intent(in) :: L_sub(n_points)
real(kind=real_cvprec), intent(in) :: L_fus(n_points)
! Total heat capacity per unit dry-mass
real(kind=real_cvprec), intent(in) :: cp_tot(n_points)

! Coefficients for melting of each ice species
real(kind=real_cvprec), intent(in) :: coefs_melt                               &
                                         ( n_points, n_coefs )

! Coefficients for deposition / sublimation
real(kind=real_cvprec), intent(in out) :: coefs_cond                           &
                                         ( n_points, n_coefs )
! Coefficients for deposition rate if ice is melting
real(kind=real_cvprec), intent(in out) :: coefs_cond_m                         &
                                         ( n_points, n_coefs )
! Note: these are intent inout because their values get swapped
! at points where melting is switched on or off, such that
! coefs_cond always stores the coefficients actually in use

! Super-arrays containing coefficients in the implicit formula
! for temperature and q_vap after phase-changes
real(kind=real_cvprec), intent(in out) :: coefs_temp                           &
                                         ( n_points, n_coefs )
real(kind=real_cvprec), intent(in out) :: coefs_q_vap                          &
                                         ( n_points, n_coefs )

! Flag for whether melting is on or off
logical, intent(in out) :: l_melt(n_points)

! Condensation rate coefficients for liquid that this ice melts into
real(kind=real_cvprec), intent(in out) :: coefs_cond_mlt                       &
                                         ( n_points, n_coefs )

! Current ice species super-array address
integer, intent(in) :: i_ice

! Flag for ice species whose melting coefficients have been used to
! modify the coefficients for liquid species they melt into, so-as to
! yield exactly zero condensate mass after melt-source is accounted for.
logical, optional, intent(in out) :: l_trunc_mlt                               &
         ( n_points, n_cond_species_liq+1 : n_cond_species )

! Store for L_fus / cp_tot
real(kind=real_cvprec) :: lrcp

! Store for coefs during swap
real(kind=real_cvprec) :: coefs(n_coefs)

! Loop counters
integer :: ic, ic2


! If switching melting on
if ( l_switch_on ) then

  ! Add contribution of melting to the implicit
  ! coefficients for the T,q solve
  do ic2 = 1, nc
    ic = index_ic(ic2)
    lrcp = L_fus(ic) / cp_tot(ic)
    coefs_temp(ic,i_q) = coefs_temp(ic,i_q)                                    &
                       - coefs_melt(ic,i_q) * lrcp
    coefs_temp(ic,i_t) = coefs_temp(ic,i_t)                                    &
                       - coefs_melt(ic,i_t) * lrcp
    coefs_temp(ic,i_0) = coefs_temp(ic,i_0)                                    &
                       - coefs_melt(ic,i_0) * lrcp
  end do

  ! If switching melting off
else

  ! Subtract contribution of melting to the implicit
  ! coefficients for the T,q solve
  do ic2 = 1, nc
    ic = index_ic(ic2)
    lrcp = L_fus(ic) / cp_tot(ic)
    coefs_temp(ic,i_q) = coefs_temp(ic,i_q)                                    &
                       + coefs_melt(ic,i_q) * lrcp
    coefs_temp(ic,i_t) = coefs_temp(ic,i_t)                                    &
                       + coefs_melt(ic,i_t) * lrcp
    coefs_temp(ic,i_0) = coefs_temp(ic,i_0)                                    &
                       + coefs_melt(ic,i_0) * lrcp
  end do

end if

! Subtract contribution from existing deposition / sublimation
! coefficients
do ic2 = 1, nc
  ic = index_ic(ic2)
  coefs_q_vap(ic,i_q) = coefs_q_vap(ic,i_q)                                    &
                      - coefs_cond(ic,i_q)
  coefs_q_vap(ic,i_t) = coefs_q_vap(ic,i_t)                                    &
                      - coefs_cond(ic,i_t)
  coefs_q_vap(ic,i_0) = coefs_q_vap(ic,i_0)                                    &
                      - coefs_cond(ic,i_0)
  lrcp = L_sub(ic) / cp_tot(ic)
  coefs_temp(ic,i_q) = coefs_temp(ic,i_q)                                      &
                     - coefs_cond(ic,i_q) * lrcp
  coefs_temp(ic,i_t) = coefs_temp(ic,i_t)                                      &
                     - coefs_cond(ic,i_t) * lrcp
  coefs_temp(ic,i_0) = coefs_temp(ic,i_0)                                      &
                     - coefs_cond(ic,i_0) * lrcp
end do

! Swap to use the alternative deposition / sublimation
! coefficients
do ic2 = 1, nc
  ic = index_ic(ic2)
  coefs(i_q) = coefs_cond(ic,i_q)
  coefs(i_t) = coefs_cond(ic,i_t)
  coefs(i_0) = coefs_cond(ic,i_0)
  coefs_cond(ic,i_q) = coefs_cond_m(ic,i_q)
  coefs_cond(ic,i_t) = coefs_cond_m(ic,i_t)
  coefs_cond(ic,i_0) = coefs_cond_m(ic,i_0)
  coefs_cond_m(ic,i_q) = coefs(i_q)
  coefs_cond_m(ic,i_t) = coefs(i_t)
  coefs_cond_m(ic,i_0) = coefs(i_0)
end do

! Add contribution from alternative deposition / sublimation
! coefficients
do ic2 = 1, nc
  ic = index_ic(ic2)
  coefs_q_vap(ic,i_q) = coefs_q_vap(ic,i_q)                                    &
                      + coefs_cond(ic,i_q)
  coefs_q_vap(ic,i_t) = coefs_q_vap(ic,i_t)                                    &
                      + coefs_cond(ic,i_t)
  coefs_q_vap(ic,i_0) = coefs_q_vap(ic,i_0)                                    &
                      + coefs_cond(ic,i_0)
  lrcp = L_sub(ic) / cp_tot(ic)
  coefs_temp(ic,i_q) = coefs_temp(ic,i_q)                                      &
                     + coefs_cond(ic,i_q) * lrcp
  coefs_temp(ic,i_t) = coefs_temp(ic,i_t)                                      &
                     + coefs_cond(ic,i_t) * lrcp
  coefs_temp(ic,i_0) = coefs_temp(ic,i_0)                                      &
                     + coefs_cond(ic,i_0) * lrcp
end do

! Set flag
do ic2 = 1, nc
  ic = index_ic(ic2)
  l_melt(ic) = l_switch_on
end do

! Turning off melting when the liquid species receiving the melt-water
! has had its process-rates truncated to avoid creating negative condensate...
if ( present(l_trunc_mlt) .and. (.not. l_switch_on) ) then
  do ic2 = 1, nc
    ic = index_ic(ic2)
    if ( l_trunc_mlt(ic,i_ice) ) then

      ! Add the melting coefficients back onto the liquid species'
      ! condensation rate coefficients
      ! (undo the subtraction in modify_coefs_liq)
      coefs_cond_mlt(ic,i_q) = coefs_cond_mlt(ic,i_q) + coefs_melt(ic,i_q)
      coefs_cond_mlt(ic,i_t) = coefs_cond_mlt(ic,i_t) + coefs_melt(ic,i_t)
      coefs_cond_mlt(ic,i_0) = coefs_cond_mlt(ic,i_0) + coefs_melt(ic,i_0)

      ! Correct the coefficients for the T,q solve consistently
      coefs_q_vap(ic,i_q) = coefs_q_vap(ic,i_q) + coefs_melt(ic,i_q)
      coefs_q_vap(ic,i_t) = coefs_q_vap(ic,i_t) + coefs_melt(ic,i_t)
      coefs_q_vap(ic,i_0) = coefs_q_vap(ic,i_0) + coefs_melt(ic,i_0)
      lrcp = (L_sub(ic)-L_fus(ic)) / cp_tot(ic)
      coefs_temp(ic,i_q) = coefs_temp(ic,i_q) + coefs_melt(ic,i_q) * lrcp
      coefs_temp(ic,i_t) = coefs_temp(ic,i_t) + coefs_melt(ic,i_t) * lrcp
      coefs_temp(ic,i_0) = coefs_temp(ic,i_0) + coefs_melt(ic,i_0) * lrcp

      ! Reset flag to false
      l_trunc_mlt(ic,i_ice) = .false.

    end if
  end do
end if

return
end subroutine toggle_melt

end module toggle_melt_mod
