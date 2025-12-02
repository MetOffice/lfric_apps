! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module modify_coefs_liq_mod

implicit none

contains


! Subroutine to compute the condensation/evaporation rate for
! liquid species, check for negative mixing ratio, and correct
! the evaporation rate and implicit solve coefficients to avoid
! any found negative values.

subroutine modify_coefs_liq( n_points, nc_cor, index_ic_cor,                   &
                             l_any_melt, i_liq, L_con, cp_tot,                 &
                             q_cond, l_melt, l_trunc_mlt,                      &
                             coefs_melt, coefs_cond,                           &
                             coefs_temp, coefs_q_vap )

use comorph_constants_mod, only: real_cvprec, zero, n_cond_species,            &
                     n_cond_species_liq, cond_params
use phase_change_coefs_mod, only: n_coefs, i_q, i_t, i_0

implicit none

! Number of points
integer, intent(in) :: n_points

! Number of points where coefficients need to be modified,
! and their indices
integer, intent(in) :: nc_cor
integer, intent(in) :: index_ic_cor(nc_cor)

! Flag for whether any melt source has been found for the
! current liquid species
logical, intent(in) :: l_any_melt

! Cond super-array address of the current liquid species
integer, intent(in) :: i_liq

! Latent heat of condensation
real(kind=real_cvprec), intent(in) :: L_con(n_points)

! Total heat capacity per unit dry-mass
real(kind=real_cvprec), intent(in) :: cp_tot(n_points)

! Mixing ratio of liquid species before implicit phase-changes
real(kind=real_cvprec), intent(in) :: q_cond(n_points)

! Melting flag for all the ice species
logical, intent(in) :: l_melt                                                  &
           ( n_points, n_cond_species_liq+1 : n_cond_species )

! Flag for ice species whose melting coefficients have been used to
! modify the coefficients for liquid species they melt into, so-as to
! yield exactly zero condensate mass after melt-source is accounted for.
logical, intent(in out) :: l_trunc_mlt                                         &
           ( n_points, n_cond_species_liq+1 : n_cond_species )


! Melting coefficients for all the ice species
real(kind=real_cvprec), intent(in) :: coefs_melt                               &
  ( n_points, n_coefs, n_cond_species_liq+1 : n_cond_species )

! Coefficients for condensation / evaporation
! (corrected by this routine to avoid negative q_cond)
real(kind=real_cvprec), intent(in out) :: coefs_cond                           &
                                         ( n_points, n_coefs )

! Super-arrays containing coefficients in the implicit formula
! for temperature and q_vap after phase-changes
real(kind=real_cvprec), intent(in out) :: coefs_temp                           &
                                         ( n_points, n_coefs )
real(kind=real_cvprec), intent(in out) :: coefs_q_vap                          &
                                         ( n_points, n_coefs )

! Sum of melting coefficients over all ice species which melt
! into the current liquid species
real(kind=real_cvprec) :: coefs_melt_source( nc_cor, n_coefs )

! Points where melting is occuring
integer :: nc_melt
integer :: index_ic_melt(nc_cor)

! Store for L/cp
real(kind=real_cvprec) :: lrcp

! loop counters
integer :: ic, ic2, ic3, i_ice


! Subtract the condensation contributions for
! this species from the implicit solve coefficients for T,q
do ic2 = 1, nc_cor
  ic = index_ic_cor(ic2)
  coefs_q_vap(ic,i_q) = coefs_q_vap(ic,i_q) - coefs_cond(ic,i_q)
  coefs_q_vap(ic,i_t) = coefs_q_vap(ic,i_t) - coefs_cond(ic,i_t)
  coefs_q_vap(ic,i_0) = coefs_q_vap(ic,i_0) - coefs_cond(ic,i_0)
  lrcp =  L_con(ic) / cp_tot(ic)
  coefs_temp(ic,i_q) = coefs_temp(ic,i_q) - coefs_cond(ic,i_q)                 &
                                            * lrcp
  coefs_temp(ic,i_t) = coefs_temp(ic,i_t) - coefs_cond(ic,i_t)                 &
                                            * lrcp
  coefs_temp(ic,i_0) = coefs_temp(ic,i_0) - coefs_cond(ic,i_0)                 &
                                            * lrcp
end do


! If any melt source found
if ( l_any_melt ) then

  ! Calculate sum of the melting coefficients which are
  ! contributing to melt_source...

  ! Initialise to zero
  do ic2 = 1, nc_cor
    coefs_melt_source(ic2,i_q) = zero
    coefs_melt_source(ic2,i_t) = zero
    coefs_melt_source(ic2,i_0) = zero
  end do

  ! Loop over ice species
  do i_ice = n_cond_species_liq+1, n_cond_species
    ! If the current ice species may melt to form the
    ! current liquid species
    if ( cond_params(i_ice)%pt % i_cond_frzmlt == i_liq ) then
      ! Find points where the current ice species is melting
      nc_melt = 0
      do ic2 = 1, nc_cor
        ic = index_ic_cor(ic2)
        if ( l_melt(ic,i_ice) ) then
          nc_melt = nc_melt + 1
          index_ic_melt(nc_melt) = ic2
          ! Set flag to indicate the current ice species is melting
          ! into a liquid species whose process rates are being truncated.
          ! This needs to be accounted for in melt_ctl, in the event that
          ! melting of this ice species gets switched off...
          l_trunc_mlt(ic,i_ice) = .true.
        end if
      end do
      ! Add contributions to sum of melt coefs
      if ( nc_melt > 0 ) then
        do ic3 = 1, nc_melt
          ic2 = index_ic_melt(ic3)
          ic = index_ic_cor(ic2)
          coefs_melt_source(ic2,i_q) = coefs_melt_source(ic2,i_q)              &
                                     + coefs_melt(ic,i_q,i_ice)
          coefs_melt_source(ic2,i_t) = coefs_melt_source(ic2,i_t)              &
                                     + coefs_melt(ic,i_t,i_ice)
          coefs_melt_source(ic2,i_0) = coefs_melt_source(ic2,i_0)              &
                                     + coefs_melt(ic,i_0,i_ice)
        end do
      end if
    end if
  end do

  ! Set new condensation coefficients with exactly the right
  ! values to balance the source from melting and end up with
  ! exactly zero q_cond
  do ic2 = 1, nc_cor
    ic = index_ic_cor(ic2)
    coefs_cond(ic,i_q) = -coefs_melt_source(ic2,i_q)
    coefs_cond(ic,i_t) = -coefs_melt_source(ic2,i_t)
    coefs_cond(ic,i_0) = -coefs_melt_source(ic2,i_0) - q_cond(ic)
  end do

  ! Add modified contributions to the coefs for the T,q solve
  do ic2 = 1, nc_cor
    ic = index_ic_cor(ic2)
    coefs_q_vap(ic,i_q) = coefs_q_vap(ic,i_q) +coefs_cond(ic,i_q)
    coefs_q_vap(ic,i_t) = coefs_q_vap(ic,i_t) +coefs_cond(ic,i_t)
    coefs_q_vap(ic,i_0) = coefs_q_vap(ic,i_0) +coefs_cond(ic,i_0)
    lrcp =  L_con(ic) / cp_tot(ic)
    coefs_temp(ic,i_q) = coefs_temp(ic,i_q) + coefs_cond(ic,i_q)               &
                                              * lrcp
    coefs_temp(ic,i_t) = coefs_temp(ic,i_t) + coefs_cond(ic,i_t)               &
                                              * lrcp
    coefs_temp(ic,i_0) = coefs_temp(ic,i_0) + coefs_cond(ic,i_0)               &
                                              * lrcp
  end do

  ! No melt sources...
else  ! ( l_any_melt )

  ! Set condensation coefficients with just an explicit term
  ! to remove exactly the remaining condensate
  do ic2 = 1, nc_cor
    ic = index_ic_cor(ic2)
    coefs_cond(ic,i_q) = zero
    coefs_cond(ic,i_t) = zero
    coefs_cond(ic,i_0) = -q_cond(ic)
  end do

  ! Add modified contributions to the coefs for the T,q solve
  do ic2 = 1, nc_cor
    ic = index_ic_cor(ic2)
    coefs_q_vap(ic,i_0) = coefs_q_vap(ic,i_0) + coefs_cond(ic,i_0)
    lrcp =  L_con(ic) / cp_tot(ic)
    coefs_temp(ic,i_0) = coefs_temp(ic,i_0) + coefs_cond(ic,i_0)               &
                                              * lrcp
  end do

end if  ! ( l_any_melt )


return
end subroutine modify_coefs_liq

end module modify_coefs_liq_mod
