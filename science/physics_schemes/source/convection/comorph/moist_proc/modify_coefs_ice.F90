! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph
module modify_coefs_ice_mod

implicit none

contains


! Subroutine to correct the evaporation rate and/or melting rate
! and implicit solve coefficients to prevent negative condensate
! mixing ratios from occuring.  The list of points where
! negative values will occur if uncorrected is passed in.
! Coefficients are modified so-as to yield exactly zero
! condensate at those points.

subroutine modify_coefs_ice( n_points, nc_cor, index_ic_cor,                   &
                             L_sub, L_fus, cp_tot,                             &
                             q_cond, dq_cond, dq_melt, l_melt,                 &
                             coefs_cond, coefs_melt,                           &
                             coefs_cond_m,                                     &
                             coefs_temp, coefs_q_vap )

use comorph_constants_mod, only: real_cvprec, zero
use phase_change_coefs_mod, only: n_coefs, i_q, i_t, i_0

implicit none

! Number of points
integer, intent(in) :: n_points

! Number of points where coefficients need to be modified,
! and their indices
integer, intent(in) :: nc_cor
integer, intent(in) :: index_ic_cor(nc_cor)

! Latent heats of sublimation and fusion
real(kind=real_cvprec), intent(in) :: L_sub(n_points)
real(kind=real_cvprec), intent(in) :: L_fus(n_points)

! Total heat capacity per unit dry-mass
real(kind=real_cvprec), intent(in) :: cp_tot(n_points)

! Mixing ratio of ice species before implicit phase-changes
real(kind=real_cvprec), intent(in) :: q_cond(n_points)

! Increments due to condensation / evaporation and melting
real(kind=real_cvprec), intent(in) :: dq_cond(n_points)
real(kind=real_cvprec), intent(in) :: dq_melt(n_points)

! Flag for points where melting of this ice species is occuring
logical, intent(in) :: l_melt(n_points)

! Coefficients for condensation / evaporation and melting
! (corrected by this routine to avoid negative q_cond)
real(kind=real_cvprec), intent(in out) :: coefs_cond                           &
                                         ( n_points, n_coefs )
real(kind=real_cvprec), intent(in out) :: coefs_melt                           &
                                         ( n_points, n_coefs )
real(kind=real_cvprec), intent(in out) :: coefs_cond_m                         &
                                         ( n_points, n_coefs )

! Super-arrays containing coefficients in the implicit formula
! for temperature and q_vap after phase-changes
real(kind=real_cvprec), intent(in out) :: coefs_temp                           &
                                         ( n_points, n_coefs )
real(kind=real_cvprec), intent(in out) :: coefs_q_vap                          &
                                         ( n_points, n_coefs )

! Points where melting is occuring
integer :: nc_melt
integer :: index_ic_melt(nc_cor)

! Store for L/cp
real(kind=real_cvprec) :: lrcp

! Loop counters
integer :: ic, ic2


! Find points to be modified where melting is occuring...
nc_melt = 0
do ic2 = 1, nc_cor
  ic = index_ic_cor(ic2)
  if ( l_melt(ic) ) then
    nc_melt = nc_melt + 1
    index_ic_melt(nc_melt) = ic
  end if
end do


! Subtract the condensation contributions for
! this species from the implicit solve coefficients for T,q
do ic2 = 1, nc_cor
  ic = index_ic_cor(ic2)
  coefs_q_vap(ic,i_q) = coefs_q_vap(ic,i_q) - coefs_cond(ic,i_q)
  coefs_q_vap(ic,i_t) = coefs_q_vap(ic,i_t) - coefs_cond(ic,i_t)
  coefs_q_vap(ic,i_0) = coefs_q_vap(ic,i_0) - coefs_cond(ic,i_0)
  lrcp =  L_sub(ic) / cp_tot(ic)
  coefs_temp(ic,i_q) = coefs_temp(ic,i_q) - coefs_cond(ic,i_q)                 &
                                            * lrcp
  coefs_temp(ic,i_t) = coefs_temp(ic,i_t) - coefs_cond(ic,i_t)                 &
                                            * lrcp
  coefs_temp(ic,i_0) = coefs_temp(ic,i_0) - coefs_cond(ic,i_0)                 &
                                            * lrcp
end do

! If melting occuring at any points which need modifying...
if ( nc_melt > 0 ) then
  ! Subtract the contribution from melting of this species from
  ! the implicit solve coefficients for T,q
  do ic2 = 1, nc_melt
    ic = index_ic_melt(ic2)
    lrcp = L_fus(ic) / cp_tot(ic)
    coefs_temp(ic,i_q) = coefs_temp(ic,i_q) + coefs_melt(ic,i_q)               &
                                              * lrcp
    coefs_temp(ic,i_t) = coefs_temp(ic,i_t) + coefs_melt(ic,i_t)               &
                                              * lrcp
    coefs_temp(ic,i_0) = coefs_temp(ic,i_0) + coefs_melt(ic,i_0)               &
                                              * lrcp
  end do
end if

! Reset condensation / evaporation and melting coefs to zero
do ic2 = 1, nc_cor
  ic = index_ic_cor(ic2)
  coefs_cond(ic,i_q) = zero
  coefs_cond(ic,i_t) = zero
  coefs_cond(ic,i_0) = zero
  coefs_melt(ic,i_q) = zero
  coefs_melt(ic,i_t) = zero
  coefs_melt(ic,i_0) = zero
  coefs_cond_m(ic,i_q) = zero
  coefs_cond_m(ic,i_t) = zero
  coefs_cond_m(ic,i_0) = zero
end do

! At melting points
if ( nc_melt > 0 ) then
  do ic2 = 1, nc_melt
    ic = index_ic_melt(ic2)

    ! Compute new melting coefficient, set purely explicitly
    ! to scale the melting rate down by q_cond / increment
    coefs_melt(ic,i_0) = dq_melt(ic)                                           &
          * ( q_cond(ic) / ( dq_melt(ic) - dq_cond(ic) ) )

    ! Add modified melting contribution to the T,q solve coefs
    lrcp = L_fus(ic) / cp_tot(ic)
    coefs_temp(ic,i_0) = coefs_temp(ic,i_0) - coefs_melt(ic,i_0)               &
                                              * lrcp
  end do
end if

! At all modified points
do ic2 = 1, nc_cor
  ic = index_ic_cor(ic2)

  ! Compute new condensation coefficient, set purely explicitly
  ! to remove exactly the condensate remaining after melting
  coefs_cond(ic,i_0) = -( q_cond(ic) - coefs_melt(ic,i_0) )

  ! Copy into alternative condensation coefficient, to ensure
  ! we get same behaviour if melting switched on/off
  coefs_cond_m(ic,i_0) = coefs_cond(ic,i_0)

  ! Add modified contributions to the coefs for the T,q solve
  coefs_q_vap(ic,i_0) = coefs_q_vap(ic,i_0) + coefs_cond(ic,i_0)
  lrcp = L_sub(ic) / cp_tot(ic)
  coefs_temp(ic,i_0) = coefs_temp(ic,i_0) + coefs_cond(ic,i_0)                 &
                                            * lrcp
end do


return
end subroutine modify_coefs_ice

end module modify_coefs_ice_mod
