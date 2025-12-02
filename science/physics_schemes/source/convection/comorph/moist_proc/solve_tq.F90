! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module solve_tq_mod

implicit none

contains

! Subroutine to compute estimates of temperature and water vapour
! mixing ratio after phase changes, using imput coefficients
! for their budgets from various implicitly-solved process-rates
subroutine solve_tq( n_points, l_full_do,                                      &
                     ref_temp, qsat_liq_ref,                                   &
                     coefs_temp, coefs_q_vap,                                  &
                     temperature, q_vap, imp_temp, imp_q_vap,                  &
                     kq_cond, kt_cond,                                         &
                     coefs_cond, coefs_cond_m, coefs_melt,                     &
                     nc, index_ic )

use comorph_constants_mod, only: real_cvprec, one, two, four,                  &
                                 sqrt_min_delta,                               &
                                 n_cond_species, n_cond_species_liq
use phase_change_coefs_mod, only: n_coefs, i_q, i_t, i_0

implicit none

! Number of points
integer, intent(in) :: n_points

! Flag for whether to use full do-loops instead of indirect indexing
logical, intent(in) :: l_full_do

! Reference temperature used for linearised qsat calculations
real(kind=real_cvprec), intent(in) :: ref_temp(n_points)

! Saturation water-vapour mixing-ratio at ref_temp w.r.t. liquid
real(kind=real_cvprec), intent(in) :: qsat_liq_ref(n_points)

! Coefficients in the budget equations for water-vapour
! mixing-ratio and temperature
real(kind=real_cvprec), intent(in out) :: coefs_temp                           &
                                      ( n_points, n_coefs )
real(kind=real_cvprec), intent(in out) :: coefs_q_vap                          &
                                      ( n_points, n_coefs )

! Temperature and water-vapour mixing-ratio with only
! explicit increments added on so far
real(kind=real_cvprec), intent(in) :: temperature(n_points)
real(kind=real_cvprec), intent(in) :: q_vap(n_points)

! Output implicitly solve temperature / water-vapour mixing-ratio
! (minus their reference values)
real(kind=real_cvprec), intent(in out) :: imp_temp(n_points)
real(kind=real_cvprec), intent(in out) :: imp_q_vap(n_points)
! (need intent inout to preserve the input values at points where
!  we don't calculate anything, when using indirect indexing)

! Vapour and heat exchange and process rate coefficients;
! only passed in-out so that they can be scaled down where
! required to avoid a rare numerical problem
real(kind=real_cvprec), intent(in out) :: kq_cond                              &
                                    ( n_points, n_cond_species )
real(kind=real_cvprec), intent(in out) :: kt_cond                              &
                                    ( n_points, n_cond_species )
real(kind=real_cvprec), intent(in out) :: coefs_cond                           &
                           ( n_points, n_coefs, n_cond_species )
real(kind=real_cvprec), intent(in out) :: coefs_cond_m                         &
    ( n_points, n_coefs, n_cond_species_liq+1 : n_cond_species )
real(kind=real_cvprec), intent(in out) :: coefs_melt                           &
    ( n_points, n_coefs, n_cond_species_liq+1 : n_cond_species )


! Optional input list of points to do indirect indexing calculations
! at only certain points
integer, optional, intent(in) :: nc
integer, optional, intent(in) :: index_ic(n_points)

! Points where numerical problem (negative denominator) found
integer :: nc_neg
integer :: index_ic_neg(n_points)

! Work variables for solving quadratic equation when scaling
! down the coefficients at points with numerical problem
real(kind=real_cvprec) :: alpha
real(kind=real_cvprec) :: a, b, c
real(kind=real_cvprec) :: work(n_points)

! Loop counters
integer :: ic, ic2, i_coef, i_cond


! Safety check to avoid coefs_T^T >= 1
! (solution yields nonsense in this limit.  Should only
!  happen very rarely in ludicrously extreme circumstances;
!  coefs_T^T is usually expected to be negative!
!  but check just in case as could lead to div-by-zero).

if ( l_full_do ) then
  ! Full-field calculation
  nc_neg = 0
  do ic = 1, n_points
    if ( coefs_temp(ic,i_t) > one - sqrt_min_delta ) then
      nc_neg = nc_neg + 1
      index_ic_neg(nc_neg) = ic
    end if
  end do
else
  ! Indirect indexing version of the above calculations
  nc_neg = 0
  do ic2 = 1, nc
    ic = index_ic(ic2)
    if ( coefs_temp(ic,i_t) > one - sqrt_min_delta ) then
      nc_neg = nc_neg + 1
      index_ic_neg(nc_neg) = ic
    end if
  end do
end if

! If coefs_T^T >= 1 found anywhere...
if ( nc_neg > 0 ) then
  do ic2 = 1, nc_neg
    ic = index_ic_neg(ic2)
    ! Scale all coefficients down consistently to avoid
    ! coefs_T^T too close to one
    alpha = ( one - sqrt_min_delta ) / coefs_temp(ic,i_t)
    ! This is equivalent to reducing the timestep, i.e.
    ! imposing a CFL limit.
    do i_coef = 1, n_coefs
      coefs_q_vap(ic,i_coef) = coefs_q_vap(ic,i_coef) * alpha
      coefs_temp(ic,i_coef)  = coefs_temp(ic,i_coef) * alpha
    end do
    do i_cond = 1, n_cond_species
      kq_cond(ic,i_cond) = kq_cond(ic,i_cond) * alpha
      kt_cond(ic,i_cond) = kt_cond(ic,i_cond) * alpha
      do i_coef = 1, n_coefs
        coefs_cond(ic,i_coef,i_cond)                                           &
          = coefs_cond(ic,i_coef,i_cond) * alpha
      end do
    end do
    do i_cond = n_cond_species_liq+1, n_cond_species
      do i_coef = 1, n_coefs
        coefs_cond_m(ic,i_coef,i_cond)                                         &
          = coefs_cond_m(ic,i_coef,i_cond) * alpha
        coefs_melt(ic,i_coef,i_cond)                                           &
          = coefs_melt(ic,i_coef,i_cond) * alpha
      end do
    end do
  end do
end if

if ( l_full_do ) then
  ! Full-field calculation

  ! Compute denominator in the formula for q_vap - qsat_liq_ref
  do ic = 1, n_points
    work(ic) = one + coefs_q_vap(ic,i_q) - coefs_temp(ic,i_t)
    ! Note: the term stored in work is guaranteed to be positive,
    ! since we always have coefs_q^q >= 0,
    ! and the check earlier in this routine ensures coefs_T^T < 1
  end do
  do ic = 1, n_points
    imp_q_vap(ic) = work(ic)                                                   &
                  + ( coefs_q_vap(ic,i_t) * coefs_temp(ic,i_q)                 &
                    - coefs_q_vap(ic,i_q) * coefs_temp(ic,i_t) )
    ! imp_q_vap stores the denominator
  end do

  ! Check to avoid denominator near or below zero
  ! (can occasionally happen with mixed-phase cloud at very low
  !  temperatures; yields a spurious solution which is an
  !  artefact of the linearisation of qsat, in which the assumed
  !  qsat_liq and qsat_ice cross over again well below the
  !  melting point, which should never happen).
  nc_neg = 0
  do ic = 1, n_points
    if ( imp_q_vap(ic) < sqrt_min_delta * work(ic) ) then
      nc_neg = nc_neg + 1
      index_ic_neg(nc_neg) = ic
    end if
  end do

else
  ! Indirect indexing version of the above calculations

  do ic2 = 1, nc
    ic = index_ic(ic2)
    work(ic) = one + coefs_q_vap(ic,i_q) - coefs_temp(ic,i_t)
    imp_q_vap(ic) = work(ic)                                                   &
                  + ( coefs_q_vap(ic,i_t) * coefs_temp(ic,i_q)                 &
                    - coefs_q_vap(ic,i_q) * coefs_temp(ic,i_t) )
  end do

  nc_neg = 0
  do ic2 = 1, nc
    ic = index_ic(ic2)
    if ( imp_q_vap(ic) < sqrt_min_delta * work(ic) ) then
      nc_neg = nc_neg + 1
      index_ic_neg(nc_neg) = ic
    end if
  end do

end if


! If negative value of denominator found...
if ( nc_neg > 0 ) then
  do ic2 = 1, nc_neg
    ic = index_ic_neg(ic2)

    ! Find scaling factor to reduce timestep and coefficients,
    ! to avoid getting nonsense (root of quadratic)
    a = ( coefs_q_vap(ic,i_t) * coefs_temp(ic,i_q)                             &
        - coefs_q_vap(ic,i_q) * coefs_temp(ic,i_t) )
    b = ( coefs_q_vap(ic,i_q) - coefs_temp(ic,i_t) )                           &
        * ( one - sqrt_min_delta )
    c = one - sqrt_min_delta
    ! To get here, I think we must have b > 0 and a < 0.
    ! Therefore the positive solution will be the negative root:
    alpha = ( -b - sqrt( b*b - four*a*c ) ) / ( two * a )

    ! Scale down all the coefficients at this point.
    ! This is equivalent to reducing the timestep, i.e.
    ! imposing a CFL limit.
    do i_coef = 1, n_coefs
      coefs_q_vap(ic,i_coef) = coefs_q_vap(ic,i_coef) * alpha
      coefs_temp(ic,i_coef)  = coefs_temp(ic,i_coef) * alpha
    end do
    do i_cond = 1, n_cond_species
      kq_cond(ic,i_cond) = kq_cond(ic,i_cond) * alpha
      kt_cond(ic,i_cond) = kt_cond(ic,i_cond) * alpha
      do i_coef = 1, n_coefs
        coefs_cond(ic,i_coef,i_cond)                                           &
          = coefs_cond(ic,i_coef,i_cond) * alpha
      end do
    end do
    do i_cond = n_cond_species_liq+1, n_cond_species
      do i_coef = 1, n_coefs
        coefs_cond_m(ic,i_coef,i_cond)                                         &
          = coefs_cond_m(ic,i_coef,i_cond) * alpha
        coefs_melt(ic,i_coef,i_cond)                                           &
          = coefs_melt(ic,i_coef,i_cond) * alpha
      end do
    end do

    ! Recalculate the denominator
    imp_q_vap(ic)                                                              &
            = one + coefs_q_vap(ic,i_q) - coefs_temp(ic,i_t)                   &
            + ( coefs_q_vap(ic,i_t) * coefs_temp(ic,i_q)                       &
              - coefs_q_vap(ic,i_q) * coefs_temp(ic,i_t) )

  end do
end if


if ( l_full_do ) then
  ! Full-field calculation

  ! Compute solution for q_vap - qsat_liq_ref
  do ic = 1, n_points
    imp_q_vap(ic)                                                              &
      = (                                                                      &
          ( q_vap(ic) - qsat_liq_ref(ic) - coefs_q_vap(ic,i_0) )               &
          * ( one - coefs_temp(ic,i_t) )                                       &
        - ( temperature(ic) - ref_temp(ic) + coefs_temp(ic,i_0) )              &
          * coefs_q_vap(ic,i_t)                                                &
        ) / imp_q_vap(ic)  ! Use denominator calculated earlier
  end do

  ! Substitute solution for q_vap back into temperature budget
  ! equation to compute implicit solution for
  ! temperature - ref_temp
  do ic = 1, n_points
    imp_temp(ic)                                                               &
      = (                                                                      &
          temperature(ic) - ref_temp(ic)                                       &
        + coefs_temp(ic,i_0)                                                   &
        + coefs_temp(ic,i_q) * imp_q_vap(ic)                                   &
        ) / (                                                                  &
          one - coefs_temp(ic,i_t)                                             &
        )
  end do

else  ! ( l_full_do )
  ! Indirect indexing version of the above

  do ic2 = 1, nc
    ic = index_ic(ic2)
    imp_q_vap(ic)                                                              &
      = (                                                                      &
          ( q_vap(ic) - qsat_liq_ref(ic) - coefs_q_vap(ic,i_0) )               &
          * ( one - coefs_temp(ic,i_t) )                                       &
        - ( temperature(ic) - ref_temp(ic) + coefs_temp(ic,i_0) )              &
          * coefs_q_vap(ic,i_t)                                                &
        ) / imp_q_vap(ic)
    imp_temp(ic)                                                               &
      = (                                                                      &
          temperature(ic) - ref_temp(ic)                                       &
        + coefs_temp(ic,i_0)                                                   &
        + coefs_temp(ic,i_q) * imp_q_vap(ic)                                   &
        ) / (                                                                  &
          one - coefs_temp(ic,i_t)                                             &
        )
  end do

end if  ! ( l_full_do )


return
end subroutine solve_tq

end module solve_tq_mod
