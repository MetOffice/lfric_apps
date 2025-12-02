! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module collision_rate_mod

implicit none

contains


! Subroutine to calculate rate of collection of hydrometeor
! species 1 by hydrometeor species 2.
subroutine collision_rate( n_points,                                           &
                           q_cond_1, n_cond_2,                                 &
                           r_cond_1, r_cond_2,                                 &
                           wf_cond_1, wf_cond_2,                               &
                           delta_t, rho_dry, rho_wet,                          &
                           rho_cond_1, rho_cond_2,                             &
                           area_coef_1, area_coef_2,                           &
                           dq_col )

use comorph_constants_mod, only: pi, coef_wf_spread, col_eff_coef,             &
                     kin_visc, drag_coef_cond,                                 &
                     real_cvprec, six, half, four_thirds, min_float,           &
                     sqrt_min_float

implicit none

! Number of points
integer, intent(in) :: n_points

! Mixing ratio of the collected hydrometeor species
real(kind=real_cvprec), intent(in) :: q_cond_1(n_points)
! Number concentration of the collecting species
real(kind=real_cvprec), intent(in) :: n_cond_2(n_points)

! Representative radii of the 2 hydrometeor species
real(kind=real_cvprec), intent(in) :: r_cond_1(n_points)
real(kind=real_cvprec), intent(in) :: r_cond_2(n_points)

! Representative fall-speeds of the 2 hydrometeor species
real(kind=real_cvprec), intent(in) :: wf_cond_1(n_points)
real(kind=real_cvprec), intent(in) :: wf_cond_2(n_points)

! Time-interval for converting process rate to increment
real(kind=real_cvprec), intent(in) :: delta_t(n_points)

! Dry-density of the air
! (needed to get number per unit volume from n_cond_2)
real(kind=real_cvprec), intent(in) :: rho_dry(n_points)
! Wet-density of the air
! (needed to calculate deflection of hydrometeors by drag)
real(kind=real_cvprec), intent(in) :: rho_wet(n_points)

! Densities of hydrometeors (either liquid or ice)
real(kind=real_cvprec), intent(in) :: rho_cond_1
real(kind=real_cvprec), intent(in) :: rho_cond_2

! Ratio of hydrometeor actual cross-section areas to what you
! get by assuming a sphere.
! This scales the area swept out by the particles.
real(kind=real_cvprec), intent(in) :: area_coef_1
real(kind=real_cvprec), intent(in) :: area_coef_2

! Increment of mixing ratio of species 1 collected by species 2
real(kind=real_cvprec), intent(out) :: dq_col(n_points)


! Rate (in s-1) at which the particles can deflect eachother
! out of the way as they approach, via the drag force from
! the deflection air-flow
real(kind=real_cvprec) :: def_rate(n_points)

! Collection efficiency / dimensionless
real(kind=real_cvprec) :: coll_eff(n_points)

! "Collision radius" = how close the centres of the 2 species'
! particles must come for them to collide
real(kind=real_cvprec) :: r_col(n_points)

! "Collision velocity" = mean relative velocity between the
! 2 species' particles
real(kind=real_cvprec) :: w_col(n_points)

! Square roots of input constants
real(kind=real_cvprec) :: sqrt_area_coef_1
real(kind=real_cvprec) :: sqrt_area_coef_2

! Difference between fall-speeds
real(kind=real_cvprec) :: wf_diff

! Precalculate square of coefficient
real(kind=real_cvprec), parameter :: coef_wf_spread_sq = coef_wf_spread        &
                                                       * coef_wf_spread

! Loop counter
integer :: ic


!----------------------------------------------------------------
! 1) Calculate deflection of particles by each-other's fall-flow
!----------------------------------------------------------------

! First, we need to calculate the rate of deflection of
! species 1's particles out of the path of species 2
! by the deflection flow induced by species 2's fall.
! We assume that species 2's fall induces a deflection
! flow with horizontal speed equal to its fall-speed wf_cond_2.
! The drag force that this wind-speed induces on species 1's
! particles is then:
!
! F_d = area_coef_1 pi rho_air ( 6 nu r_1 wf_2
!                              + 1/2 r_1^2 Cd wf_2^2 )
!
! Dividing by the mass of particle 1 (4/3 pi rho_1 r_1^3) to get
! the induced deflection acceleration, and then dividing by the
! deflection speed to get a reciprocal timescale, we have:
!
! 1/tau = area_coef_1 (rho_air/rho_1) ( 6 nu
!                                     + 1/2 r_1 Cd wf_2 )
!                       / ( 4/3 r_1^2 )
!
! We also add on the vice-versa rate of deflection; of species
! 2's particles by species 1

! Calculate 1/tau as described above, storing it in def_rate
do ic = 1, n_points
  def_rate(ic) =                                                               &
    area_coef_1 * ( rho_wet(ic) / rho_cond_1 )                                 &
              * ( six*kin_visc + half*drag_coef_cond                           &
                                 *r_cond_1(ic)*wf_cond_2(ic) )                 &
              / ( four_thirds * max( r_cond_1(ic)*r_cond_1(ic), min_float ) )  &
             +                                                                 &
    area_coef_2 * ( rho_wet(ic) / rho_cond_2 )                                 &
              * ( six*kin_visc + half*drag_coef_cond                           &
                                 *r_cond_2(ic)*wf_cond_1(ic) )                 &
              / ( four_thirds * max( r_cond_2(ic)*r_cond_2(ic), min_float ) )
end do


!----------------------------------------------------------------
! 2) Calculate collection efficiency
!----------------------------------------------------------------

! Calculate distance between centres of the hydrometeors
! required for them to actually collide
sqrt_area_coef_1 = sqrt(area_coef_1)
sqrt_area_coef_2 = sqrt(area_coef_2)
do ic = 1, n_points
  r_col(ic) = sqrt_area_coef_1 * r_cond_1(ic)                                  &
            + sqrt_area_coef_2 * r_cond_2(ic)
end do

! Effective relative velocity between the 2 colliding species
do ic = 1, n_points
  wf_diff = wf_cond_2(ic) - wf_cond_1(ic)
  w_col(ic) = sqrt( wf_diff * wf_diff                                          &
                  + coef_wf_spread_sq * wf_cond_1(ic) * wf_cond_2(ic) )
  ! This formula just yields abs( wf_2 - wf_1 ) in the limit
  ! that one species' fall-speed is much greater than the other.
  ! The extra term scaled by coef_wf_spread is an approximation to
  ! prevent us getting zero relative velocity when the mean
  ! fall-speeds of the 2 species are equal.  In reality
  ! we would have a distribution of fall-speeds for both,
  ! so that collisions should still occur.
end do

! Calculate collection efficiency; parameterised as:
!
! col_eff = exp( -col_eff_coef (r_col/w_col)
!                              ( 1/tau_1 + 1/tau_2 ) )
!
! where r_col/w_col is the timescale for the collision
do ic = 1, n_points
  coll_eff(ic) = exp( -col_eff_coef                                            &
                       * ( r_col(ic) / max( w_col(ic), sqrt_min_float ) )      &
                       * def_rate(ic) )
  ! coll_eff now stores the collection efficiency / dimensionless
end do


!----------------------------------------------------------------
! 3) Compute collection of species 1's mixing ratio by species 2
!----------------------------------------------------------------

do ic = 1, n_points
  ! Increment to species 1's mixing ratio =
  dq_col(ic) =                                                                 &
          !   Effective cross-section of the collision-range
          !   area for a single hydrometeor
                pi * r_col(ic) * r_col(ic)                                     &
          ! * Effective relative velocity between the 2 species
              * w_col(ic)                                                      &
          ! * Collection efficiency:
              * coll_eff(ic)                                                   &
  ! Units of the above sweep-out rate are m3 s-1
          ! * Number of collecting hydrometeors per unit volume:
              * n_cond_2(ic) * rho_dry(ic)                                     &
          ! * Mixing ratio of the collected species
              * q_cond_1(ic)                                                   &
  ! The above process rate is now in kg kg-1 s-1
          ! * Time interval
              * delta_t(ic)
end do


return
end subroutine collision_rate



! Compression wrapper for the above
subroutine collision_rate_cmpr( n_points, nc, index_ic,                        &
                                q_cond_1, n_cond_2,                            &
                                r_cond_1, r_cond_2,                            &
                                wf_cond_1, wf_cond_2,                          &
                                delta_t, rho_dry, rho_wet,                     &
                                rho_cond_1, rho_cond_2,                        &
                                area_coef_1, area_coef_2,                      &
                                dq_col )

use comorph_constants_mod, only: real_cvprec, zero, cmpr_thresh

implicit none

! Number of points
integer, intent(in) :: n_points

! Number of points where the two condensed water species coincide
integer, intent(in) :: nc
! Indices of those points
integer, intent(in) :: index_ic(nc)

! All other argumnts the same as in collision_rate above
real(kind=real_cvprec), intent(in) :: q_cond_1(n_points)
real(kind=real_cvprec), intent(in) :: n_cond_2(n_points)
real(kind=real_cvprec), intent(in) :: r_cond_1(n_points)
real(kind=real_cvprec), intent(in) :: r_cond_2(n_points)
real(kind=real_cvprec), intent(in) :: wf_cond_1(n_points)
real(kind=real_cvprec), intent(in) :: wf_cond_2(n_points)
real(kind=real_cvprec), intent(in) :: delta_t(n_points)
real(kind=real_cvprec), intent(in) :: rho_dry(n_points)
real(kind=real_cvprec), intent(in) :: rho_wet(n_points)

real(kind=real_cvprec), intent(in) :: rho_cond_1
real(kind=real_cvprec), intent(in) :: rho_cond_2
real(kind=real_cvprec), intent(in) :: area_coef_1
real(kind=real_cvprec), intent(in) :: area_coef_2

real(kind=real_cvprec), intent(out) :: dq_col(n_points)


! Super-array for compressed inputs
real(kind=real_cvprec), allocatable :: args_cmpr(:,:)
integer, parameter :: i_q_cond_1 = 1
integer, parameter :: i_n_cond_2 = 2
integer, parameter :: i_r_cond_1 = 3
integer, parameter :: i_r_cond_2 = 4
integer, parameter :: i_wf_cond_1 = 5
integer, parameter :: i_wf_cond_2 = 6
integer, parameter :: i_delta_t = 7
integer, parameter :: i_rho_dry = 8
integer, parameter :: i_rho_wet = 9
integer, parameter :: i_dq_col = 10
integer, parameter :: n_args = 10

! Loop counters
integer :: ic, ic2


if ( real(nc,real_cvprec) > cmpr_thresh * real(n_points,real_cvprec) ) then
  ! If majority of points have both species present...

  ! Just call the main routine on the full arrays
  call collision_rate( n_points,                                               &
                       q_cond_1, n_cond_2,                                     &
                       r_cond_1, r_cond_2,                                     &
                       wf_cond_1, wf_cond_2,                                   &
                       delta_t, rho_dry, rho_wet,                              &
                       rho_cond_1, rho_cond_2,                                 &
                       area_coef_1, area_coef_2,                               &
                       dq_col )

else  ! ( real(nc,real_cvprec) <= cmpr_thresh * real(n_points,real_cvprec) )
  ! If only a minority of points have both species present...

  ! Do a compressed call...

  ! Allocate compression array
  allocate( args_cmpr ( nc, n_args  ) )

  ! Compress inputs
  do ic2 = 1, nc
    ic = index_ic(ic2)
    args_cmpr(ic2,i_q_cond_1)     = q_cond_1(ic)
    args_cmpr(ic2,i_n_cond_2)     = n_cond_2(ic)
    args_cmpr(ic2,i_r_cond_1)     = r_cond_1(ic)
    args_cmpr(ic2,i_r_cond_2)     = r_cond_2(ic)
    args_cmpr(ic2,i_wf_cond_1)    = wf_cond_1(ic)
    args_cmpr(ic2,i_wf_cond_2)    = wf_cond_2(ic)
    args_cmpr(ic2,i_delta_t)      = delta_t(ic)
    args_cmpr(ic2,i_rho_dry)      = rho_dry(ic)
    args_cmpr(ic2,i_rho_wet)      = rho_wet(ic)
  end do

  ! Call the main routine on the compressed arrays
  call collision_rate( nc,                                                     &
         args_cmpr(:,i_q_cond_1), args_cmpr(:,i_n_cond_2),                     &
         args_cmpr(:,i_r_cond_1), args_cmpr(:,i_r_cond_2),                     &
         args_cmpr(:,i_wf_cond_1), args_cmpr(:,i_wf_cond_2),                   &
         args_cmpr(:,i_delta_t),                                               &
         args_cmpr(:,i_rho_dry), args_cmpr(:,i_rho_wet),                       &
         rho_cond_1, rho_cond_2,                                               &
         area_coef_1, area_coef_2,                                             &
         args_cmpr(:,i_dq_col) )

  ! Initialise output to zero everywhere so it gets set at
  ! points which aren't in the compression list
  do ic = 1, n_points
    dq_col(ic) = zero
  end do

  ! Decompress outputs
  do ic2 = 1, nc
    ic = index_ic(ic2)
    dq_col(ic) = args_cmpr(ic2,i_dq_col)
  end do

  ! Deallocate
  deallocate( args_cmpr )

end if  ! ( real(nc,real_cvprec) <= cmpr_thresh * real(n_points,real_cvprec) )


return
end subroutine collision_rate_cmpr


end module collision_rate_mod
