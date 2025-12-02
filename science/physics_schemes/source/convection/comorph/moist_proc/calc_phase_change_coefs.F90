! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_phase_change_coefs_mod

implicit none

contains

! Subroutine to calculate the implicit process rate coefficients
! for condensation / evaporation, deposition / sublimation,
! and melting.
! NOTE: for computational efficiency, the scaling by the
! timestep to convert from process rates to increments is done
! by scaling all the coefficients calculated in here by delta_t.
subroutine calc_phase_change_coefs( n_points, nc, index_ic, l_full_do,         &
                                    ref_temp,                                  &
                                    qsat_liq_ref, qsat_ice_ref,                &
                                    dqsatdT_liq, dqsatdT_ice,                  &
                                    L_con, L_sub, L_fus, cp_tot,               &
                                    delta_t, delta_z, wind_w,                  &
                                    prev_temp,                                 &
                                    kq_cond, kt_cond, wf_cond,                 &
                                    q_loc_cond, dq_frz_cond,                   &
                                    q_vap, coefs_cond,                         &
                                    coefs_melt, coefs_cond_m,                  &
                                    coefs_temp, coefs_q_vap )

use comorph_constants_mod, only: real_cvprec, zero, n_cond_species,            &
                                 n_cond_species_liq, one,                      &
                                 cp_dry, cp_vap, cp_liq, cp_ice, melt_temp,    &
                                 sqrt_min_float
use phase_change_coefs_mod, only: n_coefs, i_q, i_t, i_0

implicit none

! Number of points
integer, intent(in) :: n_points

! Number of points where each condensed water species is non-zero
integer, intent(in) :: nc( n_cond_species )
! Indices of those points
integer, intent(in) :: index_ic( n_points, n_cond_species )
! Flags for whether to do full-field do-loops instead of
! indirect indexing for each condensate species
logical :: l_full_do( n_cond_species )

! Reference temperature used for linearised qsat calculations
real(kind=real_cvprec), intent(in) :: ref_temp(n_points)
! Saturation water vapour mixing ratio qsat and dqsat/dT at the
! reference temperature, w.r.t. liquid water and ice.
real(kind=real_cvprec), intent(in) :: qsat_liq_ref(n_points)
real(kind=real_cvprec), intent(in) :: qsat_ice_ref(n_points)
real(kind=real_cvprec), intent(in) :: dqsatdT_liq(n_points)
real(kind=real_cvprec), intent(in) :: dqsatdT_ice(n_points)

! Latent heats of condensation, sublimation and fusion
real(kind=real_cvprec), intent(in) :: L_con(n_points)
real(kind=real_cvprec), intent(in) :: L_sub(n_points)
real(kind=real_cvprec), intent(in) :: L_fus(n_points)

! Total heat capacity per unit dry-mass
real(kind=real_cvprec), intent(in) :: cp_tot(n_points)

! Time interval
real(kind=real_cvprec), intent(in) :: delta_t(n_points)
! Vertical height interval back to where prev_temp is defined
real(kind=real_cvprec), intent(in) :: delta_z(n_points)
! Vertical wind-speed of the surrounding air
real(kind=real_cvprec), intent(in) :: wind_w(n_points)
! Temperature at the previous level
real(kind=real_cvprec), intent(in) :: prev_temp(n_points)

! Vapour and heat exchange coefficients for each species
! (dimensionless as already scaled by the timestep delta_t)
real(kind=real_cvprec), intent(in) :: kq_cond                                  &
                                    ( n_points, n_cond_species )
real(kind=real_cvprec), intent(in) :: kt_cond                                  &
                                    ( n_points, n_cond_species )

! Fall-speed of each condensed water species
real(kind=real_cvprec), intent(in) :: wf_cond                                  &
                                    ( n_points, n_cond_species )
! Local mixing ratio of each condensed water species
real(kind=real_cvprec), intent(in) :: q_loc_cond                               &
                                    ( n_points, n_cond_species )
! Total rate of freezing of liquid onto each ice species
real(kind=real_cvprec), intent(in) :: dq_frz_cond                              &
             ( n_points, n_cond_species_liq+1 : n_cond_species )

! Water-vapour mixing ratio
real(kind=real_cvprec), intent(in) :: q_vap(n_points)

! Super-array containing the 3 condensation / deposition rate
! coefficients for each condensed water species
real(kind=real_cvprec), intent(out) :: coefs_cond                              &
                          ( n_points, n_coefs, n_cond_species )

! Coefficients for melting rate of ice species
real(kind=real_cvprec), intent(out) :: coefs_melt                              &
   ( n_points, n_coefs, n_cond_species_liq+1 : n_cond_species )

! Coefficients for condensation on melting ice
real(kind=real_cvprec), intent(out) :: coefs_cond_m                            &
   ( n_points, n_coefs, n_cond_species_liq+1 : n_cond_species )

! Super-arrays containing coefficients in the implicit formula
! for temperature and q_vap after phase-changes
real(kind=real_cvprec), intent(out) :: coefs_temp                              &
                                       ( n_points, n_coefs )
real(kind=real_cvprec), intent(out) :: coefs_q_vap                             &
                                       ( n_points, n_coefs )

! Array to store Lc/cp
real(kind=real_cvprec) :: lrcp(n_points)

! Array to store the heat diffusion term
!   kt (cp_dry + cp_vap q_vap)
real(kind=real_cvprec) :: kt_cp(n_points)
! Array to store the precip fall heat transport term
!   ( q_cond cp_cond ( wind_w - wf_cond ) / delta_z ) / kt_cp
real(kind=real_cvprec) :: fall_term(n_points)

! Loop counters
integer :: ic, ic2, i_liq, i_ice, i_cond, i_coef


! Initialise all coefficients to zero
do i_cond = 1, n_cond_species
  do i_coef = 1, n_coefs
    do ic = 1, n_points
      coefs_cond(ic,i_coef,i_cond) = zero
    end do
  end do
end do
do i_ice = n_cond_species_liq+1, n_cond_species
  do i_coef = 1, n_coefs
    do ic = 1, n_points
      coefs_melt(ic,i_coef,i_ice) = zero
    end do
  end do
end do
do i_ice = n_cond_species_liq+1, n_cond_species
  do i_coef = 1, n_coefs
    do ic = 1, n_points
      coefs_cond_m(ic,i_coef,i_ice) = zero
    end do
  end do
end do
do i_coef = 1, n_coefs
  do ic = 1, n_points
    coefs_temp(ic,i_coef) = zero
  end do
end do
do i_coef = 1, n_coefs
  do ic = 1, n_points
    coefs_q_vap(ic,i_coef) = zero
  end do
end do


! Store L_con / cp_tot
do ic = 1, n_points
  lrcp(ic) = L_con(ic) / cp_tot(ic)
end do

! Loop over liquid species
do i_liq = 1, n_cond_species_liq
  ! If any points active for this species
  if ( nc(i_liq) > 0 ) then
    ! If calculations required at majority of points
    if ( l_full_do(i_liq) ) then
      ! Full-field calculation...

      ! Store the term  kt_cp =
      ! kt (cp_dry + cp_vap q_vap)
      do ic = 1, n_points
        kt_cp(ic) = kt_cond(ic,i_liq)                                          &
                    * ( cp_dry + cp_vap * q_vap(ic) )
      end do
      ! Don't let kt_cp go to zero (gives NaNs due to div-by-zero
      ! at points where there is no condensate)
      do ic = 1, n_points
        kt_cp(ic) = max( kt_cp(ic), sqrt_min_float )
      end do
      ! Store the term  B/kt_cp =
      ! ( q_cond cp_cond ( wind_w - wf_cond ) / delta_z ) / kt_cp
      do ic = 1, n_points
        fall_term(ic) = ( q_loc_cond(ic,i_liq) * cp_liq                        &
            * ( wind_w(ic) - wf_cond(ic,i_liq) ) / delta_z(ic) )               &
                      / ( kt_cp(ic) / delta_t(ic) )
        ! Note: removing factor of delta_t included in kt_cond
      end do
      ! Don't let the fall term exceed one (this causes a
      ! change of sign of the temperature-dependence coefficient
      ! below, leading to extreme solutions which are
      ! inconsistent with the thermal eqm assumptions posed).
      do ic = 1, n_points
        fall_term(ic) = min( fall_term(ic), one )
      end do

      ! Calculate condensation coefficients for current species:
      ! Set the coefficient for the vapour-dependence
      ! C_q = kq / ( 1 + kq L_con dqsat/dT / ( kt_cp )  )
      do ic = 1, n_points
        coefs_cond(ic,i_q,i_liq) = kq_cond(ic,i_liq)                           &
                   / ( one + kq_cond(ic,i_liq) * L_con(ic)                     &
                             * dqsatdT_liq(ic) / kt_cp(ic) )
      end do
      ! Set the coefficient for temperature dependence
      ! C_T = dqsat/dT C_q ( -1 + B / kt_cp  )
      do ic = 1, n_points
        coefs_cond(ic,i_t,i_liq) = dqsatdT_liq(ic)                             &
          * coefs_cond(ic,i_q,i_liq) * ( -one + fall_term(ic) )
      end do
      ! Set the explicit coefficient
      ! C_0 = dqsat/dT C_q  B/kt_cp  ( ref_temp - prev_temp )
      do ic = 1, n_points
        coefs_cond(ic,i_0,i_liq) = dqsatdT_liq(ic)                             &
          * coefs_cond(ic,i_q,i_liq)                                           &
          * fall_term(ic) * ( ref_temp(ic) - prev_temp(ic) )
      end do

      ! Add contribution to vapour and temperature coefficients
      do i_coef = 1, n_coefs
        do ic = 1, n_points
          coefs_q_vap(ic,i_coef) = coefs_q_vap(ic,i_coef)                      &
                                 + coefs_cond(ic,i_coef,i_liq)
        end do
      end do
      do i_coef = 1, n_coefs
        do ic = 1, n_points
          coefs_temp(ic,i_coef)  = coefs_temp(ic,i_coef)                       &
                                 + coefs_cond(ic,i_coef,i_liq)                 &
                                   * lrcp(ic)
        end do
      end do

      ! If calculation only needed at a minority of points
    else
      ! Indirect indexing version of the above calculations

      do ic2 = 1, nc(i_liq)
        ic = index_ic(ic2,i_liq)

        ! Calculate kt_cp and fall heat transport terms
        kt_cp(ic) = kt_cond(ic,i_liq)                                          &
                    * ( cp_dry + cp_vap * q_vap(ic) )
        kt_cp(ic) = max( kt_cp(ic), sqrt_min_float )
        fall_term(ic) = ( q_loc_cond(ic,i_liq) * cp_liq                        &
            * ( wind_w(ic) - wf_cond(ic,i_liq) ) / delta_z(ic) )               &
                      / ( kt_cp(ic) / delta_t(ic) )
        fall_term(ic) = min( fall_term(ic), one )

        ! Calculate condensation coefficients for current species
        coefs_cond(ic,i_q,i_liq) = kq_cond(ic,i_liq)                           &
                   / ( one + kq_cond(ic,i_liq) * L_con(ic)                     &
                             * dqsatdT_liq(ic) / kt_cp(ic) )
        coefs_cond(ic,i_t,i_liq) = dqsatdT_liq(ic)                             &
          * coefs_cond(ic,i_q,i_liq) * ( -one + fall_term(ic) )
        coefs_cond(ic,i_0,i_liq) = dqsatdT_liq(ic)                             &
          * coefs_cond(ic,i_q,i_liq)                                           &
          * fall_term(ic) * ( ref_temp(ic) - prev_temp(ic) )
      end do

      ! Add contribution to vapour and temperature coefficients
      do i_coef = 1, n_coefs
        do ic2 = 1, nc(i_liq)
          ic = index_ic(ic2,i_liq)
          coefs_q_vap(ic,i_coef) = coefs_q_vap(ic,i_coef)                      &
                                 + coefs_cond(ic,i_coef,i_liq)
          coefs_temp(ic,i_coef)  = coefs_temp(ic,i_coef)                       &
                                 + coefs_cond(ic,i_coef,i_liq)                 &
                                   * lrcp(ic)
        end do
      end do

    end if
  end if  ! ( nc(i_liq) > 0 )
end do  ! i_liq = 1, n_cond_species_liq



! Store L_sub / cp_tot and L_fus / cp_tot
do ic = 1, n_points
  lrcp(ic) = L_sub(ic) / cp_tot(ic)
end do

! Loop over ice species
do i_ice = n_cond_species_liq+1, n_cond_species
  ! If any points active for this species
  if ( nc(i_ice) > 0 ) then

    ! If calculations required at majority of points
    if ( l_full_do(i_ice) ) then
      ! Full-field calculation...

      ! Store the term  kt_cp =
      ! kt (cp_dry + cp_vap q_vap)
      do ic = 1, n_points
        kt_cp(ic) = kt_cond(ic,i_ice)                                          &
                    * ( cp_dry + cp_vap * q_vap(ic) )
      end do
      ! Don't let kt_cp go to zero (gives NaNs due to div-by-zero
      ! at points where there is no condensate)
      do ic = 1, n_points
        kt_cp(ic) = max( kt_cp(ic), sqrt_min_float )
      end do
      ! Store the term  B/kt_cp =
      ! ( q_cond cp_cond ( wind_w - wf_cond ) / delta_z ) / kt_cp
      do ic = 1, n_points
        fall_term(ic) = ( q_loc_cond(ic,i_ice) * cp_ice                        &
            * ( wind_w(ic) - wf_cond(ic,i_ice) ) / delta_z(ic) )               &
                      / ( kt_cp(ic) / delta_t(ic) )
        ! Note: removing factor of delta_t included in kt_cond
      end do
      ! Don't let the fall term exceed one
      do ic = 1, n_points
        fall_term(ic) = min( fall_term(ic), one )
      end do

      ! Calculate deposition coefficients for current species:
      ! Set the coefficient for the vapour-dependence
      ! C_q = kq / ( 1 + kq L_con dqsat/dT / ( kt_cp )  )
      do ic = 1, n_points
        coefs_cond(ic,i_q,i_ice) = kq_cond(ic,i_ice)                           &
                   / ( one + kq_cond(ic,i_ice) * L_sub(ic)                     &
                             * dqsatdT_ice(ic) / kt_cp(ic) )
      end do
      ! Set the coefficient for temperature dependence
      ! C_T = dqsat/dT C_q ( -1 + B / kt_cp  )
      do ic = 1, n_points
        coefs_cond(ic,i_t,i_ice) = dqsatdT_ice(ic)                             &
          * coefs_cond(ic,i_q,i_ice) * ( -one + fall_term(ic) )
      end do
      ! Set the explicit coefficient
      ! C_0 = dqsat/dT C_q ( B/kt_cp  ( ref_temp - prev_temp )
      !                    - L_fus dqdt_frz / kt_cp
      !                    - ( qsat_ice - qsat_liq ) / dqsat/dT )
      do ic = 1, n_points
        coefs_cond(ic,i_0,i_ice) = dqsatdT_ice(ic)                             &
          * coefs_cond(ic,i_q,i_ice) * (                                       &
              fall_term(ic) * ( ref_temp(ic) - prev_temp(ic) )                 &
            - L_fus(ic) * ( dq_frz_cond(ic,i_ice) / kt_cp(ic) )                &
            - ( qsat_ice_ref(ic) - qsat_liq_ref(ic) )                          &
              / dqsatdT_ice(ic)        )
      end do

      ! Calculate melting coefficients for current species:
      ! Vapour dependence
      do ic = 1, n_points
        coefs_melt(ic,i_q,i_ice) = kq_cond(ic,i_ice)                           &
                                   * L_sub(ic) / L_fus(ic)
      end do
      ! Temperature dependence
      do ic = 1, n_points
        coefs_melt(ic,i_t,i_ice) = kt_cond(ic,i_ice)                           &
                   * ( cp_dry + cp_vap * q_vap(ic) ) / L_fus(ic)
      end do
      ! Explicit term
      do ic = 1, n_points
        coefs_melt(ic,i_0,i_ice) = ( coefs_melt(ic,i_t,i_ice)                  &
                 + coefs_melt(ic,i_q,i_ice) * dqsatdT_ice(ic) )                &
                 * ( ref_temp(ic) - melt_temp )                                &
               - coefs_melt(ic,i_q,i_ice)                                      &
                 * ( qsat_ice_ref(ic) - qsat_liq_ref(ic) )                     &
               + dq_frz_cond(ic,i_ice)
      end do

      ! Calculate alternative deposition coefficients for use
      ! when the ice is melting
      do ic = 1, n_points
        coefs_cond_m(ic,i_q,i_ice) = kq_cond(ic,i_ice)
      end do
      do ic = 1, n_points
        coefs_cond_m(ic,i_t,i_ice) = zero
      end do
      do ic = 1, n_points
        coefs_cond_m(ic,i_0,i_ice) = -kq_cond(ic,i_ice)                        &
              * ( qsat_ice_ref(ic) - qsat_liq_ref(ic)                          &
                + dqsatdT_ice(ic) * (melt_temp - ref_temp(ic)) )
      end do

      ! Add contribution to vapour and temperature coefficients,
      ! assuming that melting is not occuring
      do i_coef = 1, n_coefs
        do ic = 1, n_points
          coefs_q_vap(ic,i_coef) = coefs_q_vap(ic,i_coef)                      &
                                 + coefs_cond(ic,i_coef,i_ice)
        end do
      end do
      do i_coef = 1, n_coefs
        do ic = 1, n_points
          coefs_temp(ic,i_coef)  = coefs_temp(ic,i_coef)                       &
                                 + coefs_cond(ic,i_coef,i_ice)                 &
                                   * lrcp(ic)
        end do
      end do

      ! If calculation only needed at a minority of points
    else
      ! Indirect indexing version of the above calculations

      do ic2 = 1, nc(i_ice)
        ic = index_ic(ic2,i_ice)

        ! Calculate kt_cp and fall heat transport terms
        kt_cp(ic) = kt_cond(ic,i_ice)                                          &
                    * ( cp_dry + cp_vap * q_vap(ic) )
        kt_cp(ic) = max( kt_cp(ic), sqrt_min_float )
        fall_term(ic) = ( q_loc_cond(ic,i_ice) * cp_ice                        &
            * ( wind_w(ic) - wf_cond(ic,i_ice) ) / delta_z(ic) )               &
                      / ( kt_cp(ic) / delta_t(ic) )
        fall_term(ic) = min( fall_term(ic), one )

        ! Calculate deposition coefficients for current species:
        coefs_cond(ic,i_q,i_ice) = kq_cond(ic,i_ice)                           &
                   / ( one + kq_cond(ic,i_ice) * L_sub(ic)                     &
                             * dqsatdT_ice(ic) / kt_cp(ic) )
        coefs_cond(ic,i_t,i_ice) = dqsatdT_ice(ic)                             &
          * coefs_cond(ic,i_q,i_ice) * ( -one + fall_term(ic) )
        coefs_cond(ic,i_0,i_ice) = dqsatdT_ice(ic)                             &
          * coefs_cond(ic,i_q,i_ice) * (                                       &
              fall_term(ic) * ( ref_temp(ic) - prev_temp(ic) )                 &
            - L_fus(ic) * ( dq_frz_cond(ic,i_ice) / kt_cp(ic) )                &
            - ( qsat_ice_ref(ic) - qsat_liq_ref(ic) )                          &
              / dqsatdT_ice(ic)        )

        ! Calculate melting coefficients for current species:
        coefs_melt(ic,i_q,i_ice) = kq_cond(ic,i_ice)                           &
                                   * L_sub(ic) / L_fus(ic)
        coefs_melt(ic,i_t,i_ice) = kt_cond(ic,i_ice)                           &
                   * ( cp_dry + cp_vap * q_vap(ic) ) / L_fus(ic)
        coefs_melt(ic,i_0,i_ice) = ( coefs_melt(ic,i_t,i_ice)                  &
                 + coefs_melt(ic,i_q,i_ice) * dqsatdT_ice(ic) )                &
                 * ( ref_temp(ic) - melt_temp )                                &
               - coefs_melt(ic,i_q,i_ice)                                      &
                 * ( qsat_ice_ref(ic) - qsat_liq_ref(ic) )                     &
               + dq_frz_cond(ic,i_ice)

        ! Calculate alternative deposition coefficients for use
        ! when the ice is melting
        coefs_cond_m(ic,i_q,i_ice) = kq_cond(ic,i_ice)
        coefs_cond_m(ic,i_t,i_ice) = zero
        coefs_cond_m(ic,i_0,i_ice) = -kq_cond(ic,i_ice)                        &
              * ( qsat_ice_ref(ic) - qsat_liq_ref(ic)                          &
                + dqsatdT_ice(ic) * (melt_temp - ref_temp(ic)) )

      end do

      ! Add contribution to vapour and temperature coefficients,
      ! assuming that melting is not occuring
      do i_coef = 1, n_coefs
        do ic2 = 1, nc(i_ice)
          ic = index_ic(ic2,i_ice)
          coefs_q_vap(ic,i_coef) = coefs_q_vap(ic,i_coef)                      &
                                 + coefs_cond(ic,i_coef,i_ice)
          coefs_temp(ic,i_coef)  = coefs_temp(ic,i_coef)                       &
                                 + coefs_cond(ic,i_coef,i_ice)                 &
                                   * lrcp(ic)
        end do
      end do

    end if
  end if  ! ( nc(i_ice) > 0 )
end do  ! i_ice = n_cond_species_liq+1, n_cond_species


return
end subroutine calc_phase_change_coefs

end module calc_phase_change_coefs_mod
