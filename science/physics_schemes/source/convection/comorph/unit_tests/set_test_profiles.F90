! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

#if !defined(LFRIC)
module set_test_profiles_mod

implicit none

contains

! Subroutine calculates a set of idealised atmospheric profiles
! for running various tests.
! This sets up a 2-D grid of columns, where RH varies linearly
! along one of the horizontal dimensions, and the free tropospheric
! stability varies independently along the other.
subroutine set_test_profiles( l_tracer, grid, turb, cloudfracs,                &
                              fields )

use comorph_constants_mod, only: real_hmprec, real_cvprec,                     &
                     nx_full, ny_full, k_top_conv, k_top_init, n_tracers,      &
                     l_cv_rain, l_cv_cf, l_cv_snow, l_cv_graup,                &
                     l_turb_par_gen,                                           &
                     pi_cv=>pi,                                                &
                     gravity_cv=>gravity, cp_dry_cv=>cp_dry,                   &
                     R_dry_cv=>R_dry, R_vap_cv=>R_vap

use grid_type_mod, only: grid_type
use turb_type_mod, only: turb_type
use cloudfracs_type_mod, only: cloudfracs_type
use fields_type_mod, only: fields_type, fields_list_type
use set_qsat_mod, only: set_qsat_liq
use calc_rho_dry_mod, only: calc_rho_dry

implicit none

! Flag for whether to initialise tracers
logical, intent(in) :: l_tracer

! Structure containing pointers to model-level
! heights and pressures
type(grid_type), intent(in out) :: grid

! Structure containing pointers to turbulence fields passed in
type(turb_type), intent(in out) :: turb

! Structure containing pointers to diagnostic cloud fractions
! and rain fraction
type(cloudfracs_type), intent(in out) :: cloudfracs

! Structure containing pointers to primary fields
type(fields_type), intent(in out) :: fields


real(kind=real_hmprec), parameter :: zero = 0.0_real_hmprec
real(kind=real_hmprec), parameter :: half = 0.5_real_hmprec
real(kind=real_hmprec), parameter :: one = 1.0_real_hmprec
real(kind=real_hmprec), parameter :: two = 2.0_real_hmprec
real(kind=real_hmprec), parameter :: three = 3.0_real_hmprec
real(kind=real_hmprec), parameter :: third = one/three

! Constants converted to host-model precision
real(kind=real_hmprec) :: pi
real(kind=real_hmprec) :: gravity
real(kind=real_hmprec) :: cp_dry
real(kind=real_hmprec) :: R_dry
real(kind=real_hmprec) :: R_vap

! Top of model height
real(kind=real_hmprec), parameter :: z_top = 20000.0_real_hmprec
! Other heights used to set up basic state:
!   - Top of surface statically-unstable layer:
real(kind=real_hmprec), parameter :: z_sfc = 200.0_real_hmprec
!   - Top of the boundary-layer:
real(kind=real_hmprec), parameter :: z_pbl = 1000.0_real_hmprec
!   - Tropopause height:
real(kind=real_hmprec), parameter :: z_tpp = 12000.0_real_hmprec
!   - Height-scale used to set an idealised wind profile:
real(kind=real_hmprec), parameter :: z_scl = 500.0_real_hmprec

! Power in made-up relation between water content and
! cloud-fraction
real(kind=real_cvprec), parameter :: powr = 0.2

! Work variables
real(kind=real_hmprec) :: lapse_rate_dry
real(kind=real_hmprec) :: lapse_rate
real(kind=real_hmprec) :: interp
real(kind=real_hmprec) :: tau_bl

real(kind=real_cvprec) :: work_temp(nx_full)
real(kind=real_cvprec) :: work_pres(nx_full)
real(kind=real_cvprec) :: work_qs(nx_full)

! Virtual temperature profile
real(kind=real_hmprec) :: virt_temp                                            &
    ( nx_full, ny_full, 0:k_top_conv )

! Turbulent diffusivity
real(kind=real_hmprec) :: diffusivity                                          &
    ( nx_full, ny_full, 1:k_top_init+1 )

! Loop counters
integer :: i, j, k, i_field, n


! Convert constants to host-model precision
pi      = real(pi_cv     ,real_hmprec)
gravity = real(gravity_cv,real_hmprec)
cp_dry  = real(cp_dry_cv, real_hmprec)
R_dry   = real(R_dry_cv,  real_hmprec)
R_vap   = real(R_vap_cv,  real_hmprec)

do j = 1, ny_full
  do i = 1, nx_full
    ! Set surface radius from Earth centre
    grid % r_surf(i,j) = 6.37E6_real_hmprec
  end do
end do

! Set model-level heights:
! heights increase as a quadratic formula of model-level index.
do k = 0, k_top_conv
  do j = 1, ny_full
    do i = 1, nx_full
      ! Full levels use function of k
      grid % height_full(i,j,k) = z_top                                        &
              * ( real(k,real_hmprec)                                          &
                / real(k_top_conv,real_hmprec) )**2
    end do
  end do
end do
do k = 1, k_top_conv+1
  do j = 1, ny_full
    do i = 1, nx_full
      ! Half-levels use same function but with k-1/2
      grid % height_half(i,j,k) = z_top                                        &
              * ( (real(k,real_hmprec)-half)                                   &
                /  real(k_top_conv,real_hmprec) )**2
    end do
  end do
end do

! Set made-up wind profiles
do k = 1, k_top_conv
  do j = 1, ny_full
    do i = 1, nx_full

      ! u: increases logarithmically with height from zero at the surface:
      fields % wind_u(i,j,k) = 20.0_real_hmprec                                &
        * log((grid % height_full(i,j,k)+z_scl)/z_scl)

      ! v,w: sinusoidal profile with height within the troposphere,
      ! zero in the stratosphere:
      fields % wind_v(i,j,k) = 10.0_real_hmprec                                &
        * max( sin( (half/pi) * grid % height_full(i,j,k)/z_tpp ),             &
               zero )
      fields % wind_w(i,j,k) = 0.2_real_hmprec                                 &
        * max( sin( (half/pi) * grid % height_full(i,j,k)/z_tpp ),             &
               zero )

    end do
  end do
end do


! Virtual temperature profile...
! We set this so that the stability of the free troposphere
! varies along the j-direction.

! Set surface value of 300K at all points
do j = 1, ny_full
  do i = 1, nx_full
    virt_temp(i,j,0) = 300.0_real_hmprec
  end do
end do
! Set dry-adiabatic lapse rate
lapse_rate_dry = gravity/cp_dry
! Loop over levels...
do k = 1, k_top_conv
  do j = 1, ny_full
    do i = 1, nx_full

      ! Rescale the lapse-rate at different heights...
      if ( grid % height_full(i,j,k) < z_sfc ) then
        ! Statically unstable (super-adiabatic) near surface
        lapse_rate = lapse_rate_dry * 1.1_real_hmprec
      else if ( grid % height_full(i,j,k) >= z_sfc .and.                       &
                grid % height_full(i,j,k) <  z_pbl ) then
        ! Near-neutral (but slightly stable) within PBL
        lapse_rate = lapse_rate_dry * 0.99_real_hmprec
      else if ( grid % height_full(i,j,k) >= z_pbl .and.                       &
                grid % height_full(i,j,k) <  z_tpp ) then
        ! Moderately stable in free troposphere,
        ! varying linearly in the j-direction
        lapse_rate = lapse_rate_dry * ( 0.6_real_hmprec                        &
                                      + 0.2_real_hmprec                        &
                          * real(j,real_hmprec)                                &
                          / real(ny_full,real_hmprec) )
      else if ( grid % height_full(i,j,k) > z_tpp ) then
        ! Isothermal (very stable) in stratosphere
        lapse_rate = zero
      end if

      ! Set next level Tv based on the calculated lapse rate with height:
      virt_temp(i,j,k) = virt_temp(i,j,k-1)                                    &
                 - lapse_rate * ( grid % height_full(i,j,k)                    &
                                - grid % height_full(i,j,k-1) )

    end do
  end do
end do


! Set pressure approximately in hydrostatic balance:
! dp/dz = -rho g = -p/(R Tv) g
do j = 1, ny_full
  do i = 1, nx_full
    ! Surface pressure set to 1000 hPa at all points
    grid % pressure_full(i,j,0) = 100000.0_real_hmprec
    ! First half-level; integrate hydrostatic balance using
    ! surface value of Tv and p
    grid % pressure_half(i,j,1) = grid % pressure_full(i,j,0)                  &
      - ( grid%height_half(i,j,1) - grid%height_full(i,j,0) )                  &
        * gravity * grid % pressure_full(i,j,0)                                &
        / ( R_dry * virt_temp(i,j,0) )
  end do
end do
do k = 2, k_top_conv+1
  do j = 1, ny_full
    do i = 1, nx_full
      ! Integrate hydrostatic balance over all remaining
      ! half-levels, using Tv from the intervening full-levels
      ! (very approximate as using pressure from the previous
      !  half-level).
      grid % pressure_half(i,j,k) = grid % pressure_half(i,j,k-1)              &
        - ( grid%height_half(i,j,k) - grid%height_half(i,j,k-1) )              &
          * gravity * grid % pressure_half(i,j,k-1)                            &
          / ( R_dry * virt_temp(i,j,k-1) )
    end do
  end do
end do
do k = 1, k_top_conv
  do j = 1, ny_full
    do i = 1, nx_full
      ! Interpolate between the half-levels to get pressure on full levels.
      interp                                                                   &
        = ( grid%height_full(i,j,k) - grid%height_half(i,j,k) )                &
        / ( grid%height_half(i,j,k+1) - grid%height_half(i,j,k) )
      grid % pressure_full(i,j,k)                                              &
             = (one-interp) * grid % pressure_half(i,j,k)                      &
             +      interp  * grid % pressure_half(i,j,k+1)
    end do
  end do
end do


! Initialise temperature equal to virtual temperature
do k = 1, k_top_conv
  do j = 1, ny_full
    do i = 1, nx_full
      fields % temperature(i,j,k) = virt_temp(i,j,k)
    end do
  end do
end do

! Iterate to get T,q consistent with specified Tv and RH.
! We set the RH of the profile to vary along the i-direction.
do n = 1, 10

  ! Set q_vap to varying % of saturation
  do k = 1, k_top_conv
    do j = 1, ny_full
      ! Need to make copies to convert to native precision
      ! used inside set_qsat_liq
      do i = 1, nx_full
        work_temp(i) =real(fields%temperature(i,j,k),real_cvprec)
        work_pres(i) =real(grid%pressure_full(i,j,k),real_cvprec)
      end do
      call set_qsat_liq( nx_full, work_temp, work_pres, work_qs )
      do i = 1, nx_full
        ! Set RH between 60% and 70%, varying in the i-direction
        fields % q_vap(i,j,k) = real(work_qs(i),real_hmprec)                   &
            * ( 0.6_real_hmprec + 0.1_real_hmprec                              &
                          * real(i,real_hmprec)                                &
                          / real(nx_full,real_hmprec) )
      end do
    end do
  end do

  ! Adjust temperatures to maintain the virtual temperature
  ! equal to the initially set profile
  ! (means that the set hydrostatic pressures remain valid)
  do k = 1, k_top_conv
    do j = 1, ny_full
      do i = 1, nx_full
        ! Inverting Tv = T ( 1 + Rv/Rd qv ) / ( 1 + qt ),
        ! currently assuming qt = qv (no condensate).
        fields % temperature(i,j,k) = virt_temp(i,j,k)                         &
          * ( one + fields % q_vap(i,j,k) )                                    &
        /   ( one + (R_vap/R_dry) * fields % q_vap(i,j,k) )
      end do
    end do
  end do

end do  ! n = 1, 10


! Compute dry-density
do k = 1, k_top_conv
  do j = 1, ny_full
    do i = 1, nx_full
      ! rho_d = p / ( Rd T ( 1 + Rv/Rd qv )
      grid % rho_dry(i,j,k) = grid % pressure_full(i,j,k)                      &
          / ( R_dry * fields % temperature(i,j,k)                              &
              * ( one + (R_vap/R_dry) * fields % q_vap(i,j,k) ) )
    end do
  end do
end do


! Set condensed water species.
! Where the Tv profile is moist-unstable, elevated convection is
! expected to trigger from the layer-cloud.  Want to make sure that
! convective triggering from some points but not others works correctly,
! so setup to have layer-cloud at some points but zero q_cl, q_cf etc
! at other points.

do k = 1, k_top_conv
  do j = 1, ny_full
    do i = 1, nx_full
      ! q_cl; set so that only the moistest half of the
      ! points have any liquid-cloud (i.e. where qv(i,j) > qv(nx/2,j),
      ! since qv was set earlier to increase linearly in the i-direction).
      ! Height-variation is quadratic, with a maximum at 2* the BL-height.
      fields % q_cl(i,j,k)                                                     &
         = 0.1_real_hmprec                                                     &
         * max( fields % q_vap(i,j,k)                                          &
              - fields % q_vap(nx_full/2,j,k), zero )                          &
         * max( one - ( ( grid % height_full(i,j,k)                            &
                        - 2.0_real_hmprec*z_pbl )/z_pbl )**2, zero )
    end do
  end do
end do
if ( l_cv_rain ) then
  do k = 1, k_top_conv
    do j = 1, ny_full
      do i = 1, nx_full
        ! q_rain; set so that only the moistest third of the
        ! points have any rain.
        ! Height-variation is quadratic, with a maximum at the surface,
        ! extending up to 5* the BL-height.
        fields % q_rain(i,j,k)                                                 &
         = 0.01_real_hmprec                                                    &
         * max( fields % q_vap(i,j,k)                                          &
              - fields % q_vap(2*nx_full/3,j,k), zero )                        &
         * max( one - ( grid % height_full(i,j,k)                              &
                      / ( 5.0_real_hmprec*z_pbl ) )**2, zero )
      end do
    end do
  end do
end if
if ( l_cv_cf ) then
  do k = 1, k_top_conv
    do j = 1, ny_full
      do i = 1, nx_full
        ! q_cf; set so that only the moistest two-thirds of the
        ! points have any ice-cloud.
        ! Height-variation is quadratic, with a maximum at 8* the BL-height.
        fields % q_cf(i,j,k)                                                   &
         = 0.1_real_hmprec                                                     &
         * max( fields % q_vap(i,j,k)                                          &
              - fields % q_vap(nx_full/3,j,k), zero )                          &
         * max( one - ( ( grid % height_full(i,j,k)                            &
                        - 8.0_real_hmprec*z_pbl )                              &
                      / ( 3.0_real_hmprec*z_pbl ) )**2, zero )
      end do
    end do
  end do
end if
! Set snow and graupel to zero
if ( l_cv_snow ) then
  do k = 1, k_top_conv
    do j = 1, ny_full
      do i = 1, nx_full
        fields % q_snow(i,j,k) = zero
      end do
    end do
  end do
end if
if ( l_cv_graup ) then
  do k = 1, k_top_conv
    do j = 1, ny_full
      do i = 1, nx_full
        fields % q_graup(i,j,k) = zero
      end do
    end do
  end do
end if

! Set cloud-fractions.
! Liquid and ice cloud-fraction are set as made-up power-law functions
! of the grid-mean liquid and ice cloud condensate mixing-ratios.
! Bulk cloud fraction is set assuming random overlap between the
! liquid and ice cloud fractions.
do k = 1, k_top_conv
  do j = 1, ny_full
    do i = 1, nx_full
      fields % cf_liq(i,j,k) = fields % q_cl(i,j,k)**powr
      fields % cf_ice(i,j,k) = (10_real_hmprec                                 &
                              *fields % q_cf(i,j,k))**powr
      fields % cf_bulk(i,j,k) = fields % cf_liq(i,j,k)                         &
                              + fields % cf_ice(i,j,k)                         &
          - fields % cf_liq(i,j,k) * fields % cf_ice(i,j,k)
      cloudfracs % frac_precip(i,j,k) = ( 30_real_hmprec                       &
                                    * ( fields % q_rain(i,j,k)                 &
                                      + fields % q_graup(i,j,k) ) )**powr
    end do
  end do
end do
! Initialise convective cloud variables to zero
do k = 1, k_top_conv
  do j = 1, ny_full
    do i = 1, nx_full
      cloudfracs % frac_liq_conv(i,j,k) = zero
      cloudfracs % frac_bulk_conv(i,j,k) = zero
      cloudfracs % q_cl_conv(i,j,k) = zero
    end do
  end do
end do


! Set passive tracers if requested
if ( l_tracer ) then

  if ( n_tracers >= 1 ) then
    ! First tracer set to 1 in the boundary-layer, 0 above
    do k = 1, k_top_conv
      do j = 1, ny_full
        do i = 1, nx_full
          if ( grid % height_full(i,j,k) < z_pbl ) then
            fields % tracers(1)%pt(i,j,k) = one
          else
            fields % tracers(1)%pt(i,j,k) = zero
          end if
        end do
      end do
    end do
  end if

  if ( n_tracers >= 2 ) then
    ! 2nd tracer set as a height tracer
    do k = 1, k_top_conv
      do j = 1, ny_full
        do i = 1, nx_full
          fields % tracers(2)%pt(i,j,k)                                        &
            = grid % height_full(i,j,k)
        end do
      end do
    end do
  end if

  if ( n_tracers > 2 ) then
    ! Set any other tracers to zero for now
    do i_field = 3, n_tracers
      do k = 1, k_top_conv
        do j = 1, ny_full
          do i = 1, nx_full
            fields % tracers(i_field)%pt(i,j,k) = zero
          end do
        end do
      end do
    end do
  end if

end if  ! ( l_tracer )



! Set turbulence fields if requested
if ( l_turb_par_gen ) then

  do k = 1, k_top_init+1
    do j = 1, ny_full
      do i = 1, nx_full

        ! Momentum diffusivity: minimal value except in the
        ! boundary-layer where we have a quadratic variation
        ! with height:
        diffusivity(i,j,k) = 1.0_real_hmprec
        if ( grid % height_half(i,j,k) < z_pbl ) then
          diffusivity(i,j,k) = diffusivity(i,j,k)                              &
            + 50.0_real_hmprec *                                               &
              ( one - ( two * grid % height_half(i,j,k) / z_pbl                &
                      - one )**2 )
        end if

        ! Estimate turbulence timescale
        tau_bl = 300.0_real_hmprec                                             &
          * ( diffusivity(i,j,k) / 50.0_real_hmprec )**third
        ! Use timescale to estimate TKE
        turb % w_var(i,j,k) = diffusivity(i,j,k) / tau_bl

      end do
    end do
  end do

  do k = 2, k_top_init
    do j = 1, ny_full
      do i = 1, nx_full
        ! Fluxes based on local gradients
        turb % f_templ(i,j,k) = -diffusivity(i,j,k)                            &
          * ( ( fields % temperature(i,j,k)                                    &
              - fields % temperature(i,j,k-1) )                                &
            / ( grid % height_full(i,j,k)                                      &
              - grid % height_full(i,j,k-1) ) + lapse_rate_dry)
        turb % f_q_tot(i,j,k) = -diffusivity(i,j,k)                            &
          * ( ( fields % q_vap(i,j,k)                                          &
              - fields % q_vap(i,j,k-1) )                                      &
            / ( grid % height_full(i,j,k)                                      &
              - grid % height_full(i,j,k-1) ) )
        turb % f_wind_u(i,j,k) = -diffusivity(i,j,k)                           &
          * ( ( fields % wind_u(i,j,k)                                         &
              - fields % wind_u(i,j,k-1) )                                     &
            / ( grid % height_full(i,j,k)                                      &
              - grid % height_full(i,j,k-1) ) )
        turb % f_wind_v(i,j,k) = -diffusivity(i,j,k)                           &
          * ( ( fields % wind_v(i,j,k)                                         &
              - fields % wind_v(i,j,k-1) )                                     &
            / ( grid % height_full(i,j,k)                                      &
              - grid % height_full(i,j,k-1) ) )
      end do
    end do
  end do
  do j = 1, ny_full
    do i = 1, nx_full
      ! Copy model-level 2 to get surface fluxes
      turb % f_templ(i,j,1) = turb % f_templ(i,j,2)
      turb % f_q_tot(i,j,1) = turb % f_q_tot(i,j,2)
      turb % f_wind_u(i,j,1) = turb % f_wind_u(i,j,2)
      turb % f_wind_v(i,j,1) = turb % f_wind_v(i,j,2)
      ! Ditto the upper model-level interface of the top initiation level
      turb % f_templ(i,j,k_top_init+1) = turb % f_templ(i,j,k_top_init)
      turb % f_q_tot(i,j,k_top_init+1) = turb % f_q_tot(i,j,k_top_init)
      turb % f_wind_u(i,j,k_top_init+1) = turb % f_wind_u(i,j,k_top_init)
      turb % f_wind_v(i,j,k_top_init+1) = turb % f_wind_v(i,j,k_top_init)
    end do
  end do

  do k = 1, k_top_init
    do j = 1, ny_full
      do i = 1, nx_full
        ! Set turbulent length-scale on theta-levels
        interp                                                                 &
          = ( grid%height_full(i,j,k) - grid%height_half(i,j,k) )              &
          / ( grid%height_half(i,j,k+1) - grid%height_half(i,j,k) )
        turb % lengthscale(i,j,k)                                              &
          = (one-interp) * diffusivity(i,j,k)   / sqrt( turb % w_var(i,j,k) )  &
          +      interp  * diffusivity(i,j,k+1) / sqrt( turb % w_var(i,j,k+1) )
      end do
    end do
  end do

else  ! ( .NOT. l_turb_par_gen )

  ! Set turbulence length-scale to a constant
  do k = 1, k_top_init
    do j = 1, ny_full
      do i = 1, nx_full
        turb % lengthscale(i,j,k) = 0.5 * z_pbl
      end do
    end do
  end do

end if  ! ( .NOT. l_turb_par_gen )


do j = 1, ny_full
  do i = 1, nx_full
    ! Set boundary-layer top height
    turb % z_bl_top(i,j) = z_pbl
    ! Set parcel radius amplification factor
    turb % par_radius_amp(i,j) = 1.0_real_cvprec
  end do
end do



return
end subroutine set_test_profiles

end module set_test_profiles_mod
#endif
