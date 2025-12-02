! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_turb_len_mod

use um_types, only: real_umphys

implicit none

contains

! Routine to calculate a turbulence lengthscale and
! apply sundry optional scalings to par_radius via par_radius_amp...
subroutine calc_turb_len( zh_eff, z_theta, z_rho, rho_wet_th, m_v,             &
                          rhokm, bl_w_var, ls_rain, ls_snow, w,                &
                          delta_x, delta_y,                                    &
                          turb_len, par_radius_amp_um )

use atm_fields_bounds_mod, only: pdims, tdims, tdims_s, wdims_s
use nlsizes_namelist_mod, only: bl_levels
use cv_run_mod, only: w_cape_limit
use comorph_um_namelist_mod, only:                                             &
                      par_radius_init_method,                                  &
                      rain_dependence, qfacrain_dependence, w_dependence,      &
                      linear_qfacrain_dep,                                     &
                      par_radius_knob, par_radius_knob_max, par_radius_ppn_max,&
                      l_resdep_precipramp, dx_ref
use cv_param_mod, only: refqsat

implicit none

! Effective boundary-layer tope height
real(kind=real_umphys), intent(in) ::                                          &
                    zh_eff         ( pdims%i_start:pdims%i_end,                &
                                     pdims%j_start:pdims%j_end )

! Model-level heights above surface
real(kind=real_umphys), intent(in) ::                                          &
                    z_theta        ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:tdims%k_end )
real(kind=real_umphys), intent(in) ::                                          &
                    z_rho          ( pdims%i_start:pdims%i_end,                &
                                     pdims%j_start:pdims%j_end,                &
                                     pdims%k_start:pdims%k_end )

! Wet-density on theta-levels
real(kind=real_umphys), intent(in) ::                                          &
                    rho_wet_th     ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:tdims%k_end-1 )

! Start-of-timestep water-vapour mixing-ratio
real(kind=real_umphys), intent(in) ::                                          &
                    m_v            ( tdims_s%i_start:tdims_s%i_end,            &
                                     tdims_s%j_start:tdims_s%j_end,            &
                                     tdims_s%k_start:tdims_s%k_end )

! rho * turbulent diffusivity for momentum, on theta-grid
real(kind=real_umphys), intent(in) ::                                          &
                    rhokm    ( tdims_s%i_start:tdims_s%i_end,                  &
                               tdims_s%j_start:tdims_s%j_end,                  &
                               0:bl_levels-1 )
! Sub-grid turbulent variance in vertical velocity
real(kind=real_umphys), intent(in) ::                                          &
                    bl_w_var ( tdims%i_start:tdims%i_end,                      &
                               tdims%j_start:tdims%j_end,                      &
                               1:tdims%k_end )

! Surface precipitation rates from the microphysics scheme
real(kind=real_umphys), intent(in) ::                                          &
                    ls_rain        ( pdims%i_start:pdims%i_end,                &
                                     pdims%j_start:pdims%j_end )
real(kind=real_umphys), intent(in) ::                                          &
                    ls_snow        ( pdims%i_start:pdims%i_end,                &
                                     pdims%j_start:pdims%j_end )

! Start-of-timestep vertical velocity
real(kind=real_umphys), intent(in) ::                                          &
                    w              ( wdims_s%i_start:wdims_s%i_end,            &
                                     wdims_s%j_start:wdims_s%j_end,            &
                                     wdims_s%k_start:wdims_s%k_end )

! Horizontal grid-lengths
real(kind=real_umphys), intent(in) :: delta_x ( pdims%i_start:pdims%i_end,     &
                                                pdims%j_start:pdims%j_end )
real(kind=real_umphys), intent(in) :: delta_y ( pdims%i_start:pdims%i_end,     &
                                                pdims%j_start:pdims%j_end )

! Turbulence lengthscale
real(kind=real_umphys), intent(out) ::                                         &
                     turb_len      ( pdims%i_start:pdims%i_end,                &
                                     pdims%j_start:pdims%j_end,                &
                                     1:bl_levels )

! Dimensionless amplification factor applied to the initial
! parcel radius inside comorph
real(kind=real_umphys), intent(out) ::                                         &
                     par_radius_amp_um ( pdims%i_start:pdims%i_end,            &
                                         pdims%j_start:pdims%j_end )

! Precipitation rate normalised by mean mixing-ratio
real(kind=real_umphys) :: ppnrate

! Scaling factor from the surface precip rate
real(kind=real_umphys) ::                                                      &
        rainfac                    ( pdims%i_start:pdims%i_end,                &
                                     pdims%j_start:pdims%j_end )
! Normalisation applied to the above (1/mean m_v)
real(kind=real_umphys) ::                                                      &
        qfac                       ( pdims%i_start:pdims%i_end,                &
                                     pdims%j_start:pdims%j_end )
! Vertical mean water-vapour mixing-ratio within the boundary-layer
real(kind=real_umphys) ::                                                      &
        mv_bl_mean                 ( pdims%i_start:pdims%i_end,                &
                                     pdims%j_start:pdims%j_end )
! Resolution-dependent factor in the parcel radius as a function of precip
real(kind=real_umphys) :: dxfac    ( pdims%i_start:pdims%i_end,                &
                                     pdims%j_start:pdims%j_end )

! Highest model-level within the boundary-layer
integer :: ktop                    ( pdims%i_start:pdims%i_end,                &
                                     pdims%j_start:pdims%j_end )

! Scaling factor from the grid-mean vertical velocity
real(kind=real_umphys) ::                                                      &
        w_fac                      ( pdims%i_start:pdims%i_end,                &
                                     pdims%j_start:pdims%j_end )
! Max vertical velocity below 8km
real(kind=real_umphys) ::                                                      &
        w_max                      ( pdims%i_start:pdims%i_end,                &
                                     pdims%j_start:pdims%j_end )

! Loop counters
integer :: i, j, k


!$OMP PARALLEL DEFAULT(none) private( i, j, k, ppnrate )                       &
!$OMP SHARED( pdims, tdims, bl_levels, par_radius_init_method,                 &
!$OMP         turb_len, rhokm, rho_wet_th, bl_w_var, w_fac, qfac, rainfac,     &
!$OMP         mv_bl_mean, m_v, z_rho, zh_eff, ktop, ls_rain, ls_snow,          &
!$OMP         par_radius_knob_max, par_radius_knob, par_radius_ppn_max,        &
!$OMP         w_max, z_theta, w, w_cape_limit,                                 &
!$OMP         l_resdep_precipramp, dx_ref, delta_x, delta_y, dxfac,            &
!$OMP         par_radius_amp_um )

!$OMP do SCHEDULE(STATIC)
do k = 1, bl_levels-1
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      ! Calculate turbulence-based lengthscale on theta-levels
      ! = diffusivity (m2s-1) / w' (ms-1)
      turb_len(i,j,k) = ( rhokm(i,j,k) / rho_wet_th(i,j,k) )                   &
                      / sqrt( bl_w_var(i,j,k) )
    end do
  end do
end do
!$OMP end do NOWAIT
!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    ! rhokm isn't defined on theta-level bl_levels, so set to zero.
    turb_len(i,j,bl_levels) = 0.0
  end do
end do
!$OMP end do NOWAIT

! Initialise scaling factors to 1
!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    w_fac(i,j)   = 1.0
    qfac(i,j)    = 1.0
    rainfac(i,j) = 1.0
    dxfac(i,j)   = 1.0
  end do
end do
!$OMP end do NOWAIT

! First a precip dependence
if ( par_radius_init_method >= rain_dependence ) then

  if ( par_radius_init_method >= qfacrain_dependence ) then

    ! Vertical integral; the j-loop over rows is outermost here,
    ! so that we can OMP paralellise it (we can't parallelise the k-loop).
!$OMP do SCHEDULE(STATIC)
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        mv_bl_mean(i,j) = m_v(i,j,1)*z_rho(i,j,2)
        ktop(i,j)=2
      end do
      do k = 2, bl_levels
        do i = pdims%i_start, pdims%i_end
          ! Integrate q across the PBL
          if ( z_rho(i,j,k+1) <= zh_eff(i,j) ) then
            mv_bl_mean(i,j) = mv_bl_mean(i,j) +                                &
                              m_v(i,j,k)*(z_rho(i,j,k+1)-z_rho(i,j,k))
            ktop(i,j)=k+1
          end if
        end do
      end do
      do i = pdims%i_start, pdims%i_end
        mv_bl_mean(i,j) = mv_bl_mean(i,j) / z_rho(i,j,ktop(i,j))
      end do
    end do  ! j = pdims%j_start, pdims%j_end
!$OMP end do

    if (par_radius_init_method == linear_qfacrain_dep) then
!$OMP do SCHEDULE(STATIC)
      do j = pdims%j_start, pdims%j_end
        do i = pdims%i_start, pdims%i_end
          ! Trying something that matches 1/q at qsat=0.01 and 0.02 and is
          ! quite similar between (typical range of tropical q at 500m) but
          ! increases only linearly as q decreases, so much smaller for small q
          qfac(i,j) = 2.0+1.25*max( -1.0, 100.0*(0.5*refqsat-mv_bl_mean(i,j)) )
        end do
      end do
!$OMP end do
    else  ! par_radius_init_method not linear_qfacrain_dep
!$OMP do SCHEDULE(STATIC)
      do j = pdims%j_start, pdims%j_end
        do i = pdims%i_start, pdims%i_end
          qfac(i,j) = refqsat/mv_bl_mean(i,j)
        end do
      end do
!$OMP end do
    end if ! par_radius_init_method

  end if ! >= qfacrain_dependence

  if ( l_resdep_precipramp ) then
    ! Optionally scale the precip rate dependence by horizontal grid-length.
    ! This accounts for typical resolved precip rates getting higher at
    ! smaller grid-sizes.
!$OMP do SCHEDULE(STATIC)
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        dxfac(i,j) = max( delta_x(i,j), delta_y(i,j) ) / dx_ref
      end do
    end do
!$OMP end do
  end if

!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      ppnrate = dxfac(i,j) * qfac(i,j) * (ls_rain(i,j) + ls_snow(i,j))
      rainfac(i,j) = 1.0 + ( (par_radius_knob_max/par_radius_knob)-1.0 ) *     &
                                  min( 1.0, ppnrate/par_radius_ppn_max )
    end do
  end do
!$OMP end do

end if

! Second a w-dependence
if ( par_radius_init_method == w_dependence ) then
  ! Vertical max calculation; the j-loop over rows is outermost here,
  ! so that we can OMP paralellise it (we can't parallelise the k-loop).
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      w_max(i,j) = 0.0
    end do
    do k = 1, tdims%k_end
      do i = pdims%i_start, pdims%i_end
        ! Find max w below 8km (want to avoid stratosphere)
        if (z_theta(i,j,k) < 8000.0) w_max(i,j) = max(w_max(i,j), w(i,j,k))
      end do
    end do
    do i = pdims%i_start, pdims%i_end
      ! Multiply the turb length scale by a factor that increases
      ! linearly from 1 (no enhancement) for wmax<w_cape_limit (typ =0.4 m/s)
      ! to 4 for wmax=w_cape_limit+1m/s, 7 for wmax=w_cape_limit+2m/s, etc
      w_fac(i,j) = 1.0 + max( 0.0, w_max(i,j)-w_cape_limit ) * 3.0
    end do
  end do  ! j = pdims%j_start, pdims%j_end
!$OMP end do
end if  ! w_dependence

! Enhance par_radius_amp, if requested
!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    par_radius_amp_um(i,j) = max( w_fac(i,j), rainfac(i,j) )
  end do
end do
!$OMP end do NOWAIT

!$OMP end PARALLEL


return
end subroutine calc_turb_len

end module calc_turb_len_mod
