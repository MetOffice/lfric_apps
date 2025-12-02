! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph
module interp_turb_mod

use um_types, only: real_umphys

implicit none

contains

! Subroutine to interpolate UM turbulence fields onto the
! model-levels where CoMorph expects them to be defined.
! This entails interpolating the momentum diffusivity,
! turbulent vertical velocity variance and wind-stresses
! onto rho-levels
subroutine interp_turb(                                                        &
                  z_rho, z_theta,                                              &
                  bl_w_var, fb_surf, zh, u_s, taux_p, tauy_p,                  &
                  w_var_rh, fu_rh, fv_rh )

use nlsizes_namelist_mod, only: row_length, rows, model_levels,                &
                                bl_levels
use atm_fields_bounds_mod, only: tdims, pdims

implicit none

! Model-level heights above surface
real(kind=real_umphys), intent(in) :: z_rho   ( row_length, rows, model_levels )
real(kind=real_umphys), intent(in) :: z_theta ( row_length, rows, model_levels )

! Sub-grid turbulent variance in vertical velocity
real(kind=real_umphys), intent(in) ::                                          &
                    bl_w_var ( tdims%i_start:tdims%i_end,                      &
                               tdims%j_start:tdims%j_end,                      &
                               1:tdims%k_end )

! Surface buoyancy flux / m2 s-3
real(kind=real_umphys), intent(in) ::                                          &
                    fb_surf  ( pdims%i_start:pdims%i_end,                      &
                               pdims%j_start:pdims%j_end )
! Boundary layer height / m
real(kind=real_umphys), intent(in) ::                                          &
                    zh       ( pdims%i_start:pdims%i_end,                      &
                               pdims%j_start:pdims%j_end )
! Surface friction velocity / m s-1
real(kind=real_umphys), intent(in) ::                                          &
                    u_s      ( pdims%i_start:pdims%i_end,                      &
                               pdims%j_start:pdims%j_end )

! Turbulent wind stress on p grid (on theta-levels) N m-2
real(kind=real_umphys), intent(in) ::                                          &
                    taux_p ( pdims%i_start:pdims%i_end,                        &
                             pdims%j_start:pdims%j_end,                        &
                             0:bl_levels-1 )
real(kind=real_umphys), intent(in) ::                                          &
                    tauy_p ( pdims%i_start:pdims%i_end,                        &
                             pdims%j_start:pdims%j_end,                        &
                             0:bl_levels-1 )

! Sub-grid turbulent variance in vertical velocity on rho-levels
real(kind=real_umphys), intent(out) ::                                         &
                     w_var_rh ( pdims%i_start:pdims%i_end,                     &
                                pdims%j_start:pdims%j_end,                     &
                                1:bl_levels )

! Turbulent wind stress on p grid (on rho-levels) N m-2
real(kind=real_umphys), intent(out) ::                                         &
                     fu_rh    ( pdims%i_start:pdims%i_end,                     &
                                pdims%j_start:pdims%j_end,                     &
                                1:bl_levels )
real(kind=real_umphys), intent(out) ::                                         &
                     fv_rh    ( pdims%i_start:pdims%i_end,                     &
                                pdims%j_start:pdims%j_end,                     &
                                1:bl_levels )

! Interpolation weight
real(kind=real_umphys) :: interp

! Power used in surface w-variance scale formula
real(kind=real_umphys), parameter :: two_thirds = 2.0/3.0

! Loop counters
integer :: i, j, k


!$OMP PARALLEL DEFAULT(none) private( i, j, k, interp )                        &
!$OMP SHARED( pdims, bl_levels, fb_surf, w_var_rh, zh, u_s,                    &
!$OMP         fu_rh, taux_p, fv_rh, tauy_p, z_theta, z_rho, bl_w_var )

! Interpolate w_var and wind-stresses onto rho-levels
! "Bottom rho-level" is actually the surface
!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    ! Note: bl_w_var doesn't exist at the surface in the prognostic,
    ! so using surface value calculated exactly the same as "w_m",
    ! which is used in the calculation of tv1_sd in bdy_layr
    if ( fb_surf(i,j) > 0.0 ) then
      w_var_rh(i,j,1) = ( 0.25*zh(i,j)*fb_surf(i,j)                            &
                        + u_s(i,j)*u_s(i,j)*u_s(i,j) ) ** two_thirds
    else
      w_var_rh(i,j,1) = u_s(i,j) * u_s(i,j)
    end if
    fu_rh(i,j,1) = taux_p(i,j,0)
    fv_rh(i,j,1) = tauy_p(i,j,0)
  end do
end do
!$OMP end do NOWAIT
! Linear interpolation between theta-levels k-1 and k
!$OMP do SCHEDULE(STATIC)
do k = 2, bl_levels-1
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      interp = ( z_rho(i,j,k)   - z_theta(i,j,k-1) )                           &
             / ( z_theta(i,j,k) - z_theta(i,j,k-1) )
      w_var_rh(i,j,k) = (1.0-interp) * bl_w_var(i,j,k-1)                       &
                      +      interp  * bl_w_var(i,j,k)
      fu_rh(i,j,k) = (1.0-interp) * taux_p(i,j,k-1)                            &
                   +      interp  * taux_p(i,j,k)
      fv_rh(i,j,k) = (1.0-interp) * tauy_p(i,j,k-1)                            &
                   +      interp  * tauy_p(i,j,k)
    end do
  end do
end do
!$OMP end do NOWAIT
! Interpolate at bl_levels assuming bl_w_var(:,:,bl_levels) = 0
k = bl_levels
!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    interp = ( z_rho(i,j,k)   - z_theta(i,j,k-1) )                             &
           / ( z_theta(i,j,k) - z_theta(i,j,k-1) )
    w_var_rh(i,j,k) = (1.0-interp) * bl_w_var(i,j,k-1)
    fu_rh(i,j,k) = (1.0-interp) * taux_p(i,j,k-1)
    fv_rh(i,j,k) = (1.0-interp) * tauy_p(i,j,k-1)
  end do
end do
!$OMP end do NOWAIT

!$OMP end PARALLEL


return
end subroutine interp_turb


end module interp_turb_mod
