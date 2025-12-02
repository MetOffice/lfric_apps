! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Lifting condensation level calculation
!
module lift_cond_lev_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'LIFT_COND_LEV_MOD'
contains

subroutine lift_cond_lev (npnts, nlev, k_plume,                                &
                          pstar, q, t,                                         &
                          p_theta_lev, exner_rho, z_rho,                       &
                          T_lcl, p_lcl, z_lcl, qsat_lcl )

use planet_constants_mod, only: kappa => kappa_bl, pref => pref_bl,            &
     repsilon => repsilon_bl, recip_kappa => recip_kappa_bl

use cv_diag_param_mod, only:                                                   &
    a_bolton, b_bolton, c_bolton, d_bolton
use bl_option_mod, only: zero, one
use qsat_mod, only: qsat, qsat_mix

use gen_phys_inputs_mod, only: l_mr_physics

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! ------------------------------------------------------------------------------
! Description:
!   This routine calculates the lifting condensation level (LCL) temperature,
!   pressure and height.
!
!  Is designed to work on compressed arrays for just a selected set of points.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.1.
!
! ------------------------------------------------------------------------------
! Subroutine arguments

integer, intent(in) ::                                                         &
  npnts                & ! Number of points
 ,nlev                   ! Number of model levels for calculations

integer, intent(in) ::                                                         &
  k_plume(npnts)         ! Starting model level for plume ascent

real(kind=r_bl), intent(in) ::                                                 &
  pstar(npnts)            & ! Surface pressure (Pa)
 ,q(npnts,nlev)           & ! water vapour on model levels (kg/kg)
 ,t(npnts,nlev)           & ! Temperature on model levels (K)
 ,p_theta_lev(npnts,nlev) & ! Pressure on theta levels (Pa)
 ,exner_rho(npnts,nlev)   & ! Exner Pressure on rho levels
 ,z_rho(npnts,nlev)         ! Hieght of rho levels  (m)

real(kind=r_bl), intent(out) ::                                                &
  T_lcl(npnts)            & ! Temperature of LCL  (K)
 ,p_lcl(npnts)            & ! Pressure of LCL  (Pa)
 ,z_lcl(npnts)            & ! Height of LCL  (m)
 ,qsat_lcl(npnts)           ! qsaturation at zlcl (i.e. cloud base) kg/kg

!-------------------------------------------------------------------------------
! Local variables

integer ::                                                                     &
  i,k                      ! loop counter

real(kind=r_bl) ::                                                             &
  exner_lcl              & ! Exner pressure at LCL
 ,exner_surf               ! Exner pressure at surface

real(kind=r_bl) ::                                                             &
  factor                 & ! factor used in interpolation
 ,vap_press                ! vapour pressure

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LIFT_COND_LEV'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-------------------------------------------------------------------------------
! Calculate temperature and pressure of lifting condensation level
!     using approximations from Bolton (1980)
!-------------------------------------------------------------------------------
!
!   vapour pressure e ~ qp/epsilon       q specific humidity
!   vapour pressure e ~ qp/(epsilon+q)   q mixing ratio
!-------------------------------------------------------------------------------

if (l_mr_physics) then   ! expression for mixing ratio

!$OMP  PARALLEL do DEFAULT(none) private(i, vap_press)                         &
!$OMP  SHARED(npnts,k_plume,pstar,q,t,p_theta_lev,t_lcl,p_lcl,                 &
!$OMP         recip_kappa,repsilon)                                            &
!$OMP  SCHEDULE(STATIC)
  do i=1, npnts

    vap_press = 0.01_r_bl*q(i,k_plume(i)) * p_theta_lev(i,k_plume(i))          &
                                      / (repsilon+q(i,k_plume(i)) )
    if (vap_press  >   zero) then
      T_lcl(i) = a_bolton + b_bolton/                                          &
                              (c_bolton*log(t(i,k_plume(i)))                   &
                                         - log(vap_press) - d_bolton )

      p_lcl(i) = p_theta_lev(i,k_plume(i)) *                                   &
                     ( T_lcl(i) / t(i,k_plume(i)) )**recip_kappa

    else
      ! If no moisture present, set LCL to the top of the atmosphere
      ! so that the diagnosis parcel ascent will never reach it
      T_lcl(i) = zero
      p_lcl(i) = zero
    end if

  end do     ! i loop
!$OMP end PARALLEL do

else          ! expression for specific humidity

!$OMP  PARALLEL do DEFAULT(none) private(i, vap_press)                         &
!$OMP  SHARED(npnts,k_plume,pstar,q,t,p_theta_lev,t_lcl,p_lcl,                 &
!$OMP         recip_kappa,repsilon)                                            &
!$OMP  SCHEDULE(STATIC)
  do i=1, npnts
    vap_press = q(i,k_plume(i)) *                                              &
                       p_theta_lev(i,k_plume(i)) / ( 100.0_r_bl*repsilon )
    if (vap_press  >   zero) then
      T_lcl(i) = a_bolton + b_bolton/                                          &
                         (c_bolton*log(t(i,k_plume(i)))                        &
                                         - log(vap_press) - d_bolton )
      p_lcl(i) = p_theta_lev(i,k_plume(i)) *                                   &
                        ( T_lcl(i) / t(i,k_plume(i)) )**recip_kappa
    else
      ! If no moisture present, set LCL to the top of the atmosphere
      ! so that the diagnosis parcel ascent will never reach it
      T_lcl(i) = zero
      p_lcl(i) = zero
    end if

  end do
!$OMP end PARALLEL do

end if ! test on l_mr_physics

! work out qsat at LCL
if ( l_mr_physics ) then
  call qsat_mix(qsat_lcl,T_lcl,p_lcl,npnts)
else
  call qsat(qsat_lcl,T_lcl,p_lcl,npnts)
end if

!-----------------------------------------------------------------------
! Accurate calculation of height of LCL using exner_lcl rather than p_lcl
! as UM interpolates exner linearly in height but not pressure.
!-----------------------------------------------------------------------

!$OMP  PARALLEL do DEFAULT(none) private(i, k, factor, exner_lcl, exner_surf)  &
!$OMP  SHARED(npnts,nlev,exner_rho,z_rho,z_lcl,p_lcl,pstar,pref,kappa)         &
!$OMP  SCHEDULE(STATIC)
! The main index of the loop is "i" so that it can be
! parallelised using OpenMP.
do i=1, npnts
  exner_lcl  = (p_lcl(i)/pref)**kappa
  exner_surf = (pstar(i)/pref)**kappa

  k = 1

  if ( exner_lcl >= exner_surf) then
    z_lcl(i) = zero           ! at or below surface
  else if (exner_lcl < exner_surf                                              &
                   .and. exner_lcl > exner_rho(i,k)) then
    factor= (exner_rho(i,k) - exner_lcl)/                                      &
                      (exner_rho(i,k) - exner_surf)
    z_lcl(i) = (one-factor)*z_rho(i,k)
  end if

  do k=2,nlev
    if (exner_lcl >= exner_rho(i,k)                                            &
                        .and. exner_lcl < exner_rho(i,k-1) ) then
      factor= (exner_rho(i,k) - exner_lcl)/                                    &
                        (exner_rho(i,k) - exner_rho(i,k-1))
      z_lcl(i) = (one-factor)*z_rho(i,k)+factor*z_rho(i,k-1)
    end if
  end do         ! level loop

  ! If LCL pressure is lower than the pressure at the top of the
  ! model, set z_lcl to the model-top
  if ( exner_lcl < exner_rho(i,nlev) ) then
    z_lcl(i) = z_rho(i,nlev)
  end if

  ! Check z_lcl not less than a minimum value
  ! Note z_lcl currently used by diurnal cycle diagnosis and
  ! also by Deep turbulence scheme. The diurnal cycle diagnosis will be
  ! happy with a value of zero but this will not work for the deep
  ! turbulence scheme.
  ! Set z_lcl lowest model depth - fix may imply some model resolution dependence.

  z_lcl(i) =max(z_lcl(i),z_rho(i,2))
end do
!$OMP end PARALLEL do

!-------------------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lift_cond_lev
end module lift_cond_lev_mod
