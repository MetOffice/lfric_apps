! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Advection through falling of ice
! and rain
! Subroutine Interface:
module lsp_fall_mod

! Dr Hook Modules
use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

! lsp_sedim_eulexp is common to all subroutines
use lsp_sedim_eulexp_mod, only: lsp_sedim_eulexp

! Use in kind for large scale precip, used for compressed variables
! passed down from here
use um_types, only: real_lsprec, real_64

implicit none

! Private variables (mostly for Dr Hook)
character(len=*),   parameter, private :: ModuleName='LSP_FALL_MOD'
integer(kind=jpim), parameter, private :: zhook_in  = 0
integer(kind=jpim), parameter, private :: zhook_out = 1

contains

subroutine lsp_fall(                                                           &
  points,                                                                      &
                                          ! Number of points
  qcf_cry_0, qcf_agg_0, frac_agg,                                              &
                                          ! Ice contents
  qrain_0,qgraup_0, t,                                                         &
                                          ! Rain and graupel contents
  snow_agg,snow_cry,rainrate,grauprate,                                        &
                                             ! Sedimentation into layer
  snowt_agg,snowt_cry,rainratet,graupratet,                                    &
                                                 ! Sedim. out of layer
  vf_agg, vf_cry, vf_rain, vf_graup,                                           &
                                          ! Fall speeds of hydrometeors
  area_clear,area_ice,cfice,cficei,                                            &
                                          ! Cloud fraction information
  frac_ice_above,                                                              &
                                            ! at start of microphy. ts
  frac_ice_fall,                                                               &
                                          ! Output fraction of ice falling out
  cf, cfl, cff,                                                                &
                                          ! Current cloud fractions for
                                          ! updating
  rho, rhor, tcgi, tcgci,                                                      &
                                          ! Parametrization information
  corr, dhi, dhir, rainfrac,                                                   &
                                          ! Parametrization information
  d_qcf_cry_dt,d_qcf_agg_dt,                                                   &
                                          ! Mass transfer diagnostics
  d_qrain_dt,d_qgraup_dt,                                                      &
                                          ! Mass transfer diagnostics
  one_over_tsi,                                                                &
                                          ! 1/(timestep*iterations)
  cftransfer,cfftransfer,                                                      &
                                          ! Cloud transfer diagnostics
  uk, vk, ukp1, vkp1,                                                          &
                                          ! Winds for calc of wind-shear
  r_theta_levels_c, fv_cos_theta_latitude_c,                                   &
                                          ! Grid info
  l_use_agg_vt,                                                                &
                                          ! Vt-D branch to use
  vm_used                                                                      &
  )

  ! Microphysics and Clouds modules
use mphys_inputs_mod,     only: l_mcr_qgraup, l_mcr_qrain

implicit none

! Purpose:
!   Update cloud prognostics as a result of ice particle and
!   raindrop fall

! Method:
!   Calculate particle fall speeds and solve the advection equation
!   for mixing ratios following Rotstayns method.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.

! Both small and large ice, and raindrops, will fall, hence require
! advection downwards. We include both code for the two-ice prognostics
! and the single ice prognostic that is diagnostically split. Although
! the advection methods are very similar, they are different enough
! (in their calculation of fall-speed of the ice from above) to
! use different branches of code.

! Subroutine Arguments

integer, intent(in) :: points ! Number of points to calculate

real (kind=real_lsprec), intent(in) ::                                         &
  frac_agg(points),                                                            &
                        ! Fraction of ice mass that is aggregates
  snow_cry(points),                                                            &
                        ! Crystal flux into layer / kg m-2 s-1
  snow_agg(points),                                                            &
                        ! Aggregate flux into layer / kg m-2 s-1
  rainrate(points),                                                            &
                        ! Rain flux into layer / kg m-2 s-1
  grauprate(points),                                                           &
                        ! Graupel flux into layer / kg m-2 s-1
  area_clear(points),                                                          &
                        ! Fraction of gridbox with clear sky
  area_ice(points),                                                            &
                        ! Frac of gridbox with ice but not liquid
  cfice(points),                                                               &
                        ! Fraction of gridbox with ice cloud
  cficei(points),                                                              &
                        ! 1/Fraction of gridbox with ice cloud
  rainfrac(points),                                                            &
                        ! Rain fraction in grid-box
  frac_ice_above(points)
                        ! Ice cloud fraction falling from layer above
real (kind=real_lsprec), intent(out) ::                                        &
  frac_ice_fall(points)
                        ! Ice cloud fraction falling to layer below
real (kind=real_lsprec), intent(in) ::                                         &
  rho(points),                                                                 &
                        ! Air density / kg m-3
  rhor(points),                                                                &
                        ! 1 / Air density / m3 kg-1
  dhi(points),                                                                 &
                        ! Timestep / thickness of model layer / s m-1
  dhir(points),                                                                &
                        ! 1/dhi / m s-1
  tcgi(points),                                                                &
                        ! 1/tcg (no units)
  tcgci(points),                                                               &
                        ! 1/tcgc (no units)
  corr(points),                                                                &
                        ! Air density fall speed correction (no units)
  uk(points),                                                                  &
                        ! U wind on level k
  vk(points),                                                                  &
                        ! V wind on level k
  ukp1(points),                                                                &
                        ! U wind on level k+1
  vkp1(points),                                                                &
                        ! V wind on level k+1
  r_theta_levels_c(points),                                                    &
                        ! Distance from centre of Earth and ...
  fv_cos_theta_latitude_c(points),                                             &
                        ! ... grid info for working out gridbox size.
  one_over_tsi
                        ! 1/timestep(*iterations)
logical, intent(in) ::                                                         &
  l_use_agg_vt(points)
                        ! Determines which vt-D parameters to use

real (kind=real_lsprec), intent(in out) ::                                     &
  qcf_cry_0(points),                                                           &
                        ! Ice crystal mixing ratio / kg kg-1
  qcf_agg_0(points),                                                           &
                        ! Ice aggregate mixing ratio / kg kg-1
  qrain_0(points),                                                             &
                        ! Rain mixing ratio / kg kg-1
  qgraup_0(points),                                                            &
                        ! Graupel mixing ratio / kg kg-1
  t(points),                                                                   &
                        ! Temperature / K
  snowt_cry(points),                                                           &
                        ! Snowfall rate out of this layer / kg m-2 s-1
  snowt_agg(points),                                                           &
                          ! for crystals and aggregates
  vf_cry(points),                                                              &
                        ! On input: Fall speed of hydrometeors
  vf_agg(points),                                                              &
                                  ! entering the current layer / m s-1
  vf_rain(points),                                                             &
                        ! On output: Fall speed of hydrometeors
  vf_graup(points),                                                            &
                                  ! leaving the current layer / m s-1
  cf(points),                                                                  &
                        ! Current cloud fraction
  cfl(points),                                                                 &
                        ! Current liquid cloud fraction
  cff(points)       ! Current ice cloud fraction

real (kind=real_lsprec), intent(in out) ::                                     &
  d_qcf_cry_dt(points),                                                        &
                           ! Rate of change of crystal, aggregate,
  d_qcf_agg_dt(points),                                                        &
                             ! rain and graupel mixing ratios due
  d_qrain_dt(points),                                                          &
                             ! to sedimentation / kg kg-1 s-1
  d_qgraup_dt(points),                                                         &
  cftransfer(points),                                                          &
                         ! Cloud fraction increment this tstep
  cfftransfer(points)! Ice cloud fraction inc this tstep
real (kind=real_lsprec), intent(in out) ::                                     &
  vm_used(points)
                        ! Split fallspeed variables

real (kind=real_lsprec), intent(out) ::  graupratet(points)
                        ! Graupel rate out of this layer / kg m-2 s-1

real (kind=real_lsprec), intent(out) ::  rainratet(points)
                        ! Rain rate out of this layer / kg m-2 s-1

! Local Variables
real(kind=jprb) :: zhook_handle ! Dr Hook handle

character(len=*), parameter :: RoutineName='LSP_FALL'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! When working in 32-bit, a Cray compiler bug on power breaks PROC
! comparability.  Conditionally using NOVECTOR makes this go
! away. This directive will be valid until the end of the program unit
! or when a VECTOR directive (without any extra string) is used.
! Review upon upgrade.
#if defined(LSPREC_32B) && defined (CRAYFTN_VERSION) && (CRAYFTN_VERSION <8004000)
!DIR$ NOVECTOR
#endif

!--------------------------------------------------
! Ice Sedimentation
!--------------------------------------------------

call lsp_fall_ice( points, rho, rhor, t, qcf_agg_0, qcf_cry_0, cfice,          &
                   cficei, cfl, vf_agg, vf_cry, corr, tcgi, tcgci,             &
                   snow_agg, snow_cry, vm_used, dhi, dhir, frac_agg,           &
                   snowt_agg, snowt_cry, d_qcf_agg_dt, d_qcf_cry_dt,           &
                   frac_ice_above, frac_ice_fall,                              &
                   cftransfer, cfftransfer, cf, cff,                           &
                   area_clear, area_ice, ukp1, uk, vkp1, vk,                   &
                   one_over_tsi, r_theta_levels_c,                             &
                   fv_cos_theta_latitude_c, l_use_agg_vt )

!--------------------------------------------------
! Rain Sedimentation
!--------------------------------------------------

if (l_mcr_qrain) then

  call lsp_fall_rain(points, qrain_0, rho, rhor, corr, dhi, dhir,              &
                     rainrate, rainratet, vf_rain, one_over_tsi,               &
                     rainfrac, d_qrain_dt )

end if  ! l_mcr_qrain

!--------------------------------------------------
! Graupel Sedimentation
!--------------------------------------------------

if (l_mcr_qgraup) then

  call lsp_fall_graupel(points, qgraup_0, rho, rhor, corr, dhi, dhir,          &
                        grauprate, graupratet, vf_graup, one_over_tsi,         &
                        rainfrac, d_qgraup_dt )

end if  !  l_mcr_qgraup


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_fall

subroutine lsp_fall_graupel(points, qgraup_0, rho, rhor, corr, dhi, dhir,      &
                            grauprate, graupratet, vf_graup, one_over_tsi,     &
                            rainfrac, d_qgraup_dt)

use lsprec_mod, only: cx, constp, m0, zero
use mphys_inputs_mod, only: l_subgrid_graupel_frac

implicit none

! Subroutine arguments

integer, intent(in) :: points

real(kind=real_lsprec), intent(in) :: rho(points)  ! Air density kg m-3.
real(kind=real_lsprec), intent(in) :: rhor(points) ! 1 / Air density m3 kg-1.
real(kind=real_lsprec), intent(in) :: corr(points)
                                    ! Fall speed air density correction
real(kind=real_lsprec), intent(in) :: dhi(points)
                                    ! Timestep / thickness of model layer s m-1
real(kind=real_lsprec), intent(in) :: dhir(points)
                                    ! Inverse of dhi m s-1
real(kind=real_lsprec), intent(in) :: one_over_tsi
                                    ! 1 / (timestep * number of iterations)
real(kind=real_lsprec), intent(in) :: grauprate(points)
                                    ! Graupel flux in to this layer
                                    ! kg m-2 s-1
real(kind=real_lsprec), intent(in) :: rainfrac(points)
                                    ! Sub-grid fraction of precipitation
                                    ! (contains graupel as well as rain if
                                    !  l_subgrid_graupel_frac is true)

real(kind=real_lsprec), intent(in out) :: qgraup_0(points)
                                    ! Graupel mass mixing ratio kg kg-1

real(kind=real_lsprec), intent(in out) :: vf_graup(points)
                                    ! Fall speed of graupel m s-1
real(kind=real_lsprec), intent(in out) :: d_qgraup_dt(points)
                                    ! Mass transfer diagnostic
                                    ! kg kg-1 s-1

real(kind=real_lsprec), intent(out) :: graupratet(points)
                                    ! Graupel flux out of this layer
                                    ! kg m-2 s-1

! Local variables
integer :: npts       ! Counter for number of points
integer :: i, c       ! Loop counters
integer :: ix(points) ! Original index

real(kind=real_lsprec) :: qgraup(points)    ! Local copy of qgraup_0 kg kg-1
real(kind=real_lsprec) :: fqi_graup(points) ! Bulk fall speed of gruapel
                                            ! m s-1

real(kind=jprb) :: zhook_handle ! Dr Hook handle

character(len=*), parameter :: RoutineName='LSP_FALL_GRAUPEL'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! When working in 32-bit, a Cray compiler bug on power breaks PROC
! comparability.  Conditionally using NOVECTOR makes this go
! away. This directive will be valid until the end of the program unit
! or when a VECTOR directive (without any extra string) is used.
! Review upon upgrade.
#if defined(LSPREC_32B) && defined (CRAYFTN_VERSION) && (CRAYFTN_VERSION <8004000)
!DIR$ NOVECTOR
#endif

! Identify the points that need to be calculated.
npts = 0
do i = 1, points

  ! Make a working copy of qgraup_0
  qgraup(i) = qgraup_0(i)

  if (qgraup(i)  >   m0) then
    npts = npts + 1
    ix(npts) = i

  else
    ! Fall speed is set to zero
    fqi_graup(i) = zero
  end if  ! qgraup gt m0

end do

!--------------------------------------------------------
! Estimate bulk fall speed out of this layer (FQI_GRAUP)
!--------------------------------------------------------
if ( l_subgrid_graupel_frac ) then
  ! Calculate graupel fall-speed accounting for sub-grid fraction of graupel
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do c = 1, npts
    i = ix(c)
    fqi_graup(i) = constp(64) * corr(i) *                                      &
                   (rho(i) * qgraup(i) * constp(65) / rainfrac(i)) ** cx(63)
  end do
else
  ! Calculate fall-speed assuming graupel is homogeneous over whole grid-box
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do c = 1, npts
    i = ix(c)
    fqi_graup(i) = constp(64) * corr(i) *                                      &
                   (rho(i) * qgraup(i) * constp(65)) ** cx(63)
  end do
end if

! graupratet needs to be initialised because lsp_sedim_eulexp()
! increments graupratet inplace
do i = 1, points
  graupratet(i) = zero
end do

call lsp_sedim_eulexp( points, m0, dhi, dhir, rho, rhor,                       &
                       grauprate, fqi_graup, qgraup, vf_graup,                 &
                       graupratet )

! Update fields
do i = 1, points
  ! Form transfer diagnostic
  d_qgraup_dt(i) = d_qgraup_dt(i) + (qgraup(i) - qgraup_0(i)) * one_over_tsi

  ! Update water content to take the new value
  qgraup_0(i) = qgraup(i)
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_fall_graupel

subroutine lsp_fall_ice(points, rho, rhor, t, qcf_agg_0, qcf_cry_0, cfice,     &
                        cficei, cfl, vf_agg, vf_cry, corr, tcgi, tcgci,        &
                        snow_agg, snow_cry, vm_used, dhi, dhir, frac_agg,      &
                        snowt_agg, snowt_cry, d_qcf_agg_dt, d_qcf_cry_dt,      &
                        frac_ice_above, frac_ice_fall,                         &
                        cftransfer, cfftransfer, cf, cff,                      &
                        area_clear, area_ice, ukp1, uk, vkp1, vk,              &
                        one_over_tsi, r_theta_levels_c,                        &
                        fv_cos_theta_latitude_c, l_use_agg_vt )

use mphys_inputs_mod,  only: l_psd, l_mcr_qcf2, l_diff_icevt
use mphys_bypass_mod,  only: l_crystals
use lsprec_mod,        only: cx, constp, m0, zero, one, two, timestep_mp,      &
                             wind_shear_factor, delta_lambda, delta_phi,       &
                             cff_spread_rate

use lsp_moments_mod,   only: lsp_moments
use pc2_constants_mod, only: original_but_wrong, real_shear
use cloud_inputs_mod,  only: falliceshear_method

implicit none

! Subroutine arguments
integer, intent(in) :: points

real(kind=real_lsprec), intent(in) :: rho(points)  ! Air density [kg m-3]
real(kind=real_lsprec), intent(in) :: rhor(points) ! 1 / Air density [m3 kg-1]
real(kind=real_lsprec), intent(in) :: t(points)    ! Temperature [K]
real(kind=real_lsprec), intent(in) :: cfice(points)  ! Ice cloud fraction []
real(kind=real_lsprec), intent(in) :: cficei(points) ! 1 / Ice cloud fraction []
real(kind=real_lsprec), intent(in) :: cfl(points)    ! Liquid cloud fraction []
real(kind=real_lsprec), intent(in) :: corr(points)   ! Air density-related
                                                     ! fallspeed correction []

real(kind=real_lsprec), intent(in) :: tcgi(points)
                                     ! Inverse of temperature correction factor
                                     ! for aggregates []

real(kind=real_lsprec), intent(in) :: tcgci(points)
                                     ! Inverse temperature correction factor for
                                     ! crystals []
real(kind=real_lsprec), intent(in) :: snow_agg(points)
                                     ! Aggregate flux into the layer
                                     ! [kg m-2 s-1]
real(kind=real_lsprec), intent(in) :: snow_cry(points)
                                     ! Crystal flux into the layer [kg m-2 s-1]

real(kind=real_lsprec), intent(in) :: dhi(points)
                                     ! mphys timestep / layer thickness [m s-1]
                                     ! ( Used for CFL limit )
real(kind=real_lsprec), intent(in) :: dhir(points) ! 1 / dhi [s m-1]

real(kind=real_lsprec), intent(in) :: frac_agg(points) ! Aggregate fraction []
real(kind=real_lsprec), intent(in) :: one_over_tsi
                                     ! 1 / (timestep * iterations) [s-1]
real(kind=real_lsprec), intent(in) :: frac_ice_above(points)
                                     ! Fraction of ice falling from layer above
real(kind=real_lsprec), intent(out) :: frac_ice_fall(points)
                                     ! Fraction of ice falling into layer below
real(kind=real_lsprec), intent(in) :: area_clear(points)
                                     ! Area of grid box which is clear sky []
real(kind=real_lsprec), intent(in) :: area_ice(points)
                                     ! Area of grid box which is ice []

real(kind=real_lsprec), intent(in) :: uk(points)  ! U wind at level k [m s-1]
real(kind=real_lsprec), intent(in) :: vk(points)  ! V wind at level k [m s-1]
real(kind=real_lsprec), intent(in) :: ukp1(points) ! U wind at level k+1 [m s-1]
real(kind=real_lsprec), intent(in) :: vkp1(points) ! V wind at level k+1 [m s-1]
real(kind=real_lsprec), intent(in) :: r_theta_levels_c(points)
                                     ! Distance from centre of Earth [m]
real(kind=real_lsprec), intent(in) :: fv_cos_theta_latitude_c(points)
                                     ! Finite volume cosine of latitude []

real(kind=real_lsprec), intent(in out) :: vf_agg(points)
                                     ! Fall speed of aggregates [m s-1]
                                     ! (on input - entering the layer;
                                     !  on output - exiting the layer)
real(kind=real_lsprec), intent(in out) :: vf_cry(points)
                                     ! Fall speed of crystals [m s-1]
                                     ! (on input - entering the layer;
                                     !  on output - exiting the layer)

real(kind=real_lsprec), intent(in out) :: qcf_agg_0(points)
                                        ! Ice aggregate mixing ratio [kg kg-1]
real(kind=real_lsprec), intent(in out) :: qcf_cry_0(points)
                                        ! Ice crystal mixing ratio [kg kg-1]
real(kind=real_lsprec), intent(in out) :: vm_used(points)
                                        ! Fall speed used [m s -1]
real(kind=real_lsprec), intent(in out) :: snowt_agg(points)
                                ! Aggregate flux out of the layer [kg m-2 s-1]

real(kind=real_lsprec), intent(in out) :: snowt_cry(points)
                                ! Crystal flux out of the layer [kg m-2 s-1]

real(kind=real_lsprec), intent(in out) :: d_qcf_agg_dt(points)
                                ! Diagnostic for change in aggregate mass mixing
                                ! ratio due to sedimentation [kg kg-1 s-1]

real(kind=real_lsprec), intent(in out) :: d_qcf_cry_dt(points)
                                ! Diagnostic for change in crystal mass mixing
                                ! ratio due to sedimentation [kg kg-1 s-1]

real(kind=real_lsprec), intent(in out) :: cf(points)  ! Bulk cloud fraction []
real(kind=real_lsprec), intent(in out) :: cff(points) ! Frozen cloud fraction []
real(kind=real_lsprec), intent(in out) :: cftransfer(points)
                                    ! Bulk cloud fraction transfer diagnostic
real(kind=real_lsprec), intent(in out) :: cfftransfer(points)
                                    ! Frozen cloud fraction transfer diagnostic


logical, intent(in) ::  l_use_agg_vt(points)
                        ! Logical which retermines which vt-D parameters to use

! Local variables

integer :: i, c       ! Loop counters
integer :: npts       ! Number of points
integer :: ix(points) ! Original index of points

real(kind=real_lsprec) :: m_bi_di(points)
                      ! bi+di moment of the particle size distribution []

real(kind=real_lsprec) :: m_bic_dic(points)
                      ! bic+dic moment of the particle size distribution []

real(kind=real_lsprec) :: qcf_cry(points)
                      ! Local working copy of qcf_cry_0 [kg kg-1]

real(kind=real_lsprec) :: qcf_agg(points)
                      ! Local working copy of qcf_agg_0 [kg kg-1]

real(kind=real_lsprec) :: save_vf_agg(points)
                      ! Local saved copy of vf_agg [m s-1]

real(kind=real_lsprec) :: save_vf_cry(points)
                      ! Local saved copy of vf_cry [m s-1]

real(kind=real_lsprec) :: fqi_agg(points)
                      ! Bulk fall speed of aggregates [m s-1]

real(kind=real_lsprec) :: fqi_cry(points)
                      ! Bulk fall speed of crystals [m s-1]

real(kind=real_lsprec) :: mixratio_fromabove
                      ! Mixing ratio passed in from above [kg kg-1]

real(kind=real_lsprec) :: fqirqi2_agg
                         ! Fraction of aggregate mass that remains in the
                         ! current layer after sedimentation []

real(kind=real_lsprec) :: fqirqi2_cry
                         ! Fraction of crystal mass that remains in the
                         ! current layer after sedimentation []

real(kind=real_lsprec) :: fqirqi_agg
                        ! Flux of aggregates out of layer [m s-1]

real(kind=real_lsprec) :: fqirqi_cry
                        ! Flux of crystals out of layer [m s-1]

real(kind=real_lsprec) :: temp3    ! Temporary fraction of ice []

real(kind=real_lsprec) :: overhang ! Ice cloud fraction that overhangs
                                   ! the current layer from above []

real(kind=real_lsprec) :: deltacf(points)
                        ! Change in cloud fraction across timestep []

real(kind=real_lsprec) :: deltacff(points)
                        ! Change in ice cloud fraction across tstep []

! Variables for calculating lateral displacement of falling ice due to shear.
!------------------------------------------------------------------
real (kind=real_lsprec) :: mwfv, dudz, dvdz, shear, horiz_scale, lateral_disp, &
                           cff_perimeter

real(kind=jprb) :: zhook_handle ! Dr Hook handle

character(len=*), parameter :: RoutineName='LSP_FALL_ICE'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! When working in 32-bit, a Cray compiler bug on power breaks PROC
! comparability.  Conditionally using NOVECTOR makes this go
! away. This directive will be valid until the end of the program unit
! or when a VECTOR directive (without any extra string) is used.
! Review upon upgrade.
#if defined(LSPREC_32B) && defined (CRAYFTN_VERSION) && (CRAYFTN_VERSION <8004000)
!DIR$ NOVECTOR
#endif

!-----------------------------------------------
! Calculate moment of size distribution if appropriate
!-----------------------------------------------

if (l_psd) then

  if (.not. l_diff_icevt ) then
    ! Use only one vt-D relation
    ! Calculate the bi+di (cx(82)) moment of the
    ! ice particle size distribution

    call lsp_moments(points,rho,t,qcf_agg_0,cficei,cx(82),m_bi_di)
  else
     ! The vt-D relation which gives the
     ! least mass weighted mean fallspeed will be used so
     ! calculate both required moments

    call lsp_moments(points,rho,t,qcf_agg_0,cficei,cx(182),m_bic_dic)
      ! ice mass flux moment with crystal parameters

    call lsp_moments(points,rho,t,qcf_agg_0,cficei,cx(82),m_bi_di)
     ! ice mass flux moment with aggregate parameters
  end if ! l_diff_icevt
end if ! l_psd

do i = 1, points
      !-----------------------------------------------
      ! Copy input water contents (allows sequential updates)
      !-----------------------------------------------
  qcf_cry(i) = qcf_cry_0(i)
  qcf_agg(i) = qcf_agg_0(i)

  ! Save the fall speeds coming in from layer above before
  ! they get overwritten with fall speed passed to layer below.
  ! Will be used in updating falling ice cloud fraction
  save_vf_agg(i)=vf_agg(i)
  save_vf_cry(i)=vf_cry(i)

end do

! Identify the points that need to be calculated.
npts = 0
do i = 1, points

  if (qcf_agg(i) >  m0) then

    npts = npts + 1
    ix(npts) = i

  else
        !-----------------------------------------------
        ! Ice content is small so it is numerically best to
        ! assume there is no ice here, so the fall speed is zero.
        !-----------------------------------------------
    fqi_agg(i) = zero
  end if  ! qcf_agg gt 0

end do

! In all the following loops, the 'i' variable is used for the
! original index while the 'c' variable is used for the compressed
! index.

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
do c = 1, npts

  i = ix(c)

      !-----------------------------------------------
      ! Estimate fall speed out of this layer
      !-----------------------------------------------
  if (l_psd) then
        ! Use the generic PSD
    if (.not. l_diff_icevt) then
       ! Only one set of vt-D parameters
       ! constp(82) = ci*ai
      fqi_agg(i) = constp(82) * corr(i) * m_bi_di(i)                           &
             / (rho(i) * qcf_agg(i) * cficei(i))
    else
      if (l_use_agg_vt(i)) then
         ! Use aggregate parameters
        fqi_agg(i) = constp(82) * corr(i) * m_bi_di(i)                         &
            / (rho(i) * qcf_agg(i) * cficei(i))
      else
         ! Use crystal parameters
        fqi_agg(i) = constp(182) * corr(i) * m_bic_dic(i)                      &
            / (rho(i) * qcf_agg(i) * cficei(i))
      end if
    end if ! l_diff_icevt

  else
    fqi_agg(i) = constp(24) * corr(i) *                                        &
    (rho(i)*qcf_agg(i)*constp(25)*tcgi(i)*cficei(i))**cx(23)
  end if  ! l_psd

end do

do i = 1, points

  ! Store the selected mean fallspeed
  vm_used(i) = fqi_agg(i)

end do

! Identify the points that need to be calculated.
npts = 0
do i = 1, points

  if (qcf_cry(i) >  m0 .and. l_crystals) then
    npts = npts + 1
    ix(npts) = i
  else
        !-----------------------------------------------
        ! Ice content is small so it is numerically best to
        ! assume there is no ice here, so the fall speed is zero.
        !-----------------------------------------------
    fqi_cry(i) = zero
  end if  ! qcf_cry gt 0

end do

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
do c = 1, npts

  i = ix(c)

     !-----------------------------------------------
     ! Estimate fall speed out of this layer
     !-----------------------------------------------
  fqi_cry(i) = constp(4) * corr(i)*                                            &
  (rho(i)*qcf_cry(i)*constp(5)*tcgci(i)*cficei(i))**cx(3)

end do

if (.not. l_mcr_qcf2) then

  do i = 1, points
        !-----------------------------------------------
        ! Calculate fall speed from above as an average of the
        ! crystal and aggregate fall speeds.
        !-----------------------------------------------
    mixratio_fromabove = (snow_agg(i) + snow_cry(i)) * dhi(i)*rhor(i)
    if ((qcf_cry(i)+qcf_agg(i)+mixratio_fromabove) >  m0) then

          !-----------------------------------------------
          ! Make a linear combination of fall speeds between this
          ! layer and the layer above to aid numerical solution
          !-----------------------------------------------
      fqi_cry(i)= ( fqi_cry(i)*(qcf_cry(i)+qcf_agg(i))                         &
                    + vf_agg(i)*mixratio_fromabove )                           &
                  / (qcf_cry(i)+qcf_agg(i)+mixratio_fromabove)
      fqi_agg(i)= ( fqi_agg(i)*(qcf_cry(i)+qcf_agg(i))                         &
                    + vf_agg(i)*mixratio_fromabove )                           &
                  / (qcf_cry(i)+qcf_agg(i)+mixratio_fromabove)

          !-----------------------------------------------
          ! Fall speed of ice to pass to layer below
          !-----------------------------------------------
      vf_agg(i) = fqi_cry(i) * (one-frac_agg(i))                               &
                + fqi_agg(i) *      frac_agg(i)

    else  ! qcf gt m0
          ! Little ice so set fall speed to zero
      vf_agg(i)=zero

    end if  ! qcf gt m0

        !-----------------------------------------------
        ! Solve for fraction of ice that remains in the same layer
        !-----------------------------------------------
    fqirqi2_agg = exp(-fqi_agg(i)*dhi(i))
    fqirqi2_cry = exp(-fqi_cry(i)*dhi(i))

    if (fqi_agg(i)  >   zero) then
          !-----------------------------------------------
          ! Advect aggregates
          !-----------------------------------------------
      fqirqi_agg = snow_agg(i) + dhir(i) *                                     &
                  (rho(i)*qcf_agg(i)-snow_agg(i)/fqi_agg(i))                   &
                  * (one-fqirqi2_agg)
      qcf_agg(i) = snow_agg(i) * rhor(i) / fqi_agg(i)                          &
                  * (one-fqirqi2_agg)                                          &
                  + qcf_agg(i)*fqirqi2_agg
    else
          !-----------------------------------------------
          ! No fall of ice out of the layer
          !-----------------------------------------------
      fqirqi_agg = zero
      qcf_agg(i) = snow_agg(i)*rhor(i)*dhi(i)
    end if  ! fqi_agg gt 0

    if (fqi_cry(i)  >   zero .and. l_crystals) then
          !-----------------------------------------------
          ! Advect crystals
          !-----------------------------------------------
      fqirqi_cry = snow_cry(i) + dhir(i) *                                     &
                  (rho(i)*qcf_cry(i)-snow_cry(i)/fqi_cry(i))                   &
                  * (one-fqirqi2_cry)
      qcf_cry(i) = snow_cry(i) * rhor(i) / fqi_cry(i)                          &
                  * (one-fqirqi2_cry)                                          &
                  + qcf_cry(i)*fqirqi2_cry
    else
          !-----------------------------------------------
          ! No fall of ice out of the layer
          !-----------------------------------------------
      fqirqi_cry = zero
      qcf_cry(i) = snow_cry(i)*rhor(i)*dhi(i)
    end if  ! fqi_cry gt 0

        ! --------------------------------------------------
        ! Snow is used to save flux out of layer
        ! --------------------------------------------------
    snowt_cry(i) = snowt_cry(i) + fqirqi_cry
    snowt_agg(i) = snowt_agg(i) + fqirqi_agg

  end do  ! Points

else  ! l_mcr_qcf2
      !--------------------------------------------------
      ! Prognostic Ice Crystal and Snow Aggregate Sedimentation
      !--------------------------------------------------
      ! Ice Crystals
  if (l_crystals) then
    call lsp_sedim_eulexp(                                                     &
      points, m0, dhi, dhir, rho, rhor,                                        &
      snow_cry, fqi_cry, qcf_cry, vf_cry,                                      &
      snowt_cry)
  end if  ! l_crystals

      ! Aggregates
  call lsp_sedim_eulexp(                                                       &
    points, m0, dhi, dhir, rho, rhor,                                          &
    snow_agg, fqi_agg, qcf_agg, vf_agg,                                        &
    snowt_agg)

end if ! l_mcr_qcf2


do i = 1, points

    !-----------------------------------------------
    ! Update cloud fractions
    !-----------------------------------------------

  if ((qcf_cry(i)+qcf_agg(i))  >   zero) then

    if (falliceshear_method == original_but_wrong) then

      !----------------------------------------------------------
      ! Calculate fraction of a layer the ice has fallen
      !----------------------------------------------------------
      ! This incorrectly uses the fall velocities in THIS layer
      ! rather than those from the layer ABOVE.
      !----------------------------------------------------------
      if (l_mcr_qcf2) then
        temp3 = dhi(i) * (vf_agg(i)*qcf_agg(i)                                 &
                   + vf_cry(i)*qcf_cry(i))/(qcf_cry(i)+qcf_agg(i))
      else
        temp3 = dhi(i) * vf_agg(i)
      end if
      !----------------------------------------------------------
      ! Calculate the amount of cloud overhang between levels
      !----------------------------------------------------------
      overhang = max(frac_ice_above(i)-cff(i),zero)
      if (temp3  >   zero) then
        overhang = overhang + wind_shear_factor / temp3 * timestep_mp
      end if

      !----------------------------------------------------------
      ! Physially limit the amount of overhang. Note there is
      ! no limit on the fraction of ice in the layer above
      !----------------------------------------------------------
      overhang = min(overhang,one-cff(i))
      ! Now limit the fall out quantity
      temp3 = min(max(temp3,zero),one)

      !----------------------------------------------------------
      ! Calculate change in ice cloud fraction
      !----------------------------------------------------------

      temp3 = temp3 * overhang
      deltacff(i) = temp3

      if (temp3  <=  zero) then
        !--------------------------------------------------------
        ! Total cloud fraction will be reduced
        !--------------------------------------------------------
        deltacf(i) = temp3*area_ice(i)*cficei(i) !Random overlap
      else if (cfice(i)  <   one) then
        !--------------------------------------------------------
        ! Total cloud fraction will be increased
        !--------------------------------------------------------
        !           deltacf(i) = temp3*area_clear(i) /(1.0-cfice(i))
        ! (Random overlap )
        deltacf(i) = (min(temp3,area_clear(i))) ! Minimum overlap
      end if

    else ! (falliceshear_method /= original_but_wrong)

      !--------------------------------------------------------
      ! Calculate temp3 or "fraction ice fallen":
      ! fraction of depth of this layer the ice from the layer
      ! above has fallen. Using fall speed INTO layer.
      !--------------------------------------------------------

      if (l_mcr_qcf2) then
        ! Find mass-weighted fall velocity (mwfv) for combined
        ! ice categories.
        mwfv = (save_vf_agg(i)*qcf_agg(i)+save_vf_cry(i)*qcf_cry(i))           &
             / (qcf_cry(i)+qcf_agg(i))
      else
        ! Using single ice category
        mwfv = save_vf_agg(i)
      end if

      temp3 = mwfv * dhi(i)

      ! Ensure temp3 is positive
      ! but allow "fraction fallen" to be > 1.
      temp3 = max(temp3,zero)

      !------------------------------------------------------
      ! Calculate the amount of cloud overhang between levels
      !------------------------------------------------------
      overhang = max(frac_ice_above(i)-cff(i),zero)

      if (falliceshear_method == real_shear) then
        ! Increase the overhang depending on the vertical
        ! shear of the model wind.

        ! Magnitude of vertical shear of the horizontal wind.
        ! |dU/dz| = sqrt( dudz^2 + dvdz^2 )
        dudz = ( ukp1(i) - uk(i) )
        dvdz = ( vkp1(i) - vk(i) )
        shear = sqrt( (dudz*dudz) + (dvdz*dvdz) )

        ! The horizontal scale is taken as the square root
        ! of the area of the grid box.

        horiz_scale = sqrt (   r_theta_levels_c(i) * delta_lambda              &
                             * r_theta_levels_c(i) * delta_phi                 &
                             * fv_cos_theta_latitude_c(i)     )

        ! Calculate the horizontal distance (in metres) the ice
        ! has moved across
        lateral_disp = shear * timestep_mp

        ! Convert the lateral displacement of the falling ice
        ! cloud fraction to an increase in ice cloud fraction
        ! overhang by considering the size of the grid-box.
        overhang = overhang + ( lateral_disp / horiz_scale )

      end if ! falliceshear_method == real shear

      !-----------------------------------------------
      ! Calculate change in ice cloud fraction
      !-----------------------------------------------
      ! The overhanging cloud gets advected down a
      ! certain fraction of the depth of the layer. Now assume the
      ! cloud fills the whole depth of the layer and
      ! reduce the lateral extent while conserving cloud volume.

      deltacff(i)=min(temp3 * overhang, one-cff(i))

      ! Augment the change in ice cloud fraction to account
      ! for the lateral spreading out of ice cloud (e.g. cirrus).
      ! This will increase CFF while keeping IWC the same.
      !
      ! Cloud can only spread out from its edges, so work out the
      ! perimeter of the cloud edge as a function of cloud fraction.
      cff_perimeter=-(two*cff(i)*cff(i))+(two*cff(i))
      deltacff(i)=deltacff(i)+(cff_spread_rate*cff_perimeter*timestep_mp)
      deltacff(i)=min(deltacff(i), one-cff(i))

      if (cfice(i) < one) then
        !-----------------------------------------------
        ! Total cloud fraction will be increased
        !-----------------------------------------------
        deltacf(i) = (min(deltacff(i),area_clear(i))) ! Minimum overlap
      else if (cfice(i) >= one) then
        deltacf(i) = zero
      end if

    end if ! falliceshear_method

  else  ! qcf gt 0
    ! Set ice cloud fraction to zero and total cloud
    ! fraction to the liquid cloud fraction
    deltacff(i) = -cff(i)
    deltacf(i)  = (cfl(i) - cf(i))
  end if  ! qcf gt 0

end do

do i = 1, points

  !-----------------------------------------------
  ! Update cloud fractions
  !-----------------------------------------------

  cf(i)  = cf(i)  + deltacf(i)
  cff(i) = cff(i) + deltacff(i)

  cftransfer(i)  = cftransfer(i)  + deltacf(i)  * one_over_tsi
  cfftransfer(i) = cfftransfer(i) + deltacff(i) * one_over_tsi


  !-----------------------------------------------
  ! Form transfer diagnostics
  !-----------------------------------------------

  d_qcf_cry_dt(i) = d_qcf_cry_dt(i) + (qcf_cry(i) - qcf_cry_0(i)) * one_over_tsi
  d_qcf_agg_dt(i) = d_qcf_agg_dt(i) + (qcf_agg(i) - qcf_agg_0(i)) * one_over_tsi

  !-----------------------------------------------
  ! Update water contents
  !-----------------------------------------------

  qcf_cry_0(i) = qcf_cry(i)
  qcf_agg_0(i) = qcf_agg(i)

  ! Store the fraction consistent with ice falling down to the level below,
  ! for use as frac_ice_above at the next level down.
  frac_ice_fall(i) = cff(i)
end do  ! Points

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_fall_ice

subroutine lsp_fall_rain(points, qrain_0, rho, rhor, corr, dhi, dhir,          &
                         rainrate, rainratet, vf_rain, one_over_tsi,           &
                         rainfrac, d_qrain_dt)

use lsprec_mod,          only: cx, constp, m0, zero, one
use mphys_inputs_mod,    only: l_warm_new

implicit none

! Subroutine arguments

integer, intent(in) :: points

real(kind=real_lsprec), intent(in) :: rho(points)  ! Air density kg m-3.
real(kind=real_lsprec), intent(in) :: rhor(points) ! 1 / Air density m3 kg-1.
real(kind=real_lsprec), intent(in) :: corr(points)
                                      ! Fall speed air density correction
real(kind=real_lsprec), intent(in) :: dhi(points)
                                      ! Timestep / thickness of model layer
                                      ! s m-1
real(kind=real_lsprec), intent(in) :: dhir(points) ! Inverse of dhi m s-1
real(kind=real_lsprec), intent(in) :: one_over_tsi
                                      ! 1 / (timestep * number of iterations)
real(kind=real_lsprec), intent(in) :: rainrate(points)
                                      ! Rainrate (flux) in to this layer
                                      ! kg m-2 s-1
real(kind=real_lsprec), intent(in) :: rainfrac(points) ! Rain fraction

real(kind=real_lsprec), intent(in out) :: qrain_0(points)
                                         ! Rain mass mixing ratio kg kg-1

real(kind=real_lsprec), intent(in out) :: vf_rain(points)
                                         ! Fall speed of rain m s-1
real(kind=real_lsprec), intent(in out) :: d_qrain_dt(points)
                                         ! Mass transfer diagnostic kg kg-1 s-1

real(kind=real_lsprec), intent(out) :: rainratet(points)
                                       ! Rain flux out of this layer kg m-2 s-1

! Local variables
integer :: npts       ! Counter for number of points
integer :: i, c       ! Loop counters
integer :: ix(points) ! Original index

real(kind=real_lsprec) :: qrain(points)  ! Local rain mass mixing ratio
                                         ! kg kg-1
real(kind=real_lsprec) :: fqi_rain(points) ! Bulk fall speeds of rain m s-1
real(kind=real_lsprec) :: lamr           ! Lambda for rain.
real(kind=real_64) :: lamrh1         ! Lambda for rain + h1r
real(kind=real_64) :: lamrh2         ! Lambda for rain + h2r

real(kind=jprb) :: zhook_handle ! Dr Hook handle

character(len=*), parameter :: RoutineName='LSP_FALL_RAIN'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! When working in 32-bit, a Cray compiler bug on power breaks PROC
! comparability.  Conditionally using NOVECTOR makes this go
! away. This directive will be valid until the end of the program unit
! or when a VECTOR directive (without any extra string) is used.
! Review upon upgrade.
#if defined(LSPREC_32B) && defined (CRAYFTN_VERSION) && (CRAYFTN_VERSION <8004000)
!DIR$ NOVECTOR
#endif

! Identify the points that need to be calculated.
npts = 0
do i = 1, points

  ! Set local value of qrain
  qrain(i) = qrain_0(i)

  if (qrain(i)  >   m0) then

    npts = npts + 1
    ix(npts) = i

  else
    ! Fall speed is set to zero
    fqi_rain(i) = zero
  end if  ! qrain gt m0

end do

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
do c = 1, npts

  i = ix(c)

  !-------------------------------------------------------
  ! Estimate bulk fall speed out of this layer (FQI_RAIN)
  !-------------------------------------------------------

  ! Calculate lambda for rain:
  ! assuming rain is a mixing ratio (kg kg-1) and using
  ! Abel and Shipway (2007) rain fall speeds

  if (l_warm_new) then
    ! flux should be scaled by rain fraction
    lamr = (one/(qrain(i)*rho(i)*constp(50)/rainfrac(i)))**cx(52)
  else
    ! original (wrong) method which ignored rain fraction
    lamr = (one/(qrain(i)*rho(i)*constp(50)))**cx(52)
  end if

  ! Calculate additional properties

  !-----------------------------------------------------------------------------
  ! Start of mixed precision sub-bubble

  ! It is possible that lamrh1 and lamrh2 can occasionally have values of order
  ! 5e7. Cx(59) and cx(60) are both order 5, and therefore the resulting value
  ! is > 3e38, which is the largest single precision number. Therefore these
  ! two scalars are hard-coded to 64-bit when declared and all statements in
  ! this sub-bubble are notionally mixed precision.

  ! lamda +h1r
  lamrh1 = lamr+cx(56)

  ! lamda +h2r
  lamrh2 = lamr+cx(57)

  fqi_rain(i) = ((lamr**(cx(45)))/constp(53))*corr(i)*                         &
                ( (constp(54) / lamrh1**cx(59) ) +                             &
                  (constp(55) / lamrh2**cx(60) ) )

  ! End of mixed precision sub-bubble
  !-----------------------------------------------------------------------------


  if (l_warm_new .and. fqi_rain(i) < 0.001_real_lsprec) then
    ! minimum fall speed allowed for rain
    fqi_rain(i) = 0.001_real_lsprec
  end if

end do ! points

! rainratet needs to be initialised because lsp_sedim_eulexp()
! increments rainratet inplace
do i = 1, points
  rainratet(i) = zero
end do

call lsp_sedim_eulexp( points, m0, dhi, dhir, rho, rhor,                       &
                       rainrate, fqi_rain, qrain, vf_rain,                     &
                       rainratet )

do i = 1, points
  ! Form transfer diagnostic
  d_qrain_dt(i) = d_qrain_dt(i) + (qrain(i) - qrain_0(i)) * one_over_tsi

  ! Reset qrain_0
  qrain_0(i)   = qrain(i)
end do


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_fall_rain

end module lsp_fall_mod
