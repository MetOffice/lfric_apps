! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Purpose: Calculate explicit fluxes of TL and QT

!  Programming standard: UMDP 3

!  Documentation: UM Documentation Paper No 24.

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module ex_flux_tq_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'EX_FLUX_TQ_MOD'
contains

subroutine ex_flux_tq (                                                        &
! in levels etc
  bl_levels,nSCMDpkgs,L_SCMDiags,BL_diag,                                      &
! in fields
  tl, qw, rdz, rhokh, rhokhz, grad_t_adj, grad_q_adj, rhof2, rhofsc,           &
  ft_nt, fq_nt, ft_nt_dscb, fq_nt_dscb, tothf_zh, tothf_zhsc, totqf_zh,        &
  totqf_zhsc,  weight_1dbl, ntml, ntdsc, nbdsc,                                &
! INOUT fields
  ftl, fqw, wtrac_bl                                                           &
  )

use atm_fields_bounds_mod, only: pdims, tdims, scmrowlen
use bl_option_mod, only: flux_grad, LockWhelan2006, l_converge_ga, zero
use planet_constants_mod, only: cp => cp_bl, grcp => grcp_bl
use bl_diags_mod, only: strnewbldiag
use s_scmop_mod,   only: default_streams,                                      &
    t_avg, d_bl, scmdiag_bl
use model_domain_mod, only: model_type, mt_single_column

use free_tracers_inputs_mod, only: l_wtrac, n_wtrac
use wtrac_bl_mod,            only: bl_wtrac_type

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! ARGUMENTS WITH intent in. IE: INPUT VARIABLES.

integer, intent(in) ::                                                         &
 bl_levels,                                                                    &
                            ! in No. of atmospheric levels for which
!                                boundary layer fluxes are calculated.
   ntml(pdims%i_start:pdims%i_end),                                            &
                              ! in Number of model layers in turbulent
!                                  mixing layer.
   ntdsc(pdims%i_start:pdims%i_end),                                           &
                              ! in Top level for turb mixing in any
!                                  decoupled Sc layer
   nbdsc(pdims%i_start:pdims%i_end)
                              ! in lowest flux level in DSC layer

! Additional variables for SCM diagnostics which are dummy in full UM
integer, intent(in) ::                                                         &
  nSCMDpkgs             ! No of SCM diagnostics packages

logical, intent(in) ::                                                         &
  L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

!     Declaration of new BL diagnostics.
type (strnewbldiag), intent(in out) :: BL_diag
real(kind=r_bl), intent(in) ::                                                 &
  tl(tdims%i_start:tdims%i_end, bl_levels),                                    &
                            ! in Liquid/frozen water temperture (K)
  qw(tdims%i_start:tdims%i_end, bl_levels),                                    &
                            ! in Total water content (kg/kg)
  rhokh(tdims%i_start:tdims%i_end,                                             &
        bl_levels),                                                            &
                            ! in Exchange coeffs for scalars
  rhokhz(pdims%i_start:pdims%i_end,                                            &
         2:bl_levels),                                                         &
                            ! in Non-local turbulent mixing
                            !    coefficient for heat and moisture
  weight_1dbl(pdims%i_start:pdims%i_end,                                       &
              bl_levels),                                                      &
                            ! in Weighting applied to 1D BL scheme,
                            !    to blend with Smagorinsky scheme
  rdz(pdims%i_start:pdims%i_end,bl_levels),                                    &
                            ! in RDZ(,1) is the reciprocal
                            !     height of level 1, i.e. of the
                            !     middle of layer 1.  For K > 1,
                            !     RDZ(,K) is the reciprocal of the
                            !     vertical distance from level
                            !     K-1 to level K.
  grad_t_adj(pdims%i_start:pdims%i_end),                                       &
                            ! in Temperature gradient adjustmenent
                            !    for non-local mixing in unstable
                            !    turbulent boundary layer.
  grad_q_adj(pdims%i_start:pdims%i_end)
                            ! in Humidity gradient adjustment
!                                  for non-local mixing in unstable
!                                  turbulent boundary layer.

real(kind=r_bl), intent(in) ::                                                 &
  ft_nt(pdims%i_start:pdims%i_end,                                             &
        bl_levels+1),                                                          &
                            ! in Non-turbulent heat and moisture flux
  fq_nt(pdims%i_start:pdims%i_end,                                             &
        bl_levels+1)        !    (on rho levels, surface flux(K=1)=0)
real(kind=r_bl), intent(in) ::                                                 &
  rhof2(pdims%i_start:pdims%i_end,                                             &
        2:bl_levels),                                                          &
                            ! in f2 and fsc term shape profiles
  rhofsc(pdims%i_start:pdims%i_end,                                            &
         2:bl_levels)
real(kind=r_bl), intent(in) ::                                                 &
  tothf_zh(pdims%i_start:pdims%i_end),                                         &
                            ! in Total heat fluxes at inversions
  tothf_zhsc(pdims%i_start:pdims%i_end),                                       &

  totqf_zh(pdims%i_start:pdims%i_end),                                         &
                            ! in Total moisture fluxes at inversions
  totqf_zhsc(pdims%i_start:pdims%i_end),                                       &

  ft_nt_dscb(pdims%i_start:pdims%i_end),                                       &
                            ! in Non-turbulent heat and moisture flux
  fq_nt_dscb(pdims%i_start:pdims%i_end)
                            !      at the base of the DSC layer.

! ARGUMENTS WITH intent INOUT.
real(kind=r_bl), intent(in out) ::                                             &
  ftl(pdims%i_start:pdims%i_end, bl_levels),                                   &
                           ! INOUT FTL(,K) contains net turb
!                                   sensible heat flux into layer K
!                                   from below; so FTL(,1) is the
!                                   surface sensible heat, H. (W/m2)
    fqw(pdims%i_start:pdims%i_end, bl_levels)
                             ! INOUT Moisture flux between layers
!                                   (kg per square metre per sec).
!                                   FQW(,1) is total water flux
!                                   from surface, 'E'.

! Water tracer structure containing boundary layer fields
type(bl_wtrac_type), intent(in out) :: wtrac_bl(n_wtrac)

! LOCAL VARIABLES.

character(len=*), parameter ::  RoutineName = 'EX_FLUX_TQ'

integer ::                                                                     &
  i, k, i_wt, p1

real(kind=r_bl) :: grad_ftl
real(kind=r_bl) :: grad_fqw
real(kind=r_bl) :: non_grad_ftl
real(kind=r_bl) :: non_grad_fqw
real(kind=r_bl) :: f2_ftl
real(kind=r_bl) :: f2_fqw
real(kind=r_bl) :: fsc_ftl
real(kind=r_bl) :: fsc_fqw

logical :: scm_bl_diags

real(kind=r_bl) ::                                                             &
 grad_ftl_scm(pdims%i_start:pdims%i_start+scmrowlen-1,                         &
              bl_levels),                                                      &
                                          ! K*dth/dz
 grad_fqw_scm(pdims%i_start:pdims%i_start+scmrowlen-1,                         &
              bl_levels),                                                      &
                                          ! K*dq/dz
 non_grad_ftl_scm(pdims%i_start:pdims%i_start+scmrowlen-1,                     &
                  bl_levels),                                                  &
                                          ! Non-gradient flux
 non_grad_fqw_scm(pdims%i_start:pdims%i_start+scmrowlen-1,                     &
                  bl_levels),                                                  &
                                          ! Non-gradient flux
 f2_ftl_scm(pdims%i_start:pdims%i_start+scmrowlen-1,                           &
            bl_levels),                                                        &
                                          ! Heat flux: f2 term
 f2_fqw_scm(pdims%i_start:pdims%i_start+scmrowlen-1,                           &
            bl_levels),                                                        &
                                          ! Moisture flux: f2 term
 fsc_ftl_scm(pdims%i_start:pdims%i_start+scmrowlen-1,                          &
             bl_levels),                                                       &
                                          ! Heat flux: fsc term
 fsc_fqw_scm(pdims%i_start:pdims%i_start+scmrowlen-1,                          &
             bl_levels),                                                       &
                                          ! Moisture flux: fsc term
 ftl_entr_scm(pdims%i_start:pdims%i_start+scmrowlen-1,                         &
              bl_levels),                                                      &
                                          ! Heat flux: entrainment term
 fqw_entr_scm(pdims%i_start:pdims%i_start+scmrowlen-1,                         &
              bl_levels),                                                      &
                                          ! Moisture flux: entrainment term
 ft_nt_scm(pdims%i_start:pdims%i_start+scmrowlen-1,                            &
              bl_levels),                                                      &
                                          ! Heat flux: non-turbulent term
 expl_ftl_scm(pdims%i_start:pdims%i_start+scmrowlen-1,                         &
          bl_levels)
                                          ! Heat flux: full (explicit) in W/m2

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle
!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( l_converge_ga ) THEN
  ! This option aims to improve the accuracy of SML-top height found using
  ! buoyancy flux integration by converging the gradient adjustment
  ! consistent with the latest zh...
  ! Here we ensure the gradient adjustment is applied on all model-levels
  ! where the khsurf profile is active; existing code assumes this is only
  ! up to rho-level ntml (with an explicit entrainment flux applied at ntml+1),
  ! but when buoyancy flux integration is used the khsurf profile extends
  ! up to ntml+1 (with no explicit entrainment flux).
  ! Integer "p1" is added onto ntml in the IF test, so that we go up to
  ! ntml + 1 if l_converge_ga is on, but just ntml if not.
  p1 = 1
ELSE
  p1 = 0
END IF

! Are the SCM Boundary Layer diagnostics required?
scm_bl_diags = l_scmdiags(scmdiag_bl) .and. model_type == mt_single_column

!$OMP PARALLEL DEFAULT(SHARED) private(i,k,i_wt,grad_ftl,grad_fqw,             &
!$OMP non_grad_ftl,non_grad_fqw,f2_ftl,f2_fqw,fsc_ftl, fsc_fqw)

if (scm_bl_diags) then
  ! These arrays are only for use when the SCM BL diagnostics are required.
!$OMP do SCHEDULE(STATIC)
  do k = 1, bl_levels
    do i = pdims%i_start, pdims%i_end
      grad_ftl_scm(i,k)=zero
      grad_fqw_scm(i,k)=zero
      non_grad_ftl_scm(i,k)=zero
      non_grad_fqw_scm(i,k)=zero
      f2_ftl_scm(i,k) =zero
      f2_fqw_scm(i,k) =zero
      fsc_ftl_scm(i,k)=zero
      fsc_fqw_scm(i,k)=zero
      ftl_entr_scm(i,k)=ftl(i,k)
      fqw_entr_scm(i,k)=fqw(i,k)
      ft_nt_scm(i,k)=cp*ft_nt(i,k)
    end do
  end do
!$OMP end do

  ! Remove the surface flux from the entrainment diagostic and put
  ! it in gradient flux diagnostic instead
  k = 1
!$OMP do SCHEDULE(STATIC)
  do i = pdims%i_start, pdims%i_end
    grad_ftl_scm(i,k)=ftl(i,k)
    grad_fqw_scm(i,k)=fqw(i,k)
    ftl_entr_scm(i,k)=zero
    fqw_entr_scm(i,k)=zero
  end do
!$OMP end do

end if ! scm_bl_diags

!$OMP do SCHEDULE(STATIC)
do k = 2, bl_levels
  !-----------------------------------------------------------------------
  ! 1. "Explicit" fluxes of TL and QW, on P-grid.
  !-----------------------------------------------------------------------
  do i = pdims%i_start, pdims%i_end
    grad_ftl = - rhokh(i,k) *                                                  &
      ( ( ( tl(i,k) - tl(i,k-1) ) * rdz(i,k) ) + grcp )
    grad_fqw = - rhokh(i,k) *                                                  &
          ( qw(i,k) - qw(i,k-1) ) * rdz(i,k)

      ! Copy down gradient flux into BL diagnostics array
    if (BL_diag%l_grad_ftl) then
      BL_diag%grad_ftl(i,k) = grad_ftl*cp
    end if

      ! Copy entrainment flux into BL diagnostics array before it's updated
      ! with down gradient flux
    if (BL_diag%l_ftl_e) then
      BL_diag%ftl_e(i,k) = weight_1dbl(i,k)*ftl(i,k)*cp
    end if

      ! Copy rhokhz (for nonlocal fluxes) into BL diagnostics
    if (BL_diag%l_rhokhz_ex) then
      BL_diag%rhokhz_ex(i,k) = rhokhz(i,k)
    end if

      !-------------------------------------------------------------
      ! Entrainment fluxes were specified directly in FTL,FQW in
      ! KMKHZ so now add on gradient-dependent fluxes (note that
      ! RHOKH should be zero at these levels).
      !-------------------------------------------------------------
    ftl(i,k) = weight_1dbl(i,k)*ftl(i,k) + grad_ftl
    fqw(i,k) = weight_1dbl(i,k)*fqw(i,k) + grad_fqw

      ! Copy into SCM diagnostics array if required
    if (scm_bl_diags) then
      grad_ftl_scm(i,k) = grad_ftl
      grad_fqw_scm(i,k) = grad_fqw
    end if

      !-------------------------------------------------------------
      !  Add surface-drive gradient adjustment terms to fluxes
      !  within the surface-based mixed layer.
      !-------------------------------------------------------------
    if (k  <=  ntml(i) + p1 ) then
      non_grad_ftl = weight_1dbl(i,k) *                                        &
                                 rhokhz(i,k) * grad_t_adj(i)
      non_grad_fqw = weight_1dbl(i,k) *                                        &
                                 rhokhz(i,k) * grad_q_adj(i)
      ftl(i,k) = ftl(i,k) + non_grad_ftl
      fqw(i,k) = fqw(i,k) + non_grad_fqw

      if (BL_diag%l_non_grad_ftl) then
        BL_diag%non_grad_ftl(i,k) = non_grad_ftl*cp
      end if

      if (BL_diag%l_grad_t_adj .and. k==2) then
        BL_diag%grad_t_adj(i) = grad_t_adj(i)
      end if

      if (scm_bl_diags) then
        non_grad_ftl_scm(i,k) = non_grad_ftl
        non_grad_fqw_scm(i,k) = non_grad_fqw
      end if

    end if
    end do
end do
!$OMP end do

! Repeat for water tracers  - explicit fluxes of QW on P-grid
if (l_wtrac) then
  do i_wt = 1, n_wtrac
!$OMP do SCHEDULE(STATIC)
    do k = 2, bl_levels
      do i = pdims%i_start, pdims%i_end
        grad_fqw = - rhokh(i,k) *                                              &
          ( wtrac_bl(i_wt)%qw(i,k) - wtrac_bl(i_wt)%qw(i,k-1) )                &
            * rdz(i,k)
        !-------------------------------------------------------------
        ! Entrainment fluxes were specified directly in FQW in
        ! KMKHZ so now add on gradient-dependent fluxes (note that
        ! RHOKH should be zero at these levels).
        !-------------------------------------------------------------
          wtrac_bl(i_wt)%fqw(i,k) =                                            &
            weight_1dbl(i,k)*wtrac_bl(i_wt)%fqw(i,k) + grad_fqw

        !-------------------------------------------------------------
        !  Add surface-drive gradient adjustment terms to fluxes
        !  within the surface-based mixed layer.
        !-------------------------------------------------------------
        if (k  <=  ntml(i) + p1 ) then
          non_grad_fqw = weight_1dbl(i,k) *                                    &
                        rhokhz(i,k) * wtrac_bl(i_wt)%grad_q_adj(i)
          wtrac_bl(i_wt)%fqw(i,k) = wtrac_bl(i_wt)%fqw(i,k)                    &
                                      + non_grad_fqw
        end if
      end do
    end do
!$OMP end do
  end do
end if

!-----------------------------------------------------------------------
! 2. Lock and Whelan revised non-gradient formulation
!-----------------------------------------------------------------------
if (flux_grad  ==  LockWhelan2006) then
!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do i = pdims%i_start, pdims%i_end

      if ( k  <=  ntml(i) + p1) then
        f2_ftl  = rhof2(i,k)  * tothf_zh(i)
        fsc_ftl = rhofsc(i,k) * tothf_zh(i)
        ftl(i,k)   = ftl(i,k) + weight_1dbl(i,k) *                             &
                   ( f2_ftl + fsc_ftl - ft_nt(i,k) )

        f2_fqw  = rhof2(i,k)  * totqf_zh(i)
        fsc_fqw = rhofsc(i,k) * totqf_zh(i)
        fqw(i,k)   = fqw(i,k) + weight_1dbl(i,k) *                             &
                   ( f2_fqw + fsc_fqw - fq_nt(i,k) )

        if (scm_bl_diags) then
          f2_ftl_scm(i,k) = f2_ftl
          f2_fqw_scm(i,k) = f2_fqw
          fsc_ftl_scm(i,k) = fsc_ftl
          fsc_fqw_scm(i,k) = fsc_fqw
        end if

      end if

      if ( k  <=  ntdsc(i) .and. k >= nbdsc(i) ) then

        f2_ftl  = rhof2(i,k)  *                                                &
                          ( tothf_zhsc(i)-ft_nt_dscb(i) )
        fsc_ftl = rhofsc(i,k) *                                                &
                          ( tothf_zhsc(i)-ft_nt_dscb(i) )
        ftl(i,k)   = ftl(i,k) + weight_1dbl(i,k) *                             &
                   ( f2_ftl + fsc_ftl                                          &
                       - ( ft_nt(i,k)-ft_nt_dscb(i) ) )

        f2_fqw  = rhof2(i,k)  *                                                &
                          ( totqf_zhsc(i)-fq_nt_dscb(i) )
        fsc_fqw = rhofsc(i,k) *                                                &
                          ( totqf_zhsc(i)-fq_nt_dscb(i) )
        fqw(i,k)   = fqw(i,k) + weight_1dbl(i,k) *                             &
                   ( f2_fqw + fsc_fqw                                          &
                       - ( fq_nt(i,k)-fq_nt_dscb(i) ) )

        if (scm_bl_diags) then
          f2_ftl_scm(i,k) = f2_ftl
          f2_fqw_scm(i,k) = f2_fqw
          fsc_ftl_scm(i,k) = fsc_ftl
          fsc_fqw_scm(i,k) = fsc_fqw
        end if

      end if
    end do
  end do
!$OMP end do
end if   ! FLUX_GRAD
!$OMP end PARALLEL

!-----------------------------------------------------------------------
!     SCM Boundary Layer Diagnostics Package
!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine ex_flux_tq
end module ex_flux_tq_mod
