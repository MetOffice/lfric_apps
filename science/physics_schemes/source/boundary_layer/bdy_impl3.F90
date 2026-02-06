! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!  Subroutine BDY_IMPL3

!  Purpose: Calculate downward sweep of matrix for increments to
!           U, V, T and Q in the boundary layer for the
!           unconditionally stable and non-oscillatory numerical solver

!  Programming standard: UMDP3

!  Documentation: UMDP24

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module bdy_impl3_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'BDY_IMPL3_MOD'
contains

subroutine bdy_impl3 (                                                         &
! in levels/switches
 bl_levels, l_correct,                                                         &
! in fields
 q,qcl,qcf,q_latest,qcl_latest,qcf_latest,t,t_latest,                          &
 dtrdz_charney_grid,dtrdz_u,dtrdz_v,rhokh,rhokm_u,rhokm_v,                     &
 rdz_charney_grid,rdz_u,rdz_v,gamma1,gamma2,gamma_in,                          &
 du_nt,dv_nt, r_theta_levels, r_rho_levels,                                    &
 k_blend_tq,k_blend_u,k_blend_v,                                               &
! INOUT fields
 fqw,ftl,tau_x,tau_y,du,dv,dqw,dtl,                                            &
! out fields
 dqw_nt,dtl_nt,qw,tl,ct_ctq,cq_cm_u,cq_cm_v,                                   &
 cq_cm_u_1,cq_cm_v_1,du_1,dv_1,                                                &
 dqw1_1,dtl1_1,ctctq1_1                                                        &
 )

use atm_fields_bounds_mod, only:                                               &
 udims, vdims, udims_s, vdims_s, pdims, tdims, tdims_l
use bl_option_mod, only: one
use planet_constants_mod, only: lcrcp => lcrcp_bl, lsrcp => lsrcp_bl
use vectlib_mod, only: oneover_v => oneover_v_interface
use model_domain_mod, only: model_type, mt_single_column
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

!$ use omp_lib, only: omp_get_max_threads

implicit none

! in arrays
logical, intent(in) ::                                                         &
 l_correct

integer, intent(in) ::                                                         &
 bl_levels,                                                                    &
         ! in No. of atmospheric levels for
         !    which boundary layer fluxes are
         !    calculated.
 k_blend_tq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                                   ! in Theta level for blending height.
 k_blend_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),               &
                                   ! in u level for blending height.
 k_blend_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
                                   ! in v level for blending height.

real(kind=r_bl), intent(in) ::                                                 &
 gamma1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                  &
 gamma2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                  &
                               ! in new scheme weights.
 gamma_in(bl_levels)          ! in standard implicit scheme weights.

real(kind=r_bl), intent(in) ::                                                 &
 r_theta_levels(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,   &
                tdims_l%k_start:bl_levels),                                    &
 r_rho_levels(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,     &
              1:bl_levels),                                                    &
                                 ! in height of rho and theta levels
 q(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,                &
   tdims_l%k_start:bl_levels),                                                 &
                                 ! in specific humidity
 qcl(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,              &
     tdims_l%k_start:bl_levels),                                               &
                                 ! in Cloud liquid water
 qcf(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,              &
     tdims_l%k_start:bl_levels),                                               &
                                 ! in Cloud ice (kg per kg air)
 q_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          bl_levels),                                                          &
                                 ! in specific humidity
 qcl_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,               &
            bl_levels),                                                        &
                                 ! in Cloud liquid water
 qcf_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,               &
            bl_levels),                                                        &
                                 ! in Cloud ice (kg per kg air)
 t(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),             &
                                 ! in temperature
                                 !    Latest estimates to time
                                 !    level n+1 values
 t_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          bl_levels),                                                          &
                                 ! in temperature
 dtrdz_charney_grid(tdims%i_start:tdims%i_end,                                 &
                    tdims%j_start:tdims%j_end,bl_levels),                      &
                                 ! in dz for bottom BL_LEVELS
 dtrdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,                  &
           bl_levels),                                                         &
                                 ! in -g.dt/dp for model wind layers
 dtrdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,                  &
           bl_levels),                                                         &
                                 ! in -g.dt/dp for model wind layers
 rhokh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                    &
       bl_levels),                                                             &
                                 ! in Exchange coeff for FTL above
                                 !    surface.
 rhokm_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,                  &
         bl_levels),                                                           &
                                 ! in Exchange coefficients for
                                 !    momentum, on U-grid with
                                 !    first and last rows ignored.
                                 !    for K>=2 (from KMKH).
 rhokm_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,                  &
         bl_levels),                                                           &
                                 ! in Exchange coefficients for
                                 !    momentum, on V-grid with
                                 !    first and last rows ignored.
                                 !    for K>=2 (from KMKH).
 rdz_charney_grid(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
                  bl_levels),                                                  &
                                 ! in RDZ(,1) is the reciprocal of the
                                 ! height of level 1, i.e. of the
                                 ! middle of layer 1.  For K > 1,
                                 ! RDZ(,K) is the reciprocal
                                 ! of the vertical distance
                                 ! from level K-1 to level K.
 rdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,                    &
        2:bl_levels),                                                          &
                                 ! in Reciprocal of the vertical
                                 !    distance from level K-1 to
                                 !    level K. (K > 1) on wind levels
 rdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,                    &
        2:bl_levels),                                                          &
                                 ! in Reciprocal of the vertical
                                 !    distance from level K-1 to
                                 !    level K. (K > 1) on wind levels
 du_nt(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,            &
        bl_levels),                                                            &
                                 ! in u non-turbulent increments.
 dv_nt(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,            &
        bl_levels)
                                 ! in v non-turbulent increments.
! INOUT arrays
real(kind=r_bl), intent(in out) ::                                             &
 fqw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),           &
                                 ! INOUT Flux of QW (ie., for surface,
                                 !    total evaporation). Kg/sq m/s
 ftl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),           &
                                 ! INOUT Flux of TL (ie., for surface,
                                 !    H/Cp where H is sensible heat
                                 !    in W per sq m).
 tau_x(udims%i_start:udims%i_end,udims%j_start:udims%j_end,                    &
       bl_levels),                                                             &
                                 ! INOUT x-component of turbulent
                                 !    stress at levels k-1/2;
                                 !    eg. TAUX(,1) is surface stress.
                                 !    U-grid, 1st and last rows set
                                 !    to "missing data". (N/sq m)
                                 !    in as "explicit" fluxes from
                                 !    ex_flux_uv, out as "implicit
 tau_y(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,                    &
       bl_levels),                                                             &
                                 ! INOUT y-component of turbulent
                                 !    stress at levels k-1/2;
                                 !    eg. TAUX(,1) is surface stress.
                                 !    V-grid, 1st and last rows set
                                 !    to "missing data". (N/sq m)
                                 !    in as "explicit" fluxes from
                                 !    ex_flux_uv, out as "implicit
 du(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,               &
      bl_levels),                                                              &
                                 ! INOUT BL increment to u wind field
 dv(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,               &
      bl_levels),                                                              &
                                 ! INOUT BL increment to v wind field
 dqw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),           &
                                 ! INOUT BL increment to q field
 dtl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels)
                                 ! INOUT BL increment to T field

! out arrays which are unused in LFRic and therefore need declaring as in out
real(kind=r_bl), intent(in out) ::                                             &
 cq_cm_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,                  &
           bl_levels),                                                         &
                                 ! Coefficient in U and V
                                 !  tri-diagonal implicit matrix
 cq_cm_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,                  &
           bl_levels),                                                         &
                                 ! Coefficient in U and V
                                 !  tri-diagonal implicit matrix
 cq_cm_u_1(udims%i_start:udims%i_end,udims%j_start:udims%j_end),               &
                                 ! Coefficient of taux*
                                 !  for implicit coupling
                                 !  at level k_blend_u
 du_1(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end),            &
                                 ! Coefficient needed
                                 !  for implicit coupling
                                 !  at level k_blend_u
 cq_cm_v_1(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),               &
                                 ! Coefficient of tauy*
                                 !  for implicit coupling
                                 !  at level k_blend_v
 dv_1(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end)
                                 ! Coefficient needed
                                 !  for implicit coupling
                                 !  at level k_blend_v

! out arrays
real(kind=r_bl), intent(out) ::                                                &
 qw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels),           &
                                 ! out total water
 tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels),           &
                                 ! out liquid water temperature
 ct_ctq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        bl_levels),                                                            &
                                 ! out Coefficient in T and q
                                 !     tri-diagonal implicit matrix
 dqw_nt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        bl_levels),                                                            &
                                       ! NT incr to q field
 dtl_nt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        bl_levels),                                                            &
                                       ! NT incr to T field
 ctctq1_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
                                 ! out Coefficient of H*, E*
                                 !     for implicit coupling
                                 !     at level k_blend_tq
 dqw1_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
                                 ! out Coefficient needed
                                 !     for implicit coupling
                                 !     at level k_blend_tq
 dtl1_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                 ! out Coefficient needed
                                 !     for implicit coupling
                                 !     at level k_blend_tq

! LOCAL arrays

! The three set of arrays below are needed by the uncond stable
! BL numerical solver
real(kind=r_bl) ::                                                             &
  dqw1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),         &
  dtl1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),         &
  ctctq1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels)


real(kind=r_bl) ::                                                             &
 r_theta_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,                &
           0:bl_levels),                                                       &
                                             ! Vertical grids for U
 r_theta_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,                &
           0:bl_levels),                                                       &
                                             ! and V flux levels
 ct_prod(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
                ! Product of coefficients in T and q matrix needed to
                ! calculate the coefficients of surface fluxes for
                ! implicit coupling at level k_blend_tq
 cu_prod(udims%i_start:udims%i_end,udims%j_start:udims%j_end),                 &
                ! Product of coefficients in U matrix needed to
                ! calculate the coefficient of surface momentum flux
                ! for implicit coupling at level k_blend_u
 cv_prod(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
                ! Product of coefficients in V matrix needed to
                ! calculate the coefficient of surface momentum flux
                ! for implicit coupling at level k_blend_v

real(kind=r_bl) ::                                                             &
temp(pdims%i_end*pdims%j_end),                                                 &
temp_out(pdims%i_end*pdims%j_end),                                             &
                          ! temp for pressure grid vector division
temp_u( udims%i_len * udims%j_len ),                                           &
temp_u_out( udims%i_len *                                                      &
            udims%j_len ),                                                     &
                                 ! temp for u grid vector division
temp_v( vdims%i_len * vdims%j_len ),                                           &
temp_v_out( vdims%i_len *                                                      &
            vdims%j_len )
                                 ! temp for v grid vector division
!  Local scalars :-
real(kind=r_bl) ::                                                             &
 at,                                                                           &
              ! Matrix element in "T" row in eqn P244.79.
 rbt,                                                                          &
              ! Reciprocal of BT' (eqns P244.107, 110, 113).
 am,                                                                           &
              ! Matrix element in eqn P244.80.
 rbm,                                                                          &
              ! Reciprocal of BM(') (eqns P244.81, 85, 89).
 gamma1_uv,                                                                    &
 gamma2_uv,                                                                    &
              ! gamma1 and gamma2 shifted to u or v points
 r_sq,                                                                         &
              ! square of height variables
 rr_sq    ! 1/square of height variables

integer ::                                                                     &
 blm1,                                                                         &
              ! BL_LEVELS minus 1.
 i,j,                                                                          &
              ! Loop counter (horizontal field index).
 k,                                                                            &
              ! Loop counter (vertical index).
 tdims_omp_block,                                                              &
              ! omp block length
 tdims_seg_block,                                                              &
              ! omp segment length
 ii,                                                                           &
              ! omp block loop counter
 l,                                                                            &
              ! vector counter
 max_threads

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='BDY_IMPL3'


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

blm1 = bl_levels-1

max_threads = 1
!$ max_threads = omp_get_max_threads()
tdims_omp_block  = ceiling(real(tdims%i_len)/max_threads)
tdims_seg_block = min(tdims_omp_block, tdims%i_len)

!$OMP  PARALLEL DEFAULT(none) SHARED(tdims_seg_block,l_correct,bl_levels,      &
!$OMP  blm1,tdims, dqw_nt,dtl_nt,q_latest,qcl_latest,dtrdz_v,dtrdz_u,udims,    &
!$OMP  rdz_v,gamma1,q,qcl,qcf,t_latest,t,ftl,rhokh,dtl,rdz_charney_grid,dqw,   &
!$OMP  tau_x,rhokm_u,du,rdz_u,vdims,tau_y,dv, qcf_latest,                      &
!$OMP  qw,tl,r_theta_levels,r_theta_u,r_theta_v,r_rho_levels,fqw,              &
!$OMP  dtrdz_charney_grid,gamma2,ct_ctq,dqw1,dtl1,ctctq1,model_type,           &
!$OMP  cq_cm_u_1,cq_cm_v_1,du_1,dv_1,                                          &
!$OMP  dqw1_1,dtl1_1,ctctq1_1,                                                 &
!$OMP  ct_prod, cu_prod, cv_prod,k_blend_tq,k_blend_u,k_blend_v,               &
!$OMP  gamma_in,cq_cm_u,cq_cm_v,du_nt,dv_nt,rhokm_v,lcrcp,lsrcp)               &
!$OMP  private(k,j,i,r_sq,rbt,temp,temp_u,temp_v,l,temp_out,temp_u_out,        &
!$OMP  temp_v_out,at,am,rbm,rr_sq,ii,gamma1_uv,gamma2_uv)

if ( l_correct ) then

!$OMP do SCHEDULE(STATIC)
  do k = 1, bl_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        ! Don't use QW, TL here as these are no longer at time level n
        dqw_nt(i,j,k) = q_latest(i,j,k) + qcl_latest(i,j,k)                    &
                      + qcf_latest(i,j,k)                                      &
                      - q(i,j,k) - qcl(i,j,k) - qcf(i,j,k)
        dtl_nt(i,j,k) = t_latest(i,j,k)                                        &
             - lcrcp * qcl_latest(i,j,k)                                       &
             - lsrcp * qcf_latest(i,j,k)                                       &
             - ( t(i,j,k) - lcrcp*qcl(i,j,k) - lsrcp*qcf(i,j,k) )
      end do
    end do
  end do
!$OMP end do NOWAIT

  ! Update explicit fluxes using predictor X* value as needed by the
  ! 2nd stage of the scheme. Note that: DTL=TL*-TL, DQW=QW*-QW etc

!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        ftl(i,j,k) = ftl(i,j,k) - rhokh(i,j,k) *                               &
          ( dtl(i,j,k) - dtl(i,j,k-1) ) * rdz_charney_grid(i,j,k)
        fqw(i,j,k) = fqw(i,j,k) - rhokh(i,j,k) *                               &
          ( dqw(i,j,k) - dqw(i,j,k-1) ) * rdz_charney_grid(i,j,k)
      end do
    end do
  end do
!$OMP end do

else

!$OMP do SCHEDULE(STATIC)
  do k = 1, bl_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        qw(i,j,k) = q(i,j,k) + qcl(i,j,k) + qcf(i,j,k)
        tl(i,j,k) = t(i,j,k) - lcrcp*qcl(i,j,k) - lsrcp*qcf(i,j,k)
        dqw_nt(i,j,k) = q_latest(i,j,k) + qcl_latest(i,j,k)                    &
                        + qcf_latest(i,j,k) - qw(i,j,k)
        dtl_nt(i,j,k) = t_latest(i,j,k)                                        &
                        - lcrcp * qcl_latest(i,j,k)                            &
                        - lsrcp * qcf_latest(i,j,k)                            &
                        - tl(i,j,k)
      end do
    end do
  end do
!$OMP end do

end if  ! l_correct

!-----------------------------------------------------------------------
!  1.0 Interpolate r_theta_levels to U,V columns
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 2.0 For simulations on a sphere we must use spherical geometry for
!     vertical flux-divergences.   Thus, leaving out rho for
!     simplicity, the standard cartesian flux-divergence:
!          dQ(K)/dt = -(FQ(K+1)-FQ(K))/DZ
!     becomes:
!          dQ(K)/dt = -(r_flux(K+1)^2*FQ(K+1)-r_flux(K)^2*FQ(K))
!                      / (r_full(K)^2 * DZ)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  3.0 Calculate matrix elements
!-----------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end
    ! Include non-turbulent increments.
    r_sq = r_rho_levels(i,j,bl_levels)*r_rho_levels(i,j,bl_levels)
    dqw(i,j,bl_levels) = ( dtrdz_charney_grid(i,j,bl_levels) *                 &
                    (r_sq * fqw(i,j,bl_levels)) +                              &
                     dqw_nt(i,j,bl_levels) ) * gamma2(i,j)
    dtl(i,j,bl_levels) = ( dtrdz_charney_grid(i,j,bl_levels) *                 &
                   (r_sq * ftl(i,j,bl_levels)) + dtl_nt(i,j,bl_levels)         &
                         ) * gamma2(i,j)
    ct_ctq(i,j,bl_levels) = -dtrdz_charney_grid(i,j,bl_levels) *               &
           gamma1(i,j)*(rhokh(i,j,bl_levels)*r_sq)*                            &
            rdz_charney_grid(i,j,bl_levels)
    rbt = one / ( one - ct_ctq(i,j,bl_levels) )
    dqw(i,j,bl_levels) = rbt * dqw(i,j,bl_levels)
    dtl(i,j,bl_levels) = rbt * dtl(i,j,bl_levels)
    ct_ctq(i,j,bl_levels) = rbt * ct_ctq(i,j,bl_levels)
  end do
end do
!$OMP end do

!$OMP do SCHEDULE(STATIC)
do ii = tdims%i_start, tdims%i_end, tdims_seg_block
  do k = blm1, 2, -1
    l = 0
    do j = tdims%j_start, tdims%j_end
      do i = ii, min(ii+tdims_seg_block -1, tdims%i_end)
        r_sq = r_rho_levels(i,j,k)*r_rho_levels(i,j,k)
        rr_sq = r_rho_levels(i,j,k+1)*r_rho_levels(i,j,k+1)
        dqw(i,j,k) = ( -dtrdz_charney_grid(i,j,k)*                             &
             ((rr_sq*fqw(i,j,k+1))-(r_sq*fqw(i,j,k)))+dqw_nt(i,j,k) )          &
               *gamma2(i,j)
        dtl(i,j,k) = ( -dtrdz_charney_grid(i,j,k)*                             &
             ((rr_sq*ftl(i,j,k+1))-(r_sq*ftl(i,j,k)))+dtl_nt(i,j,k) )          &
               *gamma2(i,j)
        at = -dtrdz_charney_grid(i,j,k) *                                      &
             gamma1(i,j)*(rr_sq*rhokh(i,j,k+1))*                               &
             rdz_charney_grid(i,j,k+1)
        ct_ctq(i,j,k) = -dtrdz_charney_grid(i,j,k) *                           &
             gamma1(i,j)*(r_sq*rhokh(i,j,k))*rdz_charney_grid(i,j,k)
        l = l + 1
        temp(l) = ( one - ct_ctq(i,j,k) -                                      &
             at*( one + ct_ctq(i,j,k+1) ) )
        dqw(i,j,k) = (dqw(i,j,k) - at*dqw(i,j,k+1) )
        dtl(i,j,k) = (dtl(i,j,k) - at*dtl(i,j,k+1) )
      end do
    end do

    call oneover_v(l, temp, temp_out)

    l = 0
    do j = tdims%j_start, tdims%j_end
      do i = ii, min(ii+tdims_seg_block -1, tdims%i_end)
        l = l + 1
        dqw(i,j,k) = temp_out(l) * dqw(i,j,k)
        dtl(i,j,k) = temp_out(l) * dtl(i,j,k)
        ct_ctq(i,j,k) = temp_out(l) * ct_ctq(i,j,k)
      end do
    end do

  end do !blm1,2,-1
end do
!$OMP end do

!-----------------------------------------------------------------------
!  Bottom model layer QW row of matrix equation.
!-----------------------------------------------------------------------

if ( .not. l_correct ) then

  !-----------------------------------------------------------------------
  ! The following calculations are only done on the 1st stage (predictor).
  ! Their purpose is to compute the surface scalar (T, Q) increments which
  ! are needed by the surface scheme to compute the implicit scalar fluxes.
  ! The same implicit scalar fluxes are used as a boundary condition for
  ! discrete equations of the 2nd stage.
  ! Due to the dependency to the surface scalar fluxes, the 1st stage
  ! downward sweep remains incomplete. It is completed at the
  ! beginning of bdy_impl4, since there the surface scalar fluxes
  ! have been fully updated. This is done by sf_impl2().
  ! NOTE: The standard scheme solver is used for this calculation.
  !       Incorporation of the new scheme for the scalar surface variables
  !       would have been a more preferable choice but currently not
  !       available.
  !-----------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      ! Include non-turbulent increments.
      r_sq = r_rho_levels(i,j,bl_levels)*r_rho_levels(i,j,bl_levels)
      dqw1(i,j,bl_levels) = dtrdz_charney_grid(i,j,bl_levels)*                 &
                       (r_sq*fqw(i,j,bl_levels)) +                             &
                       dqw_nt(i,j,bl_levels)
      dtl1(i,j,bl_levels) = dtrdz_charney_grid(i,j,bl_levels)*                 &
                       (r_sq*ftl(i,j,bl_levels)) +                             &
                       dtl_nt(i,j,bl_levels)
      ctctq1(i,j,bl_levels) = -dtrdz_charney_grid(i,j,bl_levels)*              &
         gamma_in(bl_levels)*r_sq*rhokh(i,j,bl_levels)*                        &
         rdz_charney_grid(i,j,bl_levels)
      rbt = one / ( one - ctctq1(i,j,bl_levels) )
      dqw1(i,j,bl_levels) = rbt * dqw1(i,j,bl_levels)
      dtl1(i,j,bl_levels) = rbt * dtl1(i,j,bl_levels)
      ctctq1(i,j,bl_levels) = rbt * ctctq1(i,j,bl_levels)
    end do
  end do
!$OMP end do

!$OMP do SCHEDULE(STATIC)
  do ii = tdims%i_start, tdims%i_end, tdims_seg_block
    do k = blm1, 2, -1
      l = 0
      do j = tdims%j_start, tdims%j_end
        do i = ii, min(ii+tdims_seg_block -1, tdims%i_end)
          r_sq = r_rho_levels(i,j,k)*r_rho_levels(i,j,k)
          rr_sq = r_rho_levels(i,j,k+1)*r_rho_levels(i,j,k+1)
          dqw1(i,j,k) = -dtrdz_charney_grid(i,j,k) *                           &
            ((rr_sq*fqw(i,j,k+1)) - (r_sq*fqw(i,j,k))) + dqw_nt(i,j,k)
          dtl1(i,j,k) = -dtrdz_charney_grid(i,j,k) *                           &
            ((rr_sq*ftl(i,j,k+1)) - (r_sq*ftl(i,j,k))) + dtl_nt(i,j,k)
          at = -dtrdz_charney_grid(i,j,k) *                                    &
           gamma_in(k+1)*(rr_sq*rhokh(i,j,k+1))*rdz_charney_grid(i,j,k+1)
          ctctq1(i,j,k) = -dtrdz_charney_grid(i,j,k) *                         &
            gamma_in(k)*(r_sq*rhokh(i,j,k))*rdz_charney_grid(i,j,k)
          ! pack
          l = l + 1
          temp(l) = ( one - ctctq1(i,j,k) -                                    &
               at*( one + ctctq1(i,j,k+1) ) )
          dqw1(i,j,k) =  (dqw1(i,j,k) - at*dqw1(i,j,k+1) )
          dtl1(i,j,k) =  (dtl1(i,j,k) - at*dtl1(i,j,k+1) )
        end do
      end do

      call oneover_v(l, temp, temp_out)
      l = 0
      do j = tdims%j_start, tdims%j_end
        do i = ii, min(ii+tdims_seg_block -1, tdims%i_end)
          l = l + 1
          dqw1(i,j,k) = temp_out(l) * dqw1(i,j,k)
          dtl1(i,j,k) = temp_out(l) * dtl1(i,j,k)
          ctctq1(i,j,k) = temp_out(l) * ctctq1(i,j,k)
        end do
      end do
    end do !blm1,2,-1
  end do
!$OMP end do

!$OMP do SCHEDULE(STATIC)
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      r_sq = r_rho_levels(i,j,2)*r_rho_levels(i,j,2)
      dqw1(i,j,1) = -dtrdz_charney_grid(i,j,1) * (r_sq*fqw(i,j,2)) +           &
                    dqw_nt(i,j,1)
      dtl1(i,j,1) = -dtrdz_charney_grid(i,j,1) * (r_sq*ftl(i,j,2)) +           &
                    dtl_nt(i,j,1)
      at = -dtrdz_charney_grid(i,j,1) *                                        &
                 gamma_in(2)*(r_sq*rhokh(i,j,2))*rdz_charney_grid(i,j,2)
      rbt = one / ( one - at*( one + ctctq1(i,j,2) ) )
      dqw1(i,j,1) = rbt * (dqw1(i,j,1) - at*dqw1(i,j,2) )
      dtl1(i,j,1) = rbt * (dtl1(i,j,1) - at*dtl1(i,j,2) )

      ! Now set CT_CTQ(1) to be r^2 * BETA
      r_sq = r_theta_levels(i,j,0)*r_theta_levels(i,j,0)
      ctctq1(i,j,1) = - (r_sq * dtrdz_charney_grid(i,j,1)) * rbt
    end do
  end do
!$OMP end do

  !-----------------------------------------------------------------------
  !! Calculate the coefficients for calculation of implicit heat and
  !! moisture fluxes, i.e. Upward sweep of matrix (to k_blend)
  !! NOTE: Gives coeffs for coupling to the bottom model level
  !! when k_blend = 1 (i.e. when blend_height_opt switch set to blend_level1)
  !-----------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
  do j = tdims%j_start,tdims%j_end
    do i = tdims%i_start,tdims%i_end

      dtl1_1(i,j) = dtl1(i,j,k_blend_tq(i,j))
      dqw1_1(i,j) = dqw1(i,j,k_blend_tq(i,j))
      ct_prod(i,j) = ctctq1(i,j,k_blend_tq(i,j))

      do k = k_blend_tq(i,j)-1, 1, -1

        dtl1_1(i,j) = dtl1_1(i,j) + ( (-1) ** (k_blend_tq(i,j)+k) ) *          &
                     dtl1(i,j,k) * ct_prod(i,j)

        dqw1_1(i,j) = dqw1_1(i,j) + ( (-1) ** (k_blend_tq(i,j)+k) ) *          &
                     dqw1(i,j,k) * ct_prod(i,j)

        ct_prod(i,j) = ct_prod(i,j) * ctctq1(i,j,k)

      end do

      ctctq1_1(i,j) = ( (-1) ** ( k_blend_tq(i,j) + 1 ) ) *                    &
                      ct_prod(i,j)
    end do
  end do
!$OMP end do

else
  !-----------------------------------------------------------------------
  ! The following calculations complete the downward sweep for  the surface
  ! scalar variables. They apply on the 2nd stage of the scheme.
  ! The equivalent calculations for the 1st stage are done at the
  ! beginning of bdy_impl4, since at this stage (bdy_impl3) the surface
  ! scalar fluxes are not fully updated. This will be done by sf_impl2().
  !-----------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      r_sq = r_rho_levels(i,j,1)*r_rho_levels(i,j,1)
      rr_sq = r_rho_levels(i,j,2)*r_rho_levels(i,j,2)
      dqw(i,j,1) = gamma2(i,j) * ( -dtrdz_charney_grid(i,j,1) *                &
          ((rr_sq*fqw(i,j,2)) - (r_sq*fqw(i,j,1))) + dqw_nt(i,j,1) )
      dtl(i,j,1) = gamma2(i,j) * ( -dtrdz_charney_grid(i,j,1) *                &
          ((rr_sq*ftl(i,j,2)) - (r_sq*ftl(i,j,1))) + dtl_nt(i,j,1) )
      at = -dtrdz_charney_grid(i,j,1) *                                        &
            gamma1(i,j)*(rr_sq*rhokh(i,j,2))*rdz_charney_grid(i,j,2)
      rbt = one / ( one - at*( one + ct_ctq(i,j,2) ) )
      dqw(i,j,1) = rbt * (dqw(i,j,1) - at*dqw(i,j,2) )
      dtl(i,j,1) = rbt * (dtl(i,j,1) - at*dtl(i,j,2) )

      ! Now set CT_CTQ(1) to be r^2 * BETA
      r_sq = r_theta_levels(i,j,0)*r_theta_levels(i,j,0)
      ct_ctq(i,j,1) = - (r_sq * dtrdz_charney_grid(i,j,1)) * rbt
    end do
  end do
!$OMP end do

end if
!-----------------------------------------------------------------------
! 4.0 Return fluxes back to their true value by dividing by r*r
!-----------------------------------------------------------------------

!$OMP end PARALLEL

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine bdy_impl3
end module bdy_impl3_mod
