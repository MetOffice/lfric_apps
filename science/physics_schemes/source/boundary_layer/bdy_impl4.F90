! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!  Subroutine BDY_IMPL4

!  Purpose: Calculate implicit correction to boundary layer fluxes of
!           heat, moisture and momentum for the unconditionally
!           stable and non-oscillatory numerical solver.

!  Programming standard: UMDP3

!  Documentation: UMDP24

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module bdy_impl4_mod

use um_types, only: real_umphys, r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'BDY_IMPL4_MOD'
contains

subroutine bdy_impl4 (                                                         &
! in levels, switches
 bl_levels, l_correct,                                                         &
! in data :
 gamma1,gamma2,rdz_charney_grid, r_rho_levels,                                 &
 dtrdz_charney_grid,ct_ctq,dqw_nt,dtl_nt,                                      &
! INOUT data :
 qw,tl,fqw,ftl,fqw_star,ftl_star,                                              &
 dqw,dtl, rhokh, BL_diag,                                                      &
! out data, NB these are really tl and qt on exit!
 t_latest,q_latest                                                             &
 )

use atm_fields_bounds_mod, only: tdims, pdims, tdims_l
use bl_diags_mod, only: strnewbldiag
use model_domain_mod, only: model_type, mt_single_column
use planet_constants_mod, only: cp => cp_bl
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
!$ use omp_lib, only: omp_get_max_threads

implicit none

!  Inputs :-
integer, intent(in) ::                                                         &
 bl_levels                   ! in Max. no. of "boundary" levels
!                                     allowed.
logical, intent(in) ::                                                         &
 l_correct

real(kind=r_bl), intent(in) ::                                                 &
 rdz_charney_grid(tdims%i_start:tdims%i_end,                                   &
                  bl_levels),                                                  &
                                 ! in RDZ(,1) is the reciprocal of the
                                 ! height of level 1, i.e. of the
                                 ! middle of layer 1.  For K > 1,
                                 ! RDZ(,K) is the reciprocal
                                 ! of the vertical distance
                                 ! from level K-1 to level K.
 r_rho_levels(tdims_l%i_start:tdims_l%i_end,                                   &
              1:bl_levels),                                                    &
                                 ! in height of rho levels
 dtrdz_charney_grid(tdims%i_start:tdims%i_end,bl_levels),                      &
 gamma1(tdims%i_start:tdims%i_end),                                            &
 gamma2(tdims%i_start:tdims%i_end),                                            &
                                 ! in new scheme weights.
 ct_ctq(tdims%i_start:tdims%i_end,                                             &
        bl_levels),                                                            &
                                 ! in Coefficient in T and q
                                 !       tri-diagonal implicit matrix
 dqw_nt(tdims%i_start:tdims%i_end,                                             &
        bl_levels),                                                            &
                                      ! in NT incr to qw
 dtl_nt(tdims%i_start:tdims%i_end,                                             &
        bl_levels)
                                      ! in NT incr to TL

!  In/outs :-
!     Declaration of BL diagnostics.
type (strnewbldiag), intent(in) :: BL_diag

real(kind=r_bl), intent(in out) ::                                             &
 rhokh(pdims%i_start:pdims%i_end,                                              &
       bl_levels),                                                             &
                                 ! INOUT Exchange coeffs for moisture.
                                 ! shouldnt change but is scaled by
                                 ! r_sq then back
 qw(tdims%i_start:tdims%i_end,bl_levels),                                      &
                                 ! INOUT Total water content, but
                                 !       replaced by specific
                                 !       humidity in LS_CLD.
 tl(tdims%i_start:tdims%i_end,bl_levels),                                      &
                                 ! INOUT Ice/liquid water temperature,
                                 !       but replaced by T in LS_CLD.
 fqw(tdims%i_start:tdims%i_end,bl_levels),                                     &
                                 ! INOUT Moisture flux between layers
                                 !       (kg per square metre per sec)
                                 !       FQW(,1) is total water flux
                                 !       from surface, 'E'.
 ftl(tdims%i_start:tdims%i_end,bl_levels),                                     &
                                 ! INOUT FTL(,K) contains net
                                 !       turbulent sensible heat flux
                                 !       into layer K from below; so
                                 !       FTL(,1) is the surface
                                 !       sensible heat, H. (W/m2)
                                   ! INOUT temp arrays for diags
   fqw_star(tdims%i_start:tdims%i_end,                                         &
            bl_levels),                                                        &
   ftl_star(tdims%i_start:tdims%i_end,                                         &
            bl_levels),                                                        &
   dqw(tdims%i_start:tdims%i_end,bl_levels),                                   &
                                   ! INOUT BL increment to q field
   dtl(tdims%i_start:tdims%i_end,bl_levels)
                                   ! INOUT BL increment to T field

! out fields
real(kind=real_umphys), intent(out) ::                                         &
 q_latest(tdims%i_start:tdims%i_end,                                           &
          bl_levels),                                                          &
      ! out specific humidity
      ! But at this stage it is qT = qv+qcl+qcf
 t_latest(tdims%i_start:tdims%i_end,                                           &
          bl_levels)
      ! out temperature
      ! But at this stage it is tL = t - lcrcp*qcl - lsrcp*qcf

!-----------------------------------------------------------------------
!  Local scalars :-

real(kind=r_bl) :: r_sq, rbt, at ,                                             &
 gamma1_uv,                                                                    &
 gamma2_uv      ! gamma1 and gamma2 shifted to u or v points

integer ::                                                                     &
 i,                                                                            &
            ! LOCAL Loop counter (horizontal field index).
 k          ! LOCAL Loop counter (vertical level index).

integer ::                                                                     &
 ii,                                                                           &
             ! omp blocking index
 tdims_omp_block,                                                              &
             ! omp block length
 tdims_seg_block,                                                              &
             ! omp segment length
 max_threads ! store the no of threads

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='BDY_IMPL4'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

max_threads = 1
!$ max_threads = omp_get_max_threads()
tdims_omp_block  = ceiling(real(tdims%i_len)/max_threads)
tdims_seg_block = min(tdims_omp_block, tdims%i_len)

!$OMP  PARALLEL DEFAULT(SHARED) private(i,k,ii,at,rbt,gamma1_uv,               &
!$OMP  gamma2_uv,r_sq)
if ( .not. l_correct ) then
  !  1st stage: predictor
  !---------------------------------------------------------------------

  ! Complete downward sweep of matrix for increments to TL and QW in the
  ! boundary layer. It needs to be done here since the scalar fluxes
  ! (FQW(,,1),FTL(,,1)) computed previously by the surface scheme are not
  ! available when the main downward sweep subroutine bdy_impl3 is first
  ! invoked (at the 1st stage of the new solver).

  !---------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
  do i = tdims%i_start, tdims%i_end
  r_sq = r_rho_levels(i,1)*r_rho_levels(i,1)
  rhokh(i,2) = r_sq * rhokh(i,2)
  fqw(i,1)   = r_sq * fqw(i,1)
  ftl(i,1)   = r_sq * ftl(i,1)
  fqw(i,2)   = r_sq * fqw(i,2)
  ftl(i,2)   = r_sq * ftl(i,2)
  dqw(i,1) = gamma2(i) * ( -dtrdz_charney_grid(i,1) *                          &
    ( fqw(i,2) - fqw(i,1) ) + dqw_nt(i,1) )
  dtl(i,1) = gamma2(i) * ( -dtrdz_charney_grid(i,1) *                          &
    ( ftl(i,2) - ftl(i,1) ) + dtl_nt(i,1) )
  at = -dtrdz_charney_grid(i,1) *                                              &
        gamma1(i)*rhokh(i,2)*rdz_charney_grid(i,2)
  rbt = 1.0_r_bl / ( 1.0_r_bl - at*( 1.0_r_bl + ct_ctq(i,2) ) )
  dqw(i,1) = rbt*(dqw(i,1) - at*dqw(i,2) )
  dtl(i,1) = rbt*(dtl(i,1) - at*dtl(i,2) )
  rhokh(i,2) = rhokh(i,2)/r_sq
  fqw(i,1) = fqw(i,1)/r_sq
  ftl(i,1) = ftl(i,1)/r_sq
  fqw(i,2) = fqw(i,2)/r_sq
  ftl(i,2) = ftl(i,2)/r_sq
  end do
!$OMP end do

else ! L_CORRECT == true: 2nd stage of the scheme

  ! Compute 2nd stage correction (total: from tn to tn+1)
end if

! Update TL, QW and their increments
!$OMP do SCHEDULE(STATIC)
do i = tdims%i_start, tdims%i_end
  tl(i,1) = tl(i,1) + dtl(i,1)
  qw(i,1) = qw(i,1) + dqw(i,1)
end do
!$OMP end do

!$OMP do SCHEDULE(STATIC)
 do ii = tdims%i_start, tdims%i_end, tdims_seg_block
  do k = 2, bl_levels
    do i = ii, min(ii+tdims_seg_block-1,tdims%i_end)
      dtl(i,k) = dtl(i,k) - ct_ctq(i,k)*dtl(i,k-1)
      tl(i,k) = tl(i,k) + dtl(i,k)
      dqw(i,k) = dqw(i,k) - ct_ctq(i,k)*dqw(i,k-1)
      qw(i,k) = qw(i,k) + dqw(i,k)
    end do
  end do !bl_levels
end do
!$OMP end do

! Calculate stress and flux diagnostics using the new scheme equations.
! The fluxes are calculated only when requested.
if ( BL_diag%l_ftl ) then

  if ( .not. l_correct ) then

!$OMP do SCHEDULE(STATIC)
    do k = 2, bl_levels
      do i = tdims%i_start, tdims%i_end
        ftl_star(i,k) = gamma2(i)*ftl(i,k)                                     &
              - gamma1(i)*rhokh(i,k)*rdz_charney_grid(i,k)                     &
                          * (dtl(i,k)-dtl(i,k-1))
      end do
    end do
!$OMP end do NOWAIT
  else

!$OMP do SCHEDULE(STATIC)
    do k = 2, bl_levels
      do i = tdims%i_start, tdims%i_end
        ftl(i,k) = ftl_star(i,k)+gamma2(i)*ftl(i,k)                            &
              - gamma1(i)*rhokh(i,k)*rdz_charney_grid(i,k)                     &
                          * (dtl(i,k)-dtl(i,k-1))
      end do
    end do
!$OMP end do NOWAIT

  end if ! L_correct
end if ! L_ftl

if ( BL_diag%l_fqw ) then

  if ( .not. l_correct ) then

!$OMP do SCHEDULE(STATIC)
    do k = 2, bl_levels
      do i = tdims%i_start, tdims%i_end
        fqw_star(i,k) = gamma2(i)*fqw(i,k)                                     &
              - gamma1(i)*rhokh(i,k)*rdz_charney_grid(i,k)                     &
              * (dqw(i,k)-dqw(i,k-1))
      end do
    end do ! bl_levels
!$OMP end do NOWAIT

  else

!$OMP do SCHEDULE(STATIC)
    do k = 2, bl_levels
      do i = tdims%i_start, tdims%i_end
        fqw(i,k) = fqw_star(i,k)+gamma2(i)*fqw(i,k)                            &
              - gamma1(i)*rhokh(i,k)*rdz_charney_grid(i,k)                     &
                * (dqw(i,k)-dqw(i,k-1))
      end do
    end do ! bl_levels
!$OMP end do NOWAIT

  end if
end if

if ( l_correct ) then

  !-----------------------------------------------------------------------
  !     Convert FTL to sensible heat flux in Watts per square metre.
  !      Also, IMPL_CAL only updates FTL_TILE(*,1) and FQW_TILE(*,1)
  !      over sea points, so copy this to remaining tiles
  !-----------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do i = tdims%i_start, tdims%i_end
      ftl(i,k) = ftl(i,k)*cp
    end do
  end do
!$OMP end do NOWAIT

  !Copy T and Q from workspace to INOUT space.

!$OMP do SCHEDULE(STATIC)
  do k = 1, bl_levels
    do i = tdims%i_start, tdims%i_end
      t_latest(i,k)=tl(i,k)
      q_latest(i,k)=qw(i,k)
    end do
  end do
!$OMP end do NOWAIT

end if ! L_CORRECT

!$OMP end PARALLEL

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine bdy_impl4
end module bdy_impl4_mod
