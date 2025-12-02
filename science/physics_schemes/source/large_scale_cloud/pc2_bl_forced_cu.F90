! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! PC2 code for forced cumulus clouds

module pc2_bl_forced_cu_mod

implicit none
!---------------------------------------------------------------------------
! Description: Called from ni_imp_ctl to do pc2 forced cumulus clouds

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_cloud
!---------------------------------------------------------------------------

character(len=*), parameter, private :: ModuleName='PC2_BL_FORCED_CU_MOD'

contains

subroutine pc2_bl_forced_cu( zhnl, dzh, zlcl, bl_type_3, bl_type_6,            &
                             z_theta, qcl_inv_top,                             &
                             cca0, ccw0, ccb0, cct0, lcbase0,                  &
                             cfl_latest, cf_latest,                            &
                             qcl_latest, q_latest, t_latest, l_wtrac_bl)

use atm_fields_bounds_mod, only: tdims, pdims
use bl_option_mod,         only: kprof_cu, on
use cloud_inputs_mod,      only: forced_cu_fac, forced_cu
use mphys_constants_mod,   only: mprog_min
use pc2_constants_mod,     only: cbl_and_cu, forced_cu_cca
use planet_constants_mod,  only: lcrcp
use wtrac_pc2_mod,         only: wtrac_pc2
use yomhook,               only: lhook, dr_hook
use parkind1,              only: jprb, jpim

implicit none

real, intent(in) ::                                                            &
    zhnl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
                                  ! in non-local BL depth
    dzh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                  &
                                  ! in inversion thickness
    zlcl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
                                  ! in height of lifting condensation level
    bl_type_3(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                                  ! in Indicator set to 1.0 if well
!                                 !     mixed b.l. diagnosed,
!                                 !     0.0 otherwise.
    bl_type_6(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                                  ! in Indicator set to 1.0 if a
!                                 !     cumulus capped b.l. diagnosed,
!                                 !     0.0 otherwise.
  z_theta(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                 &
          pdims%k_start:pdims%k_end),                                          &
                                  ! in height of model theta levels
    qcl_inv_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                  ! in Parcel water content at inv top

! Convective cloud fraction and grid-mean water content;
! optionally forced_cu modifies these instead of directly
! altering cfl and qcl.
real, intent(in out) :: cca0(tdims%i_start:tdims%i_end,                        &
                             tdims%j_start:tdims%j_end,                        &
                             1:tdims%k_end)
real, intent(in out) :: ccw0(tdims%i_start:tdims%i_end,                        &
                             tdims%j_start:tdims%j_end,                        &
                             1:tdims%k_end)

! Convective cloud base and top levels; need to be updated if force Cu
! adds convective cloud on a level where there was none calculated
! by the convection scheme.
integer, intent(in out) :: ccb0    (tdims%i_start:tdims%i_end,                 &
                                    tdims%j_start:tdims%j_end)
integer, intent(in out) :: cct0    (tdims%i_start:tdims%i_end,                 &
                                    tdims%j_start:tdims%j_end)
integer, intent(in out) :: lcbase0 (tdims%i_start:tdims%i_end,                 &
                                    tdims%j_start:tdims%j_end)

real, intent(in out) ::                                                        &
    cfl_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,            &
               tdims%k_end),                                                   &
               ! in out liquid cloud fraction current value to update
    cf_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,             &
              tdims%k_end),                                                    &
              ! in out bulk cloud fraction current value to update
    qcl_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,            &
               tdims%k_end),                                                   &
               ! in out liquid cloud water content current value to update
    t_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,              &
             tdims%k_end),                                                     &
             ! in out temperature current value to update
    q_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,              &
             tdims%k_end)
             ! in out water vapour content current value to update

logical, intent(in) :: l_wtrac_bl    ! Controls water tracer storage

integer :: i, j, k

real ::                                                                        &
 qcl_forced,                                                                   &
 cf_forced,                                                                    &
            ! forced cloud water content and fraction
 dqcl,                                                                         &
 dcfl,                                                                         &
            ! forced cloud water content and fraction increments
 qcl_tol,                                                                      &
            ! max tolerated forced cloud water content
 cf_base,                                                                      &
            ! forced cloud fraction at cloud base
 zc_depth
            ! forced cloud depth

real, parameter :: qcl_forced_min = 0.00005
                 ! minimum water content in forced cu clouds

real, parameter :: qcl_max_factor = 0.1
                 ! maximum fraction of water vapour forced cu is
                 ! allowed to condense

real, parameter :: cf_top = 0.1
                 ! forced cloud fraction at cloud top

character(len=*), parameter :: RoutineName='PC2_BL_FORCED_CU'

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle
!----------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!----------------------------------------------------------------------

!$OMP PARALLEL do SCHEDULE(DYNAMIC) DEFAULT(none)                              &
!$OMP SHARED( tdims, zhnl, dzh, zlcl, bl_type_3, z_theta,                      &
!$OMP         cfl_latest, qcl_inv_top, cf_latest, qcl_latest,                  &
!$OMP         q_latest, forced_cu_fac, t_latest, lcrcp,                        &
!$OMP         forced_cu, cca0, ccw0, cct0, ccb0, lcbase0, l_wtrac_bl,          &
!$OMP         wtrac_pc2)                                                       &
!$OMP private ( i, j, k, zc_depth, cf_base, cf_forced, qcl_forced,             &
!$OMP           dqcl, qcl_tol, dcfl )
do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end
    zc_depth = zhnl(i,j)+dzh(i,j)-zlcl(i,j)
    if ( zc_depth > 1.0 .and. bl_type_3(i,j)>0.5 ) then

      do k = 1, tdims%k_end
        if ( z_theta(i,j,k) >= zlcl(i,j)                                       &
             .and. z_theta(i,j,k) <= zhnl(i,j)+dzh(i,j)                        &
             .and. cfl_latest(i,j,k) < 0.5                                     &
             ) then
          ! Make cloud fraction at cloud base a function of
          ! the cloud depth:
          cf_base = max( cf_top,                                               &
               0.3 * min( 1.0, zc_depth/300.0 ) )
          cf_forced = cf_base - (cf_base-cf_top) *                             &
               (z_theta(i,j,k)-zlcl(i,j)) / zc_depth
          ! Use diagnostic parcel qcl for in-cloud water content
          qcl_forced = max( qcl_forced_min,                                    &
               qcl_inv_top(i,j)*(z_theta(i,j,k)-zlcl(i,j))/zc_depth )
          qcl_forced = forced_cu_fac*qcl_forced ! tuning knob
          qcl_forced = qcl_forced*cf_forced ! GBM water content
          !-----------------------------------------------------
          ! calculate resulting changes in moisture variables
          ! and check for low values
          !-----------------------------------------------------
          if ( forced_cu == forced_cu_cca ) then
            ! Option to include forced cumulus in the diagnosed convective
            ! cloud passed to radiation, instead of the "large-scale" cloud

            ! a) Impose appropriate limits on qcl_forced
            qcl_forced = min( qcl_forced,                                      &
                              min( qcl_max_factor*q_latest(i,j,k),             &
                                   q_latest(i,j,k) - mprog_min ) )
            ! b) Combine the forced cumulus with the diagnosed convective
            !    cloud, assuming maximal overlap:
            ! (note: converting ccw to grid-mean and then back to in-cloud)
            ccw0(i,j,k) = max( ccw0(i,j,k)*cca0(i,j,k), qcl_forced )
            cca0(i,j,k) = max( cca0(i,j,k), cf_forced )
            ccw0(i,j,k) = ccw0(i,j,k) / cca0(i,j,k)
            ! c) Update the base and top model-levels if needed...
            ! When there are multiple layers of convective cloud, ccb0
            ! is the base of the highest layer, and lcbase0 is the base of the
            ! lowest layer.  Assuming forced_cu will only affect the lowest
            ! layer, so only update ccb0 if it is unset, or if only 1 layer.
            if ( ccb0(i,j)==0 .or. ( k < ccb0(i,j) .and.                       &
                                     ccb0(i,j)==lcbase0(i,j) ) ) then
              ccb0(i,j) = k
            end if
            if ( k > cct0(i,j) )  cct0(i,j) = k
            if ( lcbase0(i,j)==0 .or. k < lcbase0(i,j) )  lcbase0(i,j) = k

          else  ! ( .not. forced_cu == forced_cu_cca )
            ! Include forced_cu in the prognostic large-scale cloud fields...

            dqcl = 0.0
            if ( qcl_forced > qcl_latest(i,j,k) ) then
              dqcl = qcl_forced - qcl_latest(i,j,k)
              ! make sure we won't condense more than a fraction of
              ! the existing water vapour
              qcl_tol = qcl_max_factor*q_latest(i,j,k)
              if ( qcl_forced > qcl_tol ) then
                qcl_forced = qcl_tol
                dqcl = max( 0.0, qcl_forced - qcl_latest(i,j,k) )
                qcl_forced = qcl_latest(i,j,k) + dqcl
              end if
              if (q_latest(i,j,k)-dqcl < mprog_min) then
                ! almost no water in this grid-point so leave it!
                qcl_forced = qcl_latest(i,j,k)
                dqcl = 0.0
              end if
              qcl_latest(i,j,k) = qcl_forced
              t_latest(i,j,k)  = t_latest(i,j,k) + lcrcp*dqcl
              q_latest(i,j,k)  = q_latest(i,j,k) - dqcl
              if (l_wtrac_bl) then
                wtrac_pc2%q_cond(i,j,k) = dqcl
              end if

            end if ! test on qcl_forced

            if ( cf_forced > cfl_latest(i,j,k) .and.                           &
                 dqcl > qcl_max_factor*mprog_min ) then
              dcfl = cf_forced - cfl_latest(i,j,k)
              cfl_latest(i,j,k) = cf_forced
              cf_latest(i,j,k)  = cf_latest(i,j,k) + dcfl
            end if

          end if  ! ( .not. forced_cu == forced_cu_cca )
        end if  ! test on z
      end do  ! loop over k

    end if  ! test on zc_depth and bl_type3
  end do
end do
!$OMP end PARALLEL do

if ( kprof_cu >= on .and. ( forced_cu == cbl_and_cu                            &
                       .or. forced_cu == forced_cu_cca ) ) then

!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP SHARED( tdims, zhnl, zlcl, bl_type_6, z_theta, cfl_latest,               &
!$OMP         qcl_inv_top, forced_cu_fac, qcl_latest, q_latest,                &
!$OMP         t_latest, lcrcp, cf_latest,                                      &
!$OMP         forced_cu, cca0, ccw0, cct0, ccb0, lcbase0, l_wtrac_bl,          &
!$OMP         wtrac_pc2)                                                       &
!$OMP private( i, j, k, zc_depth, cf_base, cf_forced, qcl_forced, dqcl,        &
!$OMP          qcl_tol, dcfl )
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      zc_depth = zhnl(i,j)-zlcl(i,j)
      if ( zc_depth > 1.0 .and. bl_type_6(i,j)>0.5 ) then
        do k = 1, tdims%k_end
          !------------------------------------------------------------
          ! check cumulus cloud chunky enough above LCL
          !  - if kprof_cu used then zhnl penetrates above the LCL
          !    and so marks the top of "forced" cumulus
          !------------------------------------------------------------
          if ( z_theta(i,j,k) >= zlcl(i,j)                                     &
               .and. z_theta(i,j,k) <= zhnl(i,j)                               &
               .and. cfl_latest(i,j,k) < 0.5                                   &
               ) then
            ! Make cloud fraction at cloud base a function of
            ! the cloud depth:
            cf_base = max( cf_top,                                             &
                 0.3 * min( 1.0, zc_depth/300.0 ) )
            cf_forced = cf_base - (cf_base-cf_top) *                           &
                 (z_theta(i,j,k)-zlcl(i,j)) / zc_depth
            ! Use diagnostic parcel qcl for in-cloud water content
            ! For cumulus layers this is taken at zlcl+300m
            qcl_forced = max( qcl_forced_min,                                  &
                 qcl_inv_top(i,j)*(z_theta(i,j,k)-zlcl(i,j))/300.0 )
            qcl_forced = forced_cu_fac*qcl_forced !tuning knob
            qcl_forced = qcl_forced*cf_forced ! GBM water content
            !-----------------------------------------------------
            ! calculate resulting changes in moisture variables
            ! and check for low values
            !-----------------------------------------------------
            if ( forced_cu == forced_cu_cca ) then
              ! Option to include forced cumulus in the diagnosed convective
              ! cloud passed to radiation, instead of the "large-scale" cloud

                ! a) Impose appropriate limits on qcl_forced
              qcl_forced = min( qcl_forced,                                    &
                                min( qcl_max_factor*q_latest(i,j,k),           &
                                     q_latest(i,j,k) - mprog_min ) )
              ! b) Combine the forced cumulus with the diagnosed convective
              !    cloud, assuming maximal overlap:
              ! (note: converting ccw to grid-mean and then back to in-cloud)
              ccw0(i,j,k) = max( ccw0(i,j,k)*cca0(i,j,k), qcl_forced )
              cca0(i,j,k) = max( cca0(i,j,k), cf_forced )
              ccw0(i,j,k) = ccw0(i,j,k) / cca0(i,j,k)
              ! c) Update the base and top model-levels if needed...
              ! When there are multiple layers of convective cloud, ccb0
              ! is the base of the highest layer, and lcbase0 is the base of the
              ! lowest layer.  Assuming forced_cu will only affect the lowest
              ! layer, so only update ccb0 if it is unset, or if only 1 layer.
              if ( ccb0(i,j)==0 .or. ( k < ccb0(i,j) .and.                     &
                                       ccb0(i,j)==lcbase0(i,j) ) ) then
                ccb0(i,j) = k
              end if
              if ( k > cct0(i,j) )  cct0(i,j) = k
              if ( lcbase0(i,j)==0 .or. k < lcbase0(i,j) )  lcbase0(i,j) = k

            else  ! ( .not. forced_cu == forced_cu_cca )
              ! Include forced_cu in the prognostic large-scale cloud fields...

              dqcl = 0.0
              if ( qcl_forced > qcl_latest(i,j,k) ) then
                dqcl = qcl_forced - qcl_latest(i,j,k)
                ! make sure we won't condense more than a fraction of
                ! the existing water vapour
                qcl_tol = qcl_max_factor*q_latest(i,j,k)
                if ( qcl_forced > qcl_tol ) then
                  qcl_forced = qcl_tol
                  dqcl = max( 0.0, qcl_forced - qcl_latest(i,j,k) )
                  qcl_forced = qcl_latest(i,j,k) + dqcl
                end if
                if (q_latest(i,j,k)-dqcl < mprog_min) then
                  ! almost no water in this grid-point so leave it!
                  qcl_forced = qcl_latest(i,j,k)
                  dqcl = 0.0
                end if
                qcl_latest(i,j,k) = qcl_forced
                t_latest(i,j,k)  = t_latest(i,j,k) + lcrcp*dqcl
                q_latest(i,j,k)  = q_latest(i,j,k) - dqcl
                if (l_wtrac_bl) then
                  wtrac_pc2%q_cond(i,j,k) = wtrac_pc2%q_cond(i,j,k) + dqcl
                end if
              end if ! test on qcl_forced

              if ( cf_forced > cfl_latest(i,j,k) .and.                         &
                   dqcl > qcl_max_factor*mprog_min  ) then
                dcfl = cf_forced - cfl_latest(i,j,k)
                cfl_latest(i,j,k) = cf_forced
                cf_latest(i,j,k)  = cf_latest(i,j,k) + dcfl
              end if

            end if  ! ( .not. forced_cu == forced_cu_cca )
          end if  ! test on z

        end do  ! loop over k
      end if  ! test on zc_depth and bltype6
    end do
  end do
!$OMP end PARALLEL do

end if  ! test on forced_cu eq cbl_and_cu

!----------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!----------------------------------------------------------------------
return
end subroutine pc2_bl_forced_cu

end module pc2_bl_forced_cu_mod
