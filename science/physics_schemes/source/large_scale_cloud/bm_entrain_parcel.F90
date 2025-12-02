! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_cloud

module bm_entrain_parcel_mod

implicit none

character(len=*), parameter, private :: ModuleName = 'BM_ENTRAIN_PARCEL_MOD'

contains

! Subroutine to calculate Tl, qt and turbulence properties following
! entraining parcels going up and down, for use in setting the modes
! from above and below in the bimodal cloud scheme
subroutine bm_entrain_parcel( nlevels,                                         &
                              zh, zhsc, dzh, bl_type_7,                        &
                              z_theta, bl_w_var, tau_dec_bm, tau_hom_bm,       &
                              mix_len_bm, tl_in, qt_in,                        &
                              tl_below, qt_below, wvar_below, tau_dec_below,   &
                              tau_hom_below, dtldz_below, dqtdz_below,         &
                              tl_above, qt_above, wvar_above, tau_dec_above,   &
                              tau_hom_above, dtldz_above, dqtdz_above )

use um_types,              only: real_umphys
use atm_fields_bounds_mod, only: tdims
use yomhook,               only: lhook, dr_hook
use parkind1,              only: jprb, jpim
use planet_constants_mod,  only: g, cp
use cloud_inputs_mod,      only: ent_coef_bm

implicit none

! Number of model-levels
integer, intent(in) :: nlevels

! Boundary-layer height
real(kind=real_umphys), intent(in) :: zh                                       &
                                      ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end )
! Decoupled layer height
real(kind=real_umphys), intent(in) :: zhsc                                     &
                                      ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end )
! Inversion Thickness
real(kind=real_umphys), intent(in) :: dzh                                      &
                                      ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end )
! Shear-driven boundary layer indicator
real(kind=real_umphys), intent(in) :: bl_type_7                                &
                                      ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end )

! Model-level heights above surface
real(kind=real_umphys), intent(in) :: z_theta                                  &
                                      ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end, nlevels )

! Turbulent vertical velocity variance
real(kind=real_umphys), intent(in) :: bl_w_var                                 &
                                      ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end, nlevels )

! Turbulent decorrelation and homogenisation timescales
real(kind=real_umphys), intent(in) :: tau_dec_bm                               &
                                      ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys), intent(in) :: tau_hom_bm                               &
                                      ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end, nlevels )

! Turbulent mixing length
real(kind=real_umphys), intent(in) :: mix_len_bm                               &
                                      ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end, nlevels )

! Liquid water temperature (T - Lc/cp qcl) and total-water content (q + qcl)
real(kind=real_umphys), intent(in) :: tl_in                                    &
                                      ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys), intent(in) :: qt_in                                    &
                                      ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end, nlevels )

! Properties of entraining parcel lifted from below
real(kind=real_umphys), intent(out) :: tl_below                                &
                                        ( tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys), intent(out) :: qt_below                                &
                                        ( tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys), intent(out) :: wvar_below                              &
                                        ( tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys), intent(out) :: tau_dec_below                           &
                                        ( tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys), intent(out) :: tau_hom_below                           &
                                        ( tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys), intent(out) :: dtldz_below                             &
                                        ( tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys), intent(out) :: dqtdz_below                             &
                                        ( tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end, nlevels )

! Properties of entraining parcel descended from above
real(kind=real_umphys), intent(out) :: tl_above                                &
                                        ( tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys), intent(out) :: qt_above                                &
                                        ( tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys), intent(out) :: wvar_above                              &
                                        ( tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys), intent(out) :: tau_dec_above                           &
                                        ( tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys), intent(out) :: tau_hom_above                           &
                                        ( tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys), intent(out) :: dtldz_above                             &
                                        ( tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys), intent(out) :: dqtdz_above                             &
                                        ( tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end, nlevels )

! Height of top of surface-driven mixing
real(kind=real_umphys) :: zh_eff                                               &
                          ( tdims%i_start:tdims%i_end,                         &
                            tdims%j_start:tdims%j_end )

! Local turbulent standard deviations of Tl and qt
real(kind=real_umphys) :: sig_tl                                               &
                          ( tdims%i_start:tdims%i_end,                         &
                            tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys) :: sig_qt                                               &
                          ( tdims%i_start:tdims%i_end,                         &
                            tdims%j_start:tdims%j_end, nlevels )

! Decorrelation and homogenisation timescales converted to length-scales
real(kind=real_umphys) :: len_dec                                              &
                          ( tdims%i_start:tdims%i_end,                         &
                            tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys) :: len_hom                                              &
                          ( tdims%i_start:tdims%i_end,                         &
                            tdims%j_start:tdims%j_end, nlevels )

! Work variables
real(kind=real_umphys) :: z_rho, ent_frac, w1, w2, tmp

! Loop counters
integer :: i, j, k, km1, kp1

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real   (kind=jprb)            :: zhook_handle

character(len=*), parameter :: RoutineName='BM_ENTRAIN_PARCEL'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL                                                                 &
!$OMP DEFAULT(none)                                                            &
!$OMP SHARED( nlevels, tdims, zh, dzh, zhsc, bl_type_7, zh_eff,                &
!$OMP         bl_w_var, tau_dec_bm, tau_hom_bm, ent_coef_bm,                   &
!$OMP         z_theta, tl_in, qt_in, sig_tl, sig_qt, g, cp,                    &
!$OMP         tl_below, qt_below, wvar_below, tau_dec_below, tau_hom_below,    &
!$OMP         dtldz_below, dqtdz_below,                                        &
!$OMP         tl_above, qt_above, wvar_above, tau_dec_above, tau_hom_above,    &
!$OMP         dtldz_above, dqtdz_above,                                        &
!$OMP         len_dec, len_hom, mix_len_bm )                                   &
!$OMP private( i, j, k, km1, kp1, z_rho, ent_frac, w1, w2, tmp )

!$OMP  do SCHEDULE(STATIC)
do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end
    ! Calculate an effective PBL depth
    ! include resolved inversion thickness, dzh
    ! (which is set to rmdi when subgrid)
    zh_eff(i,j) = zh(i,j) + max( dzh(i,j), 0.0 )
    if ( bl_type_7(i,j) > 0.5_real_umphys ) then
      ! use cloud top (if higher) in decoupled PBLs that are
      ! shear-dominated
      zh_eff(i,j) = max( zh_eff(i,j), zhsc(i,j) )
    end if
    if ( zhsc(i,j) > zh_eff(i,j) ) then
      ! Also revert to cloud-top if zh is so close to cloud-top that the
      ! grid-level spanning zh will calculate local gradients which
      ! span the cloud-top.  This avoids spuriously applying large variances
      ! associated with the cloud-top inversion to ascending parcels
      ! within the cloud.
      do k = 2, nlevels-1
        if ( zh_eff(i,j) >= z_theta(i,j,k-1) .and.                             &
             zh_eff(i,j) <  z_theta(i,j,k) ) then
          if ( z_theta(i,j,k+1) > zhsc(i,j) )  zh_eff(i,j) = zhsc(i,j)
        end if
      end do
    end if
  end do
end do
!$OMP end do NOWAIT

! Pre-calculate turbulent length-scales and vertical gradients
!$OMP do SCHEDULE(STATIC)
do k = 1, nlevels
  km1 = max( k-1, 1 )
  kp1 = min( k+1, nlevels )
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end

      ! Turbulent length-scales
      tmp = sqrt(bl_w_var(i,j,k))
      len_dec(i,j,k) = tau_dec_bm(i,j,k) * tmp
      len_hom(i,j,k) = tau_hom_bm(i,j,k) * tmp

      ! Vertical gradients of Tl and qt, scaled by turbulent length-scale
      ! to scale with turbulent standard deviation
      tmp = 1.0 / ( z_theta(i,j,kp1) - z_theta(i,j,km1) )
      sig_tl(i,j,k) = len_dec(i,j,k)*( (tl_in(i,j,kp1) - tl_in(i,j,km1)) * tmp &
                                       + g/cp )
      sig_qt(i,j,k) = len_dec(i,j,k) * (qt_in(i,j,kp1) - qt_in(i,j,km1)) * tmp

    end do
  end do
end do
!$OMP end do

! Initialise entraining parcel values at first level
! Note that the tau_dec, tau_hom arrays will store equivalent length-scales
! for the entraining parcel calculations, to be converted to time-scales
! afterwards.  Similarly dtldz, dqtdz store turbulent standard deviations,
! to be converted back to gradients.
!$OMP do SCHEDULE(STATIC)
do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end
    ! Variables for ascent from below
    tl_below(i,j,1) = tl_in(i,j,1)
    qt_below(i,j,1) = qt_in(i,j,1)
    wvar_below(i,j,1) = bl_w_var(i,j,1)
    tau_dec_below(i,j,1) = len_dec(i,j,1)
    tau_hom_below(i,j,1) = len_hom(i,j,1)
    dtldz_below(i,j,1) = sig_tl(i,j,1)
    dqtdz_below(i,j,1) = sig_qt(i,j,1)
    ! Variables for descent from above
    tl_above(i,j,nlevels) = tl_in(i,j,nlevels)
    qt_above(i,j,nlevels) = qt_in(i,j,nlevels)
    wvar_above(i,j,nlevels) = bl_w_var(i,j,nlevels)
    tau_dec_above(i,j,nlevels) = len_dec(i,j,nlevels)
    tau_hom_above(i,j,nlevels) = len_hom(i,j,nlevels)
    dtldz_above(i,j,nlevels) = sig_tl(i,j,nlevels)
    dqtdz_above(i,j,nlevels) = sig_qt(i,j,nlevels)
  end do
end do
!$OMP end do

! Entraining parcel calculations...
!$OMP do SCHEDULE(STATIC)
do j = tdims%j_start, tdims%j_end
  ! j-loop outermost as k-loops are sequential and can't be parallelised

  ! Entraining parcel from below
  do k = 2, nlevels
    km1 = k-1
    do i = tdims%i_start, tdims%i_end
      z_rho = 0.5 * ( z_theta(i,j,km1) + z_theta(i,j,k) )

      ! Entrain from level k-1
      ent_frac = ent_coef_bm * ( z_rho - z_theta(i,j,km1) )                    &
                             / mix_len_bm(i,j,km1)
      w2 = ent_frac / ( 1.0 + ent_frac )
      w1 = 1.0 - w2
      tl_below(i,j,k)      = w1 * tl_below(i,j,km1)    + w2 * tl_in(i,j,km1)
      qt_below(i,j,k)      = w1 * qt_below(i,j,km1)    + w2 * qt_in(i,j,km1)
      wvar_below(i,j,k)    = w1 * wvar_below(i,j,km1)  + w2 * bl_w_var(i,j,km1)
      tau_dec_below(i,j,k) = w1*tau_dec_below(i,j,km1) + w2 * len_dec(i,j,km1)
      tau_hom_below(i,j,k) = w1*tau_hom_below(i,j,km1) + w2 * len_hom(i,j,km1)
      dtldz_below(i,j,k)   = w1 * dtldz_below(i,j,km1) + w2 * sig_tl(i,j,km1)
      dqtdz_below(i,j,k)   = w1 * dqtdz_below(i,j,km1) + w2 * sig_qt(i,j,km1)

      ! Adjust Tl following dry lapse rate
      tl_below(i,j,k) = tl_below(i,j,k)                                        &
                      - g/cp * ( z_theta(i,j,k) - z_theta(i,j,km1) )

      ! Entrain from level k
      ent_frac = ent_coef_bm  * ( z_theta(i,j,k) - z_rho )                     &
                              / mix_len_bm(i,j,k)
      w2 = ent_frac / ( 1.0 + ent_frac )
      w1 = 1.0 - w2
      tl_below(i,j,k)      = w1 * tl_below(i,j,k)    + w2 * tl_in(i,j,k)
      qt_below(i,j,k)      = w1 * qt_below(i,j,k)    + w2 * qt_in(i,j,k)
      wvar_below(i,j,k)    = w1 * wvar_below(i,j,k)  + w2 * bl_w_var(i,j,k)
      tau_dec_below(i,j,k) = w1*tau_dec_below(i,j,k) + w2 * len_dec(i,j,k)
      tau_hom_below(i,j,k) = w1*tau_hom_below(i,j,k) + w2 * len_hom(i,j,k)
      dtldz_below(i,j,k)   = w1 * dtldz_below(i,j,k) + w2 * sig_tl(i,j,k)
      dqtdz_below(i,j,k)   = w1 * dqtdz_below(i,j,k) + w2 * sig_qt(i,j,k)

    end do  ! i = tdims%i_start, tdims%i_end
    do i = tdims%i_start, tdims%i_end
      if ( zh_eff(i,j) >= z_theta(i,j,km1) .and.                               &
           zh_eff(i,j) <  z_theta(i,j,k) ) then
        ! When crossing the BL-top, recalculate parcel properties based on
        ! resetting to the properties interpolated to the BL-top...
        z_rho = 0.5 * ( z_theta(i,j,km1) + z_theta(i,j,k) )

        ! Find env properties at height zh_eff (with Tl adjusted to k-1)
        w2 = ( zh_eff(i,j) - z_theta(i,j,km1) )                                &
           / ( z_theta(i,j,k) - z_theta(i,j,km1) )
        w1 = 1.0 - w2
        tl_below(i,j,k)      = w1 * tl_in(i,j,km1)    + w2 * ( tl_in(i,j,k)    &
                              + g/cp * ( z_theta(i,j,k) - z_theta(i,j,km1) ) )
        qt_below(i,j,k)      = w1 * qt_in(i,j,km1)    + w2 * qt_in(i,j,k)
        wvar_below(i,j,k)    = w1 * bl_w_var(i,j,km1) + w2 * bl_w_var(i,j,k)
        tau_dec_below(i,j,k) = w1 * len_dec(i,j,km1)  + w2 * len_dec(i,j,k)
        tau_hom_below(i,j,k) = w1 * len_hom(i,j,km1)  + w2 * len_hom(i,j,k)
        dtldz_below(i,j,k)   = w1 * sig_tl(i,j,km1)   + w2 * sig_tl(i,j,k)
        dqtdz_below(i,j,k)   = w1 * sig_qt(i,j,km1)   + w2 * sig_qt(i,j,k)

        if ( zh_eff(i,j) < z_rho ) then
          ! If zh_eff is within theta-level k-1

          ! Entrain from level k-1
          ent_frac = ent_coef_bm * ( z_rho - zh_eff(i,j) )                     &
                                 / mix_len_bm(i,j,km1)
          w2 = ent_frac / ( 1.0 + ent_frac )
          w1 = 1.0 - w2
          tl_below(i,j,k)      = w1 * tl_below(i,j,k)    + w2 * tl_in(i,j,km1)
          qt_below(i,j,k)      = w1 * qt_below(i,j,k)    + w2 * qt_in(i,j,km1)
          wvar_below(i,j,k)    = w1 * wvar_below(i,j,k)  + w2*bl_w_var(i,j,km1)
          tau_dec_below(i,j,k) = w1*tau_dec_below(i,j,k) + w2*len_dec(i,j,km1)
          tau_hom_below(i,j,k) = w1*tau_hom_below(i,j,k) + w2*len_hom(i,j,km1)
          dtldz_below(i,j,k)   = w1 * dtldz_below(i,j,k) + w2 * sig_tl(i,j,km1)
          dqtdz_below(i,j,k)   = w1 * dqtdz_below(i,j,k) + w2 * sig_qt(i,j,km1)
        end if

        ! Adjust Tl following dry lapse rate
        tl_below(i,j,k) = tl_below(i,j,k)                                      &
                        - g/cp * ( z_theta(i,j,k) - z_theta(i,j,km1) )

        ! Entrain from level k
        ent_frac = ent_coef_bm * ( z_theta(i,j,k) - max(zh_eff(i,j), z_rho) )  &
                               / mix_len_bm(i,j,k)
        w2 = ent_frac / ( 1.0 + ent_frac )
        w1 = 1.0 - w2
        tl_below(i,j,k)      = w1 * tl_below(i,j,k)    + w2 * tl_in(i,j,k)
        qt_below(i,j,k)      = w1 * qt_below(i,j,k)    + w2 * qt_in(i,j,k)
        wvar_below(i,j,k)    = w1 * wvar_below(i,j,k)  + w2 * bl_w_var(i,j,k)
        tau_dec_below(i,j,k) = w1*tau_dec_below(i,j,k) + w2 * len_dec(i,j,k)
        tau_hom_below(i,j,k) = w1*tau_hom_below(i,j,k) + w2 * len_hom(i,j,k)
        dtldz_below(i,j,k)   = w1 * dtldz_below(i,j,k) + w2 * sig_tl(i,j,k)
        dqtdz_below(i,j,k)   = w1 * dqtdz_below(i,j,k) + w2 * sig_qt(i,j,k)

      end if  ! crossing from below to above zh_eff
    end do  ! i = tdims%i_start, tdims%i_end
  end do  ! k = 2, nlevels

  ! Entraining parcel from above
  do k = nlevels-1, 1, -1
    kp1 = k+1
    do i = tdims%i_start, tdims%i_end
      z_rho = 0.5 * ( z_theta(i,j,k) + z_theta(i,j,kp1) )

      ! Entrain from level k+1
      ent_frac = ent_coef_bm * ( z_theta(i,j,kp1) - z_rho )                    &
                             / mix_len_bm(i,j,kp1)
      w2 = ent_frac / ( 1.0 + ent_frac )
      w1 = 1.0 - w2
      tl_above(i,j,k)      = w1 * tl_above(i,j,kp1)    + w2 * tl_in(i,j,kp1)
      qt_above(i,j,k)      = w1 * qt_above(i,j,kp1)    + w2 * qt_in(i,j,kp1)
      wvar_above(i,j,k)    = w1 * wvar_above(i,j,kp1)  + w2 * bl_w_var(i,j,kp1)
      tau_dec_above(i,j,k) = w1*tau_dec_above(i,j,kp1) + w2 * len_dec(i,j,kp1)
      tau_hom_above(i,j,k) = w1*tau_hom_above(i,j,kp1) + w2 * len_hom(i,j,kp1)
      dtldz_above(i,j,k)   = w1 * dtldz_above(i,j,kp1) + w2 * sig_tl(i,j,kp1)
      dqtdz_above(i,j,k)   = w1 * dqtdz_above(i,j,kp1) + w2 * sig_qt(i,j,kp1)

      ! Adjust Tl following dry lapse rate
      tl_above(i,j,k) = tl_above(i,j,k)                                        &
                      + g/cp * ( z_theta(i,j,kp1) - z_theta(i,j,k) )

      ! Entrain from level k
      ent_frac = ent_coef_bm * ( z_rho - z_theta(i,j,k) )                      &
                             / mix_len_bm(i,j,k)
      w2 = ent_frac / ( 1.0 + ent_frac )
      w1 = 1.0 - w2
      tl_above(i,j,k)      = w1 * tl_above(i,j,k)    + w2 * tl_in(i,j,k)
      qt_above(i,j,k)      = w1 * qt_above(i,j,k)    + w2 * qt_in(i,j,k)
      wvar_above(i,j,k)    = w1 * wvar_above(i,j,k)  + w2 * bl_w_var(i,j,k)
      tau_dec_above(i,j,k) = w1*tau_dec_above(i,j,k) + w2 * len_dec(i,j,k)
      tau_hom_above(i,j,k) = w1*tau_hom_above(i,j,k) + w2 * len_hom(i,j,k)
      dtldz_above(i,j,k)   = w1 * dtldz_above(i,j,k) + w2 * sig_tl(i,j,k)
      dqtdz_above(i,j,k)   = w1 * dqtdz_above(i,j,k) + w2 * sig_qt(i,j,k)

    end do  ! i = tdims%i_start, tdims%i_end
  end do  ! k = nlevels-1, 1, -1

end do  ! j = tdims%j_start, tdims%j_end
!$OMP end do

!$OMP do SCHEDULE(STATIC)
do k = 1, nlevels
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      ! Convert entraining parcel Tl,qt sigmas back to local gradients
      tmp = 1.0 / tau_dec_below(i,j,k)
      dtldz_below(i,j,k) = dtldz_below(i,j,k) * tmp - g/cp
      dqtdz_below(i,j,k) = dqtdz_below(i,j,k) * tmp
      tmp = 1.0 / tau_dec_above(i,j,k)
      dtldz_above(i,j,k) = dtldz_above(i,j,k) * tmp - g/cp
      dqtdz_above(i,j,k) = dqtdz_above(i,j,k) * tmp
      ! Convert entraining parcel length-scales back to time-scales
      tmp = 1.0 / sqrt(wvar_below(i,j,k))
      tau_dec_below(i,j,k) = tau_dec_below(i,j,k) * tmp
      tau_hom_below(i,j,k) = tau_hom_below(i,j,k) * tmp
      tmp = 1.0 / sqrt(wvar_above(i,j,k))
      tau_dec_above(i,j,k) = tau_dec_above(i,j,k) * tmp
      tau_hom_above(i,j,k) = tau_hom_above(i,j,k) * tmp
    end do
  end do
end do
!$OMP end do NOWAIT

! Sanity-checks on the mean qw and Tl from above and below
!$OMP do SCHEDULE(STATIC)
do k = 1, nlevels
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      ! Shift tl_above or tl_below if needed in order to keep tl_in
      ! safely in-between the values from above and below
      if ( tl_above(i,j,k) > tl_below(i,j,k) ) then
        if ( tl_in(i,j,k) > 0.5 * (tl_below(i,j,k) + tl_above(i,j,k)) ) then
          !    tl_b  <  tl_a  <  tl_in
          ! => tl_b      <       tl_in  <  tl_a
          tl_above(i,j,k) = max( tl_above(i,j,k),                              &
             tl_in(i,j,k) + 0.01 * ( tl_in(i,j,k) - tl_below(i,j,k) ) )
        else
          !             tl_in  <  tl_b  <  tl_a
          ! => tl_b  <  tl_in      <       tl_a
          tl_below(i,j,k) = min( tl_below(i,j,k),                              &
             tl_in(i,j,k) - 0.01 * ( tl_above(i,j,k) - tl_in(i,j,k) ) )
        end if
      else
        if ( tl_in(i,j,k) > 0.5 * (tl_below(i,j,k) + tl_above(i,j,k)) ) then
          !   tl_a  <  tl_b  <  tl_in
          !   tl_a      <       tl_in  < tl_b
          tl_below(i,j,k) = max( tl_below(i,j,k),                              &
             tl_in(i,j,k) + 0.01 * ( tl_in(i,j,k) - tl_above(i,j,k) ) )
        else
          !             tl_in  <  tl_a  <  tl_b
          ! => tl_a  <  tl_in      <       tl_b
          tl_above(i,j,k) = min( tl_above(i,j,k),                              &
             tl_in(i,j,k) - 0.01 * ( tl_below(i,j,k) - tl_in(i,j,k) ) )
        end if
      end if
      ! Same again for qt
      if ( qt_above(i,j,k) > qt_below(i,j,k) ) then
        if ( qt_in(i,j,k) > 0.5 * (qt_below(i,j,k) + qt_above(i,j,k)) ) then
          !    qt_b  <  qt_a  <  qt_in
          ! => qt_b      <       qt_in  <  qt_a
          qt_above(i,j,k) = max( qt_above(i,j,k),                              &
             qt_in(i,j,k) + 0.01 * ( qt_in(i,j,k) - qt_below(i,j,k) ) )
        else
          !             qt_in  <  qt_b  <  qt_a
          ! => qt_b  <  qt_in      <       qt_a
          qt_below(i,j,k) = min( qt_below(i,j,k),                              &
             qt_in(i,j,k) - 0.01 * ( qt_above(i,j,k) - qt_in(i,j,k) ) )
        end if
      else
        if ( qt_in(i,j,k) > 0.5 * (qt_below(i,j,k) + qt_above(i,j,k)) ) then
          !   qt_a  <  qt_b  <  qt_in
          !   qt_a      <       qt_in  < qt_b
          qt_below(i,j,k) = max( qt_below(i,j,k),                              &
             qt_in(i,j,k) + 0.01 * ( qt_in(i,j,k) - qt_above(i,j,k) ) )
        else
          !             qt_in  <  qt_a  <  qt_b
          ! => qt_a  <  qt_in      <       qt_b
          qt_above(i,j,k) = min( qt_above(i,j,k),                              &
             qt_in(i,j,k) - 0.01 * ( qt_below(i,j,k) - qt_in(i,j,k) ) )
        end if
      end if
    end do
  end do
end do
!$OMP end do

!$OMP end PARALLEL

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine bm_entrain_parcel

end module bm_entrain_parcel_mod
