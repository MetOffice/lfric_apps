! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module bm_cld_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName='BM_CLD_MOD'

contains

! Large-scale Bimodal Cloud Scheme.
! Subroutine Interface:
! ======================================================================
! Large-scale Cloud Scheme Compression routine (Cloud points only).
! Subroutine Interface:
subroutine bm_cld( p_l, qsl_l, qsi_l, q_l, t_l, inv_thm_l, tdc_l, inv_tmp_l,   &
                   wvar_l, dtldz_l, dqtdz_l, ql_mean,                          &
                   qcl_f, qcf_f, cfl_f, cff_f, cf_f, cfl_max, q_f, t_f,        &
                   sskew, svar_turb, svar_bm, indx,points, l_mixing_ratio)

use water_constants_mod,  only: lc, lf
use planet_constants_mod, only: lcrcp, r, repsilon, g, cp
use yomhook,              only: lhook, dr_hook
use parkind1,             only: jprb, jpim
use atm_fields_bounds_mod,only: tdims, pdims
use conversions_mod,      only: pi
use umErf_mod,            only: umErf
use cloud_inputs_mod,     only: ice_fraction_method, smith_orig,               &
                                min_liq_overlap, min_liq_overlap_local,        &
                                i_cld_vn, l_bm_sigma_s_grad, l_bm_tweaks,      &
                                turb_var_fac_bm, max_sigmas,                   &
                                min_sigx_ft, min_sigx_fac
use pc2_constants_mod,    only: bm_tiny,i_cld_pc2
use qsat_mod,             only: qsat_wat, qsat_wat_mix, qsat_mix, qsat

implicit none

! Purpose: Calculates liquid cloud water amounts and cloud amounts,
!          temperature and specific humidity from cloud-conserved and
!          other model variables. This is done for one model level.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_cloud
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: UMDP No.39

!  Global Variables:----------------------------------------------------

!  Subroutine Arguments:------------------------------------------------
!  Note that arrays with underscore "_l" are input arrays typically with values
!  for each of the three modes (two modes for the bottom and top of the
!  entrainment zone and one mode for the level to calculate cloud for). Arrays
!  with underscore "_f" are input/output arrays on the level to calculate cloud
!  for.
integer ::                                                                     &
                      !, intent(in)
   points,                                                                     &
!       No. of gridpoints with non-zero cloud
   indx(tdims%j_len*tdims%i_len,2)
!       index for  points with non-zero cloud from lowest model level.

real(kind=real_umphys) ::                                                      &
                      !, intent(in)
   p_l(           pdims%i_start:pdims%i_end,                                   &
                  pdims%j_start:pdims%j_end),                                  &
!       pressure (Pa) for level of interest k
   qsl_l(         tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,3),                                &
!       liquid saturated humidity at temperature TL, and pressure P_F for three
!       modes (bottom of EZ, k-level of interest and top of EZ)
   qsi_l(         tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,3),                                &
!       ice saturated humidity at temperature TL, and pressure P_F for three
!       modes (bottom of EZ, k-level of interest and top of EZ)
   q_l(           tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,3),                                &
!       On input : Vapour + liquid water content (QW) (kg per kg air) for three
!       modes (bottom of EZ, k-level of interest and top of EZ)
!                   (kg water per kg air).
   t_l(           tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,3),                                &
!       On input : Liquid water temperature (TL) (K) for three modes
!       (bottom of EZ, k-level of interest and top of EZ)
   inv_thm_l(     tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,3),                                &
!       Inverse of the turbulence homogenisation time scale for three modes
!       (bottom of EZ, k-level of interest and top of EZ)
   inv_tmp_l(     tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,3),                                &
!       Inverse of the phase-relaxation time scale for three modes
!       (bottom of EZ, k-level of interest and top of EZ)
   tdc_l(         tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,3),                                &
!       Turbulence-decorrelation time scale for three modes
!       (bottom of EZ, k-level of interest and top of EZ)
   wvar_l(        tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,3),                                &
!       Vertical velocity variance for three modes
!       (bottom of EZ, k-level of interest and top of EZ)
   dtldz_l(       tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,3),                                &
   dqtdz_l(       tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,3),                                &
!       Local vertical gradients dTl/dz and dqt/dz
   ql_mean(       tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end)
!       Mean ql in column below, used for tapering minimum-allowed variance
!       within the surface mixed-layer.

logical ::                                                                     &
                      !, intent(in)
   l_mixing_ratio

real(kind=real_umphys) ::                                                      &
                      !, intent(in)
   qcf_f(       tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end)
!       Cloud ice water content at processed levels (kg per kg air).

real(kind=real_umphys) ::                                                      &
                      !, intent(INOUT)
   q_f(         tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end),                                    &
!       On input : Vapour + liquid water content (QW) (kg per kg air).
!       On output: Specific humidity at processed levels
!                   (kg water per kg air).
   t_f(           tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end)
!       On input : Liquid water temperature (TL) (K).
!       On output: Temperature at processed levels (K).
!       Variance of turbulence-based uni-modal s-distribution on levels
real(kind=real_umphys) ::                                                      &
                      !, intent(INOUT)
   cfl_max(       tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end)
!       Max liquid cloud fraction i column, working down

real(kind=real_umphys) ::                                                      &
                      !, intent(out)
   qcl_f(       tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end),                                    &
!       Cloud liquid water content at processed levels (kg per kg air).
   cfl_f(         tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       Liquid cloud fraction at processed levels.
   cff_f(         tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       Ice cloud fraction at processed levels.
   cf_f(          tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       Total cloud fraction at processed levels.
   sskew(         tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       Skewness of mixture of two s-distributions on levels
   svar_turb(     tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       Variance of turbulence-based uni-modal s-distribution on levels
   svar_bm(       tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end)
!       Variance of mixture of two s-distributions in the bi-modal scheme


!  Local parameters and other physical constants------------------------
real(kind=real_umphys) :: alphl,alphi      ! For liquid AlphaL calculation.
real(kind=real_umphys) :: wtn              ! Weighting for ALPHAL iteration
parameter (wtn=0.75)
real(kind=real_umphys) :: mph_overlap      ! Minimal overlap between liquid    &
                                           ! and ice phase in mixed-phase
parameter (mph_overlap=0.05)

!  Local scalars--------------------------------------------------------

real(kind=real_umphys) ::                                                      &
 al(3),                                                                        &
                       ! LOCAL AL for three modes (see eq 4 of UMDP39).
 alphal(3),                                                                    &
                       ! LOCAL ALPHAL for three modes (see eq 4 of UMDP39).
 alphai(3),                                                                    &
                       ! LOCAL ALPHAI for three modes
 sq2rp,sq2,sq6,                                                                &
                       ! precalculated constants
 sigma2,phi,qsl,qsi,mixph,cflmix,qclmix,muedq,bsf,erflay, tmp
                       ! temporary scalars

real(kind=real_umphys) ::                                                      &
 phiqcf,                                                                       &
                      ! Arc-cosine term in Cloud ice fraction calc.
 cosqcf,                                                                       &
                      ! Cosine term in Cloud ice fraction calc.
 qcfrbs
                      ! qcf/bs

real(kind=real_umphys) ::                                                      &
  wgt(3),                                                                      &
                      ! weight of the three modes (index 1 and 3 used in case of
                      ! the bimodal scheme, index 2 used in case of unimodal
                      ! distribution)
  cfl(3),                                                                      &
                      ! liquid cloud fraction for each mode (unweighted)
  cff(3),                                                                      &
                      ! ice cloud fraction for each mode (unweighted)
  cf(3),                                                                       &
                      ! total cloud fraction for each mode (unweighted)
  qcl(3),                                                                      &
                      ! liquid water content for each mode (unweighted)
  sigk(3),                                                                     &
                      ! SD variance for each mode w.r.t the grid-box mean SD
  phik(3),                                                                     &
                      ! third moment of SD for each mode w.r.t. grid-box mean SD
  sgllay(3),                                                                   &
                      ! standard deviation of liquid SD distribution for each
                      ! mode
  sgilay(3),                                                                   &
                      ! standard deviation of ice SD distribution for each
                      ! mode
  muelay(3),                                                                   &
                      ! first moment of mixed-phase liquid SD distribution
                      ! for three modes
  mullay(3),                                                                   &
                      ! first moment of liquid SD distribution for three
                      ! modes
  dqsat(3)
                      ! Difference between ice and liquid saturation specific
                      ! humidity for three modes


! local arrays for bimodal cloud scheme
! ----------------------------------------------------------------------

!  (b) Others.
integer ::   i,ii,ij,n,kk,tlog,its,idn,iup,ikk
             ! Loop counters:I,II-horizontal field index.
!                                       : N - iteration counter.

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='BM_CLD'

!  Local dynamic arrays-------------------------------------------------
!    7 blocks of real workspace are required.
real(kind=real_umphys) ::                                                      &
 p(points,3),                                                                  &
!       Pressure  (Pa).
   t(points,3),                                                                &
!       temperature.
   alphal_nm1(points,3),                                                       &
!       ALPHAL at previous iteration.
   alphai_nm1(points,3)
!       ALPHAI at previous iteration.

! Factors for modifying the Gaussian distribution so that it goes to zero
! at a finite number of standard deviations
real(kind=real_umphys) :: erf_fac
real(kind=real_umphys) :: exp_fac
! Safety limits on qcl and mixed-phase fraction
real(kind=real_umphys) :: max_qcl(3)
real(kind=real_umphys) :: min_mph

!- End of Header


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ==Main Block==--------------------------------------------------------
! Operate on INDEXed points with non-zero cloud fraction.
! ----------------------------------------------------------------------

! Set pre-caclulated constants

alphl=repsilon*lc/r
alphi=repsilon*(lc+lf)/r

sq2rp  = (2.0/pi)**(0.5)
sq2    = (2.0)**0.5
sq6    = (6.0)**0.5

if ( l_bm_tweaks ) then
  ! Precompute factors used to modify the Gaussian distribution so that
  ! it is bounded with a finite width of N standard deviations.
  ! We truncate the ends off the Gaussian e^(-1/2 x^2 )
  ! and renormalise so that it still integrates to unity
  ! pdf(x) = A e^(-1/2 x^2 )
  ! cdf(x) = 1/2 ( 1 + A erf( x/sqrt(2) ) )
  ! Require  cdf(-N) = 0  cdf( N) = 1
  ! => 0 = 1/2 ( 1 + A erf( -N/sqrt(2) ) )  =>  A = -1 / erf( -N/sqrt(2) )
  !                                               =  1 / erf(  N/sqrt(2) )
  ! Then in the qcl calculation, we have:
  ! int_x1^N x pdf(x) dx
  ! = [ -A e^(-1/2 x^2 ) ]_x1^N
  ! = A ( e^(-1/2 x1^2 ) - e^(-1/2 N^2 ) )
  !
  ! Store the terms A and e^(-1/2 N^2 ):
  erf_fac = 1.0 / umerf( max_sigmas/sq2 )
  exp_fac = exp( -0.5 * max_sigmas * max_sigmas )
else
  ! Bimodal tweaks off; no truncation of the Gaussian:
  erf_fac = 1.0
  exp_fac = 0.0
end if

its = 5

do i = 1, points
  ii = indx(i,1)
  ij = indx(i,2)

  ! ----------------------------------------------------------------------
  ! Start iteration loop to approach alpha (5 iterations)
  ! ----------------------------------------------------------------------

  do n = 1, its

    ! ----------------------------------------------------------------------
    ! 1. Calculate ALPHAL (eq P292.5 - UMDP29) and AL (P292.6 - UMDP29).
    !    CAUTION: T_F acts as TL (input value) until update in final section
    !    CAUTION: Q_F acts as QW (input value) until update in final section
    ! ----------------------------------------------------------------------

    tlog = 0


    if (n < 2) then

      ! ----------------------------------------------------------------------
      !     In first iteration, calculate alphal and alphai based for each layer
      !     if this level is part of an entrainment zone, use the three input
      !     modes (mode 1 and 3 to calculate bimodal cloud and mode 2 as a
      !     constraint). Otherwise, perform uni-modal calculations only for
      !     mode 2 (k-level of interest) and ignore mode 1 and 3 in the rest of
      !     this routine
      ! ----------------------------------------------------------------------

      if (minval(t_l(ii,ij,:)) < bm_tiny ) then
        idn = 2
        iup = 2
      else
        idn = 1
        iup = 3
      end if

      do kk = idn,iup

        if (t_l(ii,ij,kk) > bm_tiny) then
          alphal(kk) = alphl*qsl_l(ii,ij,kk) / (t_l(ii,ij,kk)*t_l(ii,ij,kk))
        else
          alphal(kk) = 1.0
        end if
        alphal_nm1(i,kk) = alphal(kk)

        if (t_l(ii,ij,kk) > bm_tiny .and. qcf_f(ii,ij) > 0.0) then
          alphai(kk) = alphi*qsi_l(ii,ij,kk) / (t_l(ii,ij,kk)*t_l(ii,ij,kk))
        else
          alphai(kk) = 1.0
        end if
        alphai_nm1(i,kk) = alphai(kk)

      end do

    else

      ! ----------------------------------------------------------------------
      !   In further iterations, only continue to approximate alphal and alphai
      !   if one of the modes hasn't converged (in which case tlog is set to 1)
      ! ----------------------------------------------------------------------

      do kk=idn,iup
        if (t(i,kk) >  t_l(ii,ij,kk)) then
          if ( l_mixing_ratio ) then
            call qsat_wat_mix(qsl,t(i,kk),p(i,kk))
          else
            call qsat_wat(qsl,t(i,kk),p(i,kk))
          end if
          alphal(kk)       = (qsl-qsl_l(ii,ij,kk))/(t(i,kk)-t_l(ii,ij,kk))
          alphal(kk)       = wtn*alphal(kk)+(1.0-wtn)*alphal_nm1(i,kk)
          alphal_nm1(i,kk) = alphal(kk)
          tlog             = 1
          if (qcf_f(ii,ij) > 0.0) then
            if ( l_mixing_ratio ) then
              call qsat_mix(qsi,t(i,kk),p(i,kk))
            else
              call qsat(qsi,t(i,kk),p(i,kk))
            end if
            alphai(kk)       = (qsi-qsi_l(ii,ij,kk))/(t(i,kk)-t_l(ii,ij,kk))
            alphai(kk)       = wtn*alphai(kk)+(1.0-wtn)*alphai_nm1(i,kk)
            alphai_nm1(i,kk) = alphai(kk)
          else
            alphai(kk)       = alphai_nm1(i,kk)
          end if

        else
          alphal(kk) = alphal_nm1(i,kk)
          alphai(kk) = alphai_nm1(i,kk)
        end if
      end do

    end if

    if ((n < 2) .or. (tlog > 0)) then

      ! ----------------------------------------------------------------------
      !   In first loop, or if alpha hasn't converged, start to calculate
      !   pdf moments and cloud fraction/water content
      ! ----------------------------------------------------------------------

      do kk=idn,iup
        ! Latent heating correction term:
        al(kk)     = 1.0 / (1.0 + (lcrcp * alphal(kk)))
        ! first moment of pure liquid SD distribution for each mode:
        mullay(kk) = al(kk)*(q_l(ii,ij,kk)-qsl_l(ii,ij,kk))
        ! first moment of mixed-phase liquid SD distribution for each mode,
        ! scaled with phase-relaxation and homogenisation time scales as in
        ! Field et al. (2014):
        muelay(kk) = al(kk)*(q_l(ii,ij,kk)-qsi_l(ii,ij,kk))*                   &
                     (inv_thm_l(ii,ij,kk)/                                     &
                     (inv_thm_l(ii,ij,kk)+inv_tmp_l(ii,ij,kk)))
        if ( l_bm_sigma_s_grad ) then
          ! Account for local gradients in calculation of sigma_s
          sgllay(kk) = al(kk) * abs( alphal(kk)*( g/cp + dtldz_l(ii,ij,kk) )   &
                                   - dqtdz_l(ii,ij,kk) ) * turb_var_fac_bm     &
                     * sqrt( 0.5*wvar_l(ii,ij,kk) * tdc_l(ii,ij,kk)            &
                           / inv_thm_l(ii,ij,kk) )
          sgilay(kk) = al(kk) * abs( alphai(kk)*( g/cp + dtldz_l(ii,ij,kk) )   &
                                   - dqtdz_l(ii,ij,kk) ) * turb_var_fac_bm     &
                     * sqrt( 0.5*wvar_l(ii,ij,kk) * tdc_l(ii,ij,kk)            &
                           / ( inv_tmp_l(ii,ij,kk) + inv_thm_l(ii,ij,kk) ) )
        else
          ! Don't account for local gradients
          ! standard deviation of liquid SD distribution for each mode (ignoring
          ! phase-relaxation time scale):
          sgllay(kk) = al(kk)*alphal(kk)*g/cp * turb_var_fac_bm                &
                       *sqrt(0.5*wvar_l(ii,ij,kk)*tdc_l(ii,ij,kk)              &
                             /inv_thm_l(ii,ij,kk))
          ! standard deviation of ice SD distribution for each mode (including
          ! phase-relaxation time scale):
          sgilay(kk) = al(kk)*alphai(kk)*g/cp * turb_var_fac_bm                &
                       *sqrt(0.5*wvar_l(ii,ij,kk)*tdc_l(ii,ij,kk)              &
                             /(inv_tmp_l(ii,ij,kk)+inv_thm_l(ii,ij,kk)))
        end if
        ! departure between liquid and ice saturation specific humidity for each
        ! mode
        dqsat(kk)  = al(kk)*(qsl_l(ii,ij,kk)-qsi_l(ii,ij,kk))

        if ( l_bm_tweaks ) then
          ! Impose minimum variance, consistent with the limits imposed in
          ! bm_ctl and pc2_bm_init.
          tmp = min_sigx_fac * max( ql_mean(ii,ij) - q_l(ii,ij,kk), 0.0 )
          sgllay(kk) = max( sgllay(kk), max(                                   &
            0.01 * al(kk)*qsl_l(ii,ij,kk) / max_sigmas,                        &
            al(kk)*qsl_l(ii,ij,kk) * min_sigx_ft*tmp/(tmp+qsl_l(ii,ij,kk)) ) )
          sgilay(kk) = max( sgilay(kk), max(                                   &
            0.01 * al(kk)*qsi_l(ii,ij,kk) / max_sigmas,                        &
            al(kk)*qsi_l(ii,ij,kk) * min_sigx_ft*tmp/(tmp+qsi_l(ii,ij,kk)) ) )
        end if
      end do

      ! -------------------------------------------------------------------
      ! ------------   Calculate weights and bimodal moments   ------------
      ! -------------------------------------------------------------------

      ! calculate the weights of the EZ-top mode and the EZ-bottom mode
      ! so that the mean saturation departure (mulk) of the level of
      ! interest is conserved and initialise arrays

      ikk = 2

      cff(:) = 0.0
      cfl(:) = 0.0
      cf(:)  = 0.0
      qcl(:) = 0.0
      sigk(:)= 0.0
      phik(:)= 0.0

      ! If this point belongs to an EZ, calculate the weights of the bottom and
      ! top mode and keep weight of central mode as zero. Else, if this point
      ! does not belong to an EZ, use unimodal scheme, by keeping the bottom and
      ! top mode weights as zero and setting the central model weight to 1.
      if (idn < iup) then
        wgt(iup) = (mullay(ikk) - mullay(idn)) / (mullay(iup) - mullay(idn))
        wgt(ikk) = 0.0
        wgt(idn) = 1.0 - wgt(iup)
      else
        wgt(1)   = 0.0
        wgt(ikk) = 1.0
        wgt(3)   = 0.0
      end if

      ! If insensible weights are returned by the previous block, force the
      ! unimodal scheme to be used by resetting the top and bottom mode weights
      ! to zero and the central mode to 1. Reset the loop indices idn and iup to
      ! ikk so no needless loop iterations are done in the remainder of the
      ! code.
      if (wgt(1) <= 0.0 .or. wgt(3) <= 0.0 .or.                                &
          wgt(1) >= 1.0 .or. wgt(3) >= 1.0) then
        idn      = ikk
        iup      = ikk
        wgt(1)   = 0.0
        wgt(ikk) = 1.0
        wgt(3)   = 0.0
      end if

      ! Set max limit on qcl of 95% of the total-water, for safety
      do kk=idn,iup
        max_qcl(kk) = q_l(ii,ij,kk) * 0.95
      end do
      ! Tweak: don't impose a negative max limit when q_l is negative,
      ! as this creates spurious negative qcl
      if ( l_bm_tweaks ) then
        do kk=idn,iup
          max_qcl(kk) = max( max_qcl(kk), 0.0 )
        end do
      end if

      do kk=idn,iup

        ! Loop over all modes to calculate cloud for. If only a single mode is
        ! present, because this grid box does not belong to an entrainment zone,
        ! or no sensible solution is possible, idn = iup = 2 and only the index
        ! 2 values will be filled below. If this grid box belongs to an
        ! entrainment zone, only indices 1 and 3 will be filled with the cloud
        ! properties for the bottom and top mode of the entrainment zone. In
        ! this case, index 2 will be skipped.

        if ((kk == 2 .and. (idn == iup .or. mullay(iup) >= mullay(ikk) .or.    &
             mullay(idn) <= mullay(ikk))) .or.                                 &
            (kk /= 2 .and. idn < iup .and. mullay(iup) < mullay(ikk) .and.     &
             mullay(idn) > mullay(ikk))) then

          if (kk == 2) then
            wgt(1)   = 0.0
            wgt(ikk) = 1.0
            wgt(3)   = 0.0
          end if

          ! Calculate liquid fraction and water content based on instantaneous
          ! condensation using turbulence variances, ignoring any ice for now

          if (mullay(kk) < -max_sigmas*sgllay(kk)) then
            qcl(kk)   = 0.0
          else if ((mullay(kk) >= -max_sigmas*sgllay(kk)) .and.                &
                   (mullay(kk) <   max_sigmas*sgllay(kk))) then
            erflay = umErf( mullay(kk)/(sgllay(kk)*sq2) )
            cfl(kk)   = 0.5*( 1.0 + erf_fac*erflay )
            qcl(kk)   = min( 0.5*( mullay(kk)*( 1.0 + erf_fac*erflay )         &
                                 + erf_fac*sq2rp*sgllay(kk)                    &
                                   *( exp( -mullay(kk)*mullay(kk)              &
                                           /(2.0*sgllay(kk)*sgllay(kk)) )      &
                                    - exp_fac ) ),                             &
                             max_qcl(kk) )
          else
            cfl(kk)   = 1.0
            qcl(kk)   = min( mullay(kk), max_qcl(kk) )
          end if

          if (cfl(kk) > cfl_max(ii,ij)) then
            ! Cloud frac has increased over peak, so use this
            ! (for ice_fraction_method=min_liq_overlap)
            cfl_max(ii,ij) = cfl(kk)
          end if

          if (qcf_f(ii,ij) > 0.0) then

            ! If ice is present, first calculate the ice cloud fraction, using
            ! the diagnostic relation between liquid water content and liquid
            ! cloud fraction, but applied to ice cloud. Use triangular relation
            ! as in the Smith scheme (UMDP29, section 3.4), but use first and
            ! second moment consistent with the turbulence-based diagnostics
            ! following Field et al. 2014


            if (i_cld_vn /= i_cld_pc2) then

              bsf = sq6*sgllay(kk)

              if (ice_fraction_method  ==  smith_orig) then
                ! method based on pure liquid relation
                qcfrbs = qcf_f(ii,ij) / bsf
              else if (ice_fraction_method == min_liq_overlap) then
                ! method described in appendix of Abel et al (2017, JAS)
                qcfrbs = max(1.0 - cfl_max(ii,ij), 0.05) * qcf_f(ii,ij) / bsf
              else if (ice_fraction_method == min_liq_overlap_local) then
                ! erroneous implementation of Abel et al, only using local cfl
                qcfrbs = max(1.0 - cfl(kk), 0.05) * qcf_f(ii,ij) / bsf
              else
                ! no ice cloud fraction method defined
                qcfrbs = 1.0
              end if ! ice_fraction_method

              if (0.0 < qcfrbs .and. (6.0*qcfrbs) <= 1.0) then
                cff(kk) = 0.5 * ((6.0 * qcfrbs)**(2.0/3.0))
              else if (1.0 < (6.0*qcfrbs) .and. qcfrbs < 1.0) then
                phiqcf = acos(sq2 * 0.75 * (1.0 - qcfrbs))
                cosqcf = cos((phiqcf + (4.0 * pi)) / 3.0)
                cff(kk) = 1.0 - (4.0 * cosqcf * cosqcf)
              else if (qcfrbs  >=  1.0) then
                cff(kk) = 1.0
              end if

            else

              cff(kk) = cff_f(ii,ij)

            end if


            ! Calculate the mixed-phase liquid cloud fraction and water
            ! content following Field et al. (2014) using Gaussian
            ! distribution with turbulence based variance, accounting for
            ! ice-phase competition for water vapour.

            muedq =  muelay(kk) - dqsat(kk)

            if (muedq < -max_sigmas*sgilay(kk)) then
              cflmix   = 0.0
              qclmix   = 0.0
            else if ((muedq >= -max_sigmas*sgilay(kk)) .and.                   &
                     (muedq <   max_sigmas*sgilay(kk))) then
              erflay = umErf(muedq/(sgilay(kk)*sq2))
              cflmix   = 0.5 * ( 1.0 + erf_fac*erflay )
              qclmix   = 0.5*( muedq*( 1.0 + erf_fac*erflay )                  &
                             + erf_fac*sq2rp*sgilay(kk)                        &
                               *( exp( -muedq*muedq                            &
                                        /(2.0*sgilay(kk)*sgilay(kk)) )         &
                                - exp_fac ) )
            else
              cflmix   = 1.0
              qclmix   = muedq
            end if


            ! Use a minimum overlap assumption (consistent with large scale
            ! precipitation) to calculate overall portion of liquid-only,
            ! mixed-phase liquid, ice and total cloud cover. If ice is
            ! present, assume a minimum 5 % overlap. Use the pre-calculated
            ! values of liquid-only (cfl) and mixed-phase liquid (cflmix)
            ! proportionally to the fraction of the grid box they are assumed
            ! to cover.

            if (cfl(kk) > 0.0) then
              ! Set min possible mixed-phase cloud fraction (minimum overlap)
              ! = cfl minus the non-ice-cloud area:
              ! Tweak; add brackets to ensure mixph <= cfl
              ! (original code can give mixph > cfl due to rounding errors when
              ! cff=1 and cfl is small, which leads to negative qcl below).
              if ( l_bm_tweaks ) then
                min_mph = cfl(kk) - (1.0 - cff(kk))
              else
                min_mph = cfl(kk) + cff(kk) - 1.0
              end if
              mixph = min(max(min(mph_overlap*cfl(kk),mph_overlap*cff(kk)),    &
                              min_mph),1.0)
              qcl(kk) = qcl(kk)*(cfl(kk)-mixph)/cfl(kk)+qclmix*mixph/cfl(kk)
              cfl(kk) = min(max(0.0,cfl(kk)-mixph+(cflmix/cfl(kk))*mixph),     &
                           1.0)
              ! Update mixed-phase fraction to reflect updated, and reduced cfl
              if ( l_bm_tweaks ) then
                min_mph = cfl(kk) - (1.0 - cff(kk))
              else
                min_mph = cfl(kk) + cff(kk) - 1.0
              end if
              mixph = min(max(min(mph_overlap*cfl(kk),mph_overlap*cff(kk)),    &
                              min_mph),1.0)
              cf(kk)  = min(max(0.0,cfl(kk)+cff(kk)-mixph),1.0)
            else
              qcl(kk) = 0.0
              cfl(kk) = 0.0
              cf(kk)  = cff(kk)
            end if

          else

            ! No ice is present

            cf(kk)  = cfl(kk)

          end if

          ! Calculate distribution moments w.r.t the grid-box mean SD

          sigk(kk) = (mullay(kk)-mullay(ikk))**2.0+sgllay(kk)**2.0
          phik(kk) = 3.0*(mullay(kk)-mullay(ikk))*sgllay(kk)**2.0+             &
                     (mullay(kk)-mullay(ikk))**3.0

        else

          ! No cloud properties to be calculated for this mode index

          sigk(kk) = 0.0
          phik(kk) = 0.0

        end if

      end do

      ! Calculate the second moment of the mixture of PDFs
      sigma2 = sum(wgt(idn:iup)*sigk(idn:iup))
      ! Calculate the third moment of the mixture of PDFs
      phi = sum(wgt(idn:iup)*phik(idn:iup))

      ! ----------------------------------------------------------------------
      ! Combine all layers with their respective weights to calculate liquid,
      ! ice, and total cloud cover and liquid water content.
      ! ----------------------------------------------------------------------

      cfl_f(ii,ij) = min(sum(wgt(idn:iup)*cfl(idn:iup)),1.0)

      ! Never convert more than 95 % of the available humidity to cloud to
      ! avoid negative water vapour
      qcl_f(ii,ij) = min(sum(wgt(idn:iup)*qcl(idn:iup)),                       &
                         max_qcl(ikk))

      if (i_cld_vn /= i_cld_pc2) then
        cff_f(ii,ij) = min(sum(wgt(idn:iup)*cff(idn:iup)),1.0)
        cf_f(ii,ij)  = min(sum(wgt(idn:iup)*cf(idn:iup)),1.0)
      end if

      ! Calculate moment diagnostics of the mixture of pdfs
      sskew(ii,ij)     = phi/((sigma2)**(1.5))
      if (qcf_f(ii,ij) > 0.0) then
        svar_turb(ii,ij) = sgilay(ikk)
      else
        svar_turb(ii,ij) = sgllay(ikk)
      end if
      svar_bm(ii,ij)   = sigma2**0.5

      ! ----------------------------------------------------------------------
      ! Calculate 1st approx. specific humidity (total minus cloud water)
      ! Calculate 1st approx. to temperature, adjusting for latent heating
      ! ----------------------------------------------------------------------

      do kk=idn,iup
        t(i,kk) = t_l(ii,ij,kk) + lcrcp*qcl(kk)
        p(i,kk) = p_l(ii,ij)
      end do

      t_f(ii,ij) = t_l(ii,ij,ikk) + lcrcp*qcl_f(ii,ij)
      q_f(ii,ij) = q_l(ii,ij,ikk) - qcl_f(ii,ij)

    end if ! T_if

  end do ! Its_do
end do ! Points_do1

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine bm_cld
! ======================================================================
end module bm_cld_mod
