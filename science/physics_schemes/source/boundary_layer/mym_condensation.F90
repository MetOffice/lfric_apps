! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Purpose: To evaluate buoyancy parameters, cloud fraction, and
!           liquid water content in the MY model.
!           Fluctuation of heat and moisture is assumed to obey the
!           bi-normal distribution function.
!           Note that the cloud fraction and liquid water content
!           diagnosed here is not used in the other processes.

!  Programming standard : UMDP 3

!  Documentation: UMDP 025

!  This code is based on the code provided by the authors who wrote
!  the following papers.
!   * Nakanishi, M. and H. Niino, 2009: Development of an improved
!      turbulence closure model for the atmospheric boundary layer.
!      J. Meteor. Soc. Japan, 87, 895-912.
!   * Nakanishi, M. and H. Niino, 2006: An improved Mellor-Yamada
!      Level-3 model: Its numerical stability and application to
!      a regional prediction of advection fog.
!      Boundary-Layer Meteor., 119, 397-407.
!   * Nakanishi, M. and H. Niino, 2004: An improved Mellor-Yamada
!      Level-3 model with condensation physics: Its design and
!      verification.
!      Boundary-Layer Meteor., 112, 1-31.
!   * Nakanishi, M., 2001: Improvement of the Mellor-Yamada
!      turbulence closure model based on large-eddy simulation data.
!      Boundary-Layer Meteor., 99, 349-378.
!   The web site publicising their code:
!    http://www.nda.ac.jp/~naka/MYNN/index.html

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module mym_condensation_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'MYM_CONDENSATION_MOD'
contains

subroutine mym_condensation(                                                   &
! IN levels/switches
      bl_levels, levflag, BL_diag,                                             &
! IN fields
      qw, tl, t, p_theta_levels, tsq, qsq, cov,                                &
! OUT fields
      vt, vq, q1, cld, ql)

use atm_fields_bounds_mod, only: tdims
use bl_diags_mod, only: strnewbldiag
use conversions_mod, only: pi
use gen_phys_inputs_mod, only: l_mr_physics
use mym_option_mod, only: tke_levels
use planet_constants_mod, only:                                                &
    cp, r, repsilon, pref, kappa, c_virtual, one_minus_epsilon, ls
use water_constants_mod, only: lc

use model_domain_mod, only: model_type, mt_single_column

use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook

use mym_errfunc_mod, only: mym_errfunc
implicit none

! Intent In Variables
integer, intent(in) ::                                                         &
   bl_levels,                                                                  &
                  ! Max. no. of "boundary" levels
   levflag
                  ! flag to indicate the level of MY
                  ! 2: MY2.5, 3:MY3

real(kind=r_bl), intent(in) ::                                          &
   qw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),          &
                  ! Total water content
   tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),          &
                  ! Ice/liquid water temperature
   t(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),           &
                  ! temperature
   p_theta_levels(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
                  0:bl_levels+1),                                              &
                  ! pressure on theta levels (Pa)
   tsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
          bl_levels),                                                          &
                  ! Self covariant of liquid potential temperature
                  ! (thetal'**2) defined on theta levels K-1
   qsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
          bl_levels),                                                          &
                  ! Self covariant of total water
                  ! (qw'**2) defined on theta levels K-1
   cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
          bl_levels)
                  ! Correlation between thetal and qw
                  ! (thetal'qw') defined on theta levels K-1

!  Declaration of BL diagnostics.
type (strnewbldiag), intent(in out) :: BL_diag

! Intent OUT Variables
real(kind=r_bl), intent(out) ::                                         &
   vt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                  ! Buoyancy parameter (coefficients of <w'thetal'>)
                  ! on theta K-1
                  ! Note that g/thetav is not included.
   vq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                  ! Buoyancy parameter (coefficients of <w'qw'>)
                  ! on theta K-1
                  ! Note that g/thetav is not included.
   q1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                  ! normalized excess water content on theta K-1
   cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       tke_levels),                                                            &
                  ! cloud fraction on theta K-1
   ql(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, tke_levels)
                  ! condensed liquid water content

! Local Variables
integer ::                                                                     &
   i, j, k
                  ! loop indexes
real(kind=r_bl) ::                                                      &
   rr2,                                                                        &
                  ! 1 / sqrt(2)
   rrp,                                                                        &
                  ! 1 / sqrt(2 pi)
   hl,                                                                         &
                  ! latent heat
   qsl,                                                                        &
                  ! saturated specific humidity
   dqsl,                                                                       &
                  ! derivative of saturated specific humidity by
                  ! temperature
   t3sq,                                                                       &
                  ! work variable for <thetal'**2>
   r3sq,                                                                       &
                  ! work variable for <qw'**2>
   c3sq,                                                                       &
                  ! work variable for <thetal'qw'>
   alp_qsl,                                                                    &
                  ! alpha * qsl
   eq1,                                                                        &
                  ! work variable
   qll,                                                                        &
                  ! work variable
   r_exner,                                                                    &
                  ! reciprocal of Exner function
   q2p,                                                                        &
                  ! work variable
   pt_tmp,                                                                     &
                  ! work variable
   qt,                                                                         &
                  ! work variable
   rac
                  ! work variable

real(kind=r_bl) ::                                                      &
   rice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        tke_levels),                                                           &
                  ! ratio of ice.
                  ! temperature > 0C   : rice =0
                  ! temperature < -36C : rice =1
                  ! Between 0C and 36C, it is linearly interpolated
                  ! with temperature
   hl_ovr_cp(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,              &
             tke_levels),                                                      &
                  ! latent heat over heat capacity
   exner(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         tke_levels),                                                          &
                  ! Exner function
   qmq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       tke_levels),                                                            &
                  ! Excess of total water from saturated one
   alp(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       tke_levels),                                                            &
                  ! coefficient related to saturation
                  ! (alpha in the paper)
   bet(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       tke_levels),                                                            &
                  ! coefficient related to saturation
                  ! (beta in the paper)
   sgm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       tke_levels),                                                            &
                  ! standard deviation of the bi-normal distribution
                  ! function
   erf_arg(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           tke_levels),                                                        &
                  ! array to store an argument for error function
   erf_val(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           tke_levels),                                                        &
                  ! array to store a result of error function
   qsw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       tke_levels),                                                            &
                  ! saturated specific ratio for water
   qsi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels)
                  ! saturated specific ratio for ice

real(kind=r_bl), parameter ::                                           &
   e0cw = 6.11e2,                                                              &
   tetn1w = 17.27,                                                             &
   tetn2w = 273.15,                                                            &
   tetn3w = 35.85,                                                             &
   e0ci = 6.11e2,                                                              &
   tetn1i = 21.875,                                                            &
   tetn2i = 273.15,                                                            &
   tetn3i = 7.65,                                                              &
                  ! coefficients in the Tetens' Formula
   ttriple = 273.16,                                                           &
                  ! a triple point of water
   temp_ice = 237.15
                  ! Below this temperature, all of condensed water
                  ! should be ice. -36C
real(kind=r_bl), parameter ::                                           &
   my_sgm_min_fct = 0.0,                                                       &
                  ! factor to set the lower limit for sgm
   my_sgm_max_fct = 1.0
                  ! factor to set the upper limit for sgm

character(len=*), parameter ::  RoutineName = 'MYM_CONDENSATION'

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

rr2 = 1.0 / sqrt(2.0)
rrp = 1.0 / sqrt(2.0 * pi)

! Here, qsw and qsi are saturated vapor pressure.
! Using the Teten's formula instead of the subroutine "qmix"
! because the saturated vapor pressure on liquid water is necessary
! even in sub-zero temperature.
do k = 2, tke_levels
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end

      qsw(i, j, k) = e0cw * exp(tetn1w *                                       &
            (tl(i, j, k - 1)  - tetn2w)                                        &
                  / (tl(i, j, k - 1) - tetn3w) )
      qsi(i, j, k) = e0ci * exp(tetn1i *                                       &
            (tl(i, j, k - 1) - tetn2i)                                         &
           / (tl(i, j, k - 1) - tetn3i) )
    end do
  end do
end do

! convert to mixing ratio or specific humidity
if (l_mr_physics) then
  do k = 2, tke_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        qsw(i, j, k) = repsilon * qsw(i, j, k)                                 &
                                  / p_theta_levels(i, j, k - 1)
        qsi(i, j, k) = repsilon * qsi(i, j, k)                                 &
                                  / p_theta_levels(i, j, k - 1)
      end do
    end do
  end do
else
  do k = 2, tke_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        qsw(i, j, k) = repsilon * qsw(i, j, k)                                 &
                        / (p_theta_levels(i, j, k - 1)                         &
                             - one_minus_epsilon * qsw(i, j, k))
        qsi(i, j, k) = repsilon * qsi(i, j, k)                                 &
                        / (p_theta_levels(i, j, k - 1)                         &
                             - one_minus_epsilon * qsi(i, j, k))
      end do
    end do
  end do
end if

! Calculate sgm
do k = 2, tke_levels
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      if (tl(i, j, k -1) >= ttriple) then
        rice(i, j, k) = 0.0
      else if (tl(i, j, k - 1) < temp_ice) then
        rice(i, j, k) = 1.0
      else
        rice(i, j, k) = (ttriple - tl(i, j, k - 1))                            &
                                       / (ttriple - temp_ice)
      end if

      hl = (1.0 - rice(i, j, k)) * lc + rice(i, j, k) * ls
      qsl = (1.0 - rice(i, j, k)) * qsw(i, j, k)                               &
              + rice(i, j, k) * qsi(i, j, k)

      hl_ovr_cp(i, j, k) = hl / cp

      dqsl = qsl * repsilon * hl / (r * tl(i, j, k - 1) **2)

      exner(i, j, k) = (p_theta_levels(i, j, k - 1) / pref) ** kappa
      qmq(i, j, k) = qw(i, j, k - 1) - qsl
      alp(i, j, k) = 1.0 /(1.0 + dqsl * hl_ovr_cp(i, j, k))
      bet(i, j, k) = dqsl * exner(i, j, k)

      t3sq = max(tsq(i, j, k), 0.0)
      r3sq = max(qsq(i, j, k), 0.0)
      c3sq = cov(i, j, k)
      c3sq = sign(min(abs(c3sq), sqrt(t3sq * r3sq)), c3sq)

      r3sq = r3sq + bet(i, j, k) ** 2 * t3sq                                   &
           -2.0 * bet(i, j, k) * c3sq
      alp_qsl = min(alp(i, j, k) * qsl, qw(i, j, k - 1))
      sgm(i, j, k) = max(                                                      &
                      min(0.5 * alp(i, j, k) * sqrt(max(r3sq, 0.0)),           &
                      my_sgm_max_fct * alp_qsl),                               &
                      my_sgm_min_fct * alp_qsl, 1.0e-10)
    end do
  end do
end do

if (levflag /= 3) then
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      sgm(i, j, 2) = sgm(i, j, 3)
    end do
  end do
end if

do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end
    erf_arg(i, j, 1) = 0.0
    cld(i, j, 1) = 0.0
    ql(i, j, 1) = 0.0
    sgm(i, j, 1) = 0.0
    q1(i, j, 1) = 0.0
  end do
end do

! Preparation to calculate values of the err function
do k = 2, tke_levels
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      q1(i, j, k) = 0.5 * alp(i, j, k)                                         &
                             * qmq(i, j, k) / sgm(i, j, k)
      erf_arg(i, j, k) = q1(i, j, k) * rr2
    end do
  end do
end do

call mym_errfunc(tdims%i_end*tdims%j_end*tke_levels, erf_arg, erf_val)

! Calculate the buoyancy parameters vt and vq
do k = 2, tke_levels
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      cld(i, j, k) = 0.5 * (1.0 + erf_val(i, j, k))
      if (abs(q1(i, j, k)) > 10.0 ) then
        eq1 = 0.0
      else
        eq1  = rrp * exp(- 0.5 * q1(i, j, k) ** 2)
      end if
      ! qll = ql / (2 * sgm)
      qll  = max(cld(i, j, k) * q1(i, j, k) + eq1, 0.0)

      if (qw(i, j, k) < 1.0e-10) then
        ql(i, j, k) = 0.0
      else
        ql(i, j, k) = max(                                                     &
             2.0 * sgm(i, j, k) * qll, 0.0)
      end if
      ! To avoid negative QV (for safety)
      ql(i, j, k) = min(ql(i, j, k), qw(i, j, k - 1) * 0.5)

      r_exner = 1.0 / exner(i, j, k)
      q2p  = hl_ovr_cp(i, j, k) * r_exner
      pt_tmp = t(i, j, k - 1) * r_exner
      qt = 1.0 + c_virtual * qw(i, j, k - 1)                                   &
               - (1.0 + c_virtual) * ql(i, j, k)
      rac = alp(i, j, k) * (cld(i, j, k)                                       &
           - qll * eq1) * (q2p * qt - (1.0 + c_virtual) * pt_tmp)

      vt (i, j, k) = qt - rac * bet(i, j, k)
      vq (i, j, k) = c_virtual * pt_tmp + rac
    end do
  end do
end do

if (BL_diag%l_cf_trb) then
  do k = 2, tke_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        BL_diag%cf_trb(i, j, k) = cld(i, j, k)
      end do
    end do
  end do
end if

if (BL_diag%l_ql_trb) then
  do k = 2, tke_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        BL_diag%ql_trb(i, j, k) = ql(i, j, k)
      end do
    end do
  end do
end if

if (BL_diag%l_sgm_trb) then
  do k = 2, tke_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        BL_diag%sgm_trb(i, j, k) = sgm(i, j, k)
      end do
    end do
  end do
end if


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine mym_condensation
end module mym_condensation_mod
