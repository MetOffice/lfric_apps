! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Cloud Scheme: Homogeneous forcing and Turbulence (non-updating)

module pc2_delta_hom_turb_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'PC2_DELTA_HOM_TURB_MOD'
contains

subroutine pc2_delta_hom_turb(                                                 &
!      Pressure related fields
 p_theta_levels,                                                               &
!      Prognostic Fields
 t, q, qcl, cf, cfl, cff,                                                      &
!      Forcing quantities for driving the homogeneous forcing
 dtin, dqin,                                                                   &
!      Output increments to the prognostic fields
 dtpc2, dqpc2, dqclpc2, dcfpc2, dcflpc2,                                       &
!      Other quantities for the turbulence
 dbsdtbs0, dbsdtbs1, l_mixing_ratio)

use water_constants_mod, only: lc
use planet_constants_mod,  only: lcrcp, r, repsilon
use timestep_mod, only: timestep
use yomhook,               only: lhook, dr_hook
use parkind1,              only: jprb, jpim
use atm_fields_bounds_mod, only: pdims, tdims
use pc2_constants_mod,     only: dbsdtbs_exp, pdf_power,                       &
                                 pdf_merge_power, cloud_rounding_tol,          &
                                 i_pc2_homog_g_cf, i_pc2_homog_g_width
use cloud_inputs_mod,      only: l_fixbug_pc2_qcl_incr,l_fixbug_pc2_mixph,     &
                                 i_pc2_homog_g_method

use qsat_mod, only: qsat_wat, qsat_wat_mix

implicit none

! Description:
!   This subroutine calculates the change in liquid content, liquid
!   cloud fraction and total cloud fraction as a result of homogeneous
!   forcing of the gridbox with temperature, pressue, vapour and liquid
!   increments and turbulent mixing effects.

! Method:
!   Uses the method in Gregory et al (2002, QJRMS 128 1485-1504) and
!   Wilson and Gregory (2003, QJRMS 129 967-986)
!   which considers a probability density distribution whose
!   properties are only influenced by a change of width due to
!   turbulence.
!   There is a check to ensure we cannot remove more condensate than was
!   there to start with.
!   There is a check to ensure we remove all liquid if all
!   fraction is removed.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_cloud
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: Annexe to UMDP29 The PC2 CLoud Scheme.

!  Subroutine Arguments:------------------------------------------------

! arguments with intent in. ie: input variables.

real(kind=real_umphys), intent(in) ::                                          &
 dbsdtbs0,                                                                     &
!       Value of dbs/dt / bs which is independent of forcing (no units)
   dbsdtbs1
!       Value of dbs/dt / bs which is proportional to the forcing
!       ( (kg kg-1 s-1)-1 )

logical, intent(in) ::                                                         &
 l_mixing_ratio        ! Use mixing ratio formulation

real(kind=real_umphys), intent(in) ::                                          &
 p_theta_levels(pdims%i_start:pdims%i_end,                                     &
                pdims%j_start:pdims%j_end,                                     &
                            1:pdims%k_end),                                    &
!       Pressure at all points (Pa)
   t(             tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                              1:tdims%k_end),                                  &
!       Temperature (K)
   q(             tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                              1:tdims%k_end),                                  &
!       Vapour content (kg water per kg air)
   qcl(           tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                              1:tdims%k_end),                                  &
!       Liquid content (kg water per kg air)
   cf(            tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                              1:tdims%k_end),                                  &
!       Total cloud fraction (no units)
   cfl(           tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                              1:tdims%k_end),                                  &
!       Liquid cloud fraction (no units)
   cff(           tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                              1:tdims%k_end),                                  &
!       Ice cloud fraction (no units)
   dtin(          tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                              1:tdims%k_end),                                  &
!       Increment of temperature from forcing mechanism (K)
   dqin(          tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                              1:tdims%k_end)
!       Increment of vapour from forcing mechanism (kg kg-1)

! arguments with intent out. ie: output variables.

real(kind=real_umphys), intent(out) ::                                         &
 dtpc2(         tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,                                     &
                            1:tdims%k_end),                                    &
!       PC2 Increment to Temperature (K)
   dqpc2(         tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                              1:tdims%k_end),                                  &
!       PC2 Increment to Vapour content (kg water per kg air)
   dqclpc2(       tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                              1:tdims%k_end),                                  &
!       PC2 Increment to Liquid content (kg water per kg air)
   dcfpc2(        tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                              1:tdims%k_end),                                  &
!       PC2 Increment to Total cloud fraction (no units)
   dcflpc2(       tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                              1:tdims%k_end)
!       PC2 Increment to Liquid cloud fraction (no units)

!  Local parameters and other physical constants------------------------

real(kind=real_umphys), parameter :: b_factor=(pdf_power+1.0)/(pdf_power+2.0)
!       Premultiplier to calculate the amplitude of the probability
!       density function at the saturation boundary (G_MQC).
real(kind=real_umphys), parameter :: smallp=1.0e-10
!       Small positive value for use in if tests

!  Local scalars--------------------------------------------------------
real(kind=real_umphys) ::                                                      &
 alpha,                                                                        &
!       Rate of change of saturation specific humidity with
!       temperature calculated at dry-bulb temperature (kg kg-1 K-1)
   al,                                                                         &
!       1 / (1 + alpha L/cp)  (no units)
   c_1,                                                                        &
!       Mid-timestep liquid cloud fraction (no units)
   dbsdtbs,                                                                    &
!       Relative rate of change of distribution width (s-1)
   dqcdt,                                                                      &
!       Forcing of QC (kg kg-1 s-1)
   deltal,                                                                     &
!       Change in liquid content (kg kg-1)
   cfl_to_m,                                                                   &
!       CFL(i,j,k)**PDF_MERGE_POWER
   sky_to_m,                                                                   &
!       (1-CFL(i,j,k))**PDF_MERGE_POWER
   g_mqc,                                                                      &
!       Amplitude of the probability density function at
!       the saturation boundary (kg kg-1)-1
   qc,                                                                         &
!       aL (q + l - qsat(TL) )  (kg kg-1)
   sd
!       Saturation deficit (= aL (q - qsat(T)) )  (kg kg-1)

!  (b)  Others.
integer :: k,i,j ! Loop counters:  K - vertical level index
!                                  I,J - horizontal position index

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real   (kind=jprb)            :: zhook_handle

character(len=*), parameter :: RoutineName='PC2_DELTA_HOM_TURB'

!  Local arrays---------------------------------------------------------
real(kind=real_umphys) ::                                                      &
  qsl_t,                                                                       &
!       Saturated specific humidity for dry bulb temperature T
  qsl_tl,                                                                      &
!       Saturated specific humidity for liquid temperature TL
  tl
!- End of Header

! ==Main Block==--------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP  PARALLEL do DEFAULT(none) SCHEDULE(DYNAMIC) private(k,                  &
!$OMP  j, i, tl, qsl_t, qsl_tl, alpha, al,                                     &
!$OMP  sd, cfl_to_m, sky_to_m, g_mqc, dqcdt, qc, dbsdtbs,                      &
!$OMP  c_1, deltal)                                                            &
!$OMP  SHARED(tdims,cfl,t,lcrcp,qcl,p_theta_levels,l_mixing_ratio,             &
!$OMP     repsilon,r,q,dqin,dtin,dbsdtbs0,dbsdtbs1,                            &
!$OMP     timestep,dcflpc2,                                                    &
!$OMP     dcfpc2,cf,cff,dqclpc2,dqpc2,dtpc2,                                   &
!$OMP     i_pc2_homog_g_method)

! Loop round levels to be processed
do k = 1, tdims%k_end
  do j = tdims%j_start, tdims%j_end

    do i = tdims%i_start, tdims%i_end

      ! There is no need to perform the total cloud fraction calculation in
      ! this subroutine if there is no, or full, liquid cloud cover.

      if (cfl(i,j,k) >        cloud_rounding_tol .and.                         &
          cfl(i,j,k) < (1.0 - cloud_rounding_tol)) then

        ! ----------------------------------------------------------------------
        ! 3. Calculate the parameters relating to the probability density func.
        ! ----------------------------------------------------------------------

        ! NOTE: The following calculations are gratuitously duplicated
        ! in the following 3 routines:
        ! - pc2_homog_plus_turb
        ! - pc2_delta_hom_turb
        ! - pc2_hom_conv
        ! If you alter them in one of these routines, please make the
        ! alterations consistently in all 3.

        ! Calculate Saturated Specific Humidity with respect to liquid water
        ! for dry bulb temperature.
        if ( l_mixing_ratio ) then
          call qsat_wat_mix(qsl_t, t(i,j,k), p_theta_levels(i,j,k))
        else
          call qsat_wat(qsl_t, t(i,j,k), p_theta_levels(i,j,k))
        end if

        ! Need to estimate the rate of change of saturated specific humidity
        ! with respect to temperature (alpha) first, then use this to calculate
        ! factor aL. Also estimate the rate of change of qsat with pressure.
        alpha = repsilon*lc*qsl_t /                                            &
                (r*t(i,j,k)**2)
        al = 1.0 / (1.0 + lcrcp*alpha)

        ! Calculate the saturation deficit SD

        sd = al*(qsl_t-q(i,j,k))

        ! Calculate the amplitude of the probability density function at the
        ! saturation boundary...
        if ( i_pc2_homog_g_method == i_pc2_homog_g_cf ) then
          ! Blend two solutions as a function of cloud-fraction

          if (qcl(i,j,k) > smallp .and. sd > smallp) then

            cfl_to_m = cfl(i,j,k)**pdf_merge_power
            sky_to_m = (1.0 - cfl(i,j,k))**pdf_merge_power

            g_mqc = b_factor * ( (1.0-cfl(i,j,k))**2 *                         &
             cfl_to_m / (sd * (cfl_to_m + sky_to_m))                           &
                  +                   cfl(i,j,k)**2 *                          &
             sky_to_m / (qcl(i,j,k) * (cfl_to_m + sky_to_m)) )

          else
            g_mqc = 0.0
          end if

        else if ( i_pc2_homog_g_method == i_pc2_homog_g_width ) then
          ! Blend two solutions as a function of their PDF-widths

          ! Solution (a) is ga = b_factor cfl^2/qcl
          ! Solution (b) is gb = b_factor (1-cf)^2/sd
          ! Weight (a) is   wa = qcl/cfl
          ! Weight (b) is   wb = sd/(1-cfl)
          ! => g_mqc = ( wa ga + wb gb ) / ( wa + wb )
          !          = b_factor ( cfl + (1-cf) ) / ( qcl/cfl + sd/(1-cfl) )
          !          = b_factor / ( qcl/cfl + sd/(1-cfl) )
          g_mqc = b_factor / ( max( qcl(i,j,k) / cfl(i,j,k),   smallp )        &
                             + max( sd / ( 1.0 - cfl(i,j,k) ), smallp ) )

        end if  ! ( i_pc2_homog_g_method )

        ! Calculate the rate of change of Qc due to the forcing

        dqcdt = al * ( dqin(i,j,k) - alpha*dtin(i,j,k) )

        ! Calculate Saturated Specific Humidity with respect to liquid water
        ! for wet bulb temperature.

        tl = t(i,j,k)-lcrcp*qcl(i,j,k)
        if ( l_mixing_ratio ) then
          call qsat_wat_mix(qsl_tl, tl, p_theta_levels(i,j,k))
        else
          call qsat_wat(qsl_tl, tl, p_theta_levels(i,j,k))
        end if

        ! Calculate Qc

        qc = al * (q(i,j,k) + qcl(i,j,k) - qsl_tl)

        ! Calculate the relative rate of change of width of the distribution
        ! dbsdtbs from the forcing rate
        if (dbsdtbs0 /= 0.0 .or. dbsdtbs1 /= 0.0) then
          dbsdtbs = (dbsdtbs0 * timestep + dqcdt * dbsdtbs1) *                 &
             exp(-dbsdtbs_exp * qc / (al * qsl_tl))
        else
          dbsdtbs = 0.0
        end if
        ! ----------------------------------------------------------------------
        ! 4. Calculate the change of liquid cloud fraction. This uses the
        ! arrival value of QC for better behaved numerics.
        ! ----------------------------------------------------------------------

        ! DQCDT is the homogeneous forcing part, (QC+DQCDT)*DBSDTBS is the
        ! width narrowing part

        dcflpc2(i,j,k) = g_mqc * ( dqcdt - (qc + dqcdt)*dbsdtbs)

        ! Calculate the condensation amount DELTAL. This uses a mid value
        ! of cloud fraction for better numerical behaviour.

        c_1 = max( 0.0, min( (cfl(i,j,k) + dcflpc2(i,j,k)), 1.0) )

        dcflpc2(i,j,k) = c_1 - cfl(i,j,k)

        c_1    = 0.5 * (c_1 + cfl(i,j,k))
        deltal = c_1 * dqcdt + (qcl(i,j,k) - qc * c_1) * dbsdtbs
        !
        ! If we have removed all fraction, remove all liquid
        if (l_fixbug_pc2_qcl_incr) then
          if (cfl(i,j,k)+dcflpc2(i,j,k) == 0.0) then
            deltal = -qcl(i,j,k)
          end if
        end if
        !

        ! ----------------------------------------------------------------------
        ! 5. Calculate change in total cloud fraction.
        ! ----------------------------------------------------------------------

        ! ----------------------------------------------------------------------
        ! The following If test is a copy of the PC2_TOTAL_CF subroutine.
        ! ----------------------------------------------------------------------
        if (dcflpc2(i,j,k)  >   0.0) then
          ! ...  .and. CFL(i,j,k)  <   1.0 already assured.
          if (l_fixbug_pc2_mixph) then
            ! minimum overlap, consistent with pc2_totalcf
            dcfpc2(i,j,k) = min(dcflpc2(i,j,k),(1.0-cf(i,j,k)))
          else
            ! random overlap, this is inconsistent with pc2_totalcf
            dcfpc2(i,j,k) = dcflpc2(i,j,k) * (1.0 - cf(i,j,k)) /               &
                                             (1.0 - cfl(i,j,k))
          end if
        else if (dcflpc2(i,j,k)  <   0.0) then
          ! ...  .and. CFL(i,j,k)  >   0.0 already assured.
          if (l_fixbug_pc2_mixph) then
            ! minimum overlap, consistent with pc2_totalcf
            dcfpc2(i,j,k) = max(dcflpc2(i,j,k),(cff(i,j,k)-cf(i,j,k)))
          else
            ! random overlap, this is inconsistent with pc2_totalcf
            dcfpc2(i,j,k) = dcflpc2(i,j,k) *                                   &
                (cf(i,j,k)- cff(i,j,k)) / cfl(i,j,k)
          end if
        else
          dcfpc2(i,j,k) = 0.0
        end if
        ! ----------------------------------------------------------------------

      else if ( ( cfl(i,j,k)  >=  (1.0 - cloud_rounding_tol) .and.             &
                  cfl(i,j,k)  <=  (1.0 + cloud_rounding_tol) ) .or.            &
        ! this if test is wrong, it should be cfl >= 1
        ! add fix on a switch to preserve bit-comparison
                          ( cfl(i,j,k)  >=  (1.0 - cloud_rounding_tol) .and.   &
                          l_fixbug_pc2_mixph ) ) then

        ! Cloud fraction is 1

        dcfpc2(i,j,k)  = 0.0
        dcflpc2(i,j,k) = 0.0

        ! Calculate Saturated Specific Humidity with respect to liquid water
        ! for dry bulb temperature.
        if ( l_mixing_ratio ) then
          call qsat_wat_mix(qsl_t, t(i,j,k), p_theta_levels(i,j,k))
        else
          call qsat_wat(qsl_t, t(i,j,k), p_theta_levels(i,j,k))
        end if

        alpha = repsilon * lc * qsl_t /                                        &
                (r * t(i,j,k)**2)
        al = 1.0 / (1.0 + lcrcp*alpha)
        deltal = al * (dqin(i,j,k) - alpha*dtin(i,j,k))

      else

        ! Cloud fraction is 0

        dcfpc2(i,j,k)  = 0.0
        dcflpc2(i,j,k) = 0.0

        deltal = 0.0

      end if

      ! Increment water contents and temperature due to latent heating
      ! This subroutine will output only the condensation increments
      ! hence we comment out updates to qcl, q and t
      !           QCL(i,j,k) = QCL(i,j,k) + DELTAL
      ! Q = input Q + Forcing - Condensation
      !           Q(i,j,k)   = Q(i,j,k) + DQIN(i,j,k) - (DELTAL - DLIN(i,j,k))
      !           T(i,j,k)   = T(i,j,k) + DTIN(i,j,k)
      !    &                            + LCRCP * (DELTAL - DLIN(i,j,k))

      ! These are the condensation increments
      dqclpc2(i,j,k) = deltal

      if (l_fixbug_pc2_qcl_incr) then
        ! Don't allow more water to be removed than was there to start with.
        if (qcl(i,j,k) + dqclpc2(i,j,k) < 0.0) then
          dqclpc2(i,j,k) = -qcl(i,j,k)
        end if
      end if

      dqpc2(i,j,k)   = - dqclpc2(i,j,k)
      dtpc2(i,j,k)   = lcrcp * dqclpc2(i,j,k)

    end do !i
  end do !j

end do !k
!$OMP end PARALLEL do

! End of the subroutine

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine pc2_delta_hom_turb
end module pc2_delta_hom_turb_mod
