! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  PC2 Cloud Scheme: Homogeneous forcing and Turbulence

module pc2_homog_plus_turb_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'PC2_HOMOG_PLUS_TURB_MOD'
contains

subroutine pc2_homog_plus_turb(                                                &
!   Pressure related fields
 p_theta_levels,                                                               &
!   Array dimensions
 nlevels,                                                                      &
!   Timestep
 timestep,                                                                     &
!   Prognostic Fields
 t, cf, cfl, cff, q, qcl,                                                      &
!   Forcing quantities for driving the homogeneous forcing
 dtdt, dqdt, dldt, dpdt,                                                       &
!   Other quantities for the turbulence
 dbsdtbs0, dbsdtbs1,                                                           &
!   Model switches
 l_mixing_ratio)

use water_constants_mod,   only: lc
use planet_constants_mod,  only: lcrcp, r, repsilon
use yomhook,               only: lhook, dr_hook
use parkind1,              only: jprb, jpim
use atm_fields_bounds_mod, only: pdims, tdims
use pc2_constants_mod,     only: dbsdtbs_exp, pdf_power,                       &
                                 pdf_merge_power, cloud_rounding_tol,          &
                                 i_pc2_homog_g_cf, i_pc2_homog_g_width
use cloud_inputs_mod,      only: l_fixbug_pc2_qcl_incr,l_fixbug_pc2_mixph,     &
                                 i_pc2_homog_g_method
use science_fixes_mod,     only: l_pc2_homog_turb_q_neg
use qsat_mod,              only: qsat_wat, qsat_wat_mix
use pc2_total_cf_mod,      only: pc2_total_cf

use free_tracers_inputs_mod, only: l_wtrac
use wtrac_pc2_mod,           only: wtrac_pc2

implicit none

! Purpose:
!   This subroutine calculates the change in liquid content, liquid
!   cloud fraction and total cloud fraction as a result of homogeneous
!   forcing of the gridbox with temperature, pressure, vapour and liquid
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
! Documentation: PC2 Cloud Scheme Documentation

!  Subroutine Arguments:------------------------------------------------
integer ::                                                                     &
                      !, intent(in)
 nlevels
!    No. of levels being processed.

real(kind=real_umphys) ::                                                      &
                      !, intent(in)
 timestep,                                                                     &
!    Model timestep (s)
   p_theta_levels(pdims%i_start:pdims%i_end,                                   &
                  pdims%j_start:pdims%j_end,                                   &
                  nlevels),                                                    &
!    pressure at all points (Pa)
   cff(           tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                  nlevels),                                                    &
!    Ice cloud fraction (no units)
   dtdt(          tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                  nlevels),                                                    &
!    Increment of temperature from forcing mechanism (K)
   dqdt(          tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                  nlevels),                                                    &
!    Increment of vapour from forcing mechanism (kg kg-1)
   dldt(          tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                  nlevels),                                                    &
!    Increment of liquid from forcing mechanism (kg kg-1)
   dpdt(          pdims%i_start:pdims%i_end,                                   &
                  pdims%j_start:pdims%j_end,                                   &
                  nlevels),                                                    &
!    Increment in pressure from forcing mechanism (Pa)
   dbsdtbs0,                                                                   &
!    Value of dbs/dt / bs which is independent of forcing (no units)
   dbsdtbs1
!    Value of dbs/dt / bs which is proportional to the forcing
!    ( (kg kg-1 s-1)-1 )

logical ::                                                                     &
                       !, intent(in)
 l_mixing_ratio        ! Use mixing ratio formulation

real(kind=real_umphys) ::                                                      &
                      !, intent(INOUT)
 t(  tdims%i_start:tdims%i_end,                                                &
     tdims%j_start:tdims%j_end,                                                &
     nlevels),                                                                 &
!    Temperature (K)
   cf( tdims%i_start:tdims%i_end,                                              &
       tdims%j_start:tdims%j_end,                                              &
       nlevels),                                                               &
!    Total cloud fraction (no units)
   cfl(tdims%i_start:tdims%i_end,                                              &
       tdims%j_start:tdims%j_end,                                              &
       nlevels),                                                               &
!    Liquid cloud fraction (no units)
   q(  tdims%i_start:tdims%i_end,                                              &
       tdims%j_start:tdims%j_end,                                              &
       nlevels),                                                               &
!    Vapour content (kg water per kg air)
   qcl(tdims%i_start:tdims%i_end,                                              &
       tdims%j_start:tdims%j_end,                                              &
       nlevels)
!    Liquid content (kg water per kg air)

!    External functions:

!    Local parameters and other physical constants-----------------------
real(kind=real_umphys), parameter :: b_factor=(pdf_power+1.0)/(pdf_power+2.0)
!    Premultiplier to calculate the amplitude of the probability
!    density function at the saturation boundary (G).
real(kind=real_umphys), parameter :: smallp=1.0e-10
!    Small positive value for use in if tests

!  Local scalars--------------------------------------------------------
real(kind=real_umphys) ::                                                      &
 alpha,                                                                        &
                ! Rate of change of saturation specific humidity with
                ! temperature calculated at dry-bulb temperature
                ! (kg kg-1 K-1)
 alpha_p,                                                                      &
                ! Rate of change of saturation specific humidity with
                ! pressure calculated at dry-bulb temperature (Pa K-1)
 al,                                                                           &
                ! 1 / (1 + alpha L/cp)  (no units)
 c_1,                                                                          &
                ! Mid-timestep liquid cloud fraction (no units)
 dbsdtbs,                                                                      &
                ! Relative rate of change of distribution width (s-1)
 dqcdt,                                                                        &
                ! Forcing of QC (kg kg-1 s-1)
 deltal,                                                                       &
                ! Change in liquid content (kg kg-1)
 g_mqc,                                                                        &
                ! Amplitude of the probability density function at
                ! the saturation boundary (kg kg-1)-1
 cfl_to_m,                                                                     &
                ! CFL(i,j,k)**PDF_MERGE_POWER
 sky_to_m,                                                                     &
                ! (1-CFL(i,j,k))**PDF_MERGE_POWER
 qc,                                                                           &
                ! aL (q + l - qsat(TL) )  (kg kg-1)
 sd             ! Saturation deficit (= aL (q - qsat(T)) )  (kg kg-1)

!  (b) Others.
integer :: k,i,j,                                                              &
                ! Loop counters: K - vertical level index
                ! I,J - horizontal position index
           npt  ! Number of point on which to perform calculations

!  Local arrays---------------------------------------------------------
real(kind=real_umphys) ::                                                      &
 qsl_t,                                                                        &
!    Saturated specific humidity for dry bulb temperature T
   qsl_tl,                                                                     &
!    Saturated specific humidity for liquid temperature TL
   tl,                                                                         &
!    Liquid temperature (= T - L/cp QCL)  (K)
   cf_c(      tdims%i_len*                                                     &
              tdims%j_len ),                                                   &
!    Total cloud fraction on compressed points
   cfl_c(     tdims%i_len*                                                     &
              tdims%j_len ),                                                   &
!    Liquid cloud fraction on compressed points
   cff_c(     tdims%i_len*                                                     &
              tdims%j_len ),                                                   &
!    Ice cloud fraction on compressed points
   deltacl_c( tdims%i_len*                                                     &
              tdims%j_len ),                                                   &
!    Change in liquid cloud fraction (no units)
   deltacf_c( tdims%i_len*                                                     &
              tdims%j_len )
!    Change in ice cloud fraction (no units)
integer ::                                                                     &
 index_npt(tdims%i_start:tdims%i_end,                                          &
           tdims%j_start:tdims%j_end)


integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real   (kind=jprb)            :: zhook_handle

character(len=*), parameter :: RoutineName='PC2_HOMOG_PLUS_TURB'

!- End of Header

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


! ==Main Block==--------------------------------------------------------


! Loop round levels to be processed

!$OMP  PARALLEL do DEFAULT(SHARED) SCHEDULE(DYNAMIC) private(cfl_c,            &
!$OMP  index_npt, npt, cf_c, cff_c, deltacf_c, qsl_t, tl,                      &
!$OMP  qsl_tl, alpha, al, alpha_p, sd, g_mqc, dqcdt, dbsdtbs, qc, deltal, i,   &
!$OMP  j, k, c_1, deltacl_c, cfl_to_m, sky_to_m)
do k = 1, nlevels

  ! Copy points into compressed arrays
  npt = 0
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      if (cfl(i,j,k) > cloud_rounding_tol ) then
        npt = npt + 1
        index_npt(i,j) = npt
        cfl_c(npt) = cfl(i,j,k)
        cf_c(npt)  = cf(i,j,k)
        cff_c(npt) = cff(i,j,k)
      end if
    end do
  end do

  do i = 1, npt
    ! Initialize deltacl_c to zero
    deltacl_c(i) = 0.0

    ! Initialisation of change in frozen cloud fraction, which is always
    ! zero from this routine.
    deltacf_c(i) = 0.0
  end do

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
        ! dry bulb temperature.
        if ( l_mixing_ratio ) then
          call qsat_wat_mix(qsl_t, t(i,j,k), p_theta_levels(i,j,k))
        else
          call qsat_wat(qsl_t, t(i,j,k), p_theta_levels(i,j,k))
        end if

        ! Need to estimate the rate of change of saturated specific humidity
        ! with respect to temperature (alpha) first, then use this to calculate
        ! factor aL. Also estimate the rate of change of qsat with pressure.
        alpha=repsilon*lc*qsl_t/(r*t(i,j,k)**2)
        al=1.0/(1.0+lcrcp*alpha)
        alpha_p = -qsl_t/p_theta_levels(i,j,k)

        ! Calculate the saturation deficit SD

        sd=al*(qsl_t-q(i,j,k))

        ! Calculate the amplitude of the probability density function at the
        ! saturation boundary...
        if ( i_pc2_homog_g_method == i_pc2_homog_g_cf ) then
          ! Blend two solutions as a function of cloud-fraction

          if (qcl(i,j,k) > smallp .and. sd > smallp) then

            cfl_to_m = cfl(i,j,k)**pdf_merge_power
            sky_to_m = (1.0 - cfl(i,j,k))**pdf_merge_power

            g_mqc=b_factor*(   cfl_to_m                                        &
             *(1.0-cfl(i,j,k))**2/sd                                           &
             +                 sky_to_m                                        &
             *cfl(i,j,k)**2/qcl(i,j,k)  )                                      &
             /(   cfl_to_m + sky_to_m   )

          else
            g_mqc=0.0
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

        dqcdt=al * ( dqdt(i,j,k)-alpha*dtdt(i,j,k)                             &
                    -alpha_p*dpdt(i,j,k) ) + dldt(i,j,k)

        ! Calculate Saturated Specific Humidity with respect to liquid water
        ! wet bulb temperature.
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

        deltacl_c(index_npt(i,j)) = g_mqc * ( dqcdt - (qc + dqcdt) * dbsdtbs)

        ! Calculate the condensation amount DELTAL. This uses a mid value
        ! of cloud fraction for better numerical behaviour.

        c_1 = cfl(i,j,k) + deltacl_c(index_npt(i,j))
        if (c_1  >   1.0) then
          deltacl_c(index_npt(i,j)) = 1.0 - cfl(i,j,k)
          c_1=1.0
        else if (c_1  <   0.0) then
          c_1=0.0
          deltacl_c(index_npt(i,j)) = (- cfl(i,j,k) )
        end if
        c_1 = 0.5 * (c_1 + cfl(i,j,k))
        deltal = c_1 * dqcdt + (qcl(i,j,k) - qc * c_1) * dbsdtbs
        !
        ! If we have removed all fraction, remove all liquid
        if (l_fixbug_pc2_qcl_incr) then
          if (cfl(i,j,k)+deltacl_c(index_npt(i,j)) == 0.0) then
            deltal = -qcl(i,j,k)
          end if
        end if
        !
      else if ( ( cfl(i,j,k)  >=  (1.0 - cloud_rounding_tol) .and.             &
                  cfl(i,j,k)  <=  (1.0 + cloud_rounding_tol) ) .or.            &
        ! this if test is wrong, it should be cfl >= 1
        ! add fix on a switch to preserve bit-comparison
                          ( cfl(i,j,k)  >=  (1.0 - cloud_rounding_tol) .and.   &
                          l_fixbug_pc2_mixph ) ) then

        ! Cloud fraction is 1

        ! Calculate Saturated Specific Humidity with respect to liquid water
        ! dry/wet bulb temperature.
        if ( l_mixing_ratio ) then
          call qsat_wat_mix(qsl_t, t(i,j,k), p_theta_levels(i,j,k))
        else
          call qsat_wat(qsl_t, t(i,j,k), p_theta_levels(i,j,k))
        end if

        alpha=repsilon*lc*qsl_t/(r*t(i,j,k)**2)
        al=1.0/(1.0+lcrcp*alpha)
        alpha_p = -qsl_t/p_theta_levels(i,j,k)
        deltal=al * (dqdt(i,j,k)-alpha*dtdt(i,j,k)                             &
                     -alpha_p*dpdt(i,j,k)) + dldt(i,j,k)

        deltacl_c(index_npt(i,j)) = 0.0

      else

        ! Cloud fraction is 0 (actually CFL < cloud_rounding_tol)

        deltal = 0.0

      end if

      ! Update water contents
      !
      ! Don't allow more liquid water to be created than there is q available,
      ! making sure q doesn't go below small positive value
      ! (this is okay with the limits to deltal vs qcl below, as this is
      ! creating cloud, not destroying it):
      if ( l_pc2_homog_turb_q_neg .and.                                        &
          (q(i,j,k) + dqdt(i,j,k) - (deltal - dldt(i,j,k)) < smallp) ) then
        deltal = q(i,j,k) + dqdt(i,j,k) + dldt(i,j,k) - smallp
      end if

      if (l_fixbug_pc2_qcl_incr) then
        ! Don't allow more cloud liquid to be removed than was there to start
        ! with.
        if (qcl(i,j,k) + deltal < 0.0) then
          deltal = -qcl(i,j,k)
          ! Set qcl to be exactly zero in this case.
          qcl(i,j,k) = 0.0
        else
          qcl(i,j,k) = qcl(i,j,k) + deltal
        end if

      else
        qcl(i,j,k) = qcl(i,j,k) + deltal
      end if

      ! Update vapour content
      ! Q = input Q + Forcing - Condensation
      q(i,j,k)   = q(i,j,k)   + dqdt(i,j,k)                                    &
                  - (deltal - dldt(i,j,k))

      ! Update temperature due to latent heating
      t(i,j,k)   = t(i,j,k)   + dtdt(i,j,k)                                    &
                    + lcrcp * (deltal - dldt(i,j,k))

      if (l_wtrac) wtrac_pc2%q_cond(i,j,k) = deltal
    end do !i

  end do !j

  ! ----------------------------------------------------------------------
  ! 5. Now update cloud fractions.
  ! ----------------------------------------------------------------------

  ! Calculate change in total cloud fraction.

  if (npt > 0) then
    call pc2_total_cf(                                                         &
          npt,cfl_c,cff_c,deltacl_c,deltacf_c,cf_c)
  end if

  do j = tdims%j_start, tdims%j_end

    do i = tdims%i_start, tdims%i_end

      if (cfl(i,j,k) >        cloud_rounding_tol .and.                         &
          cfl(i,j,k) < (1.0 - cloud_rounding_tol)) then
        ! Update cloud fractions
        cf(i,j,k)  = cf_c(index_npt(i,j))
        cfl(i,j,k) = cfl(i,j,k) + deltacl_c(index_npt(i,j))
      end if ! End if for CFL gt 0 and CFL lt 1

    end do !i

  end do !j

end do
!$OMP end PARALLEL do

! End of the subroutine

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine pc2_homog_plus_turb
end module pc2_homog_plus_turb_mod
