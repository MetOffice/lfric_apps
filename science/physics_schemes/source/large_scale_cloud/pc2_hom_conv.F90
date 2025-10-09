! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module pc2_hom_conv_mod

use um_types, only: real_umphys


implicit none

character(len=*), parameter, private :: ModuleName='PC2_HOM_CONV_MOD'

contains

!  Cloud Scheme: Homogeneous forcing and Turbulence (non-updating)
! Subroutine Interface:
subroutine pc2_hom_conv(                                                       &
!      Pressure related fields
 p_theta_levels,                                                               &
!      Timestep
 timestep,                                                                     &
!      Prognostic Fields
 t, q, qcl, cf, cfl, cff,                                                      &
!      Forcing quantities for driving the homogeneous forcing
 dtin, dqin, dqclin, dpdt, dcflin,                                             &
!      Output increments to the prognostic fields
 dtpc2, dqpc2, dqclpc2, dcfpc2, dcflpc2,                                       &
!      Other quantities for the turbulence
 pc2mixingrate, dbsdtbs1 )

use water_constants_mod,   only: lc
use planet_constants_mod,  only: lcrcp, r, repsilon
use yomhook,               only: lhook, dr_hook
use parkind1,              only: jprb, jpim
use atm_fields_bounds_mod, only: pdims, tdims
use cloud_inputs_mod,      only: i_pc2_erosion_method, i_pc2_erosion_numerics, &
     l_fixbug_pc2_qcl_incr,l_fixbug_pc2_mixph, i_pc2_homog_g_method
use pc2_constants_mod,     only: pc2eros_exp_rh,                               &
     pc2eros_hybrid_sidesonly,                                                 &
     i_pc2_erosion_explicit, i_pc2_erosion_implicit, i_pc2_erosion_analytic,   &
     pdf_power, pdf_merge_power, dbsdtbs_exp, cloud_rounding_tol,              &
     i_pc2_homog_g_cf, i_pc2_homog_g_width
use ereport_mod,           only: ereport
use gen_phys_inputs_mod,   only: l_mr_physics
use qsat_mod,              only: qsat_wat, qsat_wat_mix

use errormessagelength_mod, only: errormessagelength

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
! Documentation: PC2 cloud scheme documentation

! Subroutine Arguments:------------------------------------------------

! Arguments with intent in. ie: input variables.

real(kind=real_umphys), intent(in) ::                                          &
 timestep,                                                                     &
!       Model timestep (s)
   pc2mixingrate,                                                              &
!       If i_pc2_erosion_method==pc2eros_exp_rh
!       Value of dbs/dt / bs which is independent of forcing (no units)
!       If i_pc2_erosion_method==pc2eros_hybrid_sidesonly:
!       Erosion rate to use in hybrid erosion method.
   dbsdtbs1
!       Used if i_pc2_erosion_method==pc2eros_exp_rh
!       Value of dbs/dt / bs which is proportional to the forcing
!       ( (kg kg-1 s-1)-1 )

real(kind=real_umphys), intent(in) ::                                          &
 p_theta_levels(pdims%i_start:pdims%i_end,                                     &
                pdims%j_start:pdims%j_end),                                    &
!       Pressure at all points (Pa)
   t(             tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       Temperature (K)
   q(             tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       Vapour content (kg water per kg air)
   qcl(           tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       Liquid content (kg water per kg air)
   cf(            tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       Total cloud fraction (no units)
   cfl(           tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       Liquid cloud fraction (no units)
   cff(           tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       Ice cloud fraction (no units)
   dtin(          tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       Increment of temperature from forcing mechanism (K)
   dqin(          tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       Increment of vapour from forcing mechanism (kg kg-1)
   dqclin(        tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       Increment of liquid from forcing mechanism (kg kg-1)
   dpdt(          pdims%i_start:pdims%i_end,                                   &
                  pdims%j_start:pdims%j_end),                                  &
!       Increment in pressure from forcing mechanism (Pa)
   dcflin(        tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end)
!       Increment in liquid cloud fraction (no units)

! Arguments with intent out. ie: output variables.

real(kind=real_umphys), intent(out) ::                                         &
 dtpc2(         tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end),                                    &
!       PC2 Increment to Temperature (K)
   dqpc2(         tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       PC2 Increment to Vapour content (kg water per kg air)
   dqclpc2(       tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       PC2 Increment to Liquid content (kg water per kg air)
   dcfpc2(        tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       PC2 Increment to Total cloud fraction (no units)
   dcflpc2(       tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end)
!       PC2 Increment to Liquid cloud fraction (no units)

!  External subroutine calls: ------------------------------------------

!  Local parameters and other physical constants------------------------
real(kind=real_umphys), parameter ::                                           &
                   b_factor = (pdf_power+1.0) / (pdf_power+2.0)
!       Premultiplier to calculate the amplitude of the probability
!       density function at the saturation boundary (G_MQC).
real(kind=real_umphys), parameter :: smallp = 1.0e-10
!       Small positive value for use in if tests

!  Local scalars--------------------------------------------------------

real(kind=real_umphys) ::                                                      &
 alpha,                                                                        &
!       Rate of change of saturation specific humidity with
!       temperature calculated at dry-bulb temperature (kg kg-1 K-1)
   alpha_p,                                                                    &
!       Rate of change of saturation specific humidity with
!       pressure calculated at dry-bulb temperature (Pa K-1)
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
!       CFL(i,j)**PDF_MERGE_POWER
   sky_to_m,                                                                   &
!       (1-CFL(i,j))**PDF_MERGE_POWER
   g_mqc,                                                                      &
!       Amplitude of the probability density function at
!       the saturation boundary (kg kg-1)-1
   qc,                                                                         &
!       aL (q + l - qsat(TL) )  (kg kg-1)
   sd,                                                                         &
!       Saturation deficit (= aL (q - qsat(T)) )  (kg kg-1)
   dcs
!       Injected cloud fraction

! Variables for analytical erosion method
real(kind=real_umphys) :: qcl0    ! Liquid water after homog forcing
real(kind=real_umphys) :: sde0    ! Saturation defecit after homog forcing
real(kind=real_umphys) :: cfl0    ! Liquid fraction after homog forcing
real(kind=real_umphys) :: cfc0    ! Clear fraction after homog forcing
real(kind=real_umphys) :: sde1    ! Saturation defecit guess (Qc > 0)
real(kind=real_umphys) :: cfc1    ! Clear fraction guess     (Qc > 0)
real(kind=real_umphys) :: qcl2    ! Liquid water guess       (Qc < 0)
real(kind=real_umphys) :: cfl2    ! Liquid fraction guess    (Qc < 0)
real(kind=real_umphys) :: factor  ! Store for gathered terms
real(kind=real_umphys) :: c_term  ! Constant for PDF end shape
real(kind=real_umphys) :: b_term  ! Power relating qcl and cfl increments

!  (b)  Others.
integer :: i,j,n ! Loop counters: I,J - horizontal position index

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real   (kind=jprb)            :: zhook_handle

character(len=*), parameter :: RoutineName='PC2_HOM_CONV'

real(kind=real_umphys) ::                                                      &
  qsl_t,                                                                       &
!       Saturated specific humidity for dry bulb temperature T
  qsl_tl,                                                                      &
!       Saturated specific humidity for liquid temperature TL
  tl
!    Liquid temperature (= T - L/cp QCL)  (K)

real(kind=real_umphys) :: tmp
real(kind=real_umphys) :: dqcl, dcl, midpoint_qcl
real(kind=real_umphys) :: exposed_area, satdiff

character(len=errormessagelength)       :: message
integer                  :: errorstatus

!- End of Header

! ==Main Block==--------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end

    ! There is no need to perform the total cloud fraction calculation in
    ! this subroutine if there is no, or full, liquid cloud cover.

    if (cfl(i,j) >        cloud_rounding_tol .and.                             &
        cfl(i,j) < (1.0 - cloud_rounding_tol)) then

      ! ---------------------------------------------------------------------
      ! 3. Calculate the parameters relating to the probability density func.
      ! ---------------------------------------------------------------------

      ! NOTE: The following calculations are gratuitously duplicated
      ! in the following 3 routines:
      ! - pc2_homog_plus_turb
      ! - pc2_delta_hom_turb
      ! - pc2_hom_conv
      ! If you alter them in one of these routines, please make the
      ! alterations consistently in all 3.

      if ( l_mr_physics ) then
        call qsat_wat_mix(qsl_t, t(i,j), p_theta_levels(i,j))
      else
        call qsat_wat(qsl_t, t(i,j), p_theta_levels(i,j))
      end if

      ! Need to estimate the rate of change of saturated specific humidity
      ! with respect to temperature (alpha) first, then use this to calculate
      ! factor aL. Also estimate the rate of change of qsat with pressure.
      alpha   = repsilon*lc*qsl_t / (r*t(i,j)**2)
      al      = 1.0 / ( 1.0 + lcrcp * alpha )
      alpha_p = -qsl_t / p_theta_levels(i,j)

      ! Calculate the saturation deficit SD

      sd      = al * ( qsl_t - q(i,j) )

      ! Calculate the amplitude of the probability density function at the
      ! saturation boundary...
      if ( i_pc2_homog_g_method == i_pc2_homog_g_cf ) then
        ! Blend two solutions as a function of cloud-fraction

        if (qcl(i,j) > smallp .and. sd > smallp) then

          cfl_to_m = cfl(i,j)**pdf_merge_power
          sky_to_m = (1.0 - cfl(i,j))**pdf_merge_power

          g_mqc = b_factor * ( (1.0-cfl(i,j))**2 *                             &
           cfl_to_m / (sd * (cfl_to_m + sky_to_m))                             &
                +                   cfl(i,j)**2 *                              &
           sky_to_m / (qcl(i,j) * (cfl_to_m + sky_to_m)) )

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
        g_mqc = b_factor / ( max( qcl(i,j) / cfl(i,j),   smallp )              &
                           + max( sd / ( 1.0 - cfl(i,j) ), smallp ) )

      end if  ! ( i_pc2_homog_g_method )

      ! Calculate the rate of change of Qc due to the forcing

      dqcdt = al * ( dqin(i,j) - alpha*dtin(i,j)                               &
             -alpha_p*dpdt(i,j) ) + dqclin(i,j)

      ! For the background homogeneous forcing from the convection there is
      ! an additional term because the detrained plume must be saturated.
      ! This can also be written as a forcing.

      !              IF (CFL(i,j)  >   0.0 .AND. CFL(i,j)  <   1.0) THEN
      ! This if test is already guaranteed
      if (dcflin(i,j) > 0.0) then
        dcs = dcflin(i,j) / (1.0 - cfl(i,j))
      else if (dcflin(i,j) < 0.0) then
        dcs = dcflin(i,j) / ( - cfl(i,j))
      else
        dcs = 0.0
      end if
      ! Limit DCS to 0 and 1
      dcs   = max( min(dcs,1.0) ,0.0)

      dqcdt = dqcdt - al * dcs * (qsl_t-q(i,j))

      ! Calculate Qc
      tl = t(i,j)-lcrcp*qcl(i,j)
      if ( l_mr_physics ) then
        call qsat_wat_mix(qsl_tl, tl, p_theta_levels(i,j))
      else
        call qsat_wat(qsl_tl, tl, p_theta_levels(i,j))
      end if

      qc    = al * ( q(i,j) + qcl(i,j) - qsl_tl )

      ! Calculate the relative rate of change of width of the distribution
      ! dbsdtbs from the forcing rate

      if (i_pc2_erosion_method == pc2eros_exp_rh) then
        ! Original Wilson et al (2008) formulation: rate of
        ! narrowing related to RH via an ad-hoc exponetial.
        dbsdtbs = (pc2mixingrate * timestep + dqcdt * dbsdtbs1) *              &
               exp(-dbsdtbs_exp * qc / (al * qsl_tl))
      else if (i_pc2_erosion_method == pc2eros_hybrid_sidesonly) then
        ! Hybrid method
        dbsdtbs = 0.0
        ! By setting this to zero, the next bit of code will not
        ! do any width-narrowing to represent erosion.
        ! So need to represent erosion in some other way.
      else
        errorstatus=10
        message='Attempting to use undefined i_pc2_erosion_method'
        call ereport('PC2_hom_conv',errorstatus,message)
      end if

      ! ---------------------------------------------------------------------
      ! 4. Calculate the change of liquid cloud fraction. This uses the
      ! arrival value of QC for better behaved numerics.
      ! ---------------------------------------------------------------------

      ! DQCDT is the homogeneous forcing part, (QC+DQCDT)*DBSDTBS is the
      ! width narrowing part

      dcflpc2(i,j) = g_mqc * ( dqcdt - (qc + dqcdt)*dbsdtbs)

      ! Calculate the condensation amount DELTAL. This uses a mid value
      ! of cloud fraction for better numerical behaviour.

      c_1 = max( 0.0, min( (cfl(i,j) + dcflpc2(i,j)), 1.0) )

      dcflpc2(i,j) = c_1 - cfl(i,j)

      c_1 = 0.5 * (c_1 + cfl(i,j))

      if (l_fixbug_pc2_qcl_incr) then
        ! Calculate increment here, without checking it for
        ! potential removal of more QCL than there is.
        ! A check is carried out at end of routine.
        deltal = (c_1 * dqcdt) +                                               &
               ! The homogeneous forcing part
                 ( ( qcl(i,j) - (qc * c_1) ) * dbsdtbs)
               ! The PDF width-narrowing part.

        ! If we have removed all fraction, remove all liquid
        if (cfl(i,j)+dcflpc2(i,j) == 0.0) then
          deltal = -qcl(i,j)
        end if
      else
        ! Original code. Note that only the width-narrowing
        ! (erosion) bit is limited to not remove more QCL than there
        ! is. The homog forcing part can remove too much QCL.
        deltal = c_1 * dqcdt +                                                 &
            max( (qcl(i,j) - qc * c_1) * dbsdtbs , (-qcl(i,j)) )
      end if

      if (i_pc2_erosion_method == pc2eros_exp_rh ) then
        ! Only calculate change in total cloud fraction here if not using
        ! hybrid method as need to add increments from hybrid erosion.
        ! Change in total cloud fraction will be done later.

        ! -------------------------------------------------------------------
        ! 5. Calculate change in total cloud fraction.
        ! -------------------------------------------------------------------

        ! -------------------------------------------------------------------
        ! The following If test is a copy of the PC2_TOTAL_CF subroutine.
        ! -------------------------------------------------------------------
        if (dcflpc2(i,j) > 0.0) then
          ! ...  .AND. CFL(i,j)  <   1.0 already assured.
          if (l_fixbug_pc2_mixph) then
            ! minimum overlap, this is consistent with pc2_totalcf
            dcfpc2(i,j) = min(dcflpc2(i,j),(1.0-cf(i,j)))
          else
            ! random overlap, this is inconsistent with pc2_totalcf
            dcfpc2(i,j) = dcflpc2(i,j) * (1.0 - cf(i,j)) /                     &
                                             (1.0 - cfl(i,j))
          end if
        else if (dcflpc2(i,j) < 0.0) then
          ! ...  .AND. CFL(i,j)  >   0.0 already assured.
          if (l_fixbug_pc2_mixph) then
            ! minimum overlap, this is consistent with pc2_totalcf
            dcfpc2(i,j) = max(dcflpc2(i,j),(cff(i,j)-cf(i,j)))
          else
            ! random overlap, this is inconsistent with pc2_totalcf
            dcfpc2(i,j) = dcflpc2(i,j) *                                       &
                (cf(i,j)- cff(i,j)) / cfl(i,j)
          end if
        else
          dcfpc2(i,j) = 0.0
        end if
        ! -------------------------------------------------------------------
      end if

    else if ( ( cfl(i,j)  >=  (1.0 - cloud_rounding_tol) .and.                 &
                cfl(i,j)  <=  (1.0 + cloud_rounding_tol) ) .or.                &
      ! this if test is wrong, it should be cfl >= 1
      ! add fix on a switch to preserve bit-comparison
                        ( cfl(i,j)  >=  (1.0 - cloud_rounding_tol) .and.       &
                        l_fixbug_pc2_mixph ) ) then

      ! Cloud fraction is 1

      dcfpc2(i,j)  = 0.0
      dcflpc2(i,j) = 0.0

      if ( l_mr_physics ) then
        call qsat_wat_mix(qsl_t, t(i,j), p_theta_levels(i,j))
      else
        call qsat_wat(qsl_t, t(i,j), p_theta_levels(i,j))
      end if

      alpha   = repsilon * lc * qsl_t / (r * t(i,j)**2)
      al      = 1.0 / (1.0 + lcrcp*alpha)
      alpha_p = -qsl_t / p_theta_levels(i,j)
      deltal  = al * (dqin(i,j) - alpha*dtin(i,j)                              &
            -alpha_p*dpdt(i,j)) + dqclin(i,j)

    else

      ! Cloud fraction is 0

      dcfpc2(i,j)  = 0.0
      dcflpc2(i,j) = 0.0

      deltal         = 0.0

    end if

    !===========================================================
    ! Hybrid PC2 erosion
    !===========================================================
    if (i_pc2_erosion_method == pc2eros_hybrid_sidesonly ) then
      ! Although this alternative method of doing erosion
      ! is being done in a separate, subsequent bit of code to the
      ! homogeneous forcing, it is effectively making its
      ! calculations in parallel since it uses the input values of
      ! cloud fields and thermodynamics.

      if (cfl(i,j) > cloud_rounding_tol .and.                                  &
          cfl(i,j) < (1.0-cloud_rounding_tol)) then
        ! Erosion of cloud is assumed to only happen from the
        ! cloud surface area exposed to clear sky.

        ! Three options are available below for the numerical method used to
        ! time-integrate the erosion.  For full details and derivation of
        ! the equations, see UMDP 030: The PC2 Cloud-Scheme, in the subsection
        ! "Numerical application of the hybrid erosion method".

        if ( i_pc2_erosion_numerics == i_pc2_erosion_explicit ) then
          ! Use simple explicit numerical method

          ! Calculate the difference from saturation
          satdiff = qsl_t - q(i,j)

          ! Calculate exposed lateral surface area. Define a function
          ! which is an upside-down U shape, going to zero at CFL=0
          ! and CFL=1 and with a peak value of 0.5 at CFL=0.5.
          exposed_area = (2.0*cfl(i,j)) - (2.0*cfl(i,j)*cfl(i,j))
          dqcl = - exposed_area * pc2mixingrate * satdiff * timestep

          ! Assume that the change in QCL is due to
          ! width-narrowing. Use that width-narrowing rate to find
          ! the consistent change in cloud fraction...

          ! Find the value of qcl half way through the erosion
          ! process for better numerical behaviour.
          midpoint_qcl = qcl(i,j) + ( 0.5 * dqcl )
          tmp = (midpoint_qcl - (qc * cfl(i,j)) )
          dcl = - g_mqc * qc * dqcl / tmp

        else if ( i_pc2_erosion_numerics == i_pc2_erosion_implicit ) then
          ! Use approximate implicit numerical method, assuming erosion
          ! rate reduces in proportion to cloud as it approaches zero.
          ! This avoids spuriously removing all of the cloud due to
          ! numerical overshoot.

          ! Calculate explicit qcl increment as above
          satdiff = qsl_t - q(i,j)
          exposed_area = (2.0*cfl(i,j)) - (2.0*cfl(i,j)*cfl(i,j))
          dqcl = - exposed_area * pc2mixingrate * satdiff * timestep

          ! Calculate cloud fraction increment purely explicitly
          ! (omit use of mid-point qcl) since we apply an implicit correction
          dcl = - g_mqc * qc * dqcl / ( qcl(i,j) - (qc * cfl(i,j)) )

          ! If erosion is removing cloud
          if ( dcl < 0.0 .and. dqcl < 0.0 ) then
            if ( qcl(i,j)+deltal <= 0.0 .or.                                   &
                 cfl(i,j)+dcflpc2(i,j) <= 0.0 ) then
              ! Set erosion increments to zero if homogeneous
              ! forcing has already removed all the cloud.
              dqcl = 0.0
              dcl = 0.0
            else
              ! Still some cloud left after homogeneous forcing,
              ! and erosion is trying to remove it; use implicit
              ! formula for the erosion terms:
              dqcl = dqcl * ( qcl(i,j) + deltal )                              &
                          / ( qcl(i,j) - dqcl )
              dcl  = dcl  * ( cfl(i,j) + dcflpc2(i,j) )                        &
                          / ( cfl(i,j) - dcl )

            end if  ! Still some cloud left after homogeneous forcing
          end if  ! ( dcl < 0.0 .AND. dqcl < 0.0 )

        else if ( i_pc2_erosion_numerics == i_pc2_erosion_analytic ) then
          ! Use approximate analytical solution which behaves smoothly
          ! and consistently whilst actually allowing erosion to reduce
          ! cloud to zero when appropriate.

          ! Compute latest values of qc, qcl, sde (after homogeneous forcing)
          qc = qc + dqcdt
          qcl0 = qcl(i,j) + deltal
          cfl0 = cfl(i,j) + dcflpc2(i,j)
          sde0 = qcl0     - qc
          cfc0 = 1.0      - cfl0

          if ( cfl0>0.0 .and. cfc0>0.0 .and. qcl0>0.0 .and. sde0>0.0 ) then
            ! Still have positive cloudy and clear fractions,
            ! water content and saturation defecit after homogeneous forcing

            if ( qc < -smallp ) then
              ! Qc is negative; as the PDF is narrowed, we will run out of
              ! liquid-cloud before the saturation defecit goes to zero.
              ! Take analytical solution which accounts for cfl diminishing
              ! as the cloud evaporates...

              ! Initialise guess qcl, cfl after erosion
              qcl2 = qcl0
              cfl2 = cfl0

              ! Precalculate terms
              factor = (pc2mixingrate/al) * 2.0 * cfl0 * timestep / qcl0
              if ( g_mqc > 0.0 ) then
                c_term = g_mqc * qcl0 / cfl0**2
              else
                ! Set to cloudy-end value in case where homog forcing aborted.
                c_term = b_factor
              end if

              do n = 1, 5
                ! Iterate to correct explicit terms

                ! Compute updated power b (taking mid-point in time)
                b_term = c_term * qc * 0.5*( 1.0/(qc - qcl0/cfl0)              &
                                           + 1.0/(qc - qcl2/max(cfl2,smallp)) )
                ! Safety limit; must be below 1
                b_term = min( b_term, 0.999 )

                ! Compute new value of qcl using mid-point cfl, qcl
                qcl2 = qcl0 * max( 1.0 - (1.0-b_term) * factor                 &
                                       * ( 1.0 - 0.5*(cfl0+cfl2) )             &
                                       * ( 0.5*(qcl0+qcl2) - qc ),             &
                                   0.0 )**(1.0/(1.0-b_term))

                ! Compute updated cfl consistent with this
                cfl2 = cfl0 * (qcl2/qcl0)**b_term

              end do

              ! Set erosion increments consistent with the above
              dqcl = qcl2 - qcl0
              dcl  = cfl2 - cfl0

            else if ( qc > smallp ) then
              ! Qc is positive; as the PDF is narrowed, the saturation
              ! defecit will approach zero, but since the erosion rate
              ! is proportional to sd, it will never quite reach zero.
              ! Take analytical solution which accounts for this...

              ! Initialise guess sde, cfc after erosion
              sde1 = sde0
              cfc1 = cfc0

              ! Precalculate terms
              factor = (pc2mixingrate/al) * 2.0 * cfc0 * timestep
              if ( g_mqc > 0.0 ) then
                c_term = g_mqc * sde0 / cfc0**2
              else
                ! Set to clear-end value in case where homog forcing aborted.
                c_term = b_factor
              end if

              do n = 1, 5
                ! Iterate to correct explicit terms

                ! Compute updated power b (taking mid-point in time)
                b_term = c_term * qc * 0.5*( 1.0/(qc + sde0/cfc0)              &
                                           + 1.0/(qc + sde1/max(cfc1,smallp)) )
                ! Safety limit; must be below 1
                b_term = min( b_term, 0.999 )

                ! Compute new value of sde using mid-point cfc
                sde1 = sde0 * ( 1.0 + b_term * factor                          &
                                      * ( 1.0 - 0.5*(cfc0+cfc1) )              &
                              )**(-1.0/b_term)

                ! Compute updated cfc consistent with this
                cfc1 = cfc0 * (sde1/sde0)**b_term

              end do

              ! Set erosion increments consistent with the above
              dqcl = sde1 - sde0
              dcl  = cfc0 - cfc1

            else  ! qc
              ! Qc is extremely close to zero.
              ! In this (presumably rare) case, the PDF is exactly centred
              ! on the saturation point.  This simplifies things, as the
              ! cloud fraction remains constant under width-narrowing.
              ! qcl and the saturation defecit both decline exponentially
              ! towards zero...

              sde1 = sde0 * exp( -(pc2mixingrate/al) * 2.0 * (1.0-cfc0) * cfc0 &
                                 * timestep )

              ! Set erosion increments consistent with the above
              dqcl = sde1 - sde0
              dcl  = 0.0

            end if

          else  ! ( cfl0>0.0 .AND. cfc0>0.0 .AND. qcl0>0.0 .AND. sde0>0.0 )
            ! Something gone to zero; PDF no longer contains a saturation
            ! boundary so no erosion possible; set increments to zero

            dqcl = 0.0
            dcl  = 0.0

          end if  ! ( cfl0>0.0 .AND. cfc0>0.0 .AND. qcl0>0.0 .AND. sde0>0.0 )

        end if  ! ( i_pc2_erosion_numerics = i_pc2_erosion_analytic )

      else if ( ( cfl(i,j)  >=  (1.0 - cloud_rounding_tol) .and.               &
                cfl(i,j)  <=  (1.0 + cloud_rounding_tol) ) .or.                &
        ! this if test is wrong, it should be cfl >= 1
        ! add fix on a switch to preserve bit-comparison
                          ( cfl(i,j)  >=  (1.0 - cloud_rounding_tol) .and.     &
                          l_fixbug_pc2_mixph ) ) then
                    ! No contribution from lateral exposed area,
                    ! just from top and bottom.

        !i_pc2_erosion_method == pc2eros_hybrid_sidesonly
        ! If we assume mixing is only happening from the sides
        ! and not the top and bottom, then when cfl=1 there
        ! will be no sides exposed and so no mixing at all.
        dqcl = 0.0
        dcl  = 0.0

      else
        ! Cloud fraction = 0, so no cloud there to remove.
        ! Set sink terms to zero.
        dcl=0.0
        dqcl=0.0
      end if

      ! Add the the erosion increments to the homogeneous forcing ones
      dcflpc2(i,j) = dcflpc2(i,j) + dcl

      deltal = deltal + dqcl

      ! Test for silly things using original value and total increments

      if (qcl(i,j)+deltal <= 0.0) then
        ! If we are about to removed all QCL remove all CFL.
        dcflpc2(i,j) = -cfl(i,j)
        deltal = -qcl(i,j)
      end if

      if ( cfl(i,j)+dcflpc2(i,j) <= 0.0 ) then
        ! If we are about to removed all CFL remove all QCL.
        dcflpc2(i,j) = -cfl(i,j)
        deltal = -qcl(i,j)
      end if

      if ( cfl(i,j)+dcflpc2(i,j) > 1.0 ) then
        ! If we are about to make CFL>1 then only increase it to 1.
        dcflpc2(i,j) = 1.0-cfl(i,j)
      end if

      ! Calculate change in total cloud fraction.

      ! The following only needs to be done if 0<CFL<1, but is it worth
      ! enforcing that? as dcflpc2 will only be non-zero if thing have
      ! been done to it in 0<cfl<1 bits

      ! ---------------------------------------------------------------------
      ! The following If test is a copy of the PC2_TOTAL_CF subroutine.
      ! ---------------------------------------------------------------------
      if (dcflpc2(i,j) > 0.0) then
        ! ...  .AND. CFL(i,j)  <   1.0 already assured.
        if (l_fixbug_pc2_mixph) then
          ! minimum overlap, this is consistent with pc2_totalcf
          dcfpc2(i,j) = min(dcflpc2(i,j),(1.0-cf(i,j)))
        else
          ! random overlap, this is inconsistent with pc2_totalcf
          dcfpc2(i,j) = dcflpc2(i,j) * (1.0 - cf(i,j)) /                       &
                                           (1.0 - cfl(i,j))
        end if
      else if (dcflpc2(i,j) < 0.0) then
        ! ...  .AND. CFL(i,j)  >   0.0 already assured.
        if (l_fixbug_pc2_mixph) then
          ! minimum overlap, this is consistent with pc2_totalcf
          dcfpc2(i,j) = max(dcflpc2(i,j),(cff(i,j)-cf(i,j)))
        else
          ! random overlap, this is inconsistent with pc2_totalcf
          dcfpc2(i,j) = dcflpc2(i,j) *                                         &
            (cf(i,j)- cff(i,j)) / cfl(i,j)
        end if
      else
        dcfpc2(i,j) = 0.0
      end if
      ! ----------------------------------------------------------------------

    end if ! (i_pc2_erosion_method == hybrid_method)

    !===========================================================

    ! Increment water contents and temperature due to latent heating.
    ! This subroutine will output only the condensation increments
    ! hence we comment out updates to qcl, q and t
    !           QCL(i,j) = QCL(i,j) + DELTAL
    ! Q = input Q + Forcing - Condensation
    !           Q(i,j)   = Q(i,j) + DQIN(i,j) - (DELTAL - DLIN(i,j))
    !           T(i,j)   = T(i,j) + DTIN(i,j)
    !    &                            + LCRCP * (DELTAL - DLIN(i,j))

    ! These are the condensation increments
    dqclpc2(i,j) = deltal - dqclin(i,j)

    if (l_fixbug_pc2_qcl_incr) then
      ! Ensure QCL increment cannot make QCL go negative.
      if ( qcl(i,j) + dqclpc2(i,j) < 0.0) then
        dqclpc2(i,j) = -qcl(i,j)
      end if
    end if

    dqpc2(i,j)   = - dqclpc2(i,j)
    dtpc2(i,j)   = lcrcp * dqclpc2(i,j)

  end do  ! i
end do  ! j

! End of the subroutine

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine pc2_hom_conv
end module pc2_hom_conv_mod
