! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  PC2 Cloud Scheme: Reset clouds for extreme relative total humidity.

module pc2_checks2_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'PC2_CHECKS2_MOD'
contains

subroutine pc2_checks2(                                                        &
!      Pressure related fields
 p_theta_levels, rhcrit,                                                       &
!      Array dimensions
 rhc_row_length,rhc_rows,                                                      &
!      Prognostic Fields
 t, cf, cfl, cff, q, qcl,                                                      &
!      Logical control
 l_mixing_ratio)

use water_constants_mod,   only: lc
use planet_constants_mod,  only: lcrcp, r, repsilon
use yomhook,               only: lhook, dr_hook
use parkind1,              only: jprb, jpim
use atm_fields_bounds_mod, only: pdims, tdims
use cloud_inputs_mod,      only: cloud_pc2_tol, cloud_pc2_tol_2
use qsat_mod,              only: qsat_wat, qsat_wat_mix

use free_tracers_inputs_mod, only: l_wtrac
use wtrac_pc2_mod,           only: wtrac_pc2

implicit none

! Purpose:
!   Check that cloud fraction is either zero or one when relative
!   total humidity is small or large.

! Method:
!   Calculate relative total humidity, compare to RHcrit and adjust
!   the cloud variables appropriately.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_cloud
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: PC2 Cloud Scheme Documentation

!  Subroutine Arguments:------------------------------------------------
integer ::                                                                     &
                      !, intent(in)
 rhc_row_length,rhc_rows
!       Dimensions of the RHCRIT variable.

real(kind=real_umphys) ::                                                      &
                      !, intent(in)
 p_theta_levels(pdims%i_start:pdims%i_end,                                     &
                pdims%j_start:pdims%j_end,                                     &
                pdims%k_start:pdims%k_end),                                    &
!       pressure at all points (Pa)
   rhcrit(        rhc_row_length,                                              &
                  rhc_rows,                                                    &
                              1:tdims%k_end),                                  &
!       Critical relative humidity.  See the the paragraph incorporating
!       eqs P292.11 to P292.14; the values need to be tuned for the give
!       set of levels.
   cff(           tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                              1:tdims%k_end)
!       Ice cloud fraction (no units)

logical ::                                                                     &
                      !, intent(in)
 l_mixing_ratio   ! Use mixing ratio formulation

real(kind=real_umphys) ::                                                      &
                      !, intent(INOUT)
 t(             tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,                                     &
                            1:tdims%k_end),                                    &
!       Temperature (K)
   cf(            tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                              1:tdims%k_end),                                  &
!       Total cloud fraction (no units)
   cfl(           tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                              1:tdims%k_end),                                  &
!       Liquid cloud fraction (no units)
   q(             tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                              1:tdims%k_end),                                  &
!       Vapour content (kg water per kg air)
   qcl(           tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                              1:tdims%k_end)
!       Liquid content (kg water per kg air)

!  External functions:

!  Local parameters and other physical constants------------------------
real(kind=real_umphys) :: c_thresh_low
real(kind=real_umphys) :: c_thresh_low_2
!       Low cloud fraction thresholds
real(kind=real_umphys) :: c_thresh_high
real(kind=real_umphys) :: c_thresh_high_2
!       High cloud fraction thresholds

!  Local scalars--------------------------------------------------------

!  (a)  Scalars effectively expanded to workspace by the Cray (using
!       vector registers).
real(kind=real_umphys) ::                                                      &
 alpha,                                                                        &
!       Rate of change of saturation specific humidity with
!       temperature calculated at dry-bulb temperature (kg kg-1 K-1)
   al,                                                                         &
!       1 / (1 + alpha L/cp)  (no units)
   rht,                                                                        &
!       Relative total humidity
   sd
!       Saturation deficit

!  (b)  Others.
integer :: k,i,j,                                                              &
!       Loop counters: K - vertical level index
!       I,J - horizontal position index
   irhi,irhj,                                                                  &
!       Indices for RHcrit array
   multrhc
!       Zero if (rhc_row_length*rhc_rows) le 1, else 1

!  Local dynamic arrays-------------------------------------------------
!    3 blocks of real workspace are required.
real(kind=real_umphys) ::                                                      &
  qsl_t,                                                                       &
!       Saturated specific humidity for temperature T
  qsl_tl,                                                                      &
!       Saturated specific humidity for liquid temperature TL
  tl
!    Liquid temperature (= T - L/cp QCL)  (K)

!- End of Header

! Set up a flag to state whether RHcrit is a single parameter or defined
! on all points.

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real   (kind=jprb)            :: zhook_handle

character(len=*), parameter :: RoutineName='PC2_CHECKS2'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

if (rhc_row_length * rhc_rows > 1) then
  multrhc=1
else
  multrhc=0
end if

! Set cloud-fraction thresholds based on namelist inputs
c_thresh_low    = cloud_pc2_tol
c_thresh_low_2  = cloud_pc2_tol_2
c_thresh_high   = 1.0 - cloud_pc2_tol
c_thresh_high_2 = 1.0 - cloud_pc2_tol_2

! ==Main Block==--------------------------------------------------------

! Loop round levels to be processed
! Levels_do1:

!$OMP  PARALLEL do DEFAULT(SHARED) SCHEDULE(STATIC) private(i, j, k,           &
!$OMP  irhi, irhj, rht, alpha, al, sd, qsl_t, qsl_tl,                          &
!$OMP  tl)
do k = 1, tdims%k_end
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end

      if ( cfl(i,j,k)  >=  c_thresh_high .or.                                  &
           ( cfl(i,j,k)  <=  c_thresh_low .and.                                &
             cfl(i,j,k) > 0.0 ) ) then

        ! Set up index pointers to critical relative humidity value

        irhi = (multrhc * (i - 1)) + 1
        irhj = (multrhc * (j - 1)) + 1

        ! Calculate Saturated Specific Humidity with respect to liquid water
        ! for liquid temperature.
        tl = t(i,j,k)-lcrcp*qcl(i,j,k)
        if ( l_mixing_ratio ) then
          call qsat_wat_mix(qsl_tl, tl, p_theta_levels(i,j,k))
        else
          call qsat_wat(qsl_tl, tl, p_theta_levels(i,j,k))
        end if

        ! Calculate relative total humidity with respect to the liquid
        ! temperature and threshold relative humidity
        rht=(q(i,j,k)+qcl(i,j,k))/qsl_tl

        ! ----------------------------------------------------------------------
        ! 3. Determine whether resetting is required and calculate appropriate
        !    increments if this is the case.
        ! ----------------------------------------------------------------------

        ! Is the relative total humidity greater than the threshold amount?
        ! If so, evaporate some liquid to set the saturation deficit to zero.

        if ( (rht  >   (2.0-rhcrit(irhi,irhj,k)) .and.                         &
          (cfl(i,j,k) >= c_thresh_high)) .or.                                  &
          (cfl(i,j,k)  >=  c_thresh_high_2) ) then
          cfl(i,j,k) = 1.0
          cf(i,j,k)  = 1.0

          ! Calculate Saturated Specific Humidity with respect to liquid water
          ! for dry bulb temperature.
          if ( l_mixing_ratio ) then
            call qsat_wat_mix(qsl_t, t(i,j,k), p_theta_levels(i,j,k))
          else
            call qsat_wat(qsl_t, t(i,j,k), p_theta_levels(i,j,k))
          end if

          ! Calculate the saturation deficit
          alpha = repsilon * lc * qsl_t /                                      &
                (r * t(i,j,k) ** 2)
          al    = 1.0 / (1.0 + lcrcp * alpha)
          sd    = al * (qsl_t - q(i,j,k))

          ! Update the water contents

          qcl(i,j,k) = qcl(i,j,k) - sd
          q(i,j,k)   = q(i,j,k)   + sd
          t(i,j,k)   = t(i,j,k)   - sd * lcrcp

          if (l_wtrac) then
            wtrac_pc2%q_cond(i,j,k) = -sd
          end if

        end if

        ! Is the relative total humidity less than the threshold amount?
        ! If so, evaporate all the liquid.

        if ( rht  <   (rhcrit(irhi,irhj,k)) .and.                              &
          (cfl(i,j,k)  <=  c_thresh_low) ) then
          cfl(i,j,k) = 0.0
          cf(i,j,k)  = cff(i,j,k)
          q(i,j,k)   = q(i,j,k) + qcl(i,j,k)

          if (l_wtrac) then
            wtrac_pc2%q_cond(i,j,k) = wtrac_pc2%q_cond(i,j,k) - qcl(i,j,k)
          end if

          t(i,j,k)   = t(i,j,k) - qcl(i,j,k) * lcrcp
          qcl(i,j,k) = 0.0
        end if

      end if ! cloud fraction

      if ((cfl(i,j,k)  <=  c_thresh_low_2)  ) then
        cfl(i,j,k) = 0.0
        cf(i,j,k)  = cff(i,j,k)
        q(i,j,k)   = q(i,j,k) + qcl(i,j,k)

        if (l_wtrac) then
          wtrac_pc2%q_cond(i,j,k) = wtrac_pc2%q_cond(i,j,k) - qcl(i,j,k)
        end if

        t(i,j,k)   = t(i,j,k) - qcl(i,j,k) * lcrcp
        qcl(i,j,k) = 0.0
      end if

    end do

  end do

end do
!$OMP end PARALLEL do

! End of the subroutine

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine pc2_checks2
end module pc2_checks2_mod
