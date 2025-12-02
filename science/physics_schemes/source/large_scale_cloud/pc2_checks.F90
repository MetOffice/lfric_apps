! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! PC2 Cloud Scheme: Checking cloud parameters
module pc2_checks_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'PC2_CHECKS_MOD'
contains

subroutine pc2_checks(                                                         &
!   Pressure related fields
 p_theta_levels,                                                               &
!   Prognostic Fields
 t, cf, cfl, cff, q, qcl, qcf,                                                 &
!   Logical control
 l_mixing_ratio,                                                               &
!   Sizes of input arrays
 row_length, rows, model_levels, offx, offy, halo_i, halo_j, qcf2, wtrac)

use mphys_inputs_mod,   only: l_casim, l_mcr_qcf2

use mphys_ice_mod,         only: thomo
use planet_constants_mod,  only: lcrcp, lfrcp, lsrcp, r, repsilon
use conversions_mod,       only: zerodegc
use water_constants_mod,   only: lc
use pc2_constants_mod,     only: cloud_rounding_tol,                           &
                                 one_over_qcf0,                                &
                                 min_in_cloud_qcf,                             &
                                 one_over_min_in_cloud_qcf,                    &
                                 condensate_limit,                             &
                                 wcgrow

use yomhook,               only: lhook, dr_hook
use parkind1,              only: jprb, jpim
use cloud_inputs_mod,      only: i_pc2_checks_cld_frac_method,                 &
                                 l_ensure_min_in_cloud_qcf
use science_fixes_mod,     only: l_pc2_checks_sdfix

use qsat_mod, only: qsat_wat, qsat_wat_mix

use free_tracers_inputs_mod, only: n_wtrac
use water_tracers_mod,       only: wtrac_type
use free_tracers_inputs_mod, only: l_wtrac
use wtrac_calc_ratio_mod,    only: wtrac_calc_ratio_fn
use pc2_checks_wtrac_mod,    only: pc2_checks_wtrac

use ereport_mod,             only: ereport
use errormessagelength_mod,  only: errormessagelength

implicit none

! Purpose:
!   This subroutine checks that cloud fractions, liquid water and
!   vapour contents take on physically reasonable values.
!
! Method:
!   Apply checks sequentially to the input values. It is more important
!   to ensure that liquid does not go negative than to ensure that
!   saturation deficit is correct.
!   (Note, water tracers are also updated here.)
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

! Sizes of the input arrays.
! These variables facilitate handling fields with or without halos.
integer, intent(in) ::                                                         &
 row_length, rows, model_levels, offx, offy, halo_i, halo_j

real(kind=real_umphys), intent(in) ::                                          &
 p_theta_levels(1-offx:row_length+offx,                                        &
                1-offy:rows+offy,                                              &
                     1:model_levels)
!    pressure at all points (Pa)

logical, intent(in) ::                                                         &
 l_mixing_ratio
!    Use mixing ratio formulation

real(kind=real_umphys), intent(in out) ::                                      &
 t(  1:row_length,                                                             &
     1:rows,                                                                   &
     1:model_levels),                                                          &
!    Temperature (K)
   cf( 1-halo_i:row_length+halo_i,                                             &
       1-halo_j:rows+halo_j,                                                   &
              1:model_levels),                                                 &
!    Total cloud fraction (no units)
   cfl(1-halo_i:row_length+halo_i,                                             &
       1-halo_j:rows+halo_j,                                                   &
              1:model_levels),                                                 &
!    Liquid cloud fraction (no units)
   cff(1-halo_i:row_length+halo_i,                                             &
       1-halo_j:rows+halo_j,                                                   &
              1:model_levels),                                                 &
!    Ice cloud fraction (no units)
   q(  1-halo_i:row_length+halo_i,                                             &
       1-halo_j:rows+halo_j,                                                   &
              1:model_levels),                                                 &
!    Vapour content (kg water per kg air)
   qcl(1-halo_i:row_length+halo_i,                                             &
       1-halo_j:rows+halo_j,                                                   &
              1:model_levels),                                                 &
!    Liquid content (kg water per kg air)
   qcf(1-halo_i:row_length+halo_i,                                             &
       1-halo_j:rows+halo_j,                                                   &
              1:model_levels),                                                 &
!    Ice content (kg water per kg air)
   qcf2(1-halo_i:row_length+halo_i,                                            &
       1-halo_j:rows+halo_j,                                                   &
              1:model_levels)
!    Ice content (crystals) (kg water per kg air)

! Water tracers (This field is not present when called from pc2_assim or from
! initial_pc2_checks if l_casim=T)
type(wtrac_type), optional, intent(in out) :: wtrac(n_wtrac)

!  External functions:

!  Local parameters and other physical constants------------------------
!
! Options for i_pc2_checks_cld_frac_method
integer, parameter :: original = 0
! Original. No change to liquid cloud fraction (CFL) is made when QCL is
! increased when qv>qsat.
integer, parameter :: force_cfl_cf_unity = 1
! Force CFL (and CF) to unity when creating some extra QCL when qv>qsat.
integer, parameter :: set_cfl_via_gi = 3
! Set CFL based on QCL value using Gultepe and Isaac (gi).

!  Local scalars--------------------------------------------------------

!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
real(kind=real_umphys) ::                                                      &
 alpha,                                                                        &
!      Rate of change of saturation specific humidity with
!      temperature calculated at dry-bulb temperature
!      (kg kg-1 K-1)
   al,                                                                         &
!      1 / (1 + alpha L/cp)  (no units)
   sd,                                                                         &
!      Saturation deficit (= aL (q - qsat(T)) )  (kg kg-1)
   cfl_old
!      temp store for old cfl value

!  (b) Others.
integer :: k,i,j       ! Loop counters: K - vertical level index
!                           I,J - horizontal position index
integer :: i_wt        ! Water tracer loop counter

integer :: ErrorStatus

!  Local dynamic arrays-------------------------------------------------
!    1 block of real workspace is required.
real(kind=real_umphys) ::                                                      &
 qsl_t(1:row_length,                                                           &
       1:rows)
!       Saturated specific humidity for dry bulb temperature T

! Water tracer to water ratios
real(kind=real_umphys), allocatable :: ratio_q(:,:,:,:)   ! For q
real(kind=real_umphys), allocatable :: ratio_qcl(:,:,:,:) ! For qcl

character(len=errormessagelength) :: cmessage

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real   (kind=jprb)            :: zhook_handle

character(len=*), parameter :: RoutineName='PC2_CHECKS'

!- End of Header

! ==Main Block==--------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set up water tracer to water ratio if using water tracers
if (l_wtrac) then
  ! First check that the structure wtrac is present otherwise abort as
  ! there is a problem
  if (.not. present(wtrac)) then
    ErrorStatus = 100
    cmessage = 'Water tracers are turned on but the wtrac field is not present'
    call ereport(RoutineName, ErrorStatus, cmessage)
  end if

  allocate(ratio_q(row_length,rows,model_levels,n_wtrac))
  allocate(ratio_qcl(row_length,rows,model_levels,n_wtrac))
  do i_wt = 1, n_wtrac
!$OMP  PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                              &
!$OMP  SHARED(i_wt, model_levels, rows, row_length, ratio_q, ratio_qcl,        &
!$OMP  wtrac, q, qcl) private(i, j, k)
    do k = 1,model_levels
      do j = 1, rows
        do i = 1, row_length
          ratio_q(i,j,k,i_wt)   = wtrac_calc_ratio_fn(i_wt,                    &
                                         wtrac(i_wt)%q(i,j,k), q(i,j,k))
          ratio_qcl(i,j,k,i_wt) = wtrac_calc_ratio_fn(i_wt,                    &
                                         wtrac(i_wt)%qcl(i,j,k), qcl(i,j,k))
        end do
      end do
    end do
!$OMP end PARALLEL do
  end do
end if

! Loop round levels to be processed
! Levels_do1:

!$OMP  PARALLEL do DEFAULT(SHARED) SCHEDULE(STATIC) private(i, j, k, al,       &
!$OMP  alpha, sd, cfl_old, qsl_t, i_wt)
do k = 1,model_levels

  ! ----------------------------------------------------------------------
  ! 1. Calculate Saturated Specific Humidity with respect to liquid water
  !    for dry bulb temperatures.
  ! ----------------------------------------------------------------------

  if ( l_mixing_ratio ) then
    call qsat_wat_mix(qsl_t,t(:,:,k),                                          &
          p_theta_levels(1:row_length,1:rows,k),                               &
          row_length,rows)
  else
    call qsat_wat(qsl_t,t(:,:,k),                                              &
          p_theta_levels(1:row_length,1:rows,k),                               &
          row_length,rows)
  end if

  ! Rows_do1:
  do j = 1, rows
    ! Row_length_do1:
    do i = 1, row_length

      !----------------------------------------------------------------------
      ! 2. Calculate the saturation deficit.
      ! ----------------------------------------------------------------------

      ! Need to estimate the rate of change of saturated specific humidity
      ! with respect to temperature (alpha) first, then use this to calculate
      ! factor aL.
      alpha=repsilon*lc*qsl_t(i,j)/(r*t(i,j,k)**2)
      al=1.0/(1.0+lcrcp*alpha)

      ! Calculate the saturation deficit SD

      sd=al*(qsl_t(i,j)-q(i,j,k))

      ! ----------------------------------------------------------------------
      !  3. Checks are applied here for liquid cloud
      ! ----------------------------------------------------------------------

      ! Earlier versions checked whether saturation deficit is zero (or less
      ! than zero). If so, then the liquid cloud fraction was forced to one.
      ! This check has been suspended for numerical reasons.
      !            if (SD  <=  0.0 .or. CFL(i,j,k)  >   1.0) then

      ! Instead, check simply whether input values of liquid cloud fraction
      ! are, or are between, zero and one. If not, adjust them to zero or one.
      ! Adjust the total cloud fractions accordingly.

      if (cfl(i,j,k) > (1.0 - cloud_rounding_tol)) then
        cfl(i,j,k)=1.0
        cf(i,j,k) =1.0
      end if

      ! Check also whether the liquid water content is less than zero, and
      ! set liquid cloud fraction to zero if it is.

      if (qcl(i,j,k) < condensate_limit .or.                                   &
          cfl(i,j,k) < cloud_rounding_tol) then
        cfl(i,j,k)=0.0
        cf(i,j,k) =cff(i,j,k)
      end if

      ! Check whether the saturation deficit is less than zero. If it is
      ! then condense some liquid to bring the saturation deficit to zero.
      ! Adjust the temperature for the latent heating.

      if (sd < 0.0) then
        q(i,j,k)   = q(i,j,k)   + sd
        qcl(i,j,k) = qcl(i,j,k) - sd
        t(i,j,k)   = t(i,j,k)   - sd * lcrcp

        if (l_wtrac) then    ! Update water tracers for phase change
           ! Note, sd < 0, so vapour is the source
          do i_wt = 1, n_wtrac
            wtrac(i_wt)%q(i,j,k)   =                                           &
                    wtrac(i_wt)%q(i,j,k)   + ratio_q(i,j,k,i_wt)*sd
            wtrac(i_wt)%qcl(i,j,k) =                                           &
                    wtrac(i_wt)%qcl(i,j,k) - ratio_q(i,j,k,i_wt)*sd
          end do
        end if

        ! In original PC2 (i_pc2_checks_cld_frac_method=original)
        ! there is no change to the cloud fraction
        ! as a result of the above QCL change.
        if (i_pc2_checks_cld_frac_method /= original) then
          ! Various options for adjusting the cloud fraction.
          if (i_pc2_checks_cld_frac_method == force_cfl_cf_unity) then
            ! Force CFL to 1 and hence also set CF to 1.
            cfl(i,j,k) = 1.0
            cf(i,j,k)  = 1.0
          else
            ! Increase the cloud fraction in order to
            ! keep the in-cloud condensate amount the same
            ! providing there is a well-defined in cloud value.
            cfl_old = cfl(i,j,k)
            if (cfl(i,j,k) > 0.0 .and. (qcl(i,j,k)+sd) > 0.0) then
               ! There is a well defined in cloud value
              cfl(i,j,k) = max(0.0,min(1.0,                                    &
                           qcl(i,j,k)*cfl(i,j,k)/(qcl(i,j,k)+sd) ))
            else

              if (i_pc2_checks_cld_frac_method >= set_cfl_via_gi) then
                ! Initiate cloud using Gultepe and Isaac formulation.
                ! NB 5.57 is valid for  10 km gridbox.
                !    4.37 is valid for 100 km gridbox.
                ! For now just use a parameter value valid for dx=10km
                ! as variation with gridbox size is not huge and the actual
                ! value of qcl is probably not perfect anyway,
                ! since just found from excess over 100% RH.
                ! Ensure that cloud cover is at least 0.5 though, since
                ! we are dealing with cases with RH>100%.
                cfl(i,j,k) = max(min(5.57*((qcl(i,j,k)*1000.0)**0.775),1.0),0.5)
              else
                 ! grow cloud at a fixed amount
                cfl(i,j,k) = max(0.0,min(1.0,                                  &
                             qcl(i,j,k)/wcgrow ))
              end if

            end if
            ! update total cloud fraction
            cf(i,j,k) = max(min(cf(i,j,k)+cfl(i,j,k)-cfl_old,1.0),0.0)
          end if
        end if ! i_pc2_checks_cld_frac_method /= original

      end if ! sd < 0.0

      ! Check whether the saturation deficit
      ! is greater than zero but the liquid cloud fraction is one. If it is
      ! then evaporate some liquid (provided there is enough
      ! liquid) to bring the saturation deficit to zero. Adjust the
      ! temperature for the latent heating.

      if (sd > 0.0 .and. cfl(i,j,k) == 1.0) then
        if (qcl(i,j,k) > sd) then

          q(i,j,k)   = q(i,j,k)   + sd
          qcl(i,j,k) = qcl(i,j,k) - sd
          t(i,j,k)   = t(i,j,k)   - sd * lcrcp

          if (l_wtrac) then   ! Update water tracers for phase change
            ! Note, sd > 0, so liquid is the source
            do i_wt = 1, n_wtrac
              wtrac(i_wt)%qcl(i,j,k) =                                         &
                   wtrac(i_wt)%qcl(i,j,k) - ratio_qcl(i,j,k,i_wt)*sd
              wtrac(i_wt)%q(i,j,k)   =                                         &
                   wtrac(i_wt)%q(i,j,k)   + ratio_qcl(i,j,k,i_wt)*sd
            end do
          end if

        else if ( l_pc2_checks_sdfix ) then
          ! Under switch for bug-fix:
          ! If total cloud-cover, subsaturated, but qcl <= sd,
          ! we should just evaporate all the remaining liquid
          ! and reset the cloud-fraction to zero.

          q(i,j,k)   = q(i,j,k)   + qcl(i,j,k)
          t(i,j,k)   = t(i,j,k)   - qcl(i,j,k) * lcrcp
          qcl(i,j,k) = 0.0
          cfl(i,j,k) = 0.0
          cf(i,j,k)  = cff(i,j,k)

          if (l_wtrac) then    ! Update water tracers
            do i_wt = 1, n_wtrac
              wtrac(i_wt)%q(i,j,k)   = wtrac(i_wt)%q(i,j,k)                    &
                                        + wtrac(i_wt)%qcl(i,j,k)
              wtrac(i_wt)%qcl(i,j,k) = 0.0
            end do
          end if
        end if

      end if

      ! Check whether the liquid content is less than zero, or whether it is
      ! greater than zero but the liquid cloud fraction is zero. If so then
      ! condense or evaporate liquid to bring the liquid water to zero. Adjust
      ! the temperature for latent heating.

      if (qcl(i,j,k) < condensate_limit .or.                                   &
         (qcl(i,j,k) >  0.0 .and. cfl(i,j,k) == 0.0) ) then

        q(i,j,k)   = q(i,j,k) + qcl(i,j,k)
        t(i,j,k)   = t(i,j,k) - qcl(i,j,k) * lcrcp
        qcl(i,j,k) = 0.0

        if (l_wtrac) then     ! Update water tracers
          do i_wt = 1, n_wtrac
            wtrac(i_wt)%q(i,j,k)   = wtrac(i_wt)%q(i,j,k)                      &
                                        + wtrac(i_wt)%qcl(i,j,k)
            wtrac(i_wt)%qcl(i,j,k) = 0.0
          end do
        end if

      end if

      ! ----------------------------------------------------------------------
      !  4. Check that ice content and ice cloud fraction are sensible.
      ! ----------------------------------------------------------------------

      ! Check whether ice content is zero (or less than zero). If so then
      ! force the ice cloud fraction to zero. Also check whether input values
      ! of ice cloud fraction are, or are between, zero and one. If not,
      ! adjust them to zero or one. Adjust the total cloud fractions
      ! accordingly.

      if (cff(i,j,k) > (1.0 - cloud_rounding_tol)) then
        cff(i,j,k) = 1.0
        cf(i,j,k)  = 1.0
      end if

      if (l_mcr_qcf2) then
         ! ---------------------------------------------
         ! Account for both ice prognostics in check
         !----------------------------------------------

        if ((qcf(i,j,k) < condensate_limit .and.                               &
             qcf2(i,j,k) < condensate_limit) .or.                              &
             cff(i,j,k) < cloud_rounding_tol) then
          cff(i,j,k) = 0.0
          cf(i,j,k)  = cfl(i,j,k)
        end if

        ! If ice content is negative then condense some vapour to remove the
        ! negative part. Adjust the temperature for the latent heat.

        if (qcf(i,j,k) < condensate_limit) then
          q(i,j,k)   = q(i,j,k) + qcf(i,j,k)
          t(i,j,k)   = t(i,j,k) - qcf(i,j,k) * lsrcp
          qcf(i,j,k) = 0.0
          if (l_wtrac) then    ! Update water tracers
            do i_wt = 1, n_wtrac
              wtrac(i_wt)%q(i,j,k)   = wtrac(i_wt)%q(i,j,k)                    &
                                        + wtrac(i_wt)%qcf(i,j,k)
              wtrac(i_wt)%qcf(i,j,k) = 0.0
            end do
          end if
        end if

        if (qcf2(i,j,k) < condensate_limit) then
          q(i,j,k)   = q(i,j,k) + qcf2(i,j,k)
          t(i,j,k)   = t(i,j,k) - qcf2(i,j,k) * lsrcp
          qcf2(i,j,k) = 0.0
          if (l_wtrac) then     ! Update water tracers
            do i_wt = 1, n_wtrac
              wtrac(i_wt)%q(i,j,k)    = wtrac(i_wt)%q(i,j,k)                   &
                                        + wtrac(i_wt)%qcf2(i,j,k)
              wtrac(i_wt)%qcf2(i,j,k) = 0.0
            end do
          end if
        end if

        if ((qcf(i,j,k) +qcf2(i,j,k))> 0.0 .and. cff(i,j,k) == 0.0) then
          cff(i,j,k) = (qcf(i,j,k)+qcf2(i,j,k)) * one_over_qcf0
          !             CF is not adjusted here but is checked below
        end if

        if (l_ensure_min_in_cloud_qcf) then
          if (cff(i,j,k) > 0.0) then
            if ((qcf(i,j,k)+qcf2(i,j,k))/cff(i,j,k) < min_in_cloud_qcf) then
              cff(i,j,k)=(qcf(i,j,k)+qcf2(i,j,k))*one_over_min_in_cloud_qcf
            end if
          end if
        end if

      else ! l_mcr_qcf2

        ! ----------------------------------------------------
        ! Account for just one ice prognostic (qcf) in check
        !-----------------------------------------------------

        if (qcf(i,j,k) < condensate_limit .or.                                 &
            cff(i,j,k) < cloud_rounding_tol) then
          cff(i,j,k) = 0.0
          cf(i,j,k)  = cfl(i,j,k)
        end if

        ! If ice content is negative then condense some vapour to remove the
        ! negative part. Adjust the temperature for the latent heat.

        if (qcf(i,j,k) < condensate_limit) then

          q(i,j,k)   = q(i,j,k) + qcf(i,j,k)
          t(i,j,k)   = t(i,j,k) - qcf(i,j,k) * lsrcp
          qcf(i,j,k) = 0.0

          if (l_wtrac) then     ! Update water tracers
            do i_wt = 1, n_wtrac
              wtrac(i_wt)%q(i,j,k)   = wtrac(i_wt)%q(i,j,k)                    &
                                        + wtrac(i_wt)%qcf(i,j,k)
              wtrac(i_wt)%qcf(i,j,k) = 0.0
            end do
          end if
        end if

        ! If ice content is positive but ice cloud fraction negative, create
        ! some ice cloud fraction.

        if (qcf(i,j,k) > 0.0 .and. cff(i,j,k) == 0.0) then
          cff(i,j,k) = qcf(i,j,k) * one_over_qcf0
          !             CF is not adjusted here but is checked below
        end if

        if (l_ensure_min_in_cloud_qcf) then
          if (cff(i,j,k) > 0.0) then
            if (qcf(i,j,k)/cff(i,j,k) < min_in_cloud_qcf) then
              cff(i,j,k)=qcf(i,j,k)*one_over_min_in_cloud_qcf
            end if
          end if
        end if

      end if ! l_mcr_qcf2

      ! ----------------------------------------------------------------------
      !  5. Check that total cloud fraction is sensible.
      ! ----------------------------------------------------------------------

      ! Total cloud fraction must be bounded by
      ! i) The maximum of the ice and liquid cloud fractions (maximum overlap)

      if (cf(i,j,k) < max(cfl(i,j,k),cff(i,j,k))) then
        cf(i,j,k) = max(cfl(i,j,k),cff(i,j,k))
      end if

      ! ii) The sum of the ice and liquid cloud fractions or one, whichever
      ! is the lower (minimum overlap)

      if (cf(i,j,k) > min( (cfl(i,j,k)+cff(i,j,k)),1.0 ) ) then
        cf(i,j,k) = min( (cfl(i,j,k)+cff(i,j,k)),1.0 )
      end if

      ! ----------------------------------------------------------------------
      !  6. Homogeneous nucleation.
      ! ----------------------------------------------------------------------

      ! This process is switched off for runs with CASIM as it is dealt with
      ! separately within CASIM.
      if (t(i,j,k) < (zerodegc+thomo) .and. (.not. l_casim)) then

        ! Turn all liquid to ice
        if (qcl(i,j,k) > 0.0) then
          cff(i,j,k) = cf(i,j,k)
          cfl(i,j,k) = 0.0
        end if

        qcf(i,j,k) = qcf(i,j,k) + qcl(i,j,k)
        t(i,j,k)   = t(i,j,k)   + qcl(i,j,k) * lfrcp
        qcl(i,j,k) = 0.0

        if (l_wtrac) then   ! Update water tracers
          do i_wt = 1, n_wtrac
            wtrac(i_wt)%qcf(i,j,k)   = wtrac(i_wt)%qcf(i,j,k)                  &
                                           + wtrac(i_wt)%qcl(i,j,k)
            wtrac(i_wt)%qcl(i,j,k)   = 0.0
          end do
        end if
      end if ! t < zerodegc+thomo / not l_casim

    end do !i

  end do !j

end do !k
!$OMP end PARALLEL do

if (l_wtrac) then

  ! Deallocate working arrays
  deallocate(ratio_qcl)
  deallocate(ratio_q)

  ! Call routine to ensure that water tracer condensate remains above zero
  call pc2_checks_wtrac(wtrac)

end if

! End of the subroutine

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine pc2_checks
end module pc2_checks_mod
