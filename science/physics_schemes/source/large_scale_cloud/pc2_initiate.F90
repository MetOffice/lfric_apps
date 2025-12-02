! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  PC2 Cloud Scheme: Initiation

module pc2_initiate_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'PC2_INITIATE_MOD'
contains

subroutine pc2_initiate(                                                       &
!      Pressure related fields
 p_theta_levels, cumulus, rhcrit,                                              &
!      Array dimensions
 nlevels,                                                                      &
 rhc_row_length,rhc_rows,zlcl_mixed,r_theta_levels,                            &
!      Prognostic Fields
   t, cf, cfl, cff, q, qcl, rhts,                                              &
!      Logical control
   l_mixing_ratio)

use conversions_mod,       only: zerodegc
use water_constants_mod,   only: lc
use planet_constants_mod,  only: lcrcp, r, repsilon
use yomhook,               only: lhook, dr_hook
use parkind1,              only: jprb, jpim
use atm_fields_bounds_mod, only: pdims, tdims, pdims_l
use pc2_constants_mod,     only: init_iterations, rhcrit_tol,                  &
                                 pdf_power, rhcpt_tke_based,                   &
                                 pc2init_logic_original,                       &
                                 pc2init_logic_simplified,                     &
                                 pc2init_logic_smooth
use cloud_inputs_mod,      only: i_rhcpt, i_pc2_init_logic, cloud_pc2_tol
use qsat_mod,              only: qsat_wat, qsat_wat_mix
use pc2_total_cf_mod,      only: pc2_total_cf

implicit none

! Purpose:
!   Initiate liquid and total cloud fraction and liquid water content

! Method:
!   Uses the method proposed in Annex C of the PC2 cloud scheme project
!   report, which considers a Smith-like probability density
!   distribution whose width is given by a prescribed value of RHcrit.
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
 nlevels,                                                                      &
!       No. of levels being processed.
    rhc_row_length,rhc_rows
!       Dimensions of the RHCRIT variable.

logical ::                                                                     &
                      !, intent(in)
 cumulus(       tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end),                                    &
!       Is this a boundary layer cumulus point
   l_mixing_ratio  ! Use mixing ratio formulation

! height of lcl in a well-mixed BL (types 3 or 4), 0 otherwise
real(kind=real_umphys) :: zlcl_mixed( pdims%i_start:pdims%i_end,               &
                                      pdims%j_start:pdims%j_end)

real(kind=real_umphys) ::                                                      &
                      !, intent(in)
 p_theta_levels(pdims%i_start:pdims%i_end,                                     &
                pdims%j_start:pdims%j_end,                                     &
                nlevels),                                                      &
!       Pressure at all points (Pa)
   r_theta_levels(pdims_l%i_start:pdims_l%i_end,                               &
                  pdims_l%j_start:pdims_l%j_end,                               &
                  0:nlevels),                                                  &
!       Pressure at all points (Pa)
   rhcrit(rhc_row_length,rhc_rows,nlevels),                                    &
!       Critical relative humidity.  See the the paragraph incorporating
!       eqs P292.11 to P292.14; the values need to be tuned for the give
!       set of levels.
   cff(           tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                  nlevels)
!       Ice cloud fraction (no units)

real(kind=real_umphys) ::                                                      &
                      !, intent(INOUT)
 t(             tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,                                     &
                nlevels),                                                      &
!       Temperature (K)
   cf(            tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                  nlevels),                                                    &
!       Total cloud fraction (no units)
   cfl(           tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                  nlevels),                                                    &
!       Liquid cloud fraction (no units)
   q(             tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                  nlevels),                                                    &
!       Vapour content (kg water per kg air)
   qcl(           tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                  nlevels),                                                    &
!       Liquid content (kg water per kg air)
   rhts(          tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                  nlevels)
!       Variable carrying initial RHT wrt TL from start of timestep

!  External functions:

!  Local scalars--------------------------------------------------------

!  (a)  Scalars effectively expanded to workspace by the Cray (using
!       vector registers).
real(kind=real_umphys) ::                                                      &
 alpha,                                                                        &
!       Rate of change of saturation specific humidity with
!       temperature calculated at dry-bulb temperature
!       (kg kg-1 K-1)
   al,                                                                         &
!       1 / (1 + alpha L/cp)  (no units)
   bs,                                                                         &
!       Width of distribution (kg kg-1)
   c_1,                                                                        &
!       New liquid cloud fraction
   deltal,                                                                     &
!       Change in liquid (kg kg-1)
   descent_factor,                                                             &
!       Rate at which to relax to new liquid water content
   l_bs,                                                                       &
!       New liquid water content divided by PDF width (no units)
   l_out,                                                                      &
!       New liquid water content (kg kg-1)
   q_out,                                                                      &
!       New vapour content (kg kg-1)
   qn,                                                                         &
!       Normalized value of QC for the probability density function.
   rht,                                                                        &
!       Total relative humidity (liquid+vapour)/qsat (kg kg-1)
   rh0,                                                                        &
!       Equivalent critical relative humidity
   qc,                                                                         &
!       Qc = al ( q + qcl - qsl(Tl) )
   frac_init
!       Fraction of final liquid water initiated this timestep

!  (b)  Others.
integer :: k,i,j,l,                                                            &
!       Loop counters: K   - vertical level index
!                      I,J - horizontal position index
!                      L   - counter for iterations
          irhi,irhj,                                                           &
!                      Indices for RHcrit array
          multrhc,                                                             &
!                      Zero if (rhc_row_length*rhc_rows) le 1, else 1
          npti, npt
!                      Number of points to iterate over

!  Local dynamic arrays-------------------------------------------------
real(kind=real_umphys) ::                                                      &
 qsl_tl,                                                                       &
!       Saturated specific humidity for liquid temperature TL
   tl_c,                                                                       &
!       Liquid temperature (= T - L/cp QCL)  (K)
   cf_c(    tdims%i_len*tdims%j_len),                                          &
!       Total cloud fraction on compressed points (no units)
   cfl_c(   tdims%i_len*tdims%j_len),                                          &
!       Liquid cloud fraction on compressed points (no units)
   cff_c(   tdims%i_len*tdims%j_len),                                          &
!       Ice cloud fraction on compressed points (no units)
   deltacl_c(                                                                  &
            tdims%i_len*tdims%j_len),                                          &
!       Change in liquid cloud fraction (no units)
   deltacf_c(                                                                  &
            tdims%i_len*tdims%j_len),                                          &
!       Change in ice cloud fraction (no units)
   q_c(     tdims%i_len*tdims%j_len),                                          &
!       Vapour content on compressed points (kg kg-1)
   qcl_c(   tdims%i_len*tdims%j_len),                                          &
!       Liquid water content on compressed points (kg kg-1)
   qn_c(    tdims%i_len*tdims%j_len),                                          &
!      QN on compressed points (no units)
   qsl_t_c,                                                                    &
!       Saturated specific humidity for dry-bulb temperature T on
!       compressed points (kg kg-1)
   qsl_tl_c(                                                                   &
            tdims%i_len*tdims%j_len),                                          &
!       Saturated specific humidity for liquid temperature TL on
!       compressed points (kg kg-1)
   rh0_c(   tdims%i_len*tdims%j_len),                                          &
!       Equivalent critical relative humidity on compressed points
!       (no units)
   t_c(     tdims%i_len*tdims%j_len)
!       Temperature on compressed points (K)

integer ::                                                                     &
 ni(      tdims%i_len*tdims%j_len),                                            &
 nj(      tdims%i_len*tdims%j_len),                                            &
!       Compressed point counters
   ind_i(   tdims%i_len*tdims%j_len),                                          &
   ind_j(   tdims%i_len*tdims%j_len)
!       Compressed point counters

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real   (kind=jprb)            :: zhook_handle

character(len=*), parameter :: RoutineName='PC2_INITIATE'

!- End of Header

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set up a flag to state whether RHcrit is a single parameter or defined
! on all points.

if (rhc_row_length*rhc_rows > 1) then
  multrhc=1
else
  multrhc=0
end if

! ==Main Block==--------------------------------------------------------

! Loop round levels to be processed
! Levels_do1:

!$OMP  PARALLEL do DEFAULT(SHARED) SCHEDULE(STATIC) private(npt,               &
!$OMP  npti, ind_i, ind_j, irhj, irhi, rht, rh0, qn, c_1, ni, nj, t_c,         &
!$OMP  qn_c, rh0_c, qsl_tl_c, cf_c, cfl_c, cff_c,                              &
!$OMP  qcl_c, q_c, deltacl_c, deltacf_c, qsl_t_c, l_out, l_bs, al,             &
!$OMP  deltal, descent_factor, q_out, bs, i, j, k, l, qsl_tl, tl_c,            &
!$OMP  alpha, qc, frac_init)
do k = 1, nlevels

  if ( i_pc2_init_logic == pc2init_logic_simplified ) then

    ! Determine points to calculate
    npt = 0
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        if ( ( cfl(i,j,k) < cloud_pc2_tol                                      &
              .and. (r_theta_levels(i,j,k)-r_theta_levels(i,j,0))              &
                    > zlcl_mixed(i,j) )                                        &
            .or. cfl(i,j,k) > 1.0 - cloud_pc2_tol ) then
          npt = npt+1
          ind_i(npt) = i
          ind_j(npt) = j
        end if
      end do
    end do

  else if ( i_pc2_init_logic == pc2init_logic_original ) then

    ! Determine points to calculate
    npt = 0
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        if ( ( (cfl(i,j,k)  ==  0.0                                            &
           .or. (cfl(i,j,k)  <   0.05 .and. t(i,j,k)  <   zerodegc))           &
           .and. (.not. cumulus(i,j))                                          &
           .and. (r_theta_levels(i,j,k)-r_theta_levels(i,j,0))                 &
                    > zlcl_mixed(i,j) ) .or.                                   &
                      cfl(i,j,k) == 1.0 ) then
          npt = npt+1
          ind_i(npt) = i
          ind_j(npt) = j
        end if
      end do
    end do

  else if ( i_pc2_init_logic == pc2init_logic_smooth ) then

    ! Determine points to calculate
    ! Smooth initiation logic; need to calculate qsat at all points,
    ! so just set the compression list to all points on this model-level

    npt = 0
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        npt = npt+1
        ind_i(npt) = i
        ind_j(npt) = j
      end do
    end do

  end if  ! ( i_pc2_init_logic )

  ! Set number of points to perform the compressed calculation over
  ! to zero
  npti=0.0

  if (npt > 0) then
    do i = 1, npt

      ! ----------------------------------------------------------------------
      ! 3. Calculate equivalent critical relative humidity values (RH0) for
      !    initiation.
      ! ----------------------------------------------------------------------

      ! Set up index pointers to critical relative humidity value

      irhi = (multrhc * (ind_i(i) - 1)) + 1
      irhj = (multrhc * (ind_j(i) - 1)) + 1

      ! ----------------------------------------------------------------------
      ! 2. Calculate Saturated Specific Humidity with respect to liquid water
      !    for liquid temperatures.
      ! ----------------------------------------------------------------------
      tl_c = t(ind_i(i),ind_j(i),k) - lcrcp * qcl(ind_i(i),ind_j(i),k)

      if ( l_mixing_ratio ) then
        call qsat_wat_mix(qsl_tl, tl_c, p_theta_levels(ind_i(i),ind_j(i),k))
      else
        call qsat_wat(qsl_tl, tl_c, p_theta_levels(ind_i(i),ind_j(i),k))
      end if

      ! Calculate relative total humidity with respect to the liquid
      ! temperature

      rht=(q(ind_i(i),ind_j(i),k)+qcl(ind_i(i),ind_j(i),k))                    &
          /qsl_tl

      ! Initialize the equivalent critical relative humidity to a dummy
      ! negative value

      rh0 = -1.0

      if ( i_pc2_init_logic == pc2init_logic_simplified ) then

        ! Do we need to force some initiation from zero liquid cloud amount?
        ! If so, set RH0 equal the critical relative humidity.

        if (      cfl(ind_i(i),ind_j(i),k) < cloud_pc2_tol                     &
            .and. (r_theta_levels(ind_i(i),ind_j(i),k)                         &
                  -r_theta_levels(ind_i(i),ind_j(i),0))                        &
                  > zlcl_mixed(ind_i(i),ind_j(i))                              &
            .and. rht > rhts(ind_i(i),ind_j(i),k)                              &
            .and. rht > (rhcrit(irhi,irhj,k)+rhcrit_tol)  ) then
          rh0 = rhcrit(irhi,irhj,k)
        end if

        ! Do we need to force some initiation from total liquid cloud cover?
        ! If so, set RH0 equal the critical relative humidity.

        if (      cfl(ind_i(i),ind_j(i),k) > 1.0 - cloud_pc2_tol               &
            .and. rht < rhts(ind_i(i),ind_j(i),k)                              &
            .and. rht < (2.0-rhcrit(irhi,irhj,k) - rhcrit_tol) ) then
          rh0 = rhcrit(irhi,irhj,k)
        end if

      else if ( i_pc2_init_logic == pc2init_logic_original ) then

        ! Do we need to force some initiation from zero liquid cloud amount?
        ! If so, set RH0 to equal the critical relative humidity

        if ( (cfl(ind_i(i),ind_j(i),k)  ==  0.0                                &
          .or. (cfl(ind_i(i),ind_j(i),k)  <   0.05                             &
                    .and. t(ind_i(i),ind_j(i),k) <   zerodegc))                &
          .and. (.not. cumulus(ind_i(i),ind_j(i)))                             &
          .and. (r_theta_levels(ind_i(i),ind_j(i),k)                           &
                -r_theta_levels(ind_i(i),ind_j(i),0))                          &
                > zlcl_mixed(ind_i(i),ind_j(i))                                &
          .and.  rht  >   rhts(ind_i(i),ind_j(i),k)                            &
          .and. rht  >   (rhcrit(irhi,irhj,k)+rhcrit_tol)  ) then
          rh0 = rhcrit(irhi,irhj,k)
        end if

        ! Do we need to force some initiation from total liquid cloud cover?
        ! If so, set RH0 to equal the critical relative humidity

        if ( cfl(ind_i(i),ind_j(i),k) == 1.0                                   &
        .and. rht < rhts(ind_i(i),ind_j(i),k)                                  &
            .and. rht < (2.0-rhcrit(irhi,irhj,k)                               &
            - rhcrit_tol) ) then
          rh0 = rhcrit(irhi,irhj,k)
        end if

      else if ( i_pc2_init_logic == pc2init_logic_smooth ) then
        ! Smooth initiation logic; its possible we might want to initiate
        ! at least a little bit of extra liquid water at any point where
        ! we expect the Smith scheme to diagnose non-zero qcl,
        ! i.e. wherever total-water RH exceeds RHcrit:

        if ( rht > rhcrit(irhi,irhj,k) ) then
          rh0 = rhcrit(irhi,irhj,k)
        end if

      end if ! ( i_pc2_init_logic )

      ! ----------------------------------------------------------------------
      ! 4. Calculate the cloud fraction to initiate to
      ! ----------------------------------------------------------------------

      if (rh0 > 0.0) then
        npti=npti+1

        ! Is the liquid cloud cover coming down from one or up from zero?
        ! If it is coming down then reverse the value of RHT and work with
        ! saturation deficit taking the place of liquid water content.

        if (cfl(ind_i(i),ind_j(i),k)  >   0.5) then
          rht=2.0-rht
        end if

        ! Calculate the new liquid cloud fraction. Start by calculating QN then
        ! use the Smith scheme relationships to convert this to a cloud fraction

        ! when using the TKE based RHcrit parametrization, force the scheme
        ! to always use a symmetric triangular PDF, otherwise use the original
        ! PC2 method

        qn=(rht-1.0)/(1.0-rh0)
        if (qn <= -1.0) then
          c_1 = 0.0
        else if (qn < 0.0 .and. i_rhcpt == rhcpt_tke_based) then
          c_1 = 0.5*(1.0+qn)**2
        else if (qn < 0.0 .and. i_rhcpt /= rhcpt_tke_based) then
          c_1 = 0.5*(1.0+qn)**(pdf_power+1.0)
        else if (qn < 1.0 .and. i_rhcpt == rhcpt_tke_based) then
          c_1 = 1.0 - 0.5*(1.0-qn)**2
        else if (qn < 1.0 .and. i_rhcpt /= rhcpt_tke_based) then
          c_1 = 1.0 - 0.5*(1.0-qn)**(pdf_power+1.0)
        else
          c_1 = 1.0
        end if

        ! Reverse the new cloud fraction back to its correct value if the cloud
        ! fraction is being decreased from a high value.

        if (cfl(ind_i(i),ind_j(i),k) > 0.5) then
          c_1=1.0-c_1
        end if

        ! Calculate change in total cloud fraction. This depends upon the sign
        ! of the change of liquid cloud fraction. Change in ice cloud fraction
        ! is zero.

        deltacl_c(npti) = c_1 - cfl(ind_i(i),ind_j(i),k)
        deltacf_c(npti) = 0.0

        ! ----------------------------------------------------------------------
        ! 5. Calculate the liquid water content to initiate to.
        ! ----------------------------------------------------------------------

        ! We need to iterate to obtain the liquid content because the
        ! rate of change of gradient of the saturation specific humidity is
        ! specified as a function of the temperature, not the liquid
        ! temperature. We will only iterate over the necessary points where
        ! the initiation works out that a change in cloud fraction is required.

        ! Gather variables

        ni              (npti) = ind_i(i)
        nj              (npti) = ind_j(i)
        t_c             (npti) = t(ind_i(i),ind_j(i),k)
        qn_c            (npti) = qn
        rh0_c           (npti) = rh0
        qsl_tl_c        (npti) = qsl_tl
        cf_c            (npti) = cf(ind_i(i),ind_j(i),k)
        cfl_c           (npti) = cfl(ind_i(i),ind_j(i),k)
        cff_c           (npti) = cff(ind_i(i),ind_j(i),k)
        qcl_c           (npti) = qcl(ind_i(i),ind_j(i),k)
        q_c             (npti) = q(ind_i(i),ind_j(i),k)

        ! End if for RH0 gt 0
      end if

    end do !i
  end if !npt > 0

  ! Calculate change in total cloud fraction.

  if (npti > 0) then
    call pc2_total_cf(                                                         &
          npti,cfl_c,cff_c,deltacl_c,deltacf_c,cf_c)
  end if

  ! Iterations_do1:
  do l = 1, init_iterations

    ! Now loop over the required points

    ! Points_do1
    do i = 1, npti

      ! Calculate saturated specific humidity with respect to the temperature
      ! (not the liquid temperature).
      if ( l_mixing_ratio ) then
        call qsat_wat_mix(qsl_t_c,t_c(i),p_theta_levels(ni(i),nj(i),k))
      else
        call qsat_wat(qsl_t_c,t_c(i),p_theta_levels(ni(i),nj(i),k))
      end if

      alpha=repsilon*lc*qsl_t_c/(r*t_c(i)**2)
      al=1.0/(1.0+lcrcp*alpha)
      bs=al*(1.0-rh0_c(i))*qsl_tl_c(i)

      ! when using the TKE based RHcrit parametrization, force the scheme
      ! to always use a symmetric triangular PDF, otherwise use the original
      ! PC2 method

      if (qn_c(i) <= -1.0) then
        l_bs = 0.0
      else if (qn_c(i) < 0.0 .and. i_rhcpt == rhcpt_tke_based) then
        l_bs = 0.5/3.0*(1.0+qn_c(i))**3
      else if (qn_c(i) < 0.0 .and. i_rhcpt /= rhcpt_tke_based) then
        l_bs = 0.5/(pdf_power+2.0)*(1.0+qn_c(i))**(pdf_power+2.0)
      else if (qn_c(i) < 1.0 .and. i_rhcpt == rhcpt_tke_based) then
        l_bs = (qn_c(i)+0.5/3.0                                                &
               *(1.0-qn_c(i))**3)
      else if (qn_c(i) < 1.0 .and. i_rhcpt /= rhcpt_tke_based) then
        l_bs = (qn_c(i)+0.5/(pdf_power+2.0)                                    &
               *(1.0-qn_c(i))**(pdf_power+2.0))
      else
        l_bs = qn_c(i)
      end if

      l_out = l_bs * bs

      ! Are we working with saturation deficit instead of liquid water?
      ! If so, convert back to be the correct way around.

      if (cfl_c(i) > 0.5) then
        q_out = qsl_t_c - l_out / al
        l_out = qcl_c(i) + q_c(i) - q_out
      end if


      ! Calculate amount of condensation. To ensure convergence move
      ! slowly towards the calculated liquid water content L_OUT

      descent_factor = al
      deltal         = descent_factor*(l_out-qcl_c(i))
      q_c(i)         = q_c(i)   - deltal
      qcl_c(i)       = qcl_c(i) + deltal
      t_c(i)         = t_c(i)   + deltal*lcrcp

    end do ! Points_do1

  end do ! Iterations_do1

  ! Now scatter back values which have been changed

  if ( i_pc2_init_logic == pc2init_logic_smooth ) then
    ! Smooth initiation logic; only use updated values
    ! where the diagnosed qcl exceeds the prognosed qcl
    do i = 1, npti
      if ( qcl_c(i) > qcl(ni(i),nj(i),k) ) then

        ! Calculate saturated specific humidity with respect to the temperature
        ! (not the liquid temperature).
        if ( l_mixing_ratio ) then
          call qsat_wat_mix(qsl_t_c,t_c(i),p_theta_levels(ni(i),nj(i),k))
        else
          call qsat_wat(qsl_t_c,t_c(i),p_theta_levels(ni(i),nj(i),k))
        end if

        ! Calculate Qc (corresponds to the qcl we would have with
        ! no sub-grid moisture variability)

        ! qc = al ( q + qcl - qsl(Tl) )
        ! sd = al ( qsl(T) - q )
        ! qc + sd = al ( qcl + qsl(T) - qsl(Tl) )
        !         = al ( qcl + qsl(T) - qsl(T) + alpha lcrcp qcl )
        !         = al qcl ( 1 + alpha lcrcp )
        !         = qcl
        ! => sd = qcl - qc

        alpha=repsilon*lc*qsl_t_c/(r*t_c(i)**2)
        al=1.0/(1.0+lcrcp*alpha)
        qc = al * ( q_c(i) + qcl_c(i) - qsl_tl_c(i) )

        if ( qc < 0.0 ) then
          ! If qc<0 (total-water subsaturation)
          ! Find fraction of final qcl that was just created by initiation
          frac_init = ( qcl_c(i) - qcl(ni(i),nj(i),k) ) /  qcl_c(i)
        else
          ! If qc>0 (total-water supersaturation)
          ! Find fraction of final sd = qcl-qc that was created by initiation
          frac_init = ( qcl_c(i) - qcl(ni(i),nj(i),k) )                        &
                    / ( qcl_c(i) - min( qc, qcl(ni(i),nj(i),k) ) )
          ! (safety check avoids getting frac > 1 if qcl < qc, i.e. sd < 0
          !  which shouldn't really be possible!)
        end if

        ! Use the updated q, qcl and T at these points
        q  (ni(i),nj(i),k) = q_c(i)
        qcl(ni(i),nj(i),k) = qcl_c(i)
        t  (ni(i),nj(i),k) = t_c(i)

        ! New cloud fraction is weighted towards the value from the diagnostic
        ! scheme, in proportion to the fraction of qcl or sd created
        cf(ni(i),nj(i),k)  =      frac_init  * cf_c(i)                         &
                           + (1.0-frac_init) * cf(ni(i),nj(i),k)
        cfl(ni(i),nj(i),k) =      frac_init  * ( cfl_c(i) + deltacl_c(i) )     &
                           + (1.0-frac_init) * cfl(ni(i),nj(i),k)

      end if
    end do
  else  ! ( i_pc2_init_logic )
    ! Other initiation logic options use all updated fields in the list
    do i = 1, npti

      q  (ni(i),nj(i),k) = q_c(i)
      qcl(ni(i),nj(i),k) = qcl_c(i)
      t  (ni(i),nj(i),k) = t_c(i)
      cf (ni(i),nj(i),k) = cf_c(i)
      cfl(ni(i),nj(i),k) = cfl_c(i) + deltacl_c(i)

    end do !i
  end if  ! ( i_pc2_init_logic )

end do ! Levels_do1:
!$OMP end PARALLEL do

! End of the subroutine

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine pc2_initiate
end module pc2_initiate_mod
