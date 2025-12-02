! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Ensure conservation of energy and optionally water

module cor_engy_6a_mod

use um_types, only: real_umphys

implicit none

! ------------------------------------------------------------------------------
!
! Description:
! Adjust the potential temperature increments and humidity increments
! to ensure the conservation of energy and water.
! By default the energy and water is corrected in a manner consistent
! with the global energy correction (see UMDP 85)
!
! It is possible to apply the correction in a manner consistent with
! the assumptions made in the convection scheme - this is useful for
! scheme development but is not recommended for general use.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

character(len=*), parameter, private :: ModuleName='COR_ENGY_6A_MOD'

contains

subroutine cor_engy_6a(npnts, nconv, nlev, n_wtrac, index1, timestep,          &
                       p_layer_boundaries, exner_layer_centres,                &
                       exner_rho, r_theta, r_rho, rho_dry,                     &
                       r2rho, rho_theta, rho_dry_theta,                        &
                       dr_across_rh,                                           &
                       dubydt, dvbydt, dqclbydt, dqcfbydt,                     &
                       rain, snow, q, qcl, qcf, u, v,                          &
                       !In/Out
                       dthbydt, dqbydt, wtrac_e)

use cv_derived_constants_mod, only: ra2
use gen_phys_inputs_mod, only: l_mr_physics
use planet_constants_mod, only: g, cp
use water_constants_mod, only: lc, lf
use cv_param_mod, only: method_en_rho, method_en_mx_rho, method_en_qx_p
use cv_dependent_switch_mod, only: cor_method
use science_fixes_mod, only: l_improve_cv_cons
use wtrac_conv_mod, only: l_wtrac_conv, conv_e_wtrac_type
use cor_engy_wtrac_mod, only: cor_engy_wtrac
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none
!
! ------------------------------------------------------------------------------
! Subroutine arguments

integer,intent(in) :: npnts         ! Full vector length
integer,intent(in) :: nconv         ! Number of convecting points
integer,intent(in) :: nlev          ! Number of model levels for calculations
integer,intent(in) :: n_wtrac       ! Number of water tracers

integer,intent(in) :: index1(npnts) ! index of points with convection

real(kind=real_umphys),intent(in) :: timestep         ! Timestep
real(kind=real_umphys),intent(in) :: p_layer_boundaries(npnts,0:nlev)
                                                      ! Pressure at
                                                      ! layer boundary (Pa)
real(kind=real_umphys),intent(in) :: exner_layer_centres(npnts,0:nlev)
                                                      ! Exner function
                                                      ! at layer centre
real(kind=real_umphys),intent(in) :: exner_rho(npnts,nlev)
                                            ! Exner on rho levels
real(kind=real_umphys),intent(in) :: r_theta(npnts,0:nlev)
                                            ! radius of theta levels (m)
real(kind=real_umphys),intent(in) :: r_rho(npnts,nlev)
                                            ! radius of rho levels (m)
real(kind=real_umphys),intent(in) :: r2rho(npnts,nlev)
                                            ! radius**2 wet density for
                                            ! rho lev (kg/m)
real(kind=real_umphys),intent(in) :: rho_dry(npnts,nlev)
                                            ! dry density on rho lev (kg/m3)
real(kind=real_umphys),intent(in) :: rho_theta(npnts,nlev)
                                            ! wet density on theta lev (kg/m3)
real(kind=real_umphys),intent(in) :: rho_dry_theta(npnts,nlev)
                                            ! dry density on theta lev (kg/m3)
real(kind=real_umphys),intent(in) :: dr_across_rh(npnts,nlev)
                                            ! thickness of rho levels (m)

real(kind=real_umphys),intent(in) :: dubydt(npnts,nlev)
                                            ! Increment to u-wind on rho levels
                                            ! (m/s/s)
real(kind=real_umphys),intent(in) :: dvbydt(npnts,nlev)
                                            ! Increment to u-wind on rho levels
                                            ! (m/s/s)
real(kind=real_umphys),intent(in) :: dqclbydt(npnts,nlev)
                                            ! Increment to specific liquid water
                                            ! (kg/kg/s)
real(kind=real_umphys),intent(in) :: dqcfbydt(npnts,nlev)
                                            ! Increment to specific frozen water
                                            ! (kg/kg/s)

real(kind=real_umphys),intent(in) :: rain(npnts)
                                            ! rain at surface (kg/m**2/s)
real(kind=real_umphys),intent(in) :: snow(npnts)
                                            ! Snow at surface (kg/m**2/s)

! Specific quantities required to convert the specific increments into
! mixing ratio increments.
real(kind=real_umphys),intent(in) :: q(npnts,nlev)
                                            ! Specific humidity on theta levels
                                            ! in kg/kg
real(kind=real_umphys),intent(in) :: qcl(npnts,nlev)
                                            ! Specific liquid condensate on
                                            ! theta levels in kg/kg
real(kind=real_umphys),intent(in) :: qcf(npnts,nlev)
                                            ! Specific frozen condensate on
                                            ! theta levels in kg/kg
real(kind=real_umphys),intent(in) :: u(npnts,nlev)            ! U wind (m/s)
real(kind=real_umphys),intent(in) :: v(npnts,nlev)            ! V wind (m/s)

!----------------------------------------------------------------------
! Variables which are input but which are also updated in this routine
!----------------------------------------------------------------------

real(kind=real_umphys),intent(in out) :: dthbydt(npnts,nlev)
                                             ! Increment to P. temperature
                                            ! due to convection (K/s)
                                            ! In:uncorrected
                                            ! Out:Energy Corrected
real(kind=real_umphys),intent(in out) :: dqbydt(npnts,nlev)
                                             ! Increment to specific
                                            ! water vapour (kg/kg/s)
                                            ! In:uncorrected
                                            ! Out:Mass corrected
type(conv_e_wtrac_type), intent(in out) :: wtrac_e(n_wtrac)
                                            ! Structure containing water
                                            ! tracer fields

!---------------------------------------------------------------------
! Local variables
!---------------------------------------------------------------------

integer :: i,j,k                  ! loop counters

real(kind=real_umphys) :: factor  ! factor used to convert from specific
                                  ! to mixing ratios
real(kind=real_umphys) :: sumdqxbydt
                                  ! Sum of increments to specific water
                                  ! quantities on theta levels (kg/kg/s)
real(kind=real_umphys) :: mv(nconv,nlev)
                                  ! mixing ratio of water vapour
                                  ! on theta levels in kg/kg
real(kind=real_umphys) :: mcl(nconv,nlev)
                                  ! mixing ratio of liquid condensate on
                                  ! theta levels in kg/kg
real(kind=real_umphys) :: mcf(nconv,nlev)
                                  ! mixing ratio of frozen condensate on
                                  ! theta levels in kg/kg
real(kind=real_umphys) :: mv_rho(nconv,nlev)
                                  ! mixing ratio of water vapour
                                  ! on theta levels in kg/kg
real(kind=real_umphys) :: mcl_rho(nconv,nlev)
                                  ! mixing ratio of liquid condensate on
                                  ! rho levels in kg/kg
real(kind=real_umphys) :: mcf_rho(nconv,nlev)
                                  ! mixing ratio of frozen condensate on
                                  ! rho levels in kg/kg
real(kind=real_umphys) :: dmvbydt(nconv,nlev)
                                  ! Increment to mixing ratio of water vapour
                                  ! on theta levels (kg/kg/s)
real(kind=real_umphys) :: dmclbydt(nconv,nlev)
                                  ! Increment to mixing ratio of liquid water
                                  ! on theta levels (kg/kg/s)
real(kind=real_umphys) :: dmvbydt_rho(nconv,nlev)
                                  ! Increment to mixing ratio of water vapour
                                  ! on rho levels (kg/kg/s)
real(kind=real_umphys) :: dmclbydt_rho(nconv,nlev)
                                  ! Increment to mixing ratio of liquid water
                                  ! on rho levels (kg/kg/s)
real(kind=real_umphys) :: dqbydt_rho(nconv,nlev)
                                  ! Increment to specific water vapour
                                  ! on rho levels (kg/kg/s)
real(kind=real_umphys) :: dqclbydt_rho(nconv,nlev)
                                  ! Increment to specific liquid water
                                  ! on rho levels (kg/kg/s)
real(kind=real_umphys) :: dqcfbydt_rho(nconv,nlev)
                                  ! Increment to specific frozen water
                                  ! on rho levels (kg/kg/s)
real(kind=real_umphys) :: dthbydt_rho(nconv,nlev)
                                  ! Increment to potential temperature
                                  ! on rho levels (kg/kg/s)
real(kind=real_umphys) :: dEKbydt(nconv)
                                  ! Vertically integrated rate of change of
                                  ! kinetic energy (W/m2)
real(kind=real_umphys) :: dEIbydt(nconv)
                                  ! Vertically integrated rate of change of
                                  ! internal energy (W/m2)
real(kind=real_umphys) :: dEMbydt(nconv)
                                  ! Vertically integrated rate of change of
                                  ! moist energy (W/m2)
real(kind=real_umphys) :: flux_in(nconv)
                                  ! Total energy flux into the column (W/m2)
real(kind=real_umphys) :: dEbydt(nconv)
                                  ! Vertically integrated rate of change of
                                  ! total energy (W/m2)
real(kind=real_umphys) :: dTCWbydt(nconv)
                                  ! Vertically integrated rate of change of
                                  ! total water (kg/m2/s)
real(kind=real_umphys) :: dTCVbydt(nconv)
                                  ! Vertically integrated rate of change of
                                  ! water vapour (kg/m2/s)
real(kind=real_umphys) :: corr(nconv)
                                  ! Energy correction required for the column
                                  ! (W/m2)
real(kind=real_umphys) :: TCWcorr(nconv)
                                  ! Total water correction required for the
                                  ! column (kg/m2/s)
real(kind=real_umphys) :: r2_rhod_dr
                                  ! r**2 rho_dry * dr on rho levels (kg)
real(kind=real_umphys) :: r2_rhow_dr
                                  ! r**2 rho_wet * dr on rho levels (kg)
real(kind=real_umphys) :: weight1 ! Weighting factor used in interpolation
                                  ! from theta to rho levels
real(kind=real_umphys) :: weight2 ! Weighting factor used in interpolation
                                  ! from theta to rho levels
real(kind=real_umphys) :: weight3 ! Weighting factor used in interpolation
                                  ! from theta to rho levels
real(kind=real_umphys) :: ww1     ! Weighting factor used in interpolation
                                  ! from theta to rho levels
real(kind=real_umphys) :: ww2     ! Weighting factor used in interpolation
                                  ! from theta to rho levels
real(kind=real_umphys) :: deltap  ! Pressure difference between layer boundaries
                                  ! Only used for pressure based correction.

real(kind=real_umphys) :: dtscale(nconv)
                       ! scaling factor used to correct the temperature
                       ! increments
real(kind=real_umphys) :: dmvscale(nconv)
                       ! scaling factor used to correct the water increments


integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='COR_ENGY_6A'

!----------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)



!-------------------------------------------------------------------------
! WATER CORRECTION
!-------------------------------------------------------------------------

if (cor_method == method_en_mx_rho) then
  !-------------------------------------------------------------------------
  ! Interpolate water increments to the rho grid for use in water correction
  !-------------------------------------------------------------------------
  ! k=1
  do j = 1, nconv
    ! arrays defined on all points, npnts, need to be indexed with i
    ! arrays defined on convecting points only, nconv, need to be indexed with j
    i                   = index1(j)
    ! assume bottom rho level value equal to bottom theta level value
    dqbydt_rho(j,1)     = dqbydt(i,1)
    dqclbydt_rho(j,1)   = dqclbydt(i,1)
    dqcfbydt_rho(j,1)   = dqcfbydt(i,1)
  end do

  do k=2,nlev
    do j = 1, nconv
      i = index1(j)
      weight1             = r_theta(i,k) - r_rho(i,k)
      weight2             = r_rho(i,k)   - r_theta(i,k-1)
      weight3             = r_theta(i,k) - r_theta(i,k-1)
      ww1                 = weight1/weight3
      ww2                 = weight2/weight3

      dqbydt_rho(j,k)     = ww2 * dqbydt(i,k)     + ww1 * dqbydt(i,k-1)
      dqclbydt_rho(j,k)   = ww2 * dqclbydt(i,k)   + ww1 * dqclbydt(i,k-1)
      dqcfbydt_rho(j,k)   = ww2 * dqcfbydt(i,k)   + ww1 * dqcfbydt(i,k-1)
    end do
  end do
end if !(cor_method == method_en_mx_rho)

if (cor_method == method_en_mx_rho .or. cor_method == method_en_qx_p) then
  !-------------------------------------------------------------------------
  ! Water correction: Initialise column integrals
  !-------------------------------------------------------------------------
  do j=1,nconv
    dmvscale(j) = 0.0
    dTCVbydt(j) = 0.0
    dTCWbydt(j) = 0.0
  end do
end if !(cor_method == method_en_mx_rho .or. cor_method == method_en_qx_p)


if (cor_method == method_en_mx_rho) then
  !-------------------------------------------------------------------------
  ! Water correction on rho levels: Calculate column integrals
  !-------------------------------------------------------------------------
  do k=1,nlev
    do j=1,nconv
      i           = index1(j)
      r2_rhow_dr  = r2rho(i,k)*dr_across_rh(i,k)

      dTCVbydt(j) = dTCVbydt(j) + ra2 * r2_rhow_dr * dqbydt_rho(j,k)

      dTCWbydt(j) = dTCWbydt(j) + ra2 * r2_rhow_dr *                           &
                    (dqbydt_rho(j,k) + dqclbydt_rho(j,k) + dqcfbydt_rho(j,k))

    end do
  end do

else if (cor_method == method_en_qx_p) then
  !-------------------------------------------------------------------------
  ! Water correction on pressure levels: Calculate column integrals
  !-------------------------------------------------------------------------
  do k=1,nlev
    do j=1,nconv
      i           = index1(j)
      deltap      = p_layer_boundaries(i,k-1) - p_layer_boundaries(i,k)

      dTCVbydt(j) = dTCVbydt(j) + deltap/g * dqbydt(i,k)

      dTCWbydt(j) = dTCWbydt(j) + deltap/g *                                   &
                    (dqbydt(i,k) + dqclbydt(i,k) + dqcfbydt(i,k))
    end do
  end do

end if !cor_method

if (cor_method == method_en_mx_rho .or. cor_method == method_en_qx_p) then
  !-------------------------------------------------------------------------
  ! Water correction: calculate scaling
  !-------------------------------------------------------------------------
  do j=1,nconv
    i = index1(j)

    TCWcorr(j)  = - rain(i) - snow(i) - dTCWbydt(j)

    if (abs(dTCVbydt(j)) > 1.0e-6*abs(TCWcorr(j))) then
      dmvscale(j)  = 1.0 + TCWcorr(j)/dTCVbydt(j)
    else
      dmvscale(j)  = 1.0
    end if

    if (dmvscale(j) > 1.25 .or. dmvscale(j) < 0.75) then
      ! If the scaling is too large then do not bother.
      dmvscale(j) = 1.0
    end if

  end do

  !-------------------------------------------------------------------------
  ! Apply correction to water tracers
  ! (Do this before correction applied to water as need original dqbydt)
  !-------------------------------------------------------------------------
  if (l_wtrac_conv) then
    call cor_engy_wtrac(npnts, nconv, nlev, n_wtrac, index1, timestep, q,      &
                        dqbydt, dmvscale, wtrac_e)
  end if

  !-------------------------------------------------------------------------
  ! Water correction: Apply correction
  !-------------------------------------------------------------------------
  do k=1,nlev
    do j=1,nconv
      i           = index1(j)
      dqbydt(i,k) = dqbydt(i,k)*dmvscale(j)
    end do  ! nconv
  end do ! nlev

end if !cor_method

!-------------------------------------------------------------------------
! ENERGY CORRECTION
!-------------------------------------------------------------------------

if (cor_method == method_en_mx_rho .or. cor_method == method_en_rho) then
  !-------------------------------------------------------------------------
  ! Energy correction on rho levels: Convert the specific quantities and
  ! specific water increments to mass mixing ratio and mass mixing ratio
  ! increments for use in the energy correction
  !-------------------------------------------------------------------------
  if (l_mr_physics .or. l_improve_cv_cons) then
    !------------------------------------------------------------------------
    ! Convert using m = q * rho_wet/rho_dry, where these are start of
    ! timestep densities.  Note this is then deliberately ignoring the
    ! change in rho_wet arising from changes in the specific moisture
    ! variables, but is consistent with global moisture and energy correction
    !------------------------------------------------------------------------
    do k = 1, nlev
      do j = 1, nconv
        i             = index1(j)
        factor        = rho_theta(i,k)/rho_dry_theta(i,k)

        mv(j,k)       = factor*q(i,k)
        mcl(j,k)      = factor*qcl(i,k)
        mcf(j,k)      = factor*qcf(i,k)
        dmvbydt(j,k)  = factor*dqbydt(i,k)
        dmclbydt(j,k) = factor*dqclbydt(i,k)
      end do
    end do
  else
    !-----------------------------------------------------------------------
    ! Here we effectively convert using latest rho_wet and so are including
    ! the change in rho_wet arising from changes in specific moisture
    ! variables.  Although strictly correct, this is inconsistent with
    ! global moisture and energy correction
    !-----------------------------------------------------------------------
    do k = 1, nlev
      do j = 1, nconv
        i             = index1(j)
        factor        = 1.0/(1.0-q(i,k)-qcl(i,k)-qcf(i,k))
        sumdqxbydt    = dqbydt(i,k)+dqclbydt(i,k)+dqcfbydt(i,k)

        mv(j,k)       = factor*q(i,k)
        mcl(j,k)      = factor*qcl(i,k)
        mcf(j,k)      = factor*qcf(i,k)
        dmvbydt(j,k)  = factor*(dqbydt(i,k)   + q(i,k)  * factor * sumdqxbydt)
        dmclbydt(j,k) = factor*(dqclbydt(i,k) + qcl(i,k)* factor * sumdqxbydt)
      end do
    end do
  end if
  !-------------------------------------------------------------------------
  ! Energy correction on rho levels: Interpolate mixing ratios,
  ! mixing ratio increments and theta increment to the rho grid
  ! to be used in the energy correction
  !-------------------------------------------------------------------------

  ! k=1
  do j = 1, nconv
    ! arrays defined on all points, npnts, need to be indexed with i
    ! arrays defined on convecting points only, nconv, need to be indexed with j
    i = index1(j)
    ! assume bottom rho level value equal to bottom theta level value
    mv_rho(j,1)         = mv(j,1)
    mcl_rho(j,1)        = mcl(j,1)
    mcf_rho(j,1)        = mcf(j,1)
    dthbydt_rho(j,1)    = dthbydt(i,1)
    dmvbydt_rho(j,1)    = dmvbydt(j,1)
    dmclbydt_rho(j,1)   = dmclbydt(j,1)
  end do

  do k=2,nlev
    do j = 1, nconv
      i                   = index1(j)
      weight1             = r_theta(i,k) - r_rho(i,k)
      weight2             = r_rho(i,k)   - r_theta(i,k-1)
      weight3             = r_theta(i,k) - r_theta(i,k-1)
      ww1                 = weight1/weight3
      ww2                 = weight2/weight3

      mv_rho(j,k)         = ww2 * mv(j,k)         + ww1 * mv(j,k-1)
      mcl_rho(j,k)        = ww2 * mcl(j,k)        + ww1 * mcl(j,k-1)
      mcf_rho(j,k)        = ww2 * mcf(j,k)        + ww1 * mcf(j,k-1)
      dthbydt_rho(j,k)    = ww2 * dthbydt(i,k)    + ww1 * dthbydt(i,k-1)
      dmvbydt_rho(j,k)    = ww2 * dmvbydt(j,k)    + ww1 * dmvbydt(j,k-1)
      dmclbydt_rho(j,k)   = ww2 * dmclbydt(j,k)   + ww1 * dmclbydt(j,k-1)
    end do
  end do

end if  !(cor_method == method_en_mx_rho .or. cor_method == method_en_rho)

!-------------------------------------------------------------------------
! Energy correction: Initialise column integrals
!-------------------------------------------------------------------------
do j=1,nconv
  dtscale(j)  = 0.0
  dEKbydt(j)  = 0.0
  dEIbydt(j)  = 0.0
  dEMbydt(j)  = 0.0
end do

if (cor_method == method_en_mx_rho .or. cor_method == method_en_rho) then
  !-------------------------------------------------------------------------
  ! Energy correction on rho levels: Calculate column integrals
  !-------------------------------------------------------------------------
  if (l_improve_cv_cons) then
    do k=1,nlev
      do j=1,nconv
        i           = index1(j)
        !----------------------------------------------------------------------
        ! This expression uses start of timestep density combined
        ! with the latest mixing ratios, consistent with the dynamics
        !----------------------------------------------------------------------
        r2_rhod_dr  = r_rho(i,k)*r_rho(i,k)*rho_dry(i,k)*dr_across_rh(i,k)

        dEKbydt(j)  = dEKbydt(j) + ra2 * r2_rhod_dr *                          &
                      ( u(i,k)*dubydt(i,k) + v(i,k)*dvbydt(i,k) )

        dEMbydt(j)  = dEMbydt(j) + ra2 * r2_rhod_dr *                          &
                      ( (lc+lf)*dmvbydt_rho(j,k) + lf*dmclbydt_rho(j,k) )

        dEIbydt(j)  = dEIbydt(j) + ra2 * r2_rhod_dr *                          &
                      cp*exner_rho(i,k)*dthbydt_rho(j,k)

      end do
    end do
  else
    do k=1,nlev
      do j=1,nconv
        i           = index1(j)
        !-----------------------------------------------------------------------
        ! This expression uses an approximation to dry rho and therefore
        ! is not quite as accurate as it ideally should be.
        ! It inconsistently uses start of timestep wet density combined with
        ! the latest mixing ratios. It also ignores prognostic rain. A more
        ! accurate expression should be considered in future.
        !-----------------------------------------------------------------------

        r2_rhod_dr  = r2rho(i,k)*dr_across_rh(i,k) /                           &
                      (1.0+mv_rho(j,k)+mcl_rho(j,k)+mcf_rho(j,k))

        dEKbydt(j)  = dEKbydt(j) + ra2 * r2_rhod_dr *                          &
                      ( u(i,k)*dubydt(i,k) + v(i,k)*dvbydt(i,k) )

        dEMbydt(j)  = dEMbydt(j) + ra2 * r2_rhod_dr *                          &
                      ( (lc+lf)*dmvbydt_rho(j,k) + lf*dmclbydt_rho(j,k) )

        dEIbydt(j)  = dEIbydt(j) + ra2 * r2_rhod_dr *                          &
                      cp*exner_rho(i,k)*dthbydt_rho(j,k)

      end do
    end do
  end if
else if (cor_method == method_en_qx_p) then
  !-------------------------------------------------------------------------
  !Energy correction on pressure levels: Calculate column integrals
  !-------------------------------------------------------------------------
  do k=1,nlev
    do j=1,nconv
      i           = index1(j)
      deltap      = p_layer_boundaries(i,k-1) - p_layer_boundaries(i,k)

      dEKbydt(j)  = dEKbydt(j) + deltap/g*                                     &
                    ( u(i,k)*dubydt(i,k) + v(i,k)*dvbydt(i,k) )

      dEMbydt(j)  = dEMbydt(j) + deltap/g*                                     &
                    ( (lc+lf)*dqbydt(i,k) + lf*dqclbydt(i,k) )

      dEIbydt(j)  = dEIbydt(j) + deltap/g*                                     &
                    cp*exner_layer_centres(i,k)*dthbydt(i,k)

    end do
  end do
end if !cor_method

!-------------------------------------------------------------------------
! Energy correction: calculate scaling
!-------------------------------------------------------------------------
do j=1,nconv
  i           = index1(j)
  flux_in(j)  = -lf*rain(i)
  dEbydt(j)   = dEKbydt(j) + dEMbydt(j) + dEIbydt(j)

  corr(j)     = flux_in(j) - dEbydt(j)

  if (abs(dEIbydt(j)) > 1.0e-6*abs(corr(j))) then
    dtscale(j)  = 1.0 + corr(j)/dEIbydt(j)
  else
    dtscale(j)  = 1.0
  end if

  if (dtscale(j) > 1.25 .or. dtscale(j) < 0.75) then
    !   If the scaling is too large then do not bother.
    dtscale(j) = 1.0
  end if

end do

!-------------------------------------------------------------------------
! Energy correction: Apply correction
!-------------------------------------------------------------------------
do k=1,nlev
  do j=1,nconv
    i             = index1(j)
    dthbydt(i,k)  = dthbydt(i,k)*dtscale(j)
  end do  ! nconv
end do ! nlev


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine cor_engy_6a

end module cor_engy_6a_mod
