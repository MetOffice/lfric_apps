! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Calculate visibility (in metres) from temperature, specific humidity,
! cloud liquid water or ThetaL and qt and Murk aerosol (if present).

! Description:
!   Process fields of temperature, specific humidity, cloud liquid
!   water or ThetaL and qt and Murk aerosol (if present) to give
!   visibility in metres.
!   Calculated at a single model level or level within surface layer
!   e.g. screen height (1.5m)

! Documentation:
!    Wright, B. J., 1997: Improvements to the Nimrod Visibility
!       Analysis/Forecast System. FR-Div. Tech. Rep., No. 217.
!    Wright, B. J., 1997: A New Visibility Analysis/Forecast System
!       for Nimrod. Met. Office FR Tech Rep., No. 222.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: atmos_service_visibility

! Code description:
!   This code is written to UMDP3 standards.

! Subroutine Interface:
module visbty_mod

implicit none

character(len=*), parameter, private :: ModuleName = 'VISBTY_MOD'
contains

subroutine visbty(                                                             &
             p_layer,t,q,qcl,                                                  &
                                                      !INPUT
            aerosol, prob, rhcrit, l_murk,                                     &
                                                      !INPUT
            p_field,                                                           &
                                                      !INPUT
            visibility)
                                                      !OUTPUT
! Modules
use um_types, only: real_umphys
use water_constants_mod, only: rho_water
use visbty_constants_mod, only: m0, aero0, n0, power, radius0, a0, b0,         &
                                onethird, fourthirds, visfactor,               &
                                recipvisair
use conversions_mod, only: pi

use qsat_mod, only: qsat_wat

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

integer, intent(in) ::                                                         &
        p_field                   ! in NO. points in field.
real(kind=real_umphys), intent(in) ::                                          &
  p_layer(p_field),                                                            &
        t(p_field),                                                            &
                                      ! in Temperature
        q(p_field),                                                            &
                                      ! in Qt
        qcl(p_field),                                                          &
                                      ! in cloud water array.
        aerosol(p_field),                                                      &
                                      ! in Aerosol mixing ratio(ug/kg)
        rhcrit                  ! in Critical RH (determines
                                    !    width of distribiution)
logical, intent(in) ::                                                         &
   l_murk                        ! in : Aerosol present
real(kind=real_umphys), intent(in) ::                                          &
        prob                     ! in Probability level ( e.g 0.5
                                 !    corresponds to median).
!---------------------------------------------------------------------
! output variables----------------------------------------------------
!---------------------------------------------------------------------
real(kind=real_umphys), intent(out) ::                                         &
      visibility(p_field)         ! out visibility array.
! --------------------------------------------------------------------
! Local variables:----------------------------------------------------
!---------------------------------------------------------------------
real(kind=real_umphys) ::                                                      &
       qt(p_field),                                                            &
                                    ! total of cloud water and vapour
       qs(p_field)              ! saturation vapour pressure
!---------------------------------------------------------------------
! Local parameter variables
real(kind=real_umphys), parameter :: rhmax = 0.99
                      ! Maximum value of relative humidity
                      !  which is allowed to feed into the
                      !  calculation of the 'fog' droplet radius
real(kind=real_umphys), parameter :: rhmin = 0.001
                      ! Minimum value of relative humidity which
                      !  is allowed to feed into the calculation
                      !   of the 'fog' droplet radius
real(kind=real_umphys), parameter :: weight = 0.75
                      ! Weighting on new value for iterative
                      !  solution of droplet radius
real(kind=real_umphys), parameter :: delta_radius_star = 0.001
                       ! Convergence required for iterative
                      !  solution of droplet radius
real(kind=real_umphys), parameter :: qt_limit = 0.0001
                      ! Smallest Qt value allowed
real(kind=real_umphys), parameter :: radius_star_min = 1.0
real(kind=real_umphys), parameter :: radius_star_max = 1000.0
real(kind=real_umphys), parameter :: radius_star_factor = 4.0

integer, parameter :: niterations = 20
                      !  Maximum number of iteration used to
                      !   estimate the water droplet radius

! Local workspace variables
real(kind=real_umphys) :: n             ! Local number density

integer ::                                                                     &
 point,                                                                        &
                       !  Loop variable for points
 iteration        !  Loop variable iterations used to estimate
                       !   the water droplet radius

real(kind=real_umphys) ::                                                      &
  m_over_m0,                                                                   &
                       !  Ratio of  aerosol mass mixing ratio and
                       !   the standard aerosol mass mixing ratio
  recipvis,                                                                    &
                       !  Recipirical of the visibility
  radius_dry,                                                                  &
                       !  Radius of dry aerosol particle (m)
  radius,                                                                      &
                       !  Radius of fog droplets (m)
  radius_star1,                                                                &
                       !  Previous estimate of water droplet radius
                       !   divided by the dry radius
  radius_star2,                                                                &
                       !  Current best estimate of water droplet
                       !   radius divided by the dry radius
  radius_act,                                                                  &
                       !  Activation droplet radius
  radius_star_act,                                                             &
                       !  Activation droplet rad divided by dry rad
  a,                                                                           &
                       !  A0 divided by the dry radius
  rh_lim,                                                                      &
                       !  Limited RH value (fractional)
  fn,                                                                          &
                       !  Value of droplet radius function
  deriv,                                                                       &
                       !  Derivative of droplet radius function
  radius_star_diff,                                                            &
                       !  Absolute value of radius_star1 minus
                       !    radius_star2
  rhterm,                                                                      &
                       !  Relative humidity term in function to be
                       !   minimised to find the droplet radius
  qlterm,                                                                      &
                       !  Liquid water term in function to be
                       !   minimised to find the droplet radius
  rhderiv,                                                                     &
                       !  Derivative of relative humidity term
  qlderiv,                                                                     &
                       !  Derivative of liquid water term
  bs,                                                                          &
                       !  Width of distribution in total water
                       !   mixing ratio space (kg/kg)
  qt_mod,                                                                      &
                       !  Modified total water value based on the
                       !   probability of the value occurring
                       !   assuming a triangular distriubtion
                       !   of width bs.
  qt_mod_factor    !  Factor to multiply bs to modify qt

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='VISBTY'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     Create factor to multiply bs by to modify qt
if ( prob  ==  0.5 ) then
  qt_mod_factor = 0.0
else if ( prob  >=  0.0 .and. prob  <   0.5 ) then
  qt_mod_factor = ( 1.0 - sqrt( 2.0 * prob ) )
else if ( prob  >=  0.5 .and. prob  <=  1.0 ) then
  qt_mod_factor = - ( 1.0 - sqrt( 2.0 * (1.0-prob) ) )
end if
! ----------------------------------------------------------------------
! For the new cloud and precipitation scheme only use the liquid content
! 1. Calculate total of water vapour and liquid water contents, P and
! limit aerosol
! ----------------------------------------------------------------------
do point = 1, p_field
  qt(point) = q(point)+qcl(point)
end do

call qsat_wat (qs,t,p_layer,p_field)

  !-------------------------------------------------------------------
  !  2. Calculate the ratio of the aerosol mass mixing ratio to the
  !     standard mass mixing ratio, m_over_m0, and the aerosol number
  !     density, N, the dry radius, radius_dry:
  !                       p
  !                   (m )
  !            r = r0 (--)
  !             d     (m0)
  !
  !
  !     And the activation radius:
  !
  !                              1/2
  !                   (       3 )
  !                   ( 3 B0 r  )
  !            r    = ( ------d-)
  !             act   (   A0    )
  !
  !     and A (A0 divided by the dry radius).
  ! N.B. AEROSOL is in ug/kg, m in kg/kg
  ! If not available, use 10 ug/kg
  !-------------------------------------------------------------------

! If l_murk=false, this calculation doesn't depend on points,
! so can be done outside of loop over points.
! Must be kept identical with code inside loop.
if (.not. l_murk) then
  m_over_m0 = max(10.0/m0*1.0e-9, 0.0001)

  n = n0 * m_over_m0**(1.0-3*power)

  radius_dry = radius0 * (m_over_m0)**power
  a = a0 / radius_dry

  radius_act = sqrt( (3 * b0 * radius_dry**3) / a0 )
  radius_star_act =  radius_act/radius_dry
end if

do point = 1, p_field

  if (l_murk) then
    ! Ensure that assumed aerosol conc. is at least Aero0:
    m_over_m0 = max(aerosol(point)/m0*1.0e-9,                                  &
                             aero0/m0*1.0e-9, 0.0001)

    n = n0 * m_over_m0**(1.0-3*power)

    radius_dry = radius0 * (m_over_m0)**power
    a = a0 / radius_dry

    radius_act = sqrt( (3 * b0 * radius_dry**3) / a0 )
    radius_star_act =  radius_act/radius_dry
  end if

  !-------------------------------------------------------------
  !  3. Calculate the width of the total water
  !     distribution and a modified value of total water, based
  !     on a probability.
  !-------------------------------------------------------------

  bs = (1.0-rhcrit) * qs(point)

  qt_mod = max( qt_limit, qt(point)+ qt_mod_factor* bs)

  !====================================================================
  !  4.  Use Newton-Raphson to iteratively improve on a first-guess
  !      droplet radius, using the droplet growth equation and the
  !      geometric relation between liquid water and droplet radius.
  !====================================================================
  !  4.1 Calculate a first guess relative humidity, qt/qs, but limit it
  !      to be in the range 0.001 -> 0.999.
  !      From this calculate a first-guess normalised radius using a
  !      simplified version of the droplet growth equation:
  !
  !                               1/3
  !                 (       B0   )
  !            r  = ( 1 - ------ )
  !             *   (     ln(RH) )
  !
  !----------------------------------------------------------------------

  rh_lim = min( max( qt_mod/qs(point), rhmin ) , rhmax )
  radius_star2 = (1.0-b0/log(rh_lim))**onethird

  !----------------------------------------------------------------------
  !  4.2 Initialise the iteration counter, the normalised radius
  !      difference, and the updated normalised radius value.
  !----------------------------------------------------------------------

  iteration = 0
  radius_star_diff = 1.0
  radius_star1 = radius_star2

  do while ( iteration  <   niterations .and.                                  &
               radius_star_diff  >   delta_radius_star )

    !----------------------------------------------------------------------
    !  4.3 Update the iteration counter and the normalised radius value.
    !----------------------------------------------------------------------

    iteration = iteration + 1
    radius_star1 = weight * radius_star2                                       &
                 + ( 1.0 - weight ) * radius_star1

    !----------------------------------------------------------------------
    !  4.4 Calculate the relative humidity term:
    !
    !                       ( A        B0   )
    !           RHterm = exp( --  -  ------ )
    !                       ( r       3     )
    !                       (  *     r  - 1 )
    !                       (         *     )
    !
    !       and its derivative with respect to the normalised radius:
    !
    !                     (                 2    )
    !                     (   A       3 B0 r     )
    !           RHderiv = ( - --  +  -------*- 2 ) * RHterm
    !                     (    2     (  3     )  )
    !                     (   r      ( r  - 1 )  )
    !                     (    *     (  *     )  )
    !
    !----------------------------------------------------------------------

    if ( radius_star1  <   radius_star_act ) then
      rhterm  = exp( a/radius_star1                                            &
                     - b0/(radius_star1**3-1.0) )* qs(point)
      rhderiv = - rhterm * ( -a/(radius_star1**2)                              &
                + (3.0*b0*radius_star1**2)                                     &
                /(radius_star1**3-1.0)**2 )
    else
      rhterm  = exp( a/radius_star_act                                         &
                     - b0/(radius_star_act**3-1.0) ) * qs(point)
      rhderiv = 0.0
    end if

    !----------------------------------------------------------------------
    !  4.5 Calculate the liquid water mixing ratio term:
    !
    !
    !                    4             3 (  3     )
    !           qLterm = - Pi rho_w N r  ( r  - 1 )
    !                    3             d (  *     )
    !
    !       and its derivative with respect to the normalised radius:
    !
    !                                   3  2
    !           qLderiv = 4 Pi rho_w N r  r
    !                                   d  *
    !
    !----------------------------------------------------------------------

    qlterm  = n * fourthirds * pi * rho_water * radius_dry**3                  &
                * ( radius_star1**3 - 1.0 )
    qlderiv  = - n * 4.0 * pi * rho_water                                      &
                * radius_dry**3 * radius_star1**2

    !----------------------------------------------------------------------
    !  4.6 Calculate the function, Fn, and its derivative, Deriv, and
    !      an improved estimate of the normalised radius,
    !      using Newton Raphson:
    !
    !           Fn = qt - RHterm - qLterm
    !
    !           Deriv = RHderiv + qLderiv
    !
    !                           Fn
    !           r      = r  -  -----
    !            * new    *    Deriv
    !
    !      The new estimate of the normalised radius is limited lie between
    !      prescribed maximum and minimum values and within a factor of the
    !      previous value to ensure that the soultion does not diverge.
    !----------------------------------------------------------------------

    fn    = qt_mod - rhterm - qlterm
    deriv = rhderiv + qlderiv

    radius_star2 = radius_star1 - fn/deriv

    if ( radius_star2  <   radius_star_min )                                   &
        radius_star2 = radius_star_min
    if ( radius_star2  >   radius_star_max )                                   &
        radius_star2 = radius_star_max
    if ( radius_star2  >   radius_star_factor * radius_star1 )                 &
        radius_star2 = radius_star_factor * radius_star1
    if ( radius_star2  <   radius_star1 / radius_star_factor )                 &
        radius_star2 = radius_star1 / radius_star_factor

    !---------------------------------------------------------------------
    !  4.7 Calculate difference between the old and the new values of the
    !      normalised radius.
    !---------------------------------------------------------------------

    radius_star_diff = abs( radius_star1 - radius_star2 )

  end do

  !---------------------------------------------------------------------
  !  5.  Calculate the radius from the final normalised radius.
  !---------------------------------------------------------------------

  radius = radius_star2 * radius_dry

  !---------------------------------------------------------------------
  !  6. Calculate the visibility, Vis, using the equation:
  !
  !                  ln(liminal contrast)
  !            Vis = -------------2------
  !                      Beta0 N r
  !
  !     (An extra term RecipVisAir is included in the recipical of
  !      visibility to limit visibilities to 100km in clean air).
  !---------------------------------------------------------------------

  recipvis = (n * radius**2) / visfactor + recipvisair
  visibility(point) = 1/recipvis

end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine visbty
end module visbty_mod
