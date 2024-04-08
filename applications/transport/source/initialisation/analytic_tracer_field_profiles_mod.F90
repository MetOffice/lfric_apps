!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Functions to compute a variety of analytic tracer field profiles.
!> @details Collection of functions to return the value of a field at a given
!!          point based upon a specified analytic formula to initialise the
!!          tracer fields for the transport miniapp.
module analytic_tracer_field_profiles_mod

use constants_mod,              only : r_def, pi
use log_mod,                    only : log_event,         &
                                       log_scratch_space, &
                                       LOG_LEVEL_ERROR
use coord_transform_mod,        only : xyz2llr, central_angle
use idealised_config_mod,       only : test_gaussian_hill,        &
                                       test_cosine_hill,          &
                                       test_cosine_bell,          &
                                       test_yz_cosine_hill,       &
                                       test_slotted_cylinder,     &
                                       test_constant_field,       &
                                       test_hadley_like_dcmip,    &
                                       test_cosine_stripe,        &
                                       test_cos_phi,              &
                                       test_cosine_bubble,        &
                                       test_eternal_fountain,     &
                                       test_div_free_reversible,  &
                                       test_curl_free_reversible, &
                                       test_rotational,           &
                                       test_translational,        &
                                       test_vertical_cylinder
use initial_tracer_field_config_mod,                                   &
                                only : r1, x1, y1, z1, r2, x2, y2, z2, &
                                       field_max, field_background
use base_mesh_config_mod,       only : geometry, &
                                       geometry_spherical
use planet_config_mod,          only : p_zero, Rd, kappa, scaled_radius
use extrusion_config_mod,       only : domain_top

implicit none

private

public :: analytic_tracer_field
public :: hadley_like_dcmip

contains


!> @brief Compute the tracer function from Kent et al. 2014 and
!!        Allen and Zerroukat 2016.
!> @details Equations below have been taken from Allen and Zerroukat, "A deep
!!          non-hydrostatic compressible atmospheric model on a Yin-Yang grid",
!!          JCP, 2016, equation (5.5), and Kent, Ullrich and Jablonowski,
!!          "Dynamical core model intercomparison project: Tracer transport test
!!          cases", QJRMS 2014, equation (42). Parameter values have been taken
!!          from the papers and are currently hard-wired.
!> @param[in] radius  Distance from the centre of the planet to the point of interest
!> @return    tracer  Value of the tracer field at this point
function hadley_like_dcmip(radius) result(tracer)
  implicit none
  real(kind=r_def), intent(in) :: radius
  real(kind=r_def)             :: tracer

  real(kind=r_def) :: z, z1, z2, z0

  z = radius - scaled_radius

  z1 = 2000.0_r_def
  z2 = 5000.0_r_def
  z0 = 0.5_r_def*(z1+z2)

  if ((z-z1)*(z2-z)>0) then
    tracer = 0.5_r_def*(1.0_r_def + cos(2.0_r_def*pi*(z-z0)/(z2-z1)))
  else
    tracer = 0.0_r_def
  end if

end function hadley_like_dcmip


!> @brief Compute an analytic tracer field.
!> @param[in] chi          Position in physical coordinates
!> @param[in] choice       Integer defining which specified formula to use
!> @param[in] domain_max_x Max. domain extent in x-direction.
!> @result tracer The result tracer field
function analytic_tracer_field(chi, choice, domain_max_x) result(tracer)

  implicit none
  real(kind=r_def), intent(in) :: chi(3)
  integer,          intent(in) :: choice
  real(kind=r_def), intent(in) :: domain_max_x
  real(kind=r_def)             :: tracer

  real(kind=r_def), parameter  :: XC = 0.0_r_def
  real(kind=r_def)             :: long, lat, radius
  real(kind=r_def)             :: l1, l2
  real(kind=r_def)             :: d1, d2
  real(kind=r_def)             :: h1, h2
  real(kind=r_def)             :: bubble_dist, bubble_zc
  real(kind=r_def)             :: bubble_radius, bubble_width, bubble_height
  real(kind=r_def)             :: slot_width, slot_length

  if ( geometry == geometry_spherical ) then
    call xyz2llr(chi(1),chi(2),chi(3),long,lat,radius)
    call central_angle(long,lat,x1,y1,l1)
    call central_angle(long,lat,x2,y2,l2)
  else
    long = chi(1)
    lat  = chi(2)
    l1 = sqrt((long-x1)**2 + (lat-y1)**2)
    l2 = sqrt((long-x2)**2 + (lat-y2)**2)
  end if

  select case( choice )

  case( test_gaussian_hill )
    h1 = field_max*exp( -(l1/r1)**2 )
    h2 = field_max*exp( -(l2/r2)**2 )
    tracer = field_background + h1 + h2

  case( test_cosine_hill )
    if ( l1 < r1 ) then
      h1 = field_background + (field_max/2.0_r_def)*(1.0_r_def+cos((l1/r1)*PI))
    else
      h1 = field_background
    end if
    if (l2 < r2) then
      h2 = field_background + (field_max/2.0_r_def)*(1.0_r_def+cos((l2/r2)*PI))
    else
      h2 = field_background
    end if
    tracer = h1+h2

  case( test_cosine_bell )
    bubble_height = domain_top/12.0_r_def
    d1 = min( 1.0_r_def, ( l1 / 0.5_r_def )**2 + &
         ( ( radius - scaled_radius - domain_top/2.0_r_def ) / bubble_height )**2 )
    d2 = min( 1.0_r_def, ( l2 / 0.5_r_def )**2 + &
         ( ( radius - scaled_radius - domain_top/2.0_r_def ) / bubble_height )**2 )
    tracer = field_background + ( (field_max - field_background) / 2.0_r_def ) * &
              ( ( 1.0_r_def + cos( pi*d1 ) ) + ( 1.0_r_def + cos( pi*d2 ) ) )

  case( test_yz_cosine_hill )

    l1 = sqrt((chi(2)-y1)**2 + (chi(3)-z1)**2)
    l2 = sqrt((chi(2)-y2)**2 + (chi(3)-z2)**2)

    if ( l1 < r1 ) then
      h1 = field_background + (field_max/2.0_r_def)*(1.0_r_def+cos((l1/r1)*PI))
    else
      h1 = field_background
    end if
    if (l2 < r2) then
      h2 = field_background + (field_max/2.0_r_def)*(1.0_r_def+cos((l2/r2)*PI))
    else
      h2 = field_background
    end if
    tracer = h1+h2

  case( test_slotted_cylinder )
    ! Cylinder 1
    if ( l1 < r1 ) then
      if (abs(long-x1) > r1/6.0_r_def) then
        h1 = field_max
      else
        if (lat < y1-r1*5.0_r_def/12.0_r_def) then
          h1 = field_max
        else
          h1 = field_background
        end if
      end if
    else
      h1 = field_background
    end if
    ! Cylinder 2
    if ( l2 < r2 ) then
      if (abs(long-x2) > r2/6.0_r_def) then
        h2 = field_max
      else
        if (lat > y2+r2*5.0_r_def/12.0_r_def) then
          h2 = field_max
        else
          h2 = field_background
        end if
      end if
    else
      h2 = field_background
    end if
    tracer = h1 + h2

  case( test_constant_field )
    tracer = field_background

  case( test_cosine_stripe )
    l1 = sqrt((long-x1)**2)
    if ( l1 < r1 ) then
      tracer = field_background + (field_max/2.0_r_def)*(1.0_r_def+cos((l1/r1)*PI))
    else
      tracer = field_background
    end if

  case( test_hadley_like_dcmip )
    tracer = hadley_like_dcmip(radius)

  case( test_cos_phi )
    tracer = field_max*cos(lat)**4

  case( test_cosine_bubble )
    l1 = sqrt( ((chi(1) - x1)/r1)**2 + ((chi(3) - y1)/r2)**2 )
    if ( l1 < 1.0_r_def ) then
      tracer = field_background + field_max*cos(0.5_r_def*l1*PI)**2
    else
      tracer = field_background
    end if

  case( test_eternal_fountain )
    bubble_width = 0.4_r_def * domain_max_x
    bubble_height = 0.1_r_def * domain_top

    if ( ( (chi(1) + bubble_width / 2.0_r_def) &
            * (bubble_width / 2.0_r_def - chi(1)) > 0.0_r_def ) &
      .and. ( chi(3) * (bubble_height - chi(3)) > 0.0_r_def ) ) then
      tracer = field_max
    else
      tracer = field_background
    end if

  case ( test_rotational, test_curl_free_reversible, &
         test_translational, test_div_free_reversible )
    bubble_zc = domain_top / 4.0_r_def
    bubble_width = domain_max_x / 5.0_r_def
    bubble_height = domain_top / 10.0_r_def
    bubble_radius = bubble_height / 2.0_r_def

    ! Elliptical distance from centre of bubble
    bubble_dist = bubble_radius &
      * sqrt( ((chi(1) - XC) / (bubble_width / 2.0_r_def) ) ** 2.0_r_def &
            + ((chi(3) - bubble_zc) / (bubble_height / 2.0_r_def)) ** 2.0_r_def)

    tracer = field_background + (field_max - field_background) &
                * exp(-(bubble_dist / bubble_radius)**2.0_r_def)

  case( test_vertical_cylinder )
    bubble_zc = domain_top / 4.0_r_def
    bubble_width = domain_max_x / 2.0_r_def
    bubble_height = domain_top / 4.0_r_def
    bubble_radius = bubble_height / 2.0_r_def

    ! Elliptical distance from centre of bubble
    bubble_dist = bubble_radius &
      * sqrt( ((chi(1) - XC) / (bubble_width / 2.0_r_def) ) ** 2.0_r_def &
            + ((chi(3) - bubble_zc) / (bubble_height / 2.0_r_def)) ** 2.0_r_def)

    slot_width = bubble_width / 12.0_r_def
    slot_length = 17.0_r_def * bubble_height / 24.0_r_def

    if ( bubble_dist < bubble_radius ) then
      if ( abs(chi(1) - XC) > slot_width / 2.0_r_def ) then
        tracer = field_max
      else
        if ( chi(3) < (bubble_zc + bubble_height / 2.0_r_def - slot_length) ) then
          tracer = field_max
        else
          tracer = field_background
        end if
      end if
    else
      tracer = field_background
    end if

  case default
    write( log_scratch_space, '(A)' )  'Invalid tracer profile choice, stopping'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )

  end select

end function analytic_tracer_field

end module analytic_tracer_field_profiles_mod
