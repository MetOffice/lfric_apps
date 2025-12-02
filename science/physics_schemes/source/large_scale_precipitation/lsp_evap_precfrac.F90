! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Prognostic precip fraction increment
!  accompanying evaporation of rain or graupel.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Subroutine Interface:
module lsp_evap_precfrac_mod

implicit none

character(len=*), parameter, private :: ModuleName='LSP_EVAP_PRECFRAC_MOD'

contains

subroutine lsp_evap_precfrac( points, npts, ix, l_qpr2,                        &
                              rainfrac, rain_clear, rain_ice,                  &
                              dqprec_clear, dqprec_ice,                        &
                              qpr, qpr2, dqpr, precfrac_k )

use um_types,            only: real_lsprec
use lsprec_mod,          only: small_number, zero, one
use mphys_inputs_mod,    only: i_update_precfrac,                              &
                               i_homog_areas, i_sg_correl

! Dr Hook Modules
use yomhook,             only: lhook, dr_hook
use parkind1,            only: jprb, jpim

implicit none

! Number of points in the arrays
integer, intent(in) :: points

! Compression indices of points where evaporation is occuring
integer, intent(in) :: npts
integer, intent(in) :: ix(points)

! Flag for whether to account for 2nd precip type
! (needed if both rain and graupel are localised in the precip fraction)
logical, intent(in) :: l_qpr2

! Start-of-timestep precip fraction
real(kind=real_lsprec) , intent(in):: rainfrac(points)

! Area-fractions of overlap of precip-fraction with clear-air and ice-cloud
real(kind=real_lsprec), intent(in) :: rain_clear(points)
real(kind=real_lsprec), intent(in) :: rain_ice(points)

! Increments to rain+graupel in the cloud-air and ice-only cloud regions
real(kind=real_lsprec), intent(in out) :: dqprec_clear(points)
real(kind=real_lsprec), intent(in out) :: dqprec_ice(points)

! Latest value of precip content, not yet updated with evaporation
real(kind=real_lsprec), intent(in) :: qpr(points)

! "Other" precip content (graupel when this call is for rain evaporation)
real(kind=real_lsprec), intent(in) :: qpr2(points)

! Evaporation increment to qpr
real(kind=real_lsprec), intent(in) :: dqpr(points)

! Prognostic precip fraction to be updated
real(kind=real_lsprec), intent(in out) :: precfrac_k(points)

! Fraction of increment occuring in the ice-only vs cloud-free region
real(kind=real_lsprec) :: frac

! Loop counters
integer :: i, c

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_EVAP_PRECFRAC'


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


if ( i_update_precfrac == i_homog_areas ) then
  ! Need to store contribution of rain evaporation to the
  ! total precip increment in the non-liquid-cloud partition.

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do c = 1, npts
    i = ix(c)
    ! Increment applies in both the rain_clear and rain_ice area partitions;
    ! divvy it up between the two in proportion to area
    frac = rain_ice(i) / max( rain_ice(i) + rain_clear(i), small_number )
    dqprec_clear(i) = dqprec_clear(i) + (one-frac) * dqpr(i)
    dqprec_ice(i)   = dqprec_ice(i)   -      frac  * dqpr(i)
    ! Note the term should be subtracted from dqprec_clear
    ! (consistent with dqprec_ice), not added; this is a bug,
    ! to be fixed in another ticket...
  end do

else if ( i_update_precfrac == i_sg_correl ) then
  ! Assume that the area fraction in which evaporation occurs declines
  ! in proportion to the sqrt of the precip mass it contains.

  ! qpr within evap fraction before = qpr / rainfrac
  ! qpr within evap fraction after  = qpr / rainfrac - dqpr/(frac rainfrac)
  ! ratio = ( qpr - dqpr/frac ) / qpr
  !       = ( qpr frac - dqpr ) / ( qpr frac )
  if ( l_qpr2 ) then
    ! 2nd precip species is included in the prognostic precip fraction;
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do c = 1, npts
      i = ix(c)
      ! Compute evaporating fraction (precip outside of the liquid-cloud)
      frac = max( min( ( rain_ice(i) + rain_clear(i) ) / rainfrac(i), one ),   &
                  small_number )
      ! Precip fraction (1-frac) is unchanged, while the frac bit
      ! is scaled down proportional to sqrt of precip mass:
      precfrac_k(i) = precfrac_k(i)                                            &
        * ( one-frac                                                           &
          +     frac * sqrt( ( max( qpr2(i), zero )*frac                       &
                             + max( qpr(i)*frac - dqpr(i), zero ) )            &
                           / ( ( max( qpr2(i), zero ) + qpr(i) )*frac )        &
                           ) )
      ! Due to limiting of dqpr, it should be impossible to have
      ! qpr(i)*frac - dqpr(i) < 0.0, but very small negative values
      ! actually do seem to happen due to floating point rounding errors.
    end do
  else  ! ( .not. l_qpr2 )
    ! Only qpr is included in the prognostic precip fraction;
    ! same as above but ignore qpr2.
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do c = 1, npts
      i = ix(c)
      frac = max( min( ( rain_ice(i) + rain_clear(i) ) / rainfrac(i), one ),   &
                  small_number )
      precfrac_k(i) = precfrac_k(i)                                            &
        * ( one-frac                                                           &
          +     frac * sqrt( max( qpr(i)*frac - dqpr(i), zero )                &
                           / ( qpr(i)*frac )                                   &
                           ) )
    end do
  end if  ! ( .not. l_qpr2 )

end if  ! ( i_update_precfrac )


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine lsp_evap_precfrac

end module lsp_evap_precfrac_mod
