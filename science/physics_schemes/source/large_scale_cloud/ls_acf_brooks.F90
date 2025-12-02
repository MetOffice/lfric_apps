! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module ls_acf_brooks_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'LS_ACF_BROOKS_MOD'
contains

subroutine ls_acf_brooks (                                                     &
! trig arrays
  xx_cos_theta_latitude,                                                       &
! in data fields
  bulk_cloud_fraction, cloud_fraction_liquid,                                  &
  cloud_fraction_frozen,                                                       &
! in logical control
  cumulus,                                                                     &
! out data fields
  area_cloud_fraction)

!     Description:
!       Calculates area_cloud_fraction from bulk_cloud_fraction

!     Method:
!       The calculation is  based on the parametrisation in
!       Brooks 2005 equations 2-3.
!       (Brooks et al, July 2005, JAS vol 62 pp 2248-2260)
!       The initial parametrisation uses the values in equations
!       4,5,7 and 8 of the paper for ice and liquid cloud without
!       wind shear.  For mixed phase clouds the maximum of the two
!       area_cloud_fractions resulting will be used.
!       Only area_cloud_fraction will be updated.
!       Grid box size is needed to be known.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_cloud
!
!     Code Description:
!       FORTRAN 77 with extensions recommended in the Met. Office
!       F77 Standard.

use atm_fields_bounds_mod, only: tdims_s, tdims
use cderived_mod, only: delta_lambda, delta_phi
use level_heights_mod, only:  r_theta_levels
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Array Arguments with intent(in)
! Co-ordinate arrays:
real(kind=real_umphys), intent(in) ::                                          &
 xx_cos_theta_latitude (tdims_s%i_start:tdims_s%i_end,                         &
                        tdims_s%j_start:tdims_s%j_end)
!       Finite volume cos(lat)

! Data arrays:
real(kind=real_umphys), intent(in) ::                                          &
 bulk_cloud_fraction (  tdims%i_start:tdims%i_end,                             &
                        tdims%j_start:tdims%j_end,                             &
                                    1:tdims%k_end),                            &
!       Cloud fraction at processed levels (decimal fraction).
   cloud_fraction_liquid (tdims%i_start:tdims%i_end,                           &
                          tdims%j_start:tdims%j_end,                           &
                                      1:tdims%k_end),                          &
!       Liquid cloud fraction at processed levels (decimal fraction).
   cloud_fraction_frozen (tdims%i_start:tdims%i_end,                           &
                          tdims%j_start:tdims%j_end,                           &
                                      1:tdims%k_end)
!       Frozen cloud fraction at processed levels (decimal fraction).

logical, intent(in)::                                                          &
 cumulus(               tdims%i_start:tdims%i_end,                             &
                        tdims%j_start:tdims%j_end)

! Arguments with intent(out):
! Data arrays:
real(kind=real_umphys), intent(out) ::                                         &
 area_cloud_fraction (  tdims%i_start:tdims%i_end,                             &
                        tdims%j_start:tdims%j_end,                             &
                                    1:tdims%k_end)
!       Cloud fraction at processed levels (decimal fraction).

! Local Parameters:

! Parameters for liquid clouds, from Brooks 2005 equations 7 and 8
real(kind=real_umphys), parameter :: power_law_gradient_liquid =  0.1635 ! A
real(kind=real_umphys), parameter :: vert_fit_liquid           =  0.6694 ! alpha
real(kind=real_umphys), parameter :: horiz_fit_liquid          = -0.1882 ! beta

! Parameters for frozen clouds, from Brooks 2005 equations 4 and 5
real(kind=real_umphys), parameter :: power_law_gradient_frozen =  0.0880 ! A
real(kind=real_umphys), parameter :: vert_fit_frozen           =  0.7679 ! alpha
real(kind=real_umphys), parameter :: horiz_fit_frozen          = -0.2254 ! beta

! Local Scalars:
! Loop counters
integer ::                                                                     &
 i, j, k

real(kind=real_umphys) ::                                                      &
 symmetric_adjustment_liquid,                                                  &
!    function f in eqn 7 in Brooks 2005
   symmetric_adjustment_frozen,                                                &
!    function f in eqn 4 in Brooks 2005
   horiz_scale,                                                                &
!    horizontal scale size of the grid box (m)
   vert_scale
!    vertical scale size of the grid box (m)

!  Local Arrays:
real(kind=real_umphys) ::                                                      &
 acf_liquid (           tdims%i_start:tdims%i_end,                             &
                        tdims%j_start:tdims%j_end,                             &
                                    1:tdims%k_end),                            &
!    area cloud fraction based on liquid parameters
   acf_frozen (           tdims%i_start:tdims%i_end,                           &
                          tdims%j_start:tdims%j_end,                           &
                                      1:tdims%k_end)
!    area cloud fraction based on frozen parameters

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LS_ACF_BROOKS'

!-    End of header
! ----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ==Main Block==--------------------------------------------------------

! Initialise arrays and local variables to zero
do k = 1, tdims%k_end
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      area_cloud_fraction(i,j,k) = 0.0
      acf_liquid(i,j,k) = 0.0
      acf_frozen(i,j,k) = 0.0
    end do
  end do
end do
horiz_scale = 0.0
vert_scale = 0.0
symmetric_adjustment_liquid = 0.0
symmetric_adjustment_frozen = 0.0

do k = 1, tdims%k_end
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end

      ! Test if bulk_cloud_fraction is within bounds of possibility

      if ( bulk_cloud_fraction(i,j,k) <= 0.0 ) then

        area_cloud_fraction(i,j,k) = 0.0

      else if ( bulk_cloud_fraction(i,j,k) >= 1.0 ) then

        area_cloud_fraction(i,j,k) = 1.0

      else if (cumulus(i,j)) then
            ! This is a convective point so do not apply the
            ! area cloud representation
        area_cloud_fraction(i,j,k) =                                           &
          bulk_cloud_fraction(i,j,k)

      else

        ! Only calculate area_cloud_fraction if the bulk_cloud_fraction
        ! is between (not equal to) 0.0 and 1.0

        ! Calculate horizontal and vertical scales.
        ! The horizontal scale is taken as the square root of the
        ! area of the grid box.
        ! The vertical scale is taken as the difference in radius
        ! from the centre of the Earth between the upper and lower
        ! boundaries of the grid box.

        horiz_scale = sqrt (                                                   &
                          r_theta_levels(i,j,k)                                &
                          * r_theta_levels(i,j,k)                              &
                          * delta_lambda * delta_phi                           &
                          * xx_cos_theta_latitude(i,j) )
        if (k  ==  tdims%k_end) then
              ! Assume top layer thickness is the same as the
              ! thickness of the layer below
          vert_scale =  r_theta_levels(i,j,k)                                  &
                       - r_theta_levels(i,j,k-1)
        else
          vert_scale =  r_theta_levels(i,j,k+1)                                &
                       - r_theta_levels(i,j,k)
        end if  ! k eq model_levels

        ! Calculate the symmetric_adjustment (f).
        ! This parameter controls the extent to which the area cloud fraction
        ! is greater than the bulk cloud fraction.  If f = 0, they are equal.

        symmetric_adjustment_liquid =                                          &
                   power_law_gradient_liquid                                   &
                   * ( vert_scale ** vert_fit_liquid )                         &
                   * ( horiz_scale ** horiz_fit_liquid )
        symmetric_adjustment_frozen =                                          &
                   power_law_gradient_frozen                                   &
                   * ( vert_scale ** vert_fit_frozen )                         &
                   * ( horiz_scale ** horiz_fit_frozen )

        ! Calculate the area cloud fractions for liquid and frozen cloud
        ! Calculate the liquid and frozen fractions separately to
        ! allow for greatest flexibility in future choice of decisions
        ! regarding mixed phase cloud.

        acf_liquid(i,j,k) = 1.0/                                               &
            ( 1.0 + ( exp(-1.0*symmetric_adjustment_liquid)                    &
                     * ( 1.0/bulk_cloud_fraction(i,j,k) - 1.0) ) )
        acf_frozen(i,j,k) = 1.0/                                               &
            ( 1.0 + ( exp(-1.0*symmetric_adjustment_frozen)                    &
                     * ( 1.0/bulk_cloud_fraction(i,j,k) - 1.0) ) )

        ! Calculate the final area cloud fraction for each grid box
        ! Currently this is based on which there is more of, ice or liquid.

        if ( cloud_fraction_frozen(i,j,k) == 0.0 ) then
          if ( cloud_fraction_liquid(i,j,k) == 0.0 ) then

            ! If there is no liquid or frozen cloud, there should be no area cloud
            area_cloud_fraction(i,j,k) = 0.0

          else

            ! If there is no frozen cloud but there is liquid cloud,
            ! then the area cloud fraction is given by the liquid parametrisation
            ! 0 no cloud, 1 either, 2 liq, 3 ice'
            area_cloud_fraction(i,j,k) = acf_liquid(i,j,k)
          end if

        else ! cloud_fraction_frozen

          if ( cloud_fraction_liquid(i,j,k) == 0.0 ) then

            ! If there is frozen cloud but there is no liquid cloud,
            ! then the area cloud fraction is given by the frozen parametrisation
            area_cloud_fraction(i,j,k) = acf_frozen(i,j,k)

          else

            ! If there is frozen cloud and there is liquid cloud,
            ! then the area cloud fraction is given by the maximum of the two
            ! parametrisations
            area_cloud_fraction(i,j,k) =                                       &
               max( acf_liquid(i,j,k),acf_frozen(i,j,k) )

          end if

        end if ! cloud_fraction_frozen

      end if ! bulk_cloud_fraction between 0.0 and 1.0

    end do
  end do
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine ls_acf_brooks
! ======================================================================
end module ls_acf_brooks_mod
