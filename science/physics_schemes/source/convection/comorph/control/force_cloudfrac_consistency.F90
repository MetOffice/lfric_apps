! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module force_cloudfrac_consistency_mod

implicit none

contains

! Subroutine to force the bulk cloud-fraction in a compressed
! fields array to be consistent with the liquid and ice
! cloud-fractions.  We must have:
!
! max(cf_liq,cf_ice) <= cf_bulk <= cf_liq+cf_ice
!
! If this law is broken, then attempts to calculate the
! mixed-phase and liquid-only cloud-fractions can yield
! spurious negative values.
! Various conversions of variables in CoMorph could potentially
! cause cf_bulk to go very slightly outside the above allowable
! range (e.g. rounding errors when converting the input
! cloud-fractions from 64-bit to 32-bit).
! This routine adjusts cf_bulk to force it to be within the
! above range.
subroutine force_cloudfrac_consistency( n_points, n_points_super,              &
                                        cloudfracs )

use comorph_constants_mod, only: real_cvprec
use cloudfracs_type_mod, only: i_frac_liq, i_frac_ice,                         &
                               i_frac_bulk

implicit none

! Number of points in the compression list
integer, intent(in) :: n_points

! Number of points in the array (may be > n_points where reusing
! an array for differing size compression lists)
integer, intent(in) :: n_points_super

! Super-array containing the 3 cloud-fraction fields
real(kind=real_cvprec), intent(in out) :: cloudfracs                           &
                      ( n_points_super, i_frac_liq:i_frac_bulk )

! Loop counter
integer :: ic

do ic = 1, n_points
  cloudfracs(ic,i_frac_bulk)                                                   &
    = max( min( cloudfracs(ic,i_frac_bulk),                                    &
                cloudfracs(ic,i_frac_liq)                                      &
              + cloudfracs(ic,i_frac_ice) ),                                   &
           max( cloudfracs(ic,i_frac_liq),                                     &
                cloudfracs(ic,i_frac_ice) )                                    &
         )
end do

return
end subroutine force_cloudfrac_consistency

end module force_cloudfrac_consistency_mod
