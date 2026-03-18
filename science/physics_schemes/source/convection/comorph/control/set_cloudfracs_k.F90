! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module set_cloudfracs_k_mod

implicit none

contains

! Subroutine to set the super-array containing environment cloud
! fractions at the current full model-level.  Either needs to copy them
! from the primary fields array (if they are treated as prognostics),
! or compress from separate diagnostic cloud arrays otherwise.
subroutine set_cloudfracs_k( max_points, n_fields_tot, k, cmpr,                &
                             fields_k, cloudfracs, cloudfracs_k )

use comorph_constants_mod, only: real_cvprec, zero,                            &
                                 l_cv_cloudfrac
use cmpr_type_mod, only: cmpr_type
use cloudfracs_type_mod, only: cloudfracs_type, n_cloudfracs,                  &
                               i_frac_liq, i_frac_ice, i_frac_bulk,            &
                               i_frac_precip
use fields_type_mod, only: i_cf_liq

use compress_mod, only: compress
use force_cloudfrac_consistency_mod, only: force_cloudfrac_consistency

implicit none

! Size of compression arrays
integer, intent(in) :: max_points

! Number of primary fields, accounting for whether doing tracer transport
integer, intent(in) :: n_fields_tot

! Index of current full model-level
integer, intent(in) :: k

! Structure storing compression indices
type(cmpr_type), intent(in) :: cmpr

! Environment primary fields at the current level (already compressed)
real(kind=real_cvprec), intent(in) :: fields_k(max_points,n_fields_tot)

! Structure containing cloud fractions, in the event that they
! are diagnostic rather than prognostic and so aren't included
! in the primary fields structure.
type(cloudfracs_type), intent(in) :: cloudfracs

! Environment cloud-fractions at the current level
real(kind=real_cvprec), intent(out) :: cloudfracs_k(max_points,n_cloudfracs)

! Array lower and upper bounds
integer :: lb(3), ub(3)

! Loop counters
integer :: ic, i_field, i_frac


if ( l_cv_cloudfrac ) then
  ! Cloud fractions included as primary fields

  ! Copy cloud fractions into their own arrays
  do i_frac = i_frac_liq, i_frac_bulk
    i_field = i_cf_liq-1 + i_frac
    do ic = 1, cmpr % n_points
      cloudfracs_k(ic,i_frac) = fields_k(ic,i_field)
    end do
  end do

else  ! ( l_cv_cloudfrac )
  ! Cloud fractions only passed in as diagnostic fields

  ! Compress cloud fractions from diagnostic input fields
  lb = lbound(cloudfracs % frac_liq)
  ub = ubound(cloudfracs % frac_liq)
  call compress( cmpr, lb(1:2), ub(1:2),                                       &
                 cloudfracs % frac_liq(:,:,k),                                 &
                 cloudfracs_k(:,i_frac_liq) )
  lb = lbound(cloudfracs % frac_ice)
  ub = ubound(cloudfracs % frac_ice)
  call compress( cmpr, lb(1:2), ub(1:2),                                       &
                 cloudfracs % frac_ice(:,:,k),                                 &
                 cloudfracs_k(:,i_frac_ice) )
  lb = lbound(cloudfracs % frac_bulk)
  ub = ubound(cloudfracs % frac_bulk)
  call compress( cmpr, lb(1:2), ub(1:2),                                       &
                 cloudfracs % frac_bulk(:,:,k),                                &
                 cloudfracs_k(:,i_frac_bulk) )

  do i_frac = i_frac_liq, i_frac_bulk
    ! Remove spurious negative values
    do ic = 1, cmpr % n_points
      cloudfracs_k(ic,i_frac)                                                  &
        = max( cloudfracs_k(ic,i_frac), zero )
    end do
  end do

  ! Rounding errors when converting the cloud-fractions to
  ! 32-bit in compress can cause them to become slightly
  ! inconsistent; correct if needed:
  call force_cloudfrac_consistency( cmpr % n_points, max_points,               &
                                    cloudfracs_k(:,i_frac_liq:i_frac_bulk) )

end if  ! ( l_cv_cloudfrac )

lb = lbound(cloudfracs % frac_precip)
ub = ubound(cloudfracs % frac_precip)
! Compress rain fractions
call compress( cmpr, lb(1:2), ub(1:2),                                         &
               cloudfracs % frac_precip(:,:,k),                                &
               cloudfracs_k(:,i_frac_precip) )
! Remove spurious negative values
do ic = 1, cmpr % n_points
  cloudfracs_k(ic,i_frac_precip)                                               &
    = max( cloudfracs_k(ic,i_frac_precip), zero )
end do


return
end subroutine set_cloudfracs_k

end module set_cloudfracs_k_mod
