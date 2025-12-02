! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module set_region_cond_fields_mod

implicit none

contains


! Subroutine to set the local condensed water species
! mixing ratios and cloud fractions in a specified subregion
! of the grid-box
subroutine set_region_cond_fields( n_points, nc, index_ic,                     &
                                   n_points_fields, i_region,                  &
                                   q_cond_loc, cf_ice, frac_r,                 &
                                   fields_reg )

use comorph_constants_mod, only: real_cvprec, zero, one, n_cond_species,       &
                                 cond_params, i_sg_homog, i_sg_frac_liq,       &
                                 i_sg_frac_prec,                               &
                                 l_cv_cloudfrac, sqrt_min_float
use fields_type_mod, only: n_fields, i_qc_first,                               &
                           i_cf_liq, i_cf_ice, i_cf_bulk
use subregion_mod, only: n_regions, i_dry, i_liq, i_mph, i_icr

implicit none

! Number of points in q_cond_loc
integer, intent(in) :: n_points

! Number of points in compression list
integer, intent(in) :: nc

! Indices of compression list points
integer, intent(in) :: index_ic(nc)

! Dimension of the output fields super-array
! (may differ from nc if array of larger size being reused)
integer, intent(in) :: n_points_fields

! Specified region indicator
integer, intent(in) :: i_region

! Local in-region values of condensed water mixing ratios
real(kind=real_cvprec), intent(in) :: q_cond_loc                               &
                                    ( n_points, n_cond_species )

! Total ice-cloud fraction
real(kind=real_cvprec), intent(in) :: cf_ice(n_points)

! Fraction of each region
real(kind=real_cvprec), intent(in) :: frac_r(n_points,n_regions)

! Local in-region fields compressed onto points where a
! single specified region has nonzero fraction.
real(kind=real_cvprec), intent(in out) :: fields_reg                           &
                                   ( n_points_fields, n_fields )

! Fraction of the icr region containing ice cloud
real(kind=real_cvprec) :: frac

! Loop counters
integer :: ic, ic2, i_cond


! Copy combination of water species and cloud-fractions,
! depending on which region we're in
select case(i_region)
case (i_dry)

  ! Clear-sky region; all condensed water species are zero
  ! except for any which are assumed homogeneous across
  ! the whole grid-box
  do i_cond = 1, n_cond_species
    if ( cond_params(i_cond)%pt % i_sg == i_sg_homog ) then
      do ic2 = 1, nc
        fields_reg(ic2,i_qc_first-1+i_cond) = q_cond_loc(index_ic(ic2),i_cond)
      end do
    else
      do ic2 = 1, nc
        fields_reg(ic2,i_qc_first-1+i_cond) = zero
      end do
    end if
  end do
  ! All cloud fractions are zero
  if ( l_cv_cloudfrac ) then
    do ic2 = 1, nc
      fields_reg(ic2,i_cf_liq) = zero
      fields_reg(ic2,i_cf_ice) = zero
      fields_reg(ic2,i_cf_bulk) = zero
    end do
  end if

case (i_liq)

  ! Liquid cloud only region; only species in the liquid-cloud and precip
  ! fractions are non-zero
  do i_cond = 1, n_cond_species
    if ( cond_params(i_cond)%pt % i_sg == i_sg_frac_liq .or.                   &
         cond_params(i_cond)%pt % i_sg == i_sg_frac_prec ) then
      do ic2 = 1, nc
        fields_reg(ic2,i_qc_first-1+i_cond) = q_cond_loc(index_ic(ic2),i_cond)
      end do
    else
      do ic2 = 1, nc
        fields_reg(ic2,i_qc_first-1+i_cond) = zero
      end do
    end if
  end do
  ! Full liquid cloud cover but no ice
  if ( l_cv_cloudfrac ) then
    do ic2 = 1, nc
      fields_reg(ic2,i_cf_liq) = one
      fields_reg(ic2,i_cf_ice) = zero
      fields_reg(ic2,i_cf_bulk) = one
    end do
  end if

case (i_mph)

  ! Mixed-phase cloud region; all species present
  do i_cond = 1, n_cond_species
    do ic2 = 1, nc
      fields_reg(ic2,i_qc_first-1+i_cond)                                      &
        = q_cond_loc(index_ic(ic2),i_cond)
    end do
  end do
  ! Full coverage from both liquid and ice cloud
  if ( l_cv_cloudfrac ) then
    do ic2 = 1, nc
      fields_reg(ic2,i_cf_liq) = one
      fields_reg(ic2,i_cf_ice) = one
      fields_reg(ic2,i_cf_bulk) = one
    end do
  end if

case (i_icr)

  ! Ice and rain region; all species except liquid cloud
  do i_cond = 1, n_cond_species
    if ( cond_params(i_cond)%pt % i_sg == i_sg_frac_liq ) then
      do ic2 = 1, nc
        fields_reg(ic2,i_qc_first-1+i_cond) = zero
      end do
    else
      do ic2 = 1, nc
        fields_reg(ic2,i_qc_first-1+i_cond) = q_cond_loc(index_ic(ic2),i_cond)
      end do
    end if
  end do
  ! No liquid cloud, but maybe some ice fraction...
  if ( l_cv_cloudfrac ) then
    do ic2 = 1, nc
      ic = index_ic(ic2)
      ! Find fraction of icr region that contains ice
      ! = (total ice frac - mixed-phase frac ) / icr frac
      frac = ( cf_ice(ic) - frac_r(ic,i_mph) )                                 &
           / max( frac_r(ic,i_icr), sqrt_min_float )
      fields_reg(ic2,i_cf_liq) = zero
      fields_reg(ic2,i_cf_ice) = frac
      fields_reg(ic2,i_cf_bulk) = frac
    end do
  end if

end select


return
end subroutine set_region_cond_fields

end module set_region_cond_fields_mod
