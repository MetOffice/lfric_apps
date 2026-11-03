! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module set_par_cloudfrac_mod

implicit none

contains


! Subroutine to set the in-parcel cloud fraction fields.
! Just sets them based on the presence of non-zero liquid and
! ice cloud mixing-ratios.
subroutine set_par_cloudfrac( n_points, n_points_cond, n_points_cf,            &
                              q_cond_super, cloudfracs )

use comorph_constants_mod, only: real_cvprec, zero, one,                       &
                                 n_cond_species, cond_params,                  &
                                 i_sg_frac_liq, i_sg_frac_ice,                 &
                                 i_par_cloudfrac,                              &
                                 i_par_cloudfrac_hom, i_par_cloudfrac_mph,     &
                                 overlap_power
use cloudfracs_type_mod, only: i_frac_liq, i_frac_ice, i_frac_bulk


implicit none

! Number of points
integer, intent(in) :: n_points
! Number of points in the condensate and cloud-fractions super-arrays
! (maybe larger than needed here, to save having to reallocate)
integer, intent(in) :: n_points_cond
integer, intent(in) :: n_points_cf

! Super-array containing all condensed water mixing ratios
real(kind=real_cvprec), intent(in) :: q_cond_super                             &
                                      ( n_points_cond, n_cond_species )

! Cloud-fractions super-array
real(kind=real_cvprec), intent(out) :: cloudfracs                              &
                                       ( n_points_cf, i_frac_liq:i_frac_bulk )

! Sum of liquid and ice cloud mixing-ratio
real(kind=real_cvprec) :: qc

! Loop counters
integer :: ic, i_cond

do ic = 1, n_points
  ! Initialise cloud fractions to zero
  cloudfracs(ic,i_frac_liq) = zero
  cloudfracs(ic,i_frac_ice) = zero
  cloudfracs(ic,i_frac_bulk) = zero
end do

select case(i_par_cloudfrac)
case (i_par_cloudfrac_hom)
  ! Condensate species assumed to be homogeneously
  ! distributed across the parcel, so that the fraction
  ! associated with each is just 0 if none present,
  ! 1 if any is present.

  do i_cond = 1, n_cond_species
    ! Set liquid cloud fraction to 1 if any liquid cloud present
    if ( cond_params(i_cond)%pt % i_sg == i_sg_frac_liq ) then
      do ic = 1, n_points
        if ( q_cond_super(ic,i_cond) > zero ) then
          cloudfracs(ic,i_frac_liq) = one
          cloudfracs(ic,i_frac_bulk) = one
        end if
      end do
    end if
    ! Set ice cloud fraction to 1 if any ice cloud present
    if ( cond_params(i_cond)%pt % i_sg == i_sg_frac_ice ) then
      do ic = 1, n_points
        if ( q_cond_super(ic,i_cond) > zero ) then
          cloudfracs(ic,i_frac_ice) = one
          cloudfracs(ic,i_frac_bulk) = one
        end if
      end do
    end if
  end do

case (i_par_cloudfrac_mph)
  ! Alternative option; when liquid and ice both present,
  ! they may not be fully overlapped.

  ! Accumulate total liquid and ice cloud mixing-ratios
  ! (temporarily storing in the cloud-fraction fields)
  do i_cond = 1, n_cond_species
    if ( cond_params(i_cond)%pt % i_sg == i_sg_frac_liq ) then
      do ic = 1, n_points
        cloudfracs(ic,i_frac_liq) = cloudfracs(ic,i_frac_liq)                  &
                                  + q_cond_super(ic,i_cond)
      end do
    end if
    if ( cond_params(i_cond)%pt % i_sg == i_sg_frac_ice ) then
      do ic = 1, n_points
        cloudfracs(ic,i_frac_ice) = cloudfracs(ic,i_frac_ice)                  &
                                  + q_cond_super(ic,i_cond)
      end do
    end if
  end do

  do ic = 1, n_points
    qc = cloudfracs(ic,i_frac_liq) + cloudfracs(ic,i_frac_ice)
    ! If any condensed water
    if ( qc > zero ) then

      ! Set bulk cloud fraction to 1
      cloudfracs(ic,i_frac_bulk) = one

      ! Set liquid and ice cloud fractions to their respective
      ! fractions of total condensate
      cloudfracs(ic,i_frac_liq) = cloudfracs(ic,i_frac_liq) / qc
      cloudfracs(ic,i_frac_ice) = cloudfracs(ic,i_frac_ice) / qc

      ! This yields zero overlap between liquid and ice.
      ! Parameterise some overlap by raising both fractions
      ! to a power between 0 and 1
      cloudfracs(ic,i_frac_liq) = cloudfracs(ic,i_frac_liq)**overlap_power
      cloudfracs(ic,i_frac_ice) = cloudfracs(ic,i_frac_ice)**overlap_power

    end if
  end do

end select  ! CASE(i_par_cloudfrac)


return
end subroutine set_par_cloudfrac


end module set_par_cloudfrac_mod
