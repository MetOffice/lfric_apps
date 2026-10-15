! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module core_combine_mod

implicit none

contains

! Subroutine to calculate combined parcel core properties of 2 different
! parcels; there are various different options for how to do this.
! This needs to be done when:
! a) Combining initiating parcels from multiple different sub-grid
!    regions on a given model-level.
! b) Combining initiating parcels from subsequent model-levels higher-up
!    into the rising plume.
subroutine core_combine( n_points_a, n_points_m, index_ic,                     &
                         n_points_super_a, n_points_super_m,                   &
                         i_field_first, i_field_last, l_down,                  &
                         core_a_fields, core_m_fields,                         &
                         edge_a_virt_temp, edge_m_virt_temp )

use comorph_constants_mod, only: real_cvprec, zero, one
use fields_type_mod, only: i_temperature, i_q_vap, i_qc_first, i_qc_last
use calc_virt_temp_mod, only: calc_virt_temp

implicit none

! Number of points in the parcel to be added "_a"
integer, intent(in) :: n_points_a

! Number of points in the existing parcel arrays "_m"
integer, intent(in) :: n_points_m

! Indices for referencing the "_m" arrays from the "_a" compression list
integer, intent(in) :: index_ic(n_points_a)

! Sizes of the fields super-arrays; maybe > n_points due to reusing arrays
integer, intent(in) :: n_points_super_a
integer, intent(in) :: n_points_super_m

! First and last primary fields to set
integer, intent(in) :: i_field_first
integer, intent(in) :: i_field_last

! Flag for downdraft versus updraft
logical, intent(in) :: l_down

! Parcel core fields to be combined:
! "_a": Input properties of the "added" parcel (one of the parcels to combine)
real(kind=real_cvprec), intent(in) :: core_a_fields                            &
                               ( n_points_super_a, i_field_first:i_field_last )
! "_m": IN:  properties of the other parcel to combine "_a" into
!       OUT: properties of the "merged" parcel after combining with "_a"
real(kind=real_cvprec), intent(in out) :: core_m_fields                        &
                               ( n_points_super_m, i_field_first:i_field_last )

! Parcel edge virtual temperatures to be combined
real(kind=real_cvprec), intent(in) :: edge_a_virt_temp(n_points_a)
real(kind=real_cvprec), intent(in out) :: edge_m_virt_temp(n_points_m)

! Weights for combining the core properties
real(kind=real_cvprec) :: weight_core_a(n_points_a)
real(kind=real_cvprec) :: weight_core_m(n_points_a)

! Virtual temperatures of the core properties
real(kind=real_cvprec) :: core_a_virt_temp(n_points_a)
real(kind=real_cvprec) :: core_m_virt_temp(n_points_m)

! Loop counters
integer :: ic, ic2, i_field


! Compute the core virtual temperature of the two parcels
! NOTE: this calculation relies on the fact that the parcel
! core properties are NOT in conserved variable form.
call calc_virt_temp( n_points_a, n_points_super_a,                             &
                     core_a_fields(:,i_temperature),                           &
                     core_a_fields(:,i_q_vap),                                 &
                     core_a_fields(:,i_qc_first:i_qc_last),                    &
                     core_a_virt_temp )
call calc_virt_temp( n_points_m, n_points_super_m,                             &
                     core_m_fields(:,i_temperature),                           &
                     core_m_fields(:,i_q_vap),                                 &
                     core_m_fields(:,i_qc_first:i_qc_last),                    &
                     core_m_virt_temp )

! Choose properties from the parcel with the more buoyant core.

! Reset the weights so that the core fields will inherit only
! the values from the most buoyant of the two.
! Combine the edge virtual temperatures by choosing the least buoyant edge
!  (i.e. we always try to make the PDF of Tv as wide as possible)
if ( l_down ) then
  do ic = 1, n_points_a
    ic2 = index_ic(ic)
    ! Choose most negatively buoyant core for downdrafts
    ! TEMPORARY CODE TO PRESERVE KGO: should really test on mass-fluxes > 0
    ! (this code can use core properties from parcel m with zero mass-flux,
    !  which is wrong; fix this soon...)
    if ( core_a_virt_temp(ic) <= core_m_virt_temp(ic2) .or.                    &
         ( .not. core_m_virt_temp(ic2) > zero ) ) then
      weight_core_a(ic) = one
      weight_core_m(ic) = zero
    else
      weight_core_a(ic) = zero
      weight_core_m(ic) = one
    end if
    ! Choose least negatively buoyant edge for downdrafts
    if ( edge_a_virt_temp(ic) > edge_m_virt_temp(ic2) .or.                     &
         ( .not. edge_m_virt_temp(ic2) > zero ) ) then
      edge_m_virt_temp(ic2) = edge_a_virt_temp(ic)
    end if
  end do
else
  do ic = 1, n_points_a
    ic2 = index_ic(ic)
    ! Choose most positively buoyant core for updrafts
    ! TEMPORARY CODE TO PRESERVE KGO: should really test on mass-fluxes > 0
    ! (this code can use core properties from parcel m with zero mass-flux,
    !  which is wrong; fix this soon...)
    if ( core_a_virt_temp(ic) >= core_m_virt_temp(ic2) .or.                    &
         ( .not. core_m_virt_temp(ic2) > zero ) ) then
      weight_core_a(ic) = one
      weight_core_m(ic) = zero
    else
      weight_core_a(ic) = zero
      weight_core_m(ic) = one
    end if
    ! Choose least positively buoyant edge for downdrafts
    if ( edge_a_virt_temp(ic) < edge_m_virt_temp(ic2) .or.                     &
         ( .not. edge_m_virt_temp(ic2) > zero ) ) then
      edge_m_virt_temp(ic2) = edge_a_virt_temp(ic)
    end if
  end do
end if

! Compute combined parcel core properties using the weights set above
do i_field = i_field_first, i_field_last
  do ic = 1, n_points_a
    ic2 = index_ic(ic)
    core_m_fields(ic2,i_field) = weight_core_m(ic) * core_m_fields(ic2,i_field)&
                               + weight_core_a(ic) * core_a_fields(ic,i_field)
  end do
end do


return
end subroutine core_combine

end module core_combine_mod
