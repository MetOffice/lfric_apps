! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module cfl_limit_sum_ent_mod

implicit none

contains


!----------------------------------------------------------------
! Subroutine to sum entrainment contributions over all types /
! layers for a primary updraft or downdraft
!----------------------------------------------------------------
subroutine cfl_limit_sum_ent( ij_first, ij_last, n_points_res, n_fields_tot,   &
                              scaling, res_source_cmpr,                        &
                              res_source_super, res_source_fields,             &
                              sum_ent_k, sum_res_source_q_k )

use comorph_constants_mod, only: real_cvprec, nx_full,                         &
                                 i_cfl_closure, i_cfl_closure_mass_q
use cmpr_type_mod, only: cmpr_type
use res_source_type_mod, only: n_res, i_ent
use fields_type_mod, only: i_q_vap, i_qc_last

implicit none

! First and last ij-indices on the current segment
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last

! Size of source term arrays
integer, intent(in) :: n_points_res

! Number of fields in the res_source fields array
integer, intent(in) :: n_fields_tot

! Closure scaling calculated so far
real(kind=real_cvprec), intent(in) :: scaling                                  &
                                      ( ij_first:ij_last )

! Structure containing source term compression indices
type(cmpr_type), intent(in) :: res_source_cmpr

! Arrays containing resolved-scale source terms to be added
real(kind=real_cvprec), intent(in) :: res_source_super                         &
                                      ( n_points_res, n_res )
real(kind=real_cvprec), intent(in) :: res_source_fields                        &
                                      ( n_points_res, n_fields_tot )

! Total mass sink from entrainment, to be incremented
real(kind=real_cvprec), intent(in out) :: sum_ent_k                            &
                                         ( ij_first:ij_last )

! Total source term for each water species, to be incremented
real(kind=real_cvprec), intent(in out) :: sum_res_source_q_k                   &
                        ( ij_first:ij_last, i_q_vap:i_qc_last )

! Loop counters
integer :: i, j, ij, ic, i_field


do ic = 1, res_source_cmpr % n_points
  i = res_source_cmpr % index_i(ic)
  j = res_source_cmpr % index_j(ic)
  ij = nx_full*(j-1)+i
  ! Increment total mass sink from entrainment
  sum_ent_k(ij) = sum_ent_k(ij)                                                &
    + res_source_super(ic,i_ent) * scaling(ij)
end do

if ( i_cfl_closure == i_cfl_closure_mass_q ) then
  do i_field = i_q_vap, i_qc_last
    do ic = 1, res_source_cmpr % n_points
      i = res_source_cmpr % index_i(ic)
      j = res_source_cmpr % index_j(ic)
      ij = nx_full*(j-1)+i
      ! Increment water species source term total
      sum_res_source_q_k(ij,i_field) = sum_res_source_q_k(ij,i_field)          &
        + res_source_fields(ic,i_field) * scaling(ij)
    end do
  end do
end if


return
end subroutine cfl_limit_sum_ent

end module cfl_limit_sum_ent_mod
