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
subroutine cfl_limit_sum_ent( n_conv_types, n_conv_layers,                     &
                              ij_first, ij_last,                               &
                              scaling, res_source_k,                           &
                              sum_ent_k, sum_res_source_q_k )

use comorph_constants_mod, only: real_cvprec, nx_full
use res_source_type_mod, only: res_source_type, i_ent
use fields_type_mod, only: i_q_vap, i_qc_last

implicit none

! Number of convection types
integer, intent(in) :: n_conv_types
! Number of distinct convecting layers
integer, intent(in) :: n_conv_layers

! First and last ij-indices on the current segment
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last

! Closure scaling calculated so far
real(kind=real_cvprec), intent(in) :: scaling                                  &
               ( ij_first:ij_last, n_conv_types, n_conv_layers )

! Structure containing resolved-scale source terms to be added
type(res_source_type), intent(in) :: res_source_k                              &
                             ( n_conv_types, n_conv_layers )

! Total mass sink from entrainment, to be incremented
real(kind=real_cvprec), intent(in out) :: sum_ent_k                            &
                                         ( ij_first:ij_last )

! Total source term for each water species, to be incremented
real(kind=real_cvprec), intent(in out) :: sum_res_source_q_k                   &
                        ( ij_first:ij_last, i_q_vap:i_qc_last )

! Loop counters
integer :: i_type, i_layr, i, j, ij, ic, i_field

! Loop over convecting layers and convection types
do i_layr = 1, n_conv_layers
  do i_type = 1, n_conv_types

    ! If there are any resolved-scale source terms
    ! for this convection layer / type at the current level
    if ( res_source_k(i_type,i_layr) % cmpr % n_points > 0 ) then

      do ic = 1, res_source_k(i_type,i_layr) % cmpr % n_points
        i = res_source_k(i_type,i_layr) % cmpr % index_i(ic)
        j = res_source_k(i_type,i_layr) % cmpr % index_j(ic)
        ij = nx_full*(j-1)+i
        ! Increment total mass sink from entrainment
        sum_ent_k(ij) = sum_ent_k(ij)                                          &
           + res_source_k(i_type,i_layr) % res_super(ic,i_ent)                 &
             * scaling(ij,i_type,i_layr)
      end do

      do i_field = i_q_vap, i_qc_last
        do ic = 1, res_source_k(i_type,i_layr) % cmpr % n_points
          i = res_source_k(i_type,i_layr) % cmpr % index_i(ic)
          j = res_source_k(i_type,i_layr) % cmpr % index_j(ic)
          ij = nx_full*(j-1)+i
          ! Increment water species source term total
          sum_res_source_q_k(ij,i_field)                                       &
           = sum_res_source_q_k(ij,i_field)                                    &
           + res_source_k(i_type,i_layr)%fields_super(ic,i_field)              &
             * scaling(ij,i_type,i_layr)
        end do
      end do

    end if

  end do
end do

return
end subroutine cfl_limit_sum_ent


!----------------------------------------------------------------
! Alternative version of the above to use for fall-back flows.
! In this case, fall-back initiation masses need to be subtracted
! since they don't count towards the CFL limit
!----------------------------------------------------------------
subroutine cfl_limit_sum_ent_fallback(                                         &
                              n_conv_types, n_conv_layers,                     &
                              ij_first, ij_last,                               &
                              scaling, res_source_k, par_gen_k,                &
                              sum_ent_k, sum_res_source_q_k )

use comorph_constants_mod, only: real_cvprec, nx_full
use res_source_type_mod, only: res_source_type, i_ent
use parcel_type_mod, only: parcel_type, i_massflux_d
use fields_type_mod, only: i_q_vap, i_qc_last

implicit none

! Number of convection types
integer, intent(in) :: n_conv_types
! Number of distinct convecting layers
integer, intent(in) :: n_conv_layers

! First and last ij-indices on the current segment
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last

! Closure scaling calculated so far
real(kind=real_cvprec), intent(in) :: scaling                                  &
               ( ij_first:ij_last, n_conv_types, n_conv_layers )

! Structure containing resolved-scale source terms to be added
type(res_source_type), intent(in) :: res_source_k                              &
                             ( n_conv_types, n_conv_layers )

! Structure containing initiation mass-sources; these need
! to be subtracted as they are included in the entrained
! mass in res_source, but shouldn't be included here.
type(parcel_type), intent(in) :: par_gen_k                                     &
                             ( n_conv_types, n_conv_layers )

! Total mass sink from entrainment, to be incremented
real(kind=real_cvprec), intent(in out) :: sum_ent_k                            &
                                         ( ij_first:ij_last )

! Total source term for each water species, to be incremented
real(kind=real_cvprec), intent(in out) :: sum_res_source_q_k                   &
                        ( ij_first:ij_last, i_q_vap:i_qc_last )

! Loop counters
integer :: i_type, i_layr, i, j, ij, ic, i_field

! Loop over convecting layers and convection types
do i_layr = 1, n_conv_layers
  do i_type = 1, n_conv_types

    ! If there are any resolved-scale source terms
    ! for this convection layer / type at the current level
    if ( res_source_k(i_type,i_layr) % cmpr % n_points > 0 ) then

      do ic = 1, res_source_k(i_type,i_layr) % cmpr % n_points
        ! Get cmpr_any compression index of this point from grid
        i = res_source_k(i_type,i_layr) % cmpr % index_i(ic)
        j = res_source_k(i_type,i_layr) % cmpr % index_j(ic)
        ij = nx_full*(j-1)+i
        ! Increment total mass sink from entrainment
        sum_ent_k(ij) = sum_ent_k(ij)                                          &
           + res_source_k(i_type,i_layr) % res_super(ic,i_ent)                 &
             * scaling(ij,i_type,i_layr)
      end do

      do i_field = i_q_vap, i_qc_last
        do ic = 1, res_source_k(i_type,i_layr) % cmpr % n_points
          i = res_source_k(i_type,i_layr) % cmpr % index_i(ic)
          j = res_source_k(i_type,i_layr) % cmpr % index_j(ic)
          ij = nx_full*(j-1)+i
          ! Increment water species source term total
          sum_res_source_q_k(ij,i_field)                                       &
           = sum_res_source_q_k(ij,i_field)                                    &
           + res_source_k(i_type,i_layr)%fields_super(ic,i_field)              &
             * scaling(ij,i_type,i_layr)
        end do
      end do

    end if

    ! If there are any initiating mass-sources, subtract them
    ! so that we only sum the actual entrainment
    if ( par_gen_k(i_type,i_layr) % cmpr % n_points > 0 ) then

      do ic = 1, par_gen_k(i_type,i_layr) % cmpr % n_points
        ! Get cmpr_any compression index of this point from grid
        i = par_gen_k(i_type,i_layr) % cmpr % index_i(ic)
        j = par_gen_k(i_type,i_layr) % cmpr % index_j(ic)
        ij = nx_full*(j-1)+i
        ! Decrement total mass sink from entrainment
        sum_ent_k(ij) = sum_ent_k(ij)                                          &
         - par_gen_k(i_type,i_layr) % par_super(ic,i_massflux_d)               &
           * scaling(ij,i_type,i_layr)
      end do

    end if

  end do
end do

return
end subroutine cfl_limit_sum_ent_fallback



end module cfl_limit_sum_ent_mod
