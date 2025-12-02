! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module cfl_limit_indep_mod

implicit none

contains



!----------------------------------------------------------------
! Subroutine to apply CFL-limit closure rescaling to each
! convective draft / type / layer independently of the others.
!----------------------------------------------------------------
subroutine cfl_limit_indep( n_conv_types, n_conv_layers,                       &
                            ij_first, ij_last, cmpr_any,                       &
                            layer_mass, fields_np1, res_source,                &
                            scaling )

use comorph_constants_mod, only: real_cvprec, real_hmprec, zero,               &
                     nx_full, ny_full, k_bot_conv, k_top_conv,                 &
                     max_cfl, comorph_timestep, sqrt_min_float
use cmpr_type_mod, only: cmpr_type
use fields_type_mod, only: fields_type, i_q_vap, i_qc_last
use res_source_type_mod, only: res_source_type, i_ent

implicit none

! Number of convection types
integer, intent(in) :: n_conv_types
! Number of distinct convecting layers
integer, intent(in) :: n_conv_layers

! Position indices ( nx*(j-1) + i ) of the first and last
! convecting point on the current segment
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last

! Structure storing compression list indices for the points where
! any kind of convection is occuring
type(cmpr_type), intent(in) :: cmpr_any( k_bot_conv:k_top_conv )

! Full 3-D field of layer-masses (used for finding CFL limit)
real(kind=real_hmprec), intent(in) :: layer_mass                               &
                     ( nx_full, ny_full, k_bot_conv:k_top_conv )

! Structure containing pointers to the primary fields just before
! being updated by convection increments
type(fields_type), intent(in) :: fields_np1

! Structure storing resolved-scale source terms.
type(res_source_type), intent(in) :: res_source                                &
          ( n_conv_types, n_conv_layers, k_bot_conv:k_top_conv )

! Closure scaling (to be reduced if needed to meet CFL)
real(kind=real_cvprec), intent(in out) :: scaling                              &
               ( ij_first:ij_last, n_conv_types, n_conv_layers )

! Compression arrays for layer-mass and water-species
real(kind=real_cvprec) :: layer_mass_cmpr(ij_first:ij_last)
real(kind=real_cvprec) :: q_cmpr

! Res source and layer terms compared against CFL limit
real(kind=real_cvprec) :: q_removed
real(kind=real_cvprec) :: q_available

! Loop counters
integer :: i, j, k, ij, ic, i_type, i_layr, i_field


! Loop over levels
do k = k_bot_conv, k_top_conv
  if ( cmpr_any(k) % n_points > 0 ) then

    ! Compress layer-mass from points where any type / layer of
    ! convection is active
    do ic = 1, cmpr_any(k) % n_points
      i = cmpr_any(k) % index_i(ic)
      j = cmpr_any(k) % index_j(ic)
      ij = nx_full*(j-1)+i
      layer_mass_cmpr(ij) = real( layer_mass(i,j,k), real_cvprec)
    end do

    ! Loop over types and layers
    do i_layr = 1, n_conv_layers
      do i_type = 1, n_conv_types

        ! If any active points on this level
        if ( res_source(i_type,i_layr,k) % cmpr % n_points > 0 ) then

          ! Apply CFL limit for entrained mass
          do ic = 1, res_source(i_type,i_layr,k) % cmpr % n_points
            ! Extract i,j coordinates of the current point
            i = res_source(i_type,i_layr,k) % cmpr % index_i(ic)
            j = res_source(i_type,i_layr,k) % cmpr % index_j(ic)
            ij = nx_full*(j-1)+i

            ! Calculate max closure scaling which meets the
            ! CFL limit, and adjust the actual closure scaling
            ! down to this if it exceeds it
            scaling(ij,i_type,i_layr)                                          &
             = min( scaling(ij,i_type,i_layr),                                 &
                    max_cfl * layer_mass_cmpr(ij)                              &
                  / ( comorph_timestep * max( res_source(i_type,i_layr,k)      &
                                              % res_super(ic,i_ent),           &
                                              sqrt_min_float ) ) )
          end do

          ! Apply limit to avoid creating negative values for
          ! water species
          do i_field = i_q_vap, i_qc_last
            do ic = 1, res_source(i_type,i_layr,k) % cmpr % n_points
              ! If resolved-scale source term is negative
              if ( res_source(i_type,i_layr,k)                                 &
                   % fields_super(ic,i_field) < zero ) then
                ! Extract i,j coordinates of the current point
                i = res_source(i_type,i_layr,k) % cmpr % index_i(ic)
                j = res_source(i_type,i_layr,k) % cmpr % index_j(ic)
                ij = nx_full*(j-1)+i

                ! Extract grid-mean water species value,
                ! filtering out any spurious negative values
                q_cmpr = max( real( fields_np1%list(i_field)%pt(i,j,k),        &
                                    real_cvprec ), zero )

                ! Calculate water-mass removed and water-mass available
                q_removed = -comorph_timestep                                  &
                   * res_source(i_type,i_layr,k) % fields_super(ic,i_field)
                q_available = max_cfl * q_cmpr * layer_mass_cmpr(ij)

                ! Apply CFL limit if removal term exceeds available term
                if ( q_removed > q_available ) then
                  scaling(ij,i_type,i_layr) = min( scaling(ij,i_type,i_layr),  &
                                                   q_available / q_removed )
                end if

              end if
            end do
          end do

        end if  ! ( res_source(i_type,i_layr,k) % cmpr % n_points > 0 )

      end do
    end do

  end if  ! ( cmpr_any(k) % n_points > 0 ) then
end do  ! k = k_bot_conv, k_top_conv


return
end subroutine cfl_limit_indep



!----------------------------------------------------------------
! Alternative version of the above to use in the case where
! fall-back flows are switched on.  In this case,
! the entrainment needs to be summed over the primary draft
! and its fall-back flow
!----------------------------------------------------------------
subroutine cfl_limit_indep_fallback( n_conv_types, n_conv_layers,              &
                                     ij_first, ij_last, cmpr_any,              &
                                     layer_mass, fields_np1,                   &
                                     res_source,                               &
                                     fallback_res_source,                      &
                                     fallback_par_gen,                         &
                                     scaling )

use comorph_constants_mod, only: real_cvprec, real_hmprec, zero,               &
                     nx_full, ny_full, k_bot_conv, k_top_conv,                 &
                     max_cfl, comorph_timestep, sqrt_min_float
use cmpr_type_mod, only: cmpr_type
use res_source_type_mod, only: res_source_type, i_ent
use parcel_type_mod, only: parcel_type, i_massflux_d
use fields_type_mod, only: fields_type, i_q_vap, i_qc_last

implicit none

! Number of convection types
integer, intent(in) :: n_conv_types
! Number of distinct convecting layers
integer, intent(in) :: n_conv_layers

! Position indices ( nx*(j-1) + i ) of the first and last
! convecting point on the current segment
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last

! Structure storing compression list indices for the points where
! any kind of convection is occuring
type(cmpr_type), intent(in) :: cmpr_any( k_bot_conv:k_top_conv )

! Full 3-D field of layer-masses (used for finding CFL limit)
real(kind=real_hmprec), intent(in) :: layer_mass                               &
                     ( nx_full, ny_full, k_bot_conv:k_top_conv )

! Structure containing pointers to the primary fields just before
! being updated by convection increments
type(fields_type), intent(in) :: fields_np1

! Structure storing resolved-scale source terms from the primary
! draft and its fall-back flow
type(res_source_type), intent(in) :: res_source                                &
          ( n_conv_types, n_conv_layers, k_bot_conv:k_top_conv )
type(res_source_type), intent(in) :: fallback_res_source                       &
          ( n_conv_types, n_conv_layers, k_bot_conv:k_top_conv )

! Structure storing innitiating mass-sources for fall-backs
type(parcel_type), intent(in) :: fallback_par_gen                              &
          ( n_conv_types, n_conv_layers, k_bot_conv:k_top_conv )

! Closure scaling (to be reduced if needed to meet CFL)
real(kind=real_cvprec), intent(in out) :: scaling                              &
               ( ij_first:ij_last, n_conv_types, n_conv_layers )

! Sum of entrainment over primary draft and the fall-back flow
real(kind=real_cvprec) :: sum_ent ( ij_first:ij_last )

! Compression arrays for layer-mass and water-species
real(kind=real_cvprec) :: layer_mass_cmpr(ij_first:ij_last)
real(kind=real_cvprec) :: q_cmpr

! Res source and layer terms compared against CFL limit
real(kind=real_cvprec) :: q_removed
real(kind=real_cvprec) :: q_available

! Loop counters
integer :: i, j, k, ij, ic, i_type, i_layr, i_field


! Initialise sum of entrainment to zero
do ij = ij_first, ij_last
  sum_ent(ij) = zero
end do


! Loop over levels
do k = k_bot_conv, k_top_conv
  if ( cmpr_any(k) % n_points > 0 ) then

    ! Compress layer-mass from points where any type / layer of
    ! convection is active
    do ic = 1, cmpr_any(k) % n_points
      i = cmpr_any(k) % index_i(ic)
      j = cmpr_any(k) % index_j(ic)
      ij = nx_full*(j-1)+i
      layer_mass_cmpr(ij) = real( layer_mass(i,j,k), real_cvprec)
    end do

    ! Loop over types and layers
    do i_layr = 1, n_conv_layers
      do i_type = 1, n_conv_types

        ! 1) Calculate CFL limiting of dry-mass removal

        ! Add entrainment contribution from primary draft
        if ( res_source(i_type,i_layr,k) % cmpr % n_points > 0 ) then
          do ic = 1, res_source(i_type,i_layr,k) % cmpr % n_points
            i = res_source(i_type,i_layr,k) % cmpr % index_i(ic)
            j = res_source(i_type,i_layr,k) % cmpr % index_j(ic)
            ij = nx_full*(j-1)+i
            sum_ent(ij) = sum_ent(ij)                                          &
                        + res_source(i_type,i_layr,k)                          &
                          % res_super(ic,i_ent)
          end do
        end if

        ! Add entrainment contribution from fall-back flow
        if ( fallback_res_source(i_type,i_layr,k) % cmpr % n_points > 0 ) then
          do ic = 1, fallback_res_source(i_type,i_layr,k) % cmpr % n_points
            i = fallback_res_source(i_type,i_layr,k) % cmpr % index_i(ic)
            j = fallback_res_source(i_type,i_layr,k) % cmpr % index_j(ic)
            ij = nx_full*(j-1)+i
            sum_ent(ij) = sum_ent(ij)                                          &
                        + fallback_res_source(i_type,i_layr,k)                 &
                          % res_super(ic,i_ent)
          end do
        end if

        ! Subtract fall-back initiating mass-source, which is
        ! included in the fall-back entrainment, but doesn't
        ! need to be included in the quota for the CFL limit
        ! because it is drawn from the primary draft's
        ! detrainment, not from the environment
        if ( fallback_par_gen(i_type,i_layr,k) % cmpr % n_points > 0 ) then
          do ic = 1, fallback_par_gen(i_type,i_layr,k) % cmpr % n_points
            i = fallback_par_gen(i_type,i_layr,k) % cmpr % index_i(ic)
            j = fallback_par_gen(i_type,i_layr,k) % cmpr % index_j(ic)
            ij = nx_full*(j-1)+i
            sum_ent(ij) = sum_ent(ij)                                          &
                        - fallback_par_gen(i_type,i_layr,k)                    &
                          % par_super(ic,i_massflux_d)
          end do
        end if

        ! Apply CFL limit based on sum of entrainment, and reset
        ! sum_ent to zero ready for next level:
        ! a) at primary draft points
        if ( res_source(i_type,i_layr,k) % cmpr % n_points > 0 ) then
          do ic = 1, res_source(i_type,i_layr,k) % cmpr % n_points
            i = res_source(i_type,i_layr,k) % cmpr % index_i(ic)
            j = res_source(i_type,i_layr,k) % cmpr % index_j(ic)
            ij = nx_full*(j-1)+i
            scaling(ij,i_type,i_layr)                                          &
             = min( scaling(ij,i_type,i_layr),                                 &
                    max_cfl * layer_mass_cmpr(ij)                              &
                  / ( comorph_timestep * max( sum_ent(ij),                     &
                                              sqrt_min_float ) ) )
            sum_ent(ij) = zero
          end do
        end if
        ! b) at fall-back flow points
        if ( fallback_res_source(i_type,i_layr,k) % cmpr % n_points > 0 ) then
          do ic = 1, fallback_res_source(i_type,i_layr,k) % cmpr % n_points
            i = fallback_res_source(i_type,i_layr,k) % cmpr % index_i(ic)
            j = fallback_res_source(i_type,i_layr,k) % cmpr % index_j(ic)
            ij = nx_full*(j-1)+i
            scaling(ij,i_type,i_layr)                                          &
             = min( scaling(ij,i_type,i_layr),                                 &
                    max_cfl * layer_mass_cmpr(ij)                              &
                  / ( comorph_timestep * max( sum_ent(ij),                     &
                                              sqrt_min_float ) ) )
            sum_ent(ij) = zero
          end do
        end if


        ! 2) Calculate limiting of removal of water species to
        !    avoid creating negative mixing ratios.
        do i_field = i_q_vap, i_qc_last

          ! Add contribution from primary draft
          if ( res_source(i_type,i_layr,k) % cmpr % n_points > 0 ) then
            do ic = 1, res_source(i_type,i_layr,k) % cmpr % n_points
              i = res_source(i_type,i_layr,k) % cmpr % index_i(ic)
              j = res_source(i_type,i_layr,k) % cmpr % index_j(ic)
              ij = nx_full*(j-1)+i
              sum_ent(ij) = sum_ent(ij)                                        &
                          + res_source(i_type,i_layr,k)                        &
                            % fields_super(ic,i_field)
            end do
          end if

          ! Add contribution from fall-back flow
          if ( fallback_res_source(i_type,i_layr,k) % cmpr %n_points > 0 ) then
            do ic = 1, fallback_res_source(i_type,i_layr,k) % cmpr % n_points
              i = fallback_res_source(i_type,i_layr,k) % cmpr % index_i(ic)
              j = fallback_res_source(i_type,i_layr,k) % cmpr % index_j(ic)
              ij = nx_full*(j-1)+i
              sum_ent(ij) = sum_ent(ij)                                        &
                          + fallback_res_source(i_type,i_layr,k)               &
                            % fields_super(ic,i_field)
            end do
          end if

          ! Apply limit using sum of water source terms, and
          ! reset sum_ent to zero ready for next level:
          ! a) at primary draft points
          if ( res_source(i_type,i_layr,k) % cmpr % n_points > 0 ) then
            do ic = 1, res_source(i_type,i_layr,k) % cmpr % n_points
              ! If net sink of water
              if ( sum_ent(ij) < 0 ) then
                i = res_source(i_type,i_layr,k) % cmpr % index_i(ic)
                j = res_source(i_type,i_layr,k) % cmpr % index_j(ic)
                ij = nx_full*(j-1)+i

                q_cmpr = max( real( fields_np1%list(i_field)%pt(i,j,k),        &
                                    real_cvprec ), zero )

                ! Calculate water-mass removed and water-mass available
                q_removed = -comorph_timestep * sum_ent(ij)
                q_available = max_cfl * q_cmpr * layer_mass_cmpr(ij)

                ! Apply CFL limit if removal term exceeds available term
                if ( q_removed > q_available ) then
                  scaling(ij,i_type,i_layr) = min( scaling(ij,i_type,i_layr),  &
                                                   q_available / q_removed )
                end if

              end if
              sum_ent(ij) = zero
            end do
          end if
          ! b) at fall-back flow points
          if ( fallback_res_source(i_type,i_layr,k) % cmpr %n_points > 0 ) then
            do ic = 1, fallback_res_source(i_type,i_layr,k) % cmpr % n_points
              ! If net sink of water
              if ( sum_ent(ij) < 0 ) then
                i = fallback_res_source(i_type,i_layr,k) % cmpr % index_i(ic)
                j = fallback_res_source(i_type,i_layr,k) % cmpr % index_j(ic)
                ij = nx_full*(j-1)+i

                q_cmpr = max( real( fields_np1%list(i_field)%pt(i,j,k),        &
                                    real_cvprec ), zero )

                ! Calculate water-mass removed and water-mass available
                q_removed = -comorph_timestep * sum_ent(ij)
                q_available = max_cfl * q_cmpr * layer_mass_cmpr(ij)

                ! Apply CFL limit if removal term exceeds available term
                if ( q_removed > q_available ) then
                  scaling(ij,i_type,i_layr) = min( scaling(ij,i_type,i_layr),  &
                                                   q_available / q_removed )
                end if

              end if
              sum_ent(ij) = zero
            end do
          end if

        end do  ! i_field = i_q_vap, i_qc_last

      end do  ! i_type = 1, n_conv_types
    end do  ! i_layr = 1, n_conv_layers

  end if  ! ( cmpr_any(k) % n_points > 0 )
end do  ! k = k_bot_conv, k_top_conv


return
end subroutine cfl_limit_indep_fallback


end module cfl_limit_indep_mod
