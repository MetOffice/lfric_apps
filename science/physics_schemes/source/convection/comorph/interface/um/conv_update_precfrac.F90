! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module conv_update_precfrac_mod

use um_types, only: real_umphys

implicit none

contains

! Subroutine to update the precipitation fraction using the convective
! increments to rain and graupel.
! This is a temporary simple precip fraction increment calculation,
! pending implementation of precip fraction increments inside CoMorph
subroutine conv_update_precfrac( i_call, qrain_star, qgraup_star,              &
                                 frac_bulk_conv, q_prec_b4, precfrac_star )

use atm_fields_bounds_mod, only: tdims
use calc_conv_incs_mod, only: i_call_save_before_conv, i_call_diff_to_get_incs
use mphys_inputs_mod, only: l_mcr_qgraup, l_subgrid_graupel_frac
use comorph_um_namelist_mod, only: rain_area_min

implicit none

! Integer switch indicating whether this is the call to save fields before
! convection, or subtract the values before from after to get increments
integer, intent(in) :: i_call

! Latest updated rain and graupel mixing-ratios
real(kind=real_umphys), intent(in) ::                                          &
                    qrain_star     ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:tdims%k_end )
real(kind=real_umphys), intent(in) ::                                          &
                    qgraup_star    ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:tdims%k_end )

! Convective bulk cloud fraction
real(kind=real_umphys), intent(in) ::                                          &
                    frac_bulk_conv ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:tdims%k_end )

! Saved precip mass mixing-ratio before convection
real(kind=real_umphys), intent(in out) ::                                      &
                        q_prec_b4     ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end,             &
                                        1:tdims%k_end )

! Prognostic precipitation fraction; to be updated on 2nd call
real(kind=real_umphys), intent(in out) ::                                      &
                        precfrac_star ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end,             &
                                        1:tdims%k_end )

! Loop counters
integer :: i, j, k


!$OMP PARALLEL DEFAULT(none) private( i, j, k )                                &
!$OMP SHARED( i_call, l_mcr_qgraup, l_subgrid_graupel_frac, tdims,             &
!$OMP         q_prec_b4, qrain_star, qgraup_star, precfrac_star,               &
!$OMP         frac_bulk_conv, rain_area_min )

! Which call to this routine are we in?
select case (i_call)

case (i_call_save_before_conv)
  ! 1st call: save fields before convection...

  ! Save the precip mixing ratio before convection
  ! (what counts as precip here depends on whether or not graupel
  !  is included in the precip fraction).
  if ( l_mcr_qgraup .and. l_subgrid_graupel_frac ) then
!$OMP do SCHEDULE(STATIC)
    do k = 1, tdims%k_end
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          q_prec_b4(i,j,k) = qrain_star(i,j,k) + qgraup_star(i,j,k)
        end do
      end do
    end do
!$OMP end do NOWAIT
  else
!$OMP do SCHEDULE(STATIC)
    do k = 1, tdims%k_end
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          q_prec_b4(i,j,k) = qrain_star(i,j,k)
        end do
      end do
    end do
!$OMP end do NOWAIT
  end if

case (i_call_diff_to_get_incs)
  ! 2nd call: subtract values before convection to get increments...

  ! New precip fraction is weighted towards the convective cloud
  ! fraction where the precip is produced.
  ! (what counts as precip here depends on whether or not graupel
  !  is included in the precip fraction).
  if ( l_mcr_qgraup .and. l_subgrid_graupel_frac ) then
!$OMP do SCHEDULE(STATIC)
    do k = 1, tdims%k_end
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          ! Where convection has produced precip...
          if ( qrain_star(i,j,k)+qgraup_star(i,j,k) > q_prec_b4(i,j,k) ) then
            ! Calculate precip-mass weighted new precip fraction:
            precfrac_star(i,j,k)                                               &
              = ( q_prec_b4(i,j,k) * precfrac_star(i,j,k)                      &
                + ( qrain_star(i,j,k)+qgraup_star(i,j,k) - q_prec_b4(i,j,k) )  &
                  * frac_bulk_conv(i,j,k) )                                    &
              / ( qrain_star(i,j,k) + qgraup_star(i,j,k) )
          end if
          if ( qrain_star(i,j,k)+qgraup_star(i,j,k) > 0.0 ) then
            precfrac_star(i,j,k) = max( precfrac_star(i,j,k), rain_area_min )
          end if
        end do
      end do
    end do
!$OMP end do NOWAIT
  else
!$OMP do SCHEDULE(STATIC)
    do k = 1, tdims%k_end
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          ! Where convection has produced precip...
          if ( qrain_star(i,j,k) > q_prec_b4(i,j,k) ) then
            ! Calculate rain-mass weighted new rain fraction:
            precfrac_star(i,j,k)                                               &
              = ( q_prec_b4(i,j,k) * precfrac_star(i,j,k)                      &
                + ( qrain_star(i,j,k) - q_prec_b4(i,j,k) )                     &
                  * frac_bulk_conv(i,j,k) )                                    &
              / qrain_star(i,j,k)
          end if
          if ( qrain_star(i,j,k) > 0.0 ) then
            precfrac_star(i,j,k) = max( precfrac_star(i,j,k), rain_area_min )
          end if
        end do
      end do
    end do
!$OMP end do NOWAIT
  end if

end select  ! case(i_call)

!$OMP end PARALLEL


return
end subroutine conv_update_precfrac

end module conv_update_precfrac_mod
