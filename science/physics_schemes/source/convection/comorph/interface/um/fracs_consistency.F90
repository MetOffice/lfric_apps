! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph
module fracs_consistency_mod

use um_types, only: real_umphys

implicit none

contains

! Routine to apply consistency-checks on the cloud and precip fractions.
subroutine fracs_consistency( qcl_star, qcf_star, qcf2_star,                   &
                              qrain_star, qgraup_star,                         &
                              cf_liquid_star, cf_frozen_star, bulk_cf_star,    &
                              precfrac_star )

use atm_fields_bounds_mod, only: tdims
use mphys_inputs_mod, only: l_mcr_qcf2, l_mcr_qgraup, l_mcr_precfrac,          &
                            l_subgrid_graupel_frac

implicit none

! Condensed water species masses that the fractions need to be consistent with
real(kind=real_umphys), intent(in) ::                                          &
                    qcl_star       ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:tdims%k_end )
real(kind=real_umphys), intent(in) ::                                          &
                    qcf_star       ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:tdims%k_end )
real(kind=real_umphys), intent(in) ::                                          &
                    qcf2_star      ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:tdims%k_end )
real(kind=real_umphys), intent(in) ::                                          &
                    qrain_star     ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:tdims%k_end )
real(kind=real_umphys), intent(in) ::                                          &
                    qgraup_star    ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end,                &
                                     1:tdims%k_end )

! Cloud and precip fractions to be updated where not consistent
real(kind=real_umphys), intent(in out) ::                                      &
                        cf_liquid_star ( tdims%i_start:tdims%i_end,            &
                                         tdims%j_start:tdims%j_end,            &
                                         1:tdims%k_end )
real(kind=real_umphys), intent(in out) ::                                      &
                        cf_frozen_star ( tdims%i_start:tdims%i_end,            &
                                         tdims%j_start:tdims%j_end,            &
                                         1:tdims%k_end )
real(kind=real_umphys), intent(in out) ::                                      &
                        bulk_cf_star   ( tdims%i_start:tdims%i_end,            &
                                         tdims%j_start:tdims%j_end,            &
                                         1:tdims%k_end )
real(kind=real_umphys), intent(in out) ::                                      &
                        precfrac_star  ( tdims%i_start:tdims%i_end,            &
                                         tdims%j_start:tdims%j_end,            &
                                         1:tdims%k_end )

! Temporary work variables for cloud water contents
real(kind=real_umphys) :: qcl, qcf, qc

! Max implied in-cloud condensate mixing-ratio.
! Used to compute a corresponding minimum allowed cloud-fraction
! imposed as a safety check to avoid advection or convection creating
! clouds with zero fraction but non-zero mixing-ratio
! (can happen on levels where convection doesn't detrain,
!  but condensate falls from the parcel to the environment
!  as precipitation).
real(kind=real_umphys), parameter :: qc_incloud_max = 5.0e-3  ! 5 g kg-1

! Loop counters
integer :: i, j, k


!$OMP PARALLEL DEFAULT(none) private( i, j, k, qcl, qcf, qc )                  &
!$OMP SHARED( l_mcr_qcf2, tdims,                                               &
!$OMP         cf_liquid_star, cf_frozen_star, bulk_cf_star,                    &
!$OMP         qcl_star, qcf_star, qcf2_star,                                   &
!$OMP         l_mcr_precfrac, l_mcr_qgraup, l_subgrid_graupel_frac,            &
!$OMP         precfrac_star, qrain_star, qgraup_star )

! Minimum allowed cloud-fraction consistent with maximum allowed
! in-cloud condensed water content.
! Note: this check also removes instances of
! frac < 0 or frac > 1 which SL advection can create.
if ( l_mcr_qcf2 ) then
  ! Version of checks with 2nd ice category

!$OMP do SCHEDULE(STATIC)
  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        qcl = max( qcl_star(i,j,k), 0.0 )
        qcf = max( qcf_star(i,j,k), 0.0 ) + max( qcf2_star(i,j,k), 0.0 )
        qc = qcl + qcf
        cf_liquid_star(i,j,k) = min( max( cf_liquid_star(i,j,k),               &
                                          qcl / qc_incloud_max ), 1.0 )
        cf_frozen_star(i,j,k) = min( max( cf_frozen_star(i,j,k),               &
                                          qcf / qc_incloud_max ), 1.0 )
        bulk_cf_star(i,j,k)   = min( max( bulk_cf_star(i,j,k),                 &
                                          qc  / qc_incloud_max ), 1.0 )
      end do
    end do
  end do
!$OMP end do

else  ! .not. l_mcr_qcf2
  ! Version of checks with no 2nd ice category

!$OMP do SCHEDULE(STATIC)
  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        qcl = max( qcl_star(i,j,k), 0.0 )
        qcf = max( qcf_star(i,j,k), 0.0 )
        qc = qcl + qcf
        cf_liquid_star(i,j,k) = min( max( cf_liquid_star(i,j,k),               &
                                          qcl / qc_incloud_max ), 1.0 )
        cf_frozen_star(i,j,k) = min( max( cf_frozen_star(i,j,k),               &
                                          qcf / qc_incloud_max ), 1.0 )
        bulk_cf_star(i,j,k)   = min( max( bulk_cf_star(i,j,k),                 &
                                          qc  / qc_incloud_max ), 1.0 )
      end do
    end do
  end do
!$OMP end do

end if  ! .not. l_mcr_qcf2

!$OMP do SCHEDULE(STATIC)
do k = 1, tdims%k_end
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      ! Also check that the bulk cloud fraction lies between
      ! its allowed limits for minimum and maximum overlap
      ! between the liquid and ice clouds:
      ! max( CFL, CFF )  <=  CF_bulk  <=  CFL + CFF.
      bulk_cf_star(i,j,k) = min( max(                                          &
        max( cf_liquid_star(i,j,k), cf_frozen_star(i,j,k) ),                   &
                                      bulk_cf_star(i,j,k)   ),                 &
             cf_liquid_star(i,j,k) + cf_frozen_star(i,j,k)    )
    end do
  end do
end do
!$OMP end do NOWAIT

! If using prognostic precip fraction, check that the in-precip-fraction
! value of qrain doesn't become implausibly large.  Increase
! the precip fraction where needed to avoid this problem.
if ( l_mcr_precfrac ) then
  if ( l_mcr_qgraup .and. l_subgrid_graupel_frac ) then
    ! Graupel included in the prognostic precip fraction as well as rain
!$OMP do SCHEDULE(STATIC)
    do k = 1, tdims%k_end
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          qc = max( qrain_star(i,j,k), 0.0 ) + max( qgraup_star(i,j,k), 0.0 )
          precfrac_star(i,j,k) = min( max( precfrac_star(i,j,k),               &
                                           qc / qc_incloud_max ), 1.0 )
        end do
      end do
    end do
!$OMP end do NOWAIT
  else
    ! Only rain included in the prognostic precip fraction
!$OMP do SCHEDULE(STATIC)
    do k = 1, tdims%k_end
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          qc = max( qrain_star(i,j,k), 0.0 )
          precfrac_star(i,j,k) = min( max( precfrac_star(i,j,k),               &
                                           qc / qc_incloud_max ), 1.0 )
        end do
      end do
    end do
!$OMP end do NOWAIT
  end if
end if

!$OMP end PARALLEL


return
end subroutine fracs_consistency

end module fracs_consistency_mod
