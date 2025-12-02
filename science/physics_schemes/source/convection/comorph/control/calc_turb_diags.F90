! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_turb_diags_mod

implicit none

contains


! Subroutine to calculate the turbulence-based parcel
! perturbations and parcel radius at all points, just
! for diagnostics.  This is basically just a compression
! wrapper around the routine calc_turb_parcel
subroutine calc_turb_diags( turb, comorph_diags )

use comorph_constants_mod, only: real_cvprec, real_hmprec,                     &
                                 nx_full, ny_full, k_bot_conv, k_top_init,     &
                                 par_gen_radius_fac
use turb_type_mod, only: turb_type, n_turb
use fields_type_mod, only: i_wind_u, i_q_vap
use cmpr_type_mod, only: cmpr_type, cmpr_alloc, cmpr_dealloc
use comorph_diags_type_mod, only: comorph_diags_type
use calc_turb_perts_mod, only: calc_turb_perts
use compress_mod, only: compress
use decompress_mod, only: decompress

implicit none

! Structure containing input turbulence fields
type(turb_type), intent(in) :: turb

! Structure containing all comorph diagnostics and meta-data
type(comorph_diags_type), intent(in out) :: comorph_diags

! Structure storing compression / decompression indices
type(cmpr_type) :: cmpr

! Compressed super-array to store turbulence fields on 1 level
real(kind=real_cvprec) :: turb_cmpr( nx_full*ny_full, n_turb )

! Compressed super-array to store perturbations on 1 level
real(kind=real_cvprec) :: pert_cmpr                                            &
                          ( nx_full*ny_full, i_wind_u:i_q_vap )

! Lower and upper bounds of full 3-D arrays
integer :: lb(3), ub(3)

! Loop counters
integer :: i, j, k, ic, i_field, i_diag, i_turb


! Setup compression indices for all points on a slice
cmpr % n_points = nx_full*ny_full
call cmpr_alloc( cmpr, cmpr % n_points )
do j = 1, ny_full
  do i = 1, nx_full
    ic = nx_full*(j-1) + i
    cmpr % index_i(ic) = i
    cmpr % index_j(ic) = j
  end do
end do


! Loop over levels
! (turbulence fields are not nessecarily available above level
!  k_top_init+1, so only loop up to there).
!
!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP SHARED(  k_bot_conv, k_top_init, cmpr, turb, comorph_diags,              &
!$OMP          i_wind_u, i_q_vap )                                             &
!$OMP private( k, lb, ub, i_turb, i_field, ic, i_diag, turb_cmpr, pert_cmpr )
do k = k_bot_conv, k_top_init+1

  ! Compress turbulence fields into super-array
  do i_turb = 1, n_turb
    lb = lbound( turb%list(i_turb)%pt )
    ub = ubound( turb%list(i_turb)%pt )
    call compress( cmpr, lb(1:2), ub(1:2), turb%list(i_turb)%pt(:,:,k),        &
                   turb_cmpr(:,i_turb) )
  end do

  ! Calculate perturbations
  call calc_turb_perts( cmpr%n_points, cmpr%n_points,                          &
                        turb_cmpr, pert_cmpr )

  ! Decompress requested fields into diagnostic arrays...

  ! Parcel perturbations:
  if ( comorph_diags % turb_fields_pert % n_diags > 0 ) then
    ! Loop over requested field perturbations
    do i_diag = 1, comorph_diags % turb_fields_pert % n_diags

      ! Extract the fields super-array address for this diag
      i_field = comorph_diags % turb_fields_pert                               &
                   % list(i_diag)%pt % i_field

      ! Only do anything if it is for u,v,w,T,q; perturbations
      ! are not yet available for the other fields.  If
      ! perturbations are requested for the other fields, they
      ! will just be left as initialised to zero.
      if ( i_field >= i_wind_u .and. i_field <= i_q_vap ) then
        ! Copy into 3-D output diag array
        lb = lbound( comorph_diags % turb_fields_pert                          &
                     % list(i_diag)%pt % field_3d )
        ub = ubound( comorph_diags % turb_fields_pert                          &
                     % list(i_diag)%pt % field_3d )
        call decompress( cmpr, pert_cmpr(:,i_field), lb(1:2), ub(1:2),         &
                         comorph_diags % turb_fields_pert                      &
                         % list(i_diag)%pt % field_3d(:,:,k) )
      end if

    end do
  end if

end do  ! k = k_bot_conv, k_top_init+1
!$OMP end PARALLEL do


! Deallocate compression indices
call cmpr_dealloc( cmpr )


! Turbulence-based parcel radius.
! This is on full-levels and is trivially set using the input
! turb % lengthscale
! So no need to do any compression / decompression.

if ( comorph_diags % turb_radius % flag ) then
!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP SHARED(  turb, comorph_diags, nx_full, ny_full, k_bot_conv, k_top_init,  &
!$OMP          par_gen_radius_fac )                                            &
!$OMP private( i, j, k )
  do k = k_bot_conv, k_top_init
    do j = 1, ny_full
      do i = 1, nx_full
        comorph_diags % turb_radius % field_3d(i,j,k)                          &
          = turb % lengthscale(i,j,k) * real( par_gen_radius_fac, real_hmprec )
      end do
    end do
  end do
!$OMP end PARALLEL do
end if


return
end subroutine calc_turb_diags


end module calc_turb_diags_mod
