! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module comorph_conv_cloud_extras_mod

use um_types, only: real_umphys

implicit none

contains


! Subroutine to calculate convective cloud base and top
! model-levels and other UM convective cloud fields, based on
! the 3-D convective cloud amount and liquid water output
! by CoMorph
subroutine comorph_conv_cloud_extras(                                          &
             n_conv_levels, rho_dry_th, rho_wet_th,                            &
             r_theta_levels, r_rho_levels,                                     &
             cca0, ccw0, frac_bulk_conv,                                       &
             cclwp0, ccb0, cct0, lcbase0, cca, ccw,                            &
             cclwp, cca_2d, lcca, ccb, cct, lcbase, lctop )

use nlsizes_namelist_mod, only: row_length, rows, model_levels
use gen_phys_inputs_mod, only: l_mr_physics
use atm_fields_bounds_mod, only: tdims_l, tdims

implicit none

! Highest model-level where convection is allowed
integer, intent(in) :: n_conv_levels

! Dry density and wet density on theta-levels, for calculating
! column-integrated convective liquid water path diagnostic
real(kind=real_umphys), intent(in) :: rho_dry_th                               &
                            ( row_length, rows, model_levels-1 )
real(kind=real_umphys), intent(in) :: rho_wet_th                               &
                            ( row_length, rows, model_levels-1 )
real(kind=real_umphys), intent(in) :: r_theta_levels                           &
     (tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,0:tdims%k_end)
real(kind=real_umphys), intent(in) :: r_rho_levels                             &
     (tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,tdims%k_end)
! Prognostics for 3-D convective cloud amount and cloud
! liquid water, output by CoMorph.
! Note: these are for the convective liquid cloud only!
real(kind=real_umphys), intent(in) :: cca0 (row_length,rows,model_levels)
real(kind=real_umphys), intent(in) :: ccw0 (row_length,rows,model_levels)

! Diagnosed convective bulk cloud fraction output by CoMorph
! (includes the ice cloud which is excluded from cca0)
real(kind=real_umphys), intent(in) ::                                          &
                    frac_bulk_conv (row_length,rows,model_levels)

! Prognostic for vertically-integrated convective liquid water path
real(kind=real_umphys), intent(out) :: cclwp0  (row_length,rows)

! Prognostics for integer level corresponding to convective
! cloud base and top for the highest convecting layer
integer, intent(out) :: ccb0    (row_length,rows)
integer, intent(out) :: cct0    (row_length,rows)
! Prognostic for model-level of base of lowest convecting layer
integer, intent(out) :: lcbase0 (row_length,rows)

! Diagnostics for 3-D convective cloud amount and cloud
! liquid water
real(kind=real_umphys), intent(out) :: cca (row_length,rows,model_levels)
real(kind=real_umphys), intent(out) :: ccw (row_length,rows,model_levels)

! Diagnostics for 2-D convective cloud amount and
! vertically-integrated liquid water path
real(kind=real_umphys), intent(out) :: cclwp  (row_length,rows)
real(kind=real_umphys), intent(out) :: cca_2d (row_length,rows)
! Diagnostic of convective cloud fraction on the lowest
! model-level containing convective cloud
real(kind=real_umphys), intent(out) :: lcca   (row_length,rows)

! Diagnostics for integer level corresponding to convective
! cloud base and top for the highest convecting layer
integer, intent(out) :: ccb    (row_length,rows)
integer, intent(out) :: cct    (row_length,rows)
! Diagnostics for model-level of base and top of lowest layer
integer, intent(out) :: lcbase (row_length,rows)
integer, intent(out) :: lctop  (row_length,rows)


! Flags for non-zero CCA at current and next level
logical :: l_cca_at_k
logical :: l_cca_at_kp1

! Loop counters
integer :: i, j, k


!$OMP PARALLEL DEFAULT(none) private( i, j, k, l_cca_at_k, l_cca_at_kp1 )      &
!$OMP SHARED( row_length, rows, n_conv_levels, model_levels, l_mr_physics,     &
!$OMP         cclwp0, ccb0, cct0, lcbase0, lctop, cca0, ccw0, frac_bulk_conv,  &
!$OMP         r_rho_levels, r_theta_levels, rho_dry_th, rho_wet_th,            &
!$OMP         cca_2d, cca, ccw, cclwp, ccb, cct, lcbase, lcca )

! Initialise extra prognostics to zero
!$OMP do SCHEDULE(STATIC)
do j = 1, rows
  do i = 1, row_length
    cclwp0(i,j)  = 0.0
    ccb0(i,j)    = 0
    cct0(i,j)    = 0
    lcbase0(i,j) = 0
    lctop(i,j)   = 0
  end do
end do
!$OMP end do

! Find the base and top heights based on where the
! bulk convective cloud fraction is non-zero

! The loop over levels k must be sequential;
! therefore we have the loop over rows j outermost so we can OMP parallelise it
!$OMP do SCHEDULE(STATIC)
do j = 1, rows
  do i = 1, row_length
    ! Check for cloud-base at model-level 1
    ! (missed out in the loop below)
    if ( frac_bulk_conv(i,j,1) > 0.0 ) then
      ccb0(i,j) = 1
      lcbase0(i,j) = 1
    end if
  end do
  do k = 1, n_conv_levels
    do i = 1, row_length

      ! Check whether CCA is nonzero at current and next level
      l_cca_at_k   = frac_bulk_conv(i,j,k)   > 0.0
      l_cca_at_kp1 = frac_bulk_conv(i,j,k+1) > 0.0

      ! If CCA is nonzero at k+1 but not at k, k+1 is cloud-base
      if ( l_cca_at_kp1 .and. (.not. l_cca_at_k ) ) then
        ccb0(i,j) = k + 1
        if ( lcbase0(i,j) == 0 ) lcbase0(i,j) = k + 1
      end if

      ! If CCA is nonzero at k but not at k+1, k is cloud-top
      if ( l_cca_at_k .and. (.not. l_cca_at_kp1 ) ) then
        cct0(i,j) = k
        if ( lctop(i,j) == 0 ) lctop(i,j) = k
      end if

    end do
  end do
end do  ! j = 1, rows
!$OMP end do

! Calculate column-integrated convective cloud liquid water path
! (in kg per m2).
! Note: this attempts to replicate the calculation of cclwp in the
! 6A convection scheme.  ccw0 is not a grid-mean quantity, but
! the in-cloud convective water content.  In the vertical integral,
! no account is taken of the variation of the convective cloud
! fraction cca with height.  The resulting quantity is not
! related to the total convective cloud-water in the grid-column;
! rather, it is the local maximum column water content in the
! convective core (assuming maximal overlap of the convective cloud
! at all heights).
! Which density to use depends on whether using mixing ratios.
! Vertical integral; the j-loop over rows is outermost here,
! so that we can OMP paralellise it (we can't parallelise the k-loop).
if ( l_mr_physics ) then
!$OMP do SCHEDULE(STATIC)
  do j = 1, rows
    do i = 1, row_length
      cclwp0(i,j)                                                              &
             = ( r_rho_levels(i,j,2) - r_theta_levels(i,j,0) )                 &
             * rho_dry_th(i,j,1) * ccw0(i,j,1)
    end do
    do k = 2, n_conv_levels
      do i = 1, row_length
        cclwp0(i,j) = cclwp0(i,j)                                              &
             + ( r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k) )                 &
             * rho_dry_th(i,j,k) * ccw0(i,j,k)
      end do
    end do
  end do  ! j = 1, rows
!$OMP end do
else  ! ( l_mr_physics )
!$OMP do SCHEDULE(STATIC)
  do j = 1, rows
    do i = 1, row_length
      cclwp0(i,j)                                                              &
             = ( r_rho_levels(i,j,2) - r_theta_levels(i,j,0) )                 &
             * rho_wet_th(i,j,1) * ccw0(i,j,1)
    end do
    do k = 2, n_conv_levels
      do i = 1, row_length
        cclwp0(i,j) = cclwp0(i,j)                                              &
             + ( r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k) )                 &
             * rho_wet_th(i,j,k) * ccw0(i,j,k)
      end do
    end do
  end do  ! j = 1, rows
!$OMP end do
end if  ! ( l_mr_physics )

! Set the 2-D CCA value to the max value in the column
! k-loop can't be safely parallelised, so j-loop is outermost.
!$OMP do SCHEDULE(STATIC)
do j = 1, rows
  do i = 1, row_length
    cca_2d(i,j) = 0.0
  end do
  do k = 1, n_conv_levels
    do i = 1, row_length
      cca_2d(i,j) = max( cca_2d(i,j), cca0(i,j,k) )
    end do
  end do
end do  ! j = 1, rows
!$OMP end do

! Set the diagnostics equal to the prognostics
!$OMP do SCHEDULE(STATIC)
do k = 1, model_levels
  do j = 1, rows
    do i = 1, row_length
      cca(i,j,k) = cca0(i,j,k)
      ccw(i,j,k) = ccw0(i,j,k)
    end do
  end do
end do
!$OMP end do
!$OMP do SCHEDULE(STATIC)
do j = 1, rows
  do i = 1, row_length
    cclwp(i,j)  = cclwp0(i,j)
    ccb(i,j)    = ccb0(i,j)
    cct(i,j)    = cct0(i,j)
    lcbase(i,j) = lcbase0(i,j)
  end do
end do
!$OMP end do

! Set value of CCA at lowest cloud-base
!$OMP do SCHEDULE(STATIC)
do j = 1, rows
  do i = 1, row_length
    if ( lcbase(i,j) > 0 ) then
      lcca(i,j) = cca(i,j,lcbase(i,j))
    else
      lcca(i,j) = 0.0
    end if
  end do
end do
!$OMP end do

!$OMP end PARALLEL


return
end subroutine comorph_conv_cloud_extras


end module comorph_conv_cloud_extras_mod
