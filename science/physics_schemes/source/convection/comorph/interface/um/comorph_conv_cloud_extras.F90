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
             cca, ccw, cca_bulk,                                               &
             cclwp, cca_2d, lcca, ccb, cct, lcbase, lctop,                     &
             cca0, ccw0, cclwp0, ccb0, cct0, lcbase0 )

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

! Model-level heights above Earth centre
real(kind=real_umphys), intent(in) :: r_theta_levels                           &
                            ( tdims_l%i_start:tdims_l%i_end,                   &
                              tdims_l%j_start:tdims_l%j_end,                   &
                              0:tdims%k_end )
real(kind=real_umphys), intent(in) :: r_rho_levels                             &
                            ( tdims_l%i_start:tdims_l%i_end,                   &
                              tdims_l%j_start:tdims_l%j_end,                   &
                              tdims%k_end )

! Diagnostics for 3-D convective cloud amount and cloud water content...

! Either liquid+ice fraction, or just liquid, depending on i_convcloud.
real(kind=real_umphys), intent(in) :: cca (row_length,rows,model_levels)

! Either liquid+ice content, or just liquid, depending on i_convcloud.
! Gets converted from grid-mean (output by comorph) to in-cloud value.
real(kind=real_umphys), intent(in out) :: ccw (row_length,rows,model_levels)

! Diagnosed convective bulk cloud fraction output by CoMorph
! (includes the ice cloud even when that is excluded from cca)
real(kind=real_umphys), intent(in) :: cca_bulk                                 &
                                      (row_length,rows,model_levels)

! Diagnostics for 2-D convective cloud amount and
! vertically-integrated cloud water path
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

! Prognostics for 3-D convective cloud amount and cloud water content
real(kind=real_umphys), intent(in out) :: cca0 (row_length,rows,model_levels)
real(kind=real_umphys), intent(in out) :: ccw0 (row_length,rows,model_levels)

! Prognostic for vertically-integrated convective cloud water path
real(kind=real_umphys), intent(in out) :: cclwp0  (row_length,rows)

! Prognostics for integer level corresponding to convective
! cloud base and top for the highest convecting layer
integer, intent(in out) :: ccb0    (row_length,rows)
integer, intent(in out) :: cct0    (row_length,rows)
! Prognostic for model-level of base of lowest convecting layer
integer, intent(in out) :: lcbase0 (row_length,rows)

! Density * delta z used in ertical integral
real(kind=real_umphys) :: rhodz

! Flags for non-zero CCA at each level
logical :: l_cc (row_length,rows,0:n_conv_levels+1)

! Loop counters
integer :: i, j, k


!$OMP PARALLEL DEFAULT(NONE) PRIVATE( i, j, k, rhodz )                         &
!$OMP SHARED( row_length, rows, n_conv_levels, l_mr_physics, l_cc,             &
!$OMP         cclwp0, ccb0, cct0, lcbase0, lctop, cca0, ccw0, cca_bulk,        &
!$OMP         r_rho_levels, r_theta_levels, rho_dry_th, rho_wet_th,            &
!$OMP         cca_2d, cca, ccw, cclwp, ccb, cct, lcbase, lcca )

! Initialise base and top levels to zero
!$OMP DO SCHEDULE(STATIC)
do j = 1, rows
  do i = 1, row_length
    ccb(i,j)    = 0
    cct(i,j)    = 0
    lcbase(i,j) = 0
    lctop(i,j)  = 0
  end do
end do
!$OMP END DO NOWAIT

! Find the base and top heights based on where the
! bulk convective cloud fraction is non-zero...

! Set flag for nonzero CCA
!$OMP DO SCHEDULE(STATIC)
do k = 1, n_conv_levels
  do j = 1, rows
    do i = 1, row_length
      l_cc(i,j,k) = cca_bulk(i,j,k) > 0.0
    end do
  end do
end do
!$OMP END DO NOWAIT
! Set false at upper and lower boundaries, so we find cloud-base at k=1
! or cloud-top at n_conv_levels
!$OMP DO SCHEDULE(STATIC)
do j = 1, rows
  do i = 1, row_length
    l_cc(i,j,0) = .false.
    l_cc(i,j,n_conv_levels+1) = .false.
  end do
end do
!$OMP END DO

! The loop over levels k must be sequential;
! therefore we have the loop over rows j outermost so we can OMP parallelise it
!$OMP DO SCHEDULE(STATIC)
do j = 1, rows
  do k = 1, n_conv_levels
    do i = 1, row_length
      if ( l_cc(i,j,k) ) then
        ! If CCA is nonzero at k but not at k-1, k is cloud-base
        if ( .not. l_cc(i,j,k-1) ) then
          ccb(i,j) = k
          if ( lcbase(i,j) == 0 ) lcbase(i,j) = k
        end if
        ! If CCA is nonzero at k but not at k+1, k is cloud-top
        if ( .not. l_cc(i,j,k+1) ) then
          cct(i,j) = k
          if ( lctop(i,j) == 0 ) lctop(i,j) = k
        end if
      end if  ! ( l_cc(i,j,k) )
    end do
  end do
end do  ! j = 1, rows
!$OMP END DO

! Set the 2-D CCA value to the max value in the column
! k-loop can't be safely parallelised, so j-loop is outermost.
!$OMP DO SCHEDULE(STATIC)
do j = 1, rows
  do i = 1, row_length
    cca_2d(i,j) = 0.0
  end do
  do k = 1, n_conv_levels
    do i = 1, row_length
      cca_2d(i,j) = max( cca_2d(i,j), cca(i,j,k) )
    end do
  end do
end do  ! j = 1, rows
!$OMP END DO NOWAIT

! Set value of CCA at lowest cloud-base
!$OMP DO SCHEDULE(STATIC)
do j = 1, rows
  do i = 1, row_length
    if ( lcbase(i,j) > 0 ) then
      lcca(i,j) = cca(i,j,lcbase(i,j))
    else
      lcca(i,j) = 0.0
    end if
  end do
end do
!$OMP END DO NOWAIT

! Update the prognostics using MAX of the new diagnostic values and the
! existing values passed in...
!$OMP DO SCHEDULE(STATIC)
do k = 1, n_conv_levels
  do j = 1, rows
    do i = 1, row_length
      ! Convert the prognostic ccw from in-cloud to grid-mean
      ccw0(i,j,k) = ccw0(i,j,k) * cca0(i,j,k)
      ! Take max of grid-mean conv cloud water contents
      ccw0(i,j,k) = max( ccw0(i,j,k), ccw(i,j,k) )
      ! Take max of conv cloud areas
      cca0(i,j,k) = max( cca0(i,j,k), cca(i,j,k) )
      ! Convert the prognostic and diagnostic ccw to in-cloud water contents,
      ! as expected by the rest of the UM
      if ( ccw0(i,j,k) > 0.0 )  ccw0(i,j,k) = ccw0(i,j,k) / cca0(i,j,k)
      if ( ccw(i,j,k)  > 0.0 )  ccw(i,j,k)  = ccw(i,j,k)  / cca(i,j,k)
    end do
  end do
end do
!$OMP END DO

! Update prognostics for base and top levels...
! First update the 3D flag for nonzero CCA to account for any pre-existing
! cca0 (so far only set to true where latest diagnosed cca is nonzero)
!$OMP DO SCHEDULE(STATIC)
do k = 1, n_conv_levels
  do j = 1, rows
    do i = 1, row_length
      ! Note if i_convcloud_liqonly, then some convective cloud points will
      ! still have zero cca0 due to there being ice but no liquid.
      ! So need to explicitly set the mask to true between the prognostic
      ! base and top levels.
      ! Note we can't do the same for the "lowest" base and top, since
      ! there is no prognostic for lctop.
      if ( cca0(i,j,k) > 0.0 .or. ( k>=ccb0(i,j) .and. k<=cct0(i,j) ) ) then
        l_cc(i,j,k) = .true.
      end if
    end do
  end do
end do
!$OMP END DO
! Reset prognostic ccb and cct based on latest mask of convective cloud points
!$OMP DO SCHEDULE(STATIC)
do j = 1, rows
  do k = 1, n_conv_levels
    do i = 1, row_length
      if ( l_cc(i,j,k) ) then
        ! If CCA is nonzero at k but not at k-1, k is cloud-base
        if ( .not. l_cc(i,j,k-1) )  ccb0(i,j) = k
        ! If CCA is nonzero at k but not at k+1, k is cloud-top
        if ( .not. l_cc(i,j,k+1) )  cct0(i,j) = k
      end if  ! ( l_cc(i,j,k) )
    end do
  end do
end do  ! j = 1, rows
!$OMP END DO NOWAIT
! Set prog lowest cloud-base to min of new and existing values if both set
!$OMP DO SCHEDULE(STATIC)
do j = 1, rows
  do i = 1, row_length
    if ( lcbase(i,j) > 0 ) then
      if ( lcbase0(i,j) > 0 ) then
        lcbase0(i,j) = min( lcbase0(i,j), lcbase(i,j) )
      else
        lcbase0(i,j) = lcbase(i,j)
      end if
    end if
  end do
end do
!$OMP END DO NOWAIT

! Calculate column-integrated convective cloud water path (in kg per m2).
! Note: this attempts to replicate the calculation of cclwp in the
! 6A convection scheme.  ccw is NOT a grid-mean quantity, but
! the in-cloud convective water content.  In the vertical integral,
! no account is taken of the variation of the convective cloud
! fraction cca with height.  The resulting quantity is NOT
! related to the total convective cloud-water in the grid-column;
! rather, it is the local maximum column water content in the
! convective core (assuming maximal overlap of the convective cloud
! at all heights).
! Which density to use depends on whether using mixing ratios.
! Vertical integral; the j-loop over rows is outermost here,
! so that we can OMP paralellise it (we can't parallelise the k-loop).
if ( l_mr_physics ) then
!$OMP DO SCHEDULE(STATIC)
  do j = 1, rows
    do i = 1, row_length
      rhodz = (r_rho_levels(i,j,2)-r_theta_levels(i,j,0)) * rho_dry_th(i,j,1)
      cclwp(i,j)  = rhodz * ccw(i,j,1)
      cclwp0(i,j) = rhodz * ccw0(i,j,1)
    end do
    do k = 2, n_conv_levels
      do i = 1, row_length
        rhodz = (r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k)) * rho_dry_th(i,j,k)
        cclwp(i,j)  = cclwp(i,j)  + rhodz * ccw(i,j,k)
        cclwp0(i,j) = cclwp0(i,j) + rhodz * ccw0(i,j,k)
      end do
    end do
  end do  ! j = 1, rows
!$OMP END DO NOWAIT
else  ! ( l_mr_physics )
!$OMP DO SCHEDULE(STATIC)
  do j = 1, rows
    do i = 1, row_length
      rhodz = (r_rho_levels(i,j,2)-r_theta_levels(i,j,0)) * rho_wet_th(i,j,1)
      cclwp(i,j)  = rhodz * ccw(i,j,1)
      cclwp0(i,j) = rhodz * ccw0(i,j,1)
    end do
    do k = 2, n_conv_levels
      do i = 1, row_length
        rhodz = (r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k)) * rho_wet_th(i,j,k)
        cclwp(i,j)  = cclwp(i,j)  + rhodz * ccw(i,j,k)
        cclwp0(i,j) = cclwp0(i,j) + rhodz * ccw0(i,j,k)
      end do
    end do
  end do  ! j = 1, rows
!$OMP END DO NOWAIT
end if  ! ( l_mr_physics )

!$OMP END PARALLEL


return
end subroutine comorph_conv_cloud_extras


end module comorph_conv_cloud_extras_mod
