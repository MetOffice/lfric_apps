! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_layer_mass_mod

implicit none

contains


! Subroutine to calculate the dry-mass per unit surface area on
! each model-level
subroutine calc_layer_mass( lb_rho, ub_rho, rho_dry, lb_rs, ub_rs, r_surf,     &
                            lb_h, ub_h, height_half,                           &
                            layer_mass )

use comorph_constants_mod, only: real_hmprec,                                  &
                     nx_full, ny_full, k_bot_conv, k_top_conv,                 &
                     l_spherical_coord

implicit none

! Note: the array arguments below may have halos, and in here
! we don't care about the halos.
! We need to pass in the bounds of each array to use in
! the declarations, to ensure the do loops pick out the same
! indices of the arrays regardless of whether or not they
! have halos.  The lower and upper bounds are passed in through the
! argument list in the lb_* and ub_* integer arrays.  Those storing the
! bounds for 3D arrays must have 3 elements; one for each dimension
! of the array (and similarly 2 elements for 2D arrays).

! Dry density on thermodynamic model-layers
integer, intent(in) :: lb_rho(3)
integer, intent(in) :: ub_rho(3)
real(kind=real_hmprec), intent(in) :: rho_dry                                  &
             ( lb_rho(1):ub_rho(1), lb_rho(2):ub_rho(2), lb_rho(3):ub_rho(3) )

! Height of surface above the centre of the Earth
integer, intent(in) :: lb_rs(2)
integer, intent(in) :: ub_rs(2)
real(kind=real_hmprec), intent(in) :: r_surf                                   &
                                      ( lb_rs(1):ub_rs(1), lb_rs(2):ub_rs(2) )

! Height of layer interfaces above surface
integer, intent(in) :: lb_h(3)
integer, intent(in) :: ub_h(3)
real(kind=real_hmprec), intent(in) :: height_half                              &
                         ( lb_h(1):ub_h(1), lb_h(2):ub_h(2), lb_h(3):ub_h(3) )

! Output grid-cell dry-mass per unit surface area
real(kind=real_hmprec), intent(out) :: layer_mass                              &
                     ( nx_full, ny_full, k_bot_conv:k_top_conv )

! Work arrays for spherical non-shallow atmosphere case
real(kind=real_hmprec) :: z1_plus_z2 ( nx_full, ny_full )
real(kind=real_hmprec) :: rrs_sq     ( nx_full, ny_full )

real(kind=real_hmprec), parameter :: one = 1.0_real_hmprec
real(kind=real_hmprec), parameter :: third = 1.0_real_hmprec                   &
                                           / 3.0_real_hmprec

! Loop counters
integer :: i, j, k


!$OMP PARALLEL DEFAULT(none) private( i, j, k, z1_plus_z2 )                    &
!$OMP SHARED( nx_full, ny_full, k_bot_conv, k_top_conv,                        &
!$OMP         layer_mass, height_half, l_spherical_coord, rrs_sq, r_surf,      &
!$OMP         rho_dry )

! Calculate dz
!$OMP do SCHEDULE(STATIC)
do k = k_bot_conv, k_top_conv
  do j = 1, ny_full
    do i = 1, nx_full
      layer_mass(i,j,k) = height_half(i,j,k+1)                                 &
                        - height_half(i,j,k)
    end do
  end do
end do
!$OMP end do


! Next part depends on whether we're doing this in
! spherical coordinates without the shallow atmosphere
! approximation...
if ( l_spherical_coord ) then

  ! Need to account for spherical geometry in order to get
  ! conservation right.
  ! What we want is the volume of the grid-cell per unit
  ! surface area.  In the Cartesian case this is just dz
  ! (which we already calculated above),
  ! but in the spherical case it isn't, because grid-cells
  ! get wider with height.
  !
  ! Volume per unit surface area is:
  ! V(k) = 1/3 ( r(k+1/2)^3 - r(k-1/2)^3 ) / r_surf^2
  !
  ! where r is height above the centre of the Earth;
  ! r(k) = r_surf + z(k)
  ! where r_surf is the height of the surface above the
  ! centre of the Earth
  !
  ! => V(k) = 1/3 ( ( r_surf + z(k+1/2) )^3
  !               - ( r_surf + z(k-1/2) )^3 ) / r_surf^2
  !
  ! This formula would likely not have sufficient precision
  ! to get the right answer if running in 32-bit.
  ! Therefore, expand out the brackets and rearrange...
  !
  ! => V(k) = ( z(k+1/2) - z(k-1/2) )
  !         + ( z(k+1/2)^2 - z(k-1/2)^2 ) / r_surf
  !     + 1/3 ( z(k+1/2)^3 - z(k-1/2)^3 ) / r_surf^2
  !
  ! This factorises to give:
  !
  ! => V(k) = ( z(k+1/2) - z(k-1/2) )
  !           ( 1 + ( z(k+1/2) + z(k-1/2) ) / r_surf
  !               + 1/3 ( ( z(k+1/2) + z(k-1/2) )^2
  !                     - z(k+1/2) z(k-1/2)         )
  !                   / r_surf^2 )
  !
  ! This gives us a scaling factor to apply to dz
  ! (which is already stored in the layer_mass array),
  ! which should avoid loss of precision when
  ! dz / r_surf < EPS.

  ! Precalculate 1 / r_surf^2 term
!$OMP do SCHEDULE(STATIC)
  do j = 1, ny_full
    do i = 1, nx_full
      rrs_sq(i,j) = one / ( r_surf(i,j) * r_surf(i,j) )
    end do
  end do
!$OMP end do

!$OMP do SCHEDULE(STATIC)
  do k = k_bot_conv, k_top_conv
    ! Precalculate z(k-1/2) + z(k+1/2)
    do j = 1, ny_full
      do i = 1, nx_full
        z1_plus_z2(i,j) = height_half(i,j,k)                                   &
                        + height_half(i,j,k+1)
      end do
    end do
    ! Compute rest of the grid volume / surface area formula
    do j = 1, ny_full
      do i = 1, nx_full
        layer_mass(i,j,k) = layer_mass(i,j,k)                                  &
          * ( one + z1_plus_z2(i,j) / r_surf(i,j)                              &
                  + third * ( ( z1_plus_z2(i,j)                                &
                              * z1_plus_z2(i,j) )                              &
                            - ( height_half(i,j,k)                             &
                              * height_half(i,j,k+1) ) )                       &
                          * rrs_sq(i,j)  )
      end do
    end do
  end do
!$OMP end do

end if  ! ( l_spherical_coord )


! Scale by density to get dry-mass on layer
!$OMP do SCHEDULE(STATIC)
do k = k_bot_conv, k_top_conv
  do j = 1, ny_full
    do i = 1, nx_full
      layer_mass(i,j,k) = layer_mass(i,j,k) * rho_dry(i,j,k)
    end do
  end do
end do
!$OMP end do

!$OMP end PARALLEL


return
end subroutine calc_layer_mass


end module calc_layer_mass_mod
