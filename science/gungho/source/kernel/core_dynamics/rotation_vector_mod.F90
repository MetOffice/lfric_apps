!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel computes the omega (rotation) for the u-equation.

!> @details The kernel computes the element rotation (omega) vector in W2 space for
!> the RHS of momentum equation on both F-PLANE and the SPHERE

module rotation_vector_mod

use constants_mod,     only: r_def, i_def
use planet_config_mod, only: scaled_omega

implicit none

private

public :: rotation_vector_fplane
public :: rotation_vector_sphere
public :: vert_vector_sphere

contains
!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
!> Subroutine Computes the element rotation (omega) vector in W2 space for the
!! RHS of momentum equation on the F-PLANE.
!! @param[in] ngp_h          Number of quadrature points in horizontal direction
!! @param[in] ngp_v          Number of quadrature points in vertical direction
!! @param[in] omegaf         Planetary rotation rate
!! @param[in] latitude       Pre-set value for the f-plane.
!! @param[out] rotation_vec  Rotation vector on quadrature points
subroutine rotation_vector_fplane(ngp_h, ngp_v, omegaf, latitude, rotation_vec)
!-------------------------------------------------------------------------------
! Compute the rotation vector Omega = (0, 2*cos(lat), 2*sin(lat)) on quadrature points
!-------------------------------------------------------------------------------

implicit none

integer(kind=i_def), intent(in)  :: ngp_h, ngp_v
real(kind=r_def),    intent(in)  :: omegaf, latitude
real(kind=r_def),    intent(out) :: rotation_vec(3,ngp_h,ngp_v)

integer(kind=i_def) :: i, j

rotation_vec(:,:,:) = 0.0_r_def

do j = 1, ngp_v
  do i = 1, ngp_h
     rotation_vec(1,i,j) = 0.0_r_def
     rotation_vec(2,i,j) = 2.0_r_def*omegaf*cos(latitude)
     rotation_vec(3,i,j) = 2.0_r_def*omegaf*sin(latitude)
  end do
end do

end subroutine rotation_vector_fplane

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
!> Subroutine Computes the element rotation (omega) vector in W2 space for the
!! RHS of momentum equation on the SPHERE. The coordinate field is in Wchi.
!! @param[in] ndf_chi        The size of the chi arrays
!! @param[in] ngp_h          The number of quadrature points in horizontal direction
!! @param[in] ngp_v          The number of quadrature points in vertical direction
!! @param[in] chi_1          Holds the chi_1 coordinate field
!! @param[in] chi_2          Holds the chi_2 coordinate field
!! @param[in] chi_3          Holds the chi_3 coordinate field
!! @param[in] panel_id       ID of mesh panel
!! @param[in] chi_basis      Holds the chi basis functions
!! @param[out] rotation_vec  Holds the values of the rotation vector on quadrature points
subroutine rotation_vector_sphere(ndf_chi, ngp_h, ngp_v, chi_1, chi_2, chi_3, &
                                  panel_id, chi_basis, rotation_vec)
!-------------------------------------------------------------------------------
! Compute the rotation vector Omega = (0, 2*cos(lat), 2*sin(lat)) on quadrature points
!-------------------------------------------------------------------------------

use sci_chi_transform_mod,   only: chi2llr
use coord_transform_mod,     only: sphere2cart_vector

implicit none

integer(kind=i_def), intent(in)  :: ndf_chi, ngp_h, ngp_v, panel_id
real(kind=r_def),    intent(in)  :: chi_1(ndf_chi), chi_2(ndf_chi), chi_3(ndf_chi)
real(kind=r_def),    intent(out) :: rotation_vec(3,ngp_h,ngp_v)
real(kind=r_def),    intent(in), dimension(1,ndf_chi,ngp_h,ngp_v) :: chi_basis


integer(kind=i_def) :: i, j, df
real(kind=r_def)    :: lat, long, r
real(kind=r_def)    :: llr(3), coords(3)

lat = 0.0_r_def
long = 0.0_r_def
rotation_vec(:,:,:) = 0.0_r_def

do j = 1, ngp_v
  do i = 1, ngp_h
    ! Calculate the position vector at this quadrature point
    coords(:) = 0.0_r_def
    do df = 1, ndf_chi
      coords(1) = coords(1) + chi_1(df)*chi_basis(1,df,i,j)
      coords(2) = coords(2) + chi_2(df)*chi_basis(1,df,i,j)
      coords(3) = coords(3) + chi_3(df)*chi_basis(1,df,i,j)
    end do

    ! Need to obtain longitude, latitude and radius from position vector
    call chi2llr(coords(1), coords(2), coords(3), panel_id, long, lat, r)

    ! Get (long,lat,r) components of planet rotation vector
    rotation_vec(1,i,j) = 0.0_r_def
    rotation_vec(2,i,j) = 2.0_r_def*scaled_omega*cos(lat)
    rotation_vec(3,i,j) = 2.0_r_def*scaled_omega*sin(lat)

    ! Obtain (X,Y,Z) components of rotation vector
    llr = (/long, lat, r/)
    rotation_vec(:,i,j) = sphere2cart_vector( rotation_vec(:,i,j), llr )

  end do
end do

end subroutine rotation_vector_sphere

!> Subroutine Computes the vertical unit vector
!! @param[in] ndf_chi    The size of the chi arrays
!! @param[in] ngp_h      The number of quadrature points in horizontal direction
!! @param[in] ngp_v      The number of quadrature points in vertical direction
!! @param[in] chi_1      Holds the chi_1 coordinate field
!! @param[in] chi_2      Holds the chi_2 coordinate field
!! @param[in] chi_3      Holds the chi_3 coordinate field
!! @param[in] panel_id   ID of mesh panel
!! @param[in] chi_basis  Holds the chi basis functions
!! @param[out] vert_vec  Holds the values of the vertical vector on quadrature points
subroutine vert_vector_sphere(ndf_chi, ngp_h, ngp_v, chi_1, chi_2, chi_3, &
                                  panel_id, chi_basis, vert_vec)
!-------------------------------------------------------------------------------
! Compute the vertical vector Omega on quadrature points
!-------------------------------------------------------------------------------

use sci_chi_transform_mod,   only: chi2llr
use coord_transform_mod,     only: sphere2cart_vector

implicit none

integer(kind=i_def), intent(in)  :: ndf_chi, ngp_h, ngp_v, panel_id
real(kind=r_def),    intent(in)  :: chi_1(ndf_chi), chi_2(ndf_chi), chi_3(ndf_chi)
real(kind=r_def),    intent(out) :: vert_vec(3,ngp_h,ngp_v)
real(kind=r_def),    intent(in), dimension(1,ndf_chi,ngp_h,ngp_v) :: chi_basis


integer(kind=i_def) :: i, j, df
real(kind=r_def)    :: lat, long, r
real(kind=r_def)    :: llr(3), coords(3)

lat = 0.0_r_def
long = 0.0_r_def
vert_vec(:,:,:) = 0.0_r_def

do j = 1, ngp_v
  do i = 1, ngp_h
    ! Calculate the position vector at this quadrature point
    coords(:) = 0.0_r_def
    do df = 1, ndf_chi
      coords(1) = coords(1) + chi_1(df)*chi_basis(1,df,i,j)
      coords(2) = coords(2) + chi_2(df)*chi_basis(1,df,i,j)
      coords(3) = coords(3) + chi_3(df)*chi_basis(1,df,i,j)
    end do

    ! Need to obtain longitude, latitude and radius from position vector
    call chi2llr(coords(1), coords(2), coords(3), panel_id, long, lat, r)

    ! Get (long,lat,r) components of planet rotation vector
    vert_vec(1,i,j) = 0.0_r_def
    vert_vec(2,i,j) = 0.0_r_def
    vert_vec(3,i,j) = 1.0_r_def

    ! Obtain (X,Y,Z) components of rotation vector
    llr = (/long, lat, r/)
    vert_vec(:,i,j) = sphere2cart_vector( vert_vec(:,i,j), llr )

  end do
end do

end subroutine vert_vector_sphere

end module rotation_vector_mod
