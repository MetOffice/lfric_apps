!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module adj_opt_apply_variable_hx_kernel_mod

  use, intrinsic :: iso_fortran_env, only: real32, real64

  use argument_mod,      only : arg_type,              &
                                GH_FIELD, GH_SCALAR,   &
                                GH_OPERATOR, GH_REAL,  &
                                GH_READ, GH_READWRITE, &
                                GH_INC,                &
                                CELL_COLUMN
  use constants_mod,     only : i_def
  use fs_continuity_mod, only : W2, W3, Wtheta
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  type, public, extends(kernel_type) :: adj_opt_apply_variable_hx_kernel_type
    private
    type(arg_type) :: meta_args(11) = (/                               &
         arg_type(GH_FIELD,    GH_REAL, GH_READWRITE, W3),             &
         arg_type(GH_FIELD,    GH_REAL, GH_INC,       W2),             &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,      Wtheta),         &
         arg_type(GH_FIELD,    GH_REAL, GH_READWRITE, W3),             &
         arg_type(GH_OPERATOR, GH_REAL, GH_READ,      W3,     W2),     &
         arg_type(GH_OPERATOR, GH_REAL, GH_READ,      W3,     Wtheta), &
         arg_type(GH_OPERATOR, GH_REAL, GH_READ,      Wtheta, W2),     &
         arg_type(GH_OPERATOR, GH_REAL, GH_READ,      W3,     W3),     &
         arg_type(GH_SCALAR,   GH_REAL, GH_READ),                      &
         arg_type(GH_SCALAR,   GH_REAL, GH_READ),                      &
         arg_type(GH_FIELD,    GH_REAL, GH_READWRITE, W3)              &
         /)
    integer :: operates_on = CELL_COLUMN
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: adj_opt_apply_variable_hx_code

  interface adj_opt_apply_variable_hx_code
    module procedure &
      adj_opt_apply_variable_hx_code_real32, &
      adj_opt_apply_variable_hx_code_real64
  end interface


contains

!> @brief Adjoint of the application of the component of the helmholtz operator that maps
!>        from velocity space to the pressure space as well as the constant in space part.
!> @param[in] cell Horizontal cell index
!> @param[in] nlayers Number of layers
!> @param[in,out] lhs Pressure field with helmholtz operator applied to it
!> @param[in,out] x Gradient of the pressure field in the velocity space
!> @param[in] mt_inv Lumped inverse mass matrix for the temperature space
!> @param[in,out] pressure Field that helmholtz operator is being applied to
!> @param[in] ncell_3d_1 Total number of cells for divergence matrix
!> @param[in] div Generalised divergence matrix
!> @param[in] ncell_3d_2 Total number of cells for p3t matrix
!> @param[in] p3t Mapping from temperature space to pressure space
!> @param[in] ncell_3d_3 Total number of cells for pt2 matrix
!> @param[in] pt2 Mapping from velocity space to temperature space
!> @param[in] ncell_3d_4 Total number of cells for m3 matrix
!> @param[in] m3 Mass matrix for the pressure space
!> @param[in] const_div Scaling constant (sgn*tau_r*dt) applied to the divergence term
!> @param[in] const_theta Scaling constant (sgn*tau_t*dt) applied to the potential
!!            temperature term
!> @param[in,out] rhs_p Pressure field in w3
!> @param[in] ndf_w3 Number of degrees of freedom per cell for the pressure space
!> @param[in] undf_w3 Unique number of degrees of freedom  for the pressure space
!> @param[in] map_w3 Dofmap for the cell at the base of the column for the pressure space
!> @param[in] ndf_w2 Number of degrees of freedom per cell for the velocity space
!> @param[in] undf_w2 Unique number of degrees of freedom  for the velocity space
!> @param[in] map_w2 Dofmap for the cell at the base of the column for the velocity space
!> @param[in] ndf_wt Number of degrees of freedom per cell for the temperature space
!> @param[in] undf_wt Unique number of degrees of freedom  for the temperature space
!> @param[in] map_wt Dofmap for the cell at the base of the column for the temperature space

! REAL64 PRECISION
! ==================
subroutine adj_opt_apply_variable_hx_code_real64(cell,                    &
                                                 nlayers,                 &
                                                 lhs, x,                  &
                                                 mt_inv,                  &
                                                 pressure,                &
                                                 ncell_3d_1,              &
                                                 div,                     &
                                                 ncell_3d_2,              &
                                                 p3t,                     &
                                                 ncell_3d_3,              &
                                                 pt2,                     &
                                                 ncell_3d_4,              &
                                                 m3,                      &
                                                 const_div,               &
                                                 const_theta,             &
                                                 rhs_p,                   &
                                                 ndf_w3, undf_w3, map_w3, &
                                                 ndf_w2, undf_w2, map_w2, &
                                                 ndf_wt, undf_wt, map_wt)

  implicit none

  ! Arguments
  integer(kind=i_def),                    intent(in) :: cell, nlayers
  integer(kind=i_def),                    intent(in) :: ncell_3d_1, ncell_3d_2
  integer(kind=i_def),                    intent(in) :: ncell_3d_3, ncell_3d_4
  integer(kind=i_def),                    intent(in) :: undf_w2, ndf_w2
  integer(kind=i_def),                    intent(in) :: undf_w3, ndf_w3
  integer(kind=i_def),                    intent(in) :: undf_wt, ndf_wt
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt

  real(kind=real64), dimension(undf_w2), intent(inout) :: x
  real(kind=real64), dimension(undf_wt), intent(in)    :: mt_inv
  real(kind=real64), dimension(undf_w3), intent(inout) :: lhs
  real(kind=real64), dimension(undf_w3), intent(inout) :: pressure
  real(kind=real64),                     intent(in)    :: const_div
  real(kind=real64),                     intent(in)    :: const_theta
  real(kind=real64), dimension(undf_w3), intent(inout) :: rhs_p

  real(kind=real64), dimension(ncell_3d_1,1,6), intent(in) :: div
  real(kind=real64), dimension(ncell_3d_2,2,6), intent(in) :: pt2
  real(kind=real64), dimension(ncell_3d_3,1,2), intent(in) :: p3t
  real(kind=real64), dimension(ncell_3d_4,1,1), intent(in) :: m3

  ! Internal variables
  integer(kind=i_def)            :: k, ik
  real(kind=real64), dimension(2) :: t_e
  real(kind=real64)               :: div_u
  real(kind=real64), dimension(0:nlayers-1) :: t_e1_vec, t_e2_vec
  integer(kind=i_def) :: map_w21
  integer(kind=i_def) :: map_w22
  integer(kind=i_def) :: map_w23
  integer(kind=i_def) :: map_w24
  integer(kind=i_def) :: map_w25
  integer(kind=i_def) :: map_w26
  integer(kind=i_def) :: map_w31
  integer(kind=i_def) :: map_wt1
  integer(kind=i_def) :: map_wt2

  map_w21 = map_w2(1)
  map_w22 = map_w2(2)
  map_w23 = map_w2(3)
  map_w24 = map_w2(4)
  map_w25 = map_w2(5)
  map_w26 = map_w2(6)
  map_w31 = map_w3(1)
  map_wt1 = map_wt(1)
  map_wt2 = map_wt(2)

  k = nlayers-1
  ik = (cell-1)*nlayers + k + 1

  rhs_p(map_w3(1)+k) = rhs_p(map_w3(1)+k) + lhs(map_w31+k)
  t_e(1) = const_theta*lhs(map_w31+k)
  t_e(2) = const_theta*lhs(map_w31+k)
  div_u = const_div*lhs(map_w31+k)
  pressure(map_w31+k) = pressure(map_w31+k) + m3(ik,1,1)*lhs(map_w31+k)
  lhs(map_w31+k) = 0.0_real64

  x(map_w26+k) = x(map_w26+k) + div(ik,1,6)*div_u
  x(map_w25+k) = x(map_w25+k) + div(ik,1,5)*div_u
  x(map_w24+k) = x(map_w24+k) + div(ik,1,4)*div_u
  x(map_w23+k) = x(map_w23+k) + div(ik,1,3)*div_u
  x(map_w22+k) = x(map_w22+k) + div(ik,1,2)*div_u
  x(map_w21+k) = x(map_w21+k) + div(ik,1,1)*div_u

  x(map_w26+k) = x(map_w26+k)                   &
               + mt_inv(map_wt2+k)*p3t(ik,1,2)* &
                 pt2(ik,2,6)*t_e(2)
  x(map_w25+k) = x(map_w25+k)                   &
               + mt_inv(map_wt2+k)*p3t(ik,1,2)* &
                 pt2(ik,2,5)*t_e(2)

  x(map_w26+k-1) = x(map_w26+k-1)                 &
                 + mt_inv(map_wt1+k)*p3t(ik,1,1)* &
                   pt2(ik-1,2,6)*t_e(1)
  x(map_w25+k-1) = x(map_w25+k-1)                 &
                 + mt_inv(map_wt1+k)*p3t(ik,1,1)* &
                   pt2(ik-1,2,5)*t_e(1)
  x(map_w26+k) = x(map_w26+k)                     &
                 + mt_inv(map_wt1+k)*p3t(ik,1,1)* &
                   pt2(ik,1,6)*t_e(1)
  x(map_w25+k) = x(map_w25+k)                     &
                 + mt_inv(map_wt1+k)*p3t(ik,1,1)* &
                   pt2(ik,1,5)*t_e(1)

  do k = nlayers-2, 1, -1
    ik = (cell-1)*nlayers + k + 1

    rhs_p(map_w3(1)+k) = rhs_p(map_w3(1)+k) + lhs(map_w31+k)
    t_e1_vec(k) = const_theta*lhs(map_w31+k)
    t_e2_vec(k) = const_theta*lhs(map_w31+k)
    div_u = const_div*lhs(map_w31+k)
    pressure(map_w31+k) = pressure(map_w31+k) + m3(ik,1,1)*lhs(map_w31+k)
    lhs(map_w31+k) = 0.0_real64

    x(map_w26+k) = x(map_w26+k) + div(ik,1,6)*div_u
    x(map_w25+k) = x(map_w25+k) + div(ik,1,5)*div_u
    x(map_w24+k) = x(map_w24+k) + div(ik,1,4)*div_u
    x(map_w23+k) = x(map_w23+k) + div(ik,1,3)*div_u
    x(map_w22+k) = x(map_w22+k) + div(ik,1,2)*div_u
    x(map_w21+k) = x(map_w21+k) + div(ik,1,1)*div_u
  end do

  do k = nlayers-2, 1, -1
    ik = (cell-1)*nlayers + k + 1

    x(map_w26+k+1) = x(map_w26+k+1)                 &
                   + mt_inv(map_wt2+k)*p3t(ik,1,2)* &
                     pt2(ik+1,1,6)*t_e2_vec(k)
    x(map_w25+k+1) = x(map_w25+k+1)                 &
                   + mt_inv(map_wt2+k)*p3t(ik,1,2)* &
                     pt2(ik+1,1,5)*t_e2_vec(k)
    x(map_w26+k) = x(map_w26+k)                     &
                   + mt_inv(map_wt2+k)*p3t(ik,1,2)* &
                     pt2(ik,2,6)*t_e2_vec(k)
    x(map_w25+k) = x(map_w25+k)                     &
                   + mt_inv(map_wt2+k)*p3t(ik,1,2)* &
                     pt2(ik,2,5)*t_e2_vec(k)

    x(map_w26+k-1) = x(map_w26+k-1)                 &
                   + mt_inv(map_wt1+k)*p3t(ik,1,1)* &
                     pt2(ik-1,2,6)*t_e1_vec(k)
    x(map_w25+k-1) = x(map_w25+k-1)                 &
                   + mt_inv(map_wt1+k)*p3t(ik,1,1)* &
                     pt2(ik-1,2,5)*t_e1_vec(k)
    x(map_w26+k) = x(map_w26+k)                     &
                   + mt_inv(map_wt1+k)*p3t(ik,1,1)* &
                     pt2(ik-1,1,6)*t_e1_vec(k)
    x(map_w25+k) = x(map_w25+k)                     &
                   + mt_inv(map_wt1+k)*p3t(ik,1,1)* &
                     pt2(ik-1,1,5)*t_e1_vec(k)
  end do

  k = 0
  ik = (cell-1)*nlayers + k + 1

  rhs_p(map_w3(1)+k) = rhs_p(map_w3(1)+k) + lhs(map_w31+k)
  t_e(1) = const_theta*lhs(map_w31+k)
  t_e(2) = const_theta*lhs(map_w31+k)
  div_u = const_div*lhs(map_w31+k)
  pressure(map_w31+k) = pressure(map_w31+k) + m3(ik,1,1)*lhs(map_w31+k)
  lhs(map_w31+k) = 0.0_real64

  x(map_w26+k) = x(map_w26+k) + div(ik,1,6)*div_u
  x(map_w25+k) = x(map_w25+k) + div(ik,1,5)*div_u
  x(map_w24+k) = x(map_w24+k) + div(ik,1,4)*div_u
  x(map_w23+k) = x(map_w23+k) + div(ik,1,3)*div_u
  x(map_w22+k) = x(map_w22+k) + div(ik,1,2)*div_u
  x(map_w21+k) = x(map_w21+k) + div(ik,1,1)*div_u

  x(map_w26+k+1) = x(map_w26+k+1)                 &
                 + mt_inv(map_wt2+k)*p3t(ik,1,2)* &
                   pt2(ik+1,1,6)*t_e(2)
  x(map_w25+k+1) = x(map_w25+k+1)                 &
                 + mt_inv(map_wt2+k)*p3t(ik,1,2)* &
                   pt2(ik+1,1,5)*t_e(2)
  x(map_w26+k) = x(map_w26+k)                     &
                 + mt_inv(map_wt2+k)*p3t(ik,1,2)* &
                   pt2(ik+1,2,6)*t_e(2)
  x(map_w25+k) = x(map_w25+k)                     &
                 + mt_inv(map_wt2+k)*p3t(ik,1,2)* &
                   pt2(ik+1,2,5)*t_e(2)

  x(map_w26+k) = x(map_w26+k)                   &
               + mt_inv(map_wt1+k)*p3t(ik,1,1)* &
                 pt2(ik,1,6)*t_e(1)
  x(map_w25+k) = x(map_w25+k)                   &
               + mt_inv(map_wt1+k)*p3t(ik,1,1)* &
                 pt2(ik,1,5)*t_e(1)

end subroutine adj_opt_apply_variable_hx_code_real64

! REAL32 PRECISION
! ==================
subroutine adj_opt_apply_variable_hx_code_real32(cell,                    &
                                               nlayers,                 &
                                               lhs, x,                  &
                                               mt_inv,                  &
                                               pressure,                &
                                               ncell_3d_1,              &
                                               div,                     &
                                               ncell_3d_2,              &
                                               p3t,                     &
                                               ncell_3d_3,              &
                                               pt2,                     &
                                               ncell_3d_4,              &
                                               m3,                      &
                                               const_div,               &
                                               const_theta,             &
                                               rhs_p,                   &
                                               ndf_w3, undf_w3, map_w3, &
                                               ndf_w2, undf_w2, map_w2, &
                                               ndf_wt, undf_wt, map_wt)

  implicit none

  ! Arguments
  integer(kind=i_def),                    intent(in) :: cell, nlayers
  integer(kind=i_def),                    intent(in) :: ncell_3d_1, ncell_3d_2
  integer(kind=i_def),                    intent(in) :: ncell_3d_3, ncell_3d_4
  integer(kind=i_def),                    intent(in) :: undf_w2, ndf_w2
  integer(kind=i_def),                    intent(in) :: undf_w3, ndf_w3
  integer(kind=i_def),                    intent(in) :: undf_wt, ndf_wt
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt

  real(kind=real32), dimension(undf_w2), intent(inout) :: x
  real(kind=real32), dimension(undf_wt), intent(in)    :: mt_inv
  real(kind=real32), dimension(undf_w3), intent(inout) :: lhs
  real(kind=real32), dimension(undf_w3), intent(inout) :: pressure
  real(kind=real32),                     intent(in)    :: const_div
  real(kind=real32),                     intent(in)    :: const_theta
  real(kind=real32), dimension(undf_w3), intent(inout) :: rhs_p

  real(kind=real32), dimension(ncell_3d_1,1,6), intent(in) :: div
  real(kind=real32), dimension(ncell_3d_2,2,6), intent(in) :: pt2
  real(kind=real32), dimension(ncell_3d_3,1,2), intent(in) :: p3t
  real(kind=real32), dimension(ncell_3d_4,1,1), intent(in) :: m3

  ! Internal variables
  integer(kind=i_def)            :: k, ik
  real(kind=real32), dimension(2) :: t_e
  real(kind=real32)               :: div_u
  real(kind=real32), dimension(0:nlayers-1) :: t_e1_vec, t_e2_vec
  integer(kind=i_def) :: map_w21
  integer(kind=i_def) :: map_w22
  integer(kind=i_def) :: map_w23
  integer(kind=i_def) :: map_w24
  integer(kind=i_def) :: map_w25
  integer(kind=i_def) :: map_w26
  integer(kind=i_def) :: map_w31
  integer(kind=i_def) :: map_wt1
  integer(kind=i_def) :: map_wt2

  map_w21 = map_w2(1)
  map_w22 = map_w2(2)
  map_w23 = map_w2(3)
  map_w24 = map_w2(4)
  map_w25 = map_w2(5)
  map_w26 = map_w2(6)
  map_w31 = map_w3(1)
  map_wt1 = map_wt(1)
  map_wt2 = map_wt(2)

  k = nlayers-1
  ik = (cell-1)*nlayers + k + 1

  rhs_p(map_w3(1)+k) = rhs_p(map_w3(1)+k) + lhs(map_w31+k)
  t_e(1) = const_theta*lhs(map_w31+k)
  t_e(2) = const_theta*lhs(map_w31+k)
  div_u = const_div*lhs(map_w31+k)
  pressure(map_w31+k) = pressure(map_w31+k) + m3(ik,1,1)*lhs(map_w31+k)
  lhs(map_w31+k) = 0.0_real32

  x(map_w26+k) = x(map_w26+k) + div(ik,1,6)*div_u
  x(map_w25+k) = x(map_w25+k) + div(ik,1,5)*div_u
  x(map_w24+k) = x(map_w24+k) + div(ik,1,4)*div_u
  x(map_w23+k) = x(map_w23+k) + div(ik,1,3)*div_u
  x(map_w22+k) = x(map_w22+k) + div(ik,1,2)*div_u
  x(map_w21+k) = x(map_w21+k) + div(ik,1,1)*div_u

  x(map_w26+k) = x(map_w26+k)                   &
               + mt_inv(map_wt2+k)*p3t(ik,1,2)* &
                 pt2(ik,2,6)*t_e(2)
  x(map_w25+k) = x(map_w25+k)                   &
               + mt_inv(map_wt2+k)*p3t(ik,1,2)* &
                 pt2(ik,2,5)*t_e(2)

  x(map_w26+k-1) = x(map_w26+k-1)                 &
                 + mt_inv(map_wt1+k)*p3t(ik,1,1)* &
                   pt2(ik-1,2,6)*t_e(1)
  x(map_w25+k-1) = x(map_w25+k-1)                 &
                 + mt_inv(map_wt1+k)*p3t(ik,1,1)* &
                   pt2(ik-1,2,5)*t_e(1)
  x(map_w26+k) = x(map_w26+k)                     &
                 + mt_inv(map_wt1+k)*p3t(ik,1,1)* &
                   pt2(ik,1,6)*t_e(1)
  x(map_w25+k) = x(map_w25+k)                     &
                 + mt_inv(map_wt1+k)*p3t(ik,1,1)* &
                   pt2(ik,1,5)*t_e(1)

  do k = nlayers-2, 1, -1
    ik = (cell-1)*nlayers + k + 1

    rhs_p(map_w3(1)+k) = rhs_p(map_w3(1)+k) + lhs(map_w31+k)
    t_e1_vec(k) = const_theta*lhs(map_w31+k)
    t_e2_vec(k) = const_theta*lhs(map_w31+k)
    div_u = const_div*lhs(map_w31+k)
    pressure(map_w31+k) = pressure(map_w31+k) + m3(ik,1,1)*lhs(map_w31+k)
    lhs(map_w31+k) = 0.0_real32

    x(map_w26+k) = x(map_w26+k) + div(ik,1,6)*div_u
    x(map_w25+k) = x(map_w25+k) + div(ik,1,5)*div_u
    x(map_w24+k) = x(map_w24+k) + div(ik,1,4)*div_u
    x(map_w23+k) = x(map_w23+k) + div(ik,1,3)*div_u
    x(map_w22+k) = x(map_w22+k) + div(ik,1,2)*div_u
    x(map_w21+k) = x(map_w21+k) + div(ik,1,1)*div_u
  end do

  do k = nlayers-2, 1, -1
    ik = (cell-1)*nlayers + k + 1

    x(map_w26+k+1) = x(map_w26+k+1)                 &
                   + mt_inv(map_wt2+k)*p3t(ik,1,2)* &
                     pt2(ik+1,1,6)*t_e2_vec(k)
    x(map_w25+k+1) = x(map_w25+k+1)                 &
                   + mt_inv(map_wt2+k)*p3t(ik,1,2)* &
                     pt2(ik+1,1,5)*t_e2_vec(k)
    x(map_w26+k) = x(map_w26+k)                     &
                   + mt_inv(map_wt2+k)*p3t(ik,1,2)* &
                     pt2(ik,2,6)*t_e2_vec(k)
    x(map_w25+k) = x(map_w25+k)                     &
                   + mt_inv(map_wt2+k)*p3t(ik,1,2)* &
                     pt2(ik,2,5)*t_e2_vec(k)

    x(map_w26+k-1) = x(map_w26+k-1)                 &
                   + mt_inv(map_wt1+k)*p3t(ik,1,1)* &
                     pt2(ik-1,2,6)*t_e1_vec(k)
    x(map_w25+k-1) = x(map_w25+k-1)                 &
                   + mt_inv(map_wt1+k)*p3t(ik,1,1)* &
                     pt2(ik-1,2,5)*t_e1_vec(k)
    x(map_w26+k) = x(map_w26+k)                     &
                   + mt_inv(map_wt1+k)*p3t(ik,1,1)* &
                     pt2(ik-1,1,6)*t_e1_vec(k)
    x(map_w25+k) = x(map_w25+k)                     &
                   + mt_inv(map_wt1+k)*p3t(ik,1,1)* &
                     pt2(ik-1,1,5)*t_e1_vec(k)
  end do

  k = 0
  ik = (cell-1)*nlayers + k + 1

  rhs_p(map_w3(1)+k) = rhs_p(map_w3(1)+k) + lhs(map_w31+k)
  t_e(1) = const_theta*lhs(map_w31+k)
  t_e(2) = const_theta*lhs(map_w31+k)
  div_u = const_div*lhs(map_w31+k)
  pressure(map_w31+k) = pressure(map_w31+k) + m3(ik,1,1)*lhs(map_w31+k)
  lhs(map_w31+k) = 0.0_real32

  x(map_w26+k) = x(map_w26+k) + div(ik,1,6)*div_u
  x(map_w25+k) = x(map_w25+k) + div(ik,1,5)*div_u
  x(map_w24+k) = x(map_w24+k) + div(ik,1,4)*div_u
  x(map_w23+k) = x(map_w23+k) + div(ik,1,3)*div_u
  x(map_w22+k) = x(map_w22+k) + div(ik,1,2)*div_u
  x(map_w21+k) = x(map_w21+k) + div(ik,1,1)*div_u

  x(map_w26+k+1) = x(map_w26+k+1)                 &
                 + mt_inv(map_wt2+k)*p3t(ik,1,2)* &
                   pt2(ik+1,1,6)*t_e(2)
  x(map_w25+k+1) = x(map_w25+k+1)                 &
                 + mt_inv(map_wt2+k)*p3t(ik,1,2)* &
                   pt2(ik+1,1,5)*t_e(2)
  x(map_w26+k) = x(map_w26+k)                     &
                 + mt_inv(map_wt2+k)*p3t(ik,1,2)* &
                   pt2(ik+1,2,6)*t_e(2)
  x(map_w25+k) = x(map_w25+k)                     &
                 + mt_inv(map_wt2+k)*p3t(ik,1,2)* &
                   pt2(ik+1,2,5)*t_e(2)

  x(map_w26+k) = x(map_w26+k)                   &
               + mt_inv(map_wt1+k)*p3t(ik,1,1)* &
                 pt2(ik,1,6)*t_e(1)
  x(map_w25+k) = x(map_w25+k)                   &
               + mt_inv(map_wt1+k)*p3t(ik,1,1)* &
                 pt2(ik,1,5)*t_e(1)

end subroutine adj_opt_apply_variable_hx_code_real32

end module adj_opt_apply_variable_hx_kernel_mod
