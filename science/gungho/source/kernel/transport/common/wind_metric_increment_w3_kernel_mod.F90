!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!
! Some of the content of this file has been produced with the assistance of
! Met Office Github Copilot Enterprise
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!> @brief Computes metric term increments to a transported wind at W3 points.
!> @details The wind components are normalised contravariant components in the
!!          native coordinates (lon-lat-r or alpha-beta-r) and live on W2H and
!!          Wtheta points. They are averaged to W3 points together with the
!!          transporting wind, and the change in the transported components
!!          due to the rotation of the coordinate frame over a time dt_metric
!!          is computed with a forward Euler step. The horizontal coordinates
!!          are taken from the lowest layer only.
!-----------------------------------------------------------------------------
module wind_metric_increment_w3_kernel_mod

  use argument_mod,      only: arg_type, func_type,       &
                               GH_FIELD, GH_REAL,         &
                               GH_SCALAR, GH_LOGICAL,     &
                               GH_READ, GH_WRITE,         &
                               CELL_COLUMN,               &
                               ANY_SPACE_9,               &
                               GH_BASIS, GH_DIFF_BASIS,   &
                               GH_EVALUATOR
  use constants_mod,     only: r_def, r_tran, i_def, l_def
  use fs_continuity_mod, only: W3, W2, W2h, Wtheta
  use kernel_mod,        only: kernel_type
  use reference_element_mod, &
                         only: W, S, E, N

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  type, public, extends(kernel_type) :: wind_metric_increment_w3_kernel_type
    private
    type(arg_type) :: meta_args(13) = (/                                      &
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W3),                       & ! d_uv1_w3
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W3),                       & ! d_uv2_w3
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W3),                       & ! d_w_w3
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2h),                      & ! uv_w2h
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                   & ! w_wt
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2),                       & ! transporting wind
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2),                       & ! dA
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_SPACE_9),              & ! chi1
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_SPACE_9),              & ! chi2
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),                       & ! height_w3
         arg_type(GH_SCALAR, GH_REAL,    GH_READ),                            & ! dt_metric
         arg_type(GH_SCALAR, GH_REAL,    GH_READ),                            & ! scaled_radius
         arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                             & ! cubed_sphere
         /)
    type(func_type) :: meta_funcs(1) = (/                                     &
         func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                      &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
    integer :: gh_evaluator_targets(1) = (/ W3 /)
  contains
    procedure, nopass :: wind_metric_increment_w3_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: wind_metric_increment_w3_code

contains

!> @brief Compute metric term increments to the wind components at W3 points.
!! @param[in]     nlayers           Number of layers
!! @param[in,out] d_uv1_w3          Increment to component along xi1 in W3
!! @param[in,out] d_uv2_w3          Increment to component along xi2 in W3
!! @param[in,out] d_w_w3            Increment to vertical component in W3
!! @param[in]     uv_w2h            Transported horizontal wind in W2H
!! @param[in]     w_wt              Transported vertical wind in Wtheta
!! @param[in]     transporting_wind Transporting wind in W2
!! @param[in]     da                Area element in W2
!! @param[in]     chi1              First coordinate field in Wchi
!! @param[in]     chi2              Second coordinate field in Wchi
!! @param[in]     height_w3         Height above the surface radius in W3
!! @param[in]     dt_metric         Time over which the metric terms act
!! @param[in]     scaled_radius     Radius of the planet
!! @param[in]     cubed_sphere      True for alpha-beta coordinates, false for
!!                                  lon-lat coordinates
!! @param[in]     ndf_w3            Num DoFs per cell for W3
!! @param[in]     undf_w3           Num DoFs in this partition for W3
!! @param[in]     map_w3            Map of base cell DoFs for W3
!! @param[in]     ndf_w2h           Num DoFs per cell for W2H
!! @param[in]     undf_w2h          Num DoFs in this partition for W2H
!! @param[in]     map_w2h           Map of base cell DoFs for W2H
!! @param[in]     ndf_wt            Num DoFs per cell for Wtheta
!! @param[in]     undf_wt           Num DoFs in this partition for Wtheta
!! @param[in]     map_wt            Map of base cell DoFs for Wtheta
!! @param[in]     ndf_w2            Num DoFs per cell for W2
!! @param[in]     undf_w2           Num DoFs in this partition for W2
!! @param[in]     map_w2            Map of base cell DoFs for W2
!! @param[in]     ndf_wx            Num DoFs per cell for Wchi
!! @param[in]     undf_wx           Num DoFs in this partition for Wchi
!! @param[in]     map_wx            Map of base cell DoFs for Wchi
!! @param[in]     basis_wx          Basis functions for Wchi at W3 nodes
!! @param[in]     diff_basis_wx     Differential basis functions for Wchi at
!!                                  W3 nodes
subroutine wind_metric_increment_w3_code(                                    &
    nlayers,                                                                 &
    d_uv1_w3, d_uv2_w3, d_w_w3,                                              &
    uv_w2h, w_wt, transporting_wind, da,                                     &
    chi1, chi2, height_w3,                                                   &
    dt_metric, scaled_radius, cubed_sphere,                                  &
    ndf_w3, undf_w3, map_w3,                                                 &
    ndf_w2h, undf_w2h, map_w2h,                                              &
    ndf_wt, undf_wt, map_wt,                                                 &
    ndf_w2, undf_w2, map_w2,                                                 &
    ndf_wx, undf_wx, map_wx, basis_wx, diff_basis_wx                         &
)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
  integer(kind=i_def), intent(in) :: ndf_w2h, undf_w2h
  integer(kind=i_def), intent(in) :: ndf_wt, undf_wt
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
  integer(kind=i_def), intent(in) :: ndf_wx, undf_wx

  real(kind=r_def),    intent(inout) :: d_uv1_w3(undf_w3)
  real(kind=r_def),    intent(inout) :: d_uv2_w3(undf_w3)
  real(kind=r_def),    intent(inout) :: d_w_w3(undf_w3)
  real(kind=r_def),    intent(in)    :: uv_w2h(undf_w2h)
  real(kind=r_def),    intent(in)    :: w_wt(undf_wt)
  real(kind=r_tran),   intent(in)    :: transporting_wind(undf_w2)
  real(kind=r_def),    intent(in)    :: da(undf_w2)
  real(kind=r_def),    intent(in)    :: chi1(undf_wx)
  real(kind=r_def),    intent(in)    :: chi2(undf_wx)
  real(kind=r_def),    intent(in)    :: height_w3(undf_w3)
  real(kind=r_def),    intent(in)    :: dt_metric
  real(kind=r_def),    intent(in)    :: scaled_radius
  logical(kind=l_def), intent(in)    :: cubed_sphere

  integer(kind=i_def), intent(in) :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)
  integer(kind=i_def), intent(in) :: map_wt(ndf_wt)
  integer(kind=i_def), intent(in) :: map_w2(ndf_w2)
  integer(kind=i_def), intent(in) :: map_wx(ndf_wx)

  real(kind=r_def), intent(in) :: basis_wx(1,ndf_wx,ndf_w3)
  real(kind=r_def), intent(in) :: diff_basis_wx(3,ndf_wx,ndf_w3)

  ! Internal variables
  integer(kind=i_def) :: k, df
  real(kind=r_def)    :: c1, c2, d11, d12, d21, d22
  real(kind=r_def)    :: tan_lat
  real(kind=r_def)    :: ca, cb, x, y, p, q, delta
  real(kind=r_def)    :: c11, c12, c21, c22, c_diag, c_off
  real(kind=r_def)    :: u1, u2, a1, a2, a3, scale

  ! Horizontal coordinates and their derivatives with respect to the
  ! reference coordinates, at the cell centre in the lowest layer
  c1  = 0.0_r_def
  c2  = 0.0_r_def
  d11 = 0.0_r_def
  d12 = 0.0_r_def
  d21 = 0.0_r_def
  d22 = 0.0_r_def
  do df = 1, ndf_wx
    c1  = c1  + chi1(map_wx(df))*basis_wx(1,df,1)
    c2  = c2  + chi2(map_wx(df))*basis_wx(1,df,1)
    d11 = d11 + chi1(map_wx(df))*diff_basis_wx(1,df,1)
    d12 = d12 + chi1(map_wx(df))*diff_basis_wx(2,df,1)
    d21 = d21 + chi2(map_wx(df))*diff_basis_wx(1,df,1)
    d22 = d22 + chi2(map_wx(df))*diff_basis_wx(2,df,1)
  end do

  if ( cubed_sphere ) then

    ! Find the coordinates varying along xi1 (ca) and xi2 (cb), allowing for
    ! cells whose axes are rotated or reversed relative to (alpha, beta).
    ! The formulae are unchanged by these symmetries of a cubed-sphere panel
    if ( abs(d11) >= abs(d12) ) then
      ca = sign(1.0_r_def, d11)*c1
      cb = sign(1.0_r_def, d22)*c2
    else
      ca = sign(1.0_r_def, d21)*c2
      cb = sign(1.0_r_def, d12)*c1
    end if

    x = tan(ca)
    y = tan(cb)
    p = sqrt(1.0_r_def + x**2)
    q = sqrt(1.0_r_def + y**2)
    delta = sqrt(1.0_r_def + x**2 + y**2)

    c11 = - x*y**2/(p*delta)
    c12 = y*q/delta
    c21 = x*p/delta
    c22 = - x**2*y/(q*delta)
    c_diag = (p*q/delta)**2
    c_off  = - x*y*p*q/delta**2

    do k = 0, nlayers - 1
      scale = dt_metric/(height_w3(map_w3(1)+k) + scaled_radius)

      ! Transporting wind in the cell frame
      u1 = 0.5_r_def*( real(transporting_wind(map_w2(W)+k), r_def)           &
                       /da(map_w2(W)+k)                                      &
                     + real(transporting_wind(map_w2(E)+k), r_def)           &
                       /da(map_w2(E)+k) )
      u2 = - 0.5_r_def*( real(transporting_wind(map_w2(S)+k), r_def)         &
                         /da(map_w2(S)+k)                                    &
                       + real(transporting_wind(map_w2(N)+k), r_def)         &
                         /da(map_w2(N)+k) )

      ! Transported wind in the cell frame
      a1 = 0.5_r_def*(uv_w2h(map_w2h(W)+k) + uv_w2h(map_w2h(E)+k))
      a2 = - 0.5_r_def*(uv_w2h(map_w2h(S)+k) + uv_w2h(map_w2h(N)+k))
      a3 = 0.5_r_def*(w_wt(map_wt(1)+k) + w_wt(map_wt(1)+k+1))

      d_uv1_w3(map_w3(1)+k) = scale*u1*(c11*a1 + c12*a2 - a3)
      d_uv2_w3(map_w3(1)+k) = scale*u2*(c21*a1 + c22*a2 - a3)
      d_w_w3(map_w3(1)+k) = scale*( c_diag*(u1*a1 + u2*a2)                   &
                                  + c_off*(u1*a2 + u2*a1) )
    end do

  else

    ! The reference coordinate xi2 is aligned with latitude
    tan_lat = sign(1.0_r_def, d22)*tan(c2)

    do k = 0, nlayers - 1
      scale = dt_metric/(height_w3(map_w3(1)+k) + scaled_radius)

      ! Transporting wind in the cell frame
      u1 = 0.5_r_def*( real(transporting_wind(map_w2(W)+k), r_def)           &
                       /da(map_w2(W)+k)                                      &
                     + real(transporting_wind(map_w2(E)+k), r_def)           &
                       /da(map_w2(E)+k) )
      u2 = - 0.5_r_def*( real(transporting_wind(map_w2(S)+k), r_def)         &
                         /da(map_w2(S)+k)                                    &
                       + real(transporting_wind(map_w2(N)+k), r_def)         &
                         /da(map_w2(N)+k) )

      ! Transported wind in the cell frame
      a1 = 0.5_r_def*(uv_w2h(map_w2h(W)+k) + uv_w2h(map_w2h(E)+k))
      a2 = - 0.5_r_def*(uv_w2h(map_w2h(S)+k) + uv_w2h(map_w2h(N)+k))
      a3 = 0.5_r_def*(w_wt(map_wt(1)+k) + w_wt(map_wt(1)+k+1))

      d_uv1_w3(map_w3(1)+k) = scale*u1*(tan_lat*a2 - a3)
      d_uv2_w3(map_w3(1)+k) = - scale*(tan_lat*u1*a1 + u2*a3)
      d_w_w3(map_w3(1)+k) = scale*(u1*a1 + u2*a2)
    end do

  end if

end subroutine wind_metric_increment_w3_code

end module wind_metric_increment_w3_kernel_mod
