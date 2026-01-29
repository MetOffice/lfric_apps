!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes u_inc, the change in TLM velocity due to TLM boundary layer processes.
module tl_bl_inc_kernel_mod

  use argument_mod,           only : arg_type,              &
                                     GH_FIELD, GH_OPERATOR, &
                                     GH_SCALAR, GH_INTEGER, &
                                     GH_READ, GH_INC,       &
                                     GH_REAL, CELL_COLUMN,  &
                                     ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,          only : r_def, i_def, r_um
  use fs_continuity_mod,      only : W1, W2, W3
  use kernel_mod,             only : kernel_type
  use reference_element_mod,  only : N

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  type, public, extends(kernel_type) :: tl_bl_inc_kernel_type
    private
    type(arg_type) :: meta_args(7) = (/                                       &
         arg_type(GH_FIELD,  GH_REAL,    GH_INC,  W2),                        & ! u_inc
         arg_type(GH_FIELD,  GH_REAL,    GH_READ, W2),                        & ! u
         arg_type(GH_FIELD,  GH_REAL,    GH_READ, W2),                        & ! Auv
         arg_type(GH_FIELD,  GH_REAL,    GH_READ, W2),                        & ! Buv_inv
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ, ANY_DISCONTINUOUS_SPACE_1), & ! face_selector_ew
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ, ANY_DISCONTINUOUS_SPACE_1), & ! face_selector_ns
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ )                            & ! Blevs_m
    /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: tl_bl_inc_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: tl_bl_inc_code

contains

!> @brief Computes u_inc, the change in TLM velocity due to TLM boundary layer processes.
!> @details The algorithm uses coefficients Auv and Buv_inv computed in tl_compute_aubu_kernel_mod.
!>          The TLM BL scheme is described in Var Scientific Documentation Paper 55,
!>          https://wwwspice/~frva/VAR/view/var-2025.03.0/doc/VSDP55_1A.pdf
!! @param[in]     nlayers           Number of layers
!! @param[in,out] u_inc             Change in TLM velocity due to TLM boundary layer processes
!! @param[in]     u                 TLM velocity
!! @param[in]     nlayers           Number of layers
!! @param[in]     Auv               Coefficient for TLM boundary layer
!! @param[in]     Buv_inv           Inverse of coefficient for TLM boundary layer
!! @param[in]     face_selector_ew  2D field indicating which W/E faces to loop over in this column
!! @param[in]     face_selector_ns  2D field indicating which N/S faces to loop over in this column
!! @param[in]     Blevs_m           Number of levels in momentum boundary layer
!! @param[in]     ndf_w2            Number of degrees of freedom per cell for w2 space
!! @param[in]     undf_w2           Number of unique degrees of freedom for w2 space
!! @param[in]     map_w2            Dofmap for the cell at the base of the column for w2
!! @param[in]     ndf_w3_2d         Number of DoFs for 2D W3 per cell
!! @param[in]     undf_w3_2d        Number of DoFs for this partition for 2D W3
!! @param[in]     map_w3_2d         Map for 2D W3
subroutine tl_bl_inc_code( nlayers,                 &
                           u_inc,                   &
                           u,                       &
                           Auv,                     &
                           Buv_inv,                 &
                           face_selector_ew,        &
                           face_selector_ns,        &
                           Blevs_m,                 &
                           ndf_w2, undf_w2, map_w2, &
                           ndf_w3_2d, undf_w3_2d, map_w3_2d )

  implicit none

  ! Arguments
  integer(kind=i_def),                        intent(in)    :: nlayers
  integer(kind=i_def),                        intent(in)    :: undf_w2, ndf_w2
  real(kind=r_def),    dimension(undf_w2),    intent(inout) :: u_inc
  integer(kind=i_def), dimension(ndf_w2),     intent(in)    :: map_w2
  real(kind=r_def),    dimension(undf_w2),    intent(in)    :: u
  real(kind=r_def),    dimension(undf_w2),    intent(in)    :: Auv
  real(kind=r_def),    dimension(undf_w2),    intent(in)    :: Buv_inv
  integer(kind=i_def),                        intent(in)    :: ndf_w3_2d, undf_w3_2d
  integer(kind=i_def), dimension(ndf_w3_2d),  intent(in)    :: map_w3_2d
  integer(kind=i_def), dimension(undf_w3_2d), intent(in)    :: face_selector_ew
  integer(kind=i_def), dimension(undf_w3_2d), intent(in)    :: face_selector_ns
  integer(kind=i_def),                        intent(in)    :: Blevs_m

  ! Internal variables
  integer(kind=i_def) :: df, k, j
  real(kind=r_def)    :: a0(1:BLevs_m) ! Coefficient
  real(kind=r_def)    :: a1(1:BLevs_m) ! Coefficient
  real(kind=r_def)    :: a2(1:BLevs_m) ! Coefficient
  real(kind=r_def)    :: u_rhs(1:BLevs_m) ! Local perturbation velocity variable
  real(kind=r_def)    :: u_out(1:BLevs_m) ! Local perturbation velocity variable
  real(kind=r_def)    :: factor_u(1:BLevs_m)

  ! Loop over horizontal W2 DoFs whilst minimising double counting.
  ! (Looping over all dofs would mean that faces are visited twice – for the cells on both sides.
  ! So here the loop is only for a specific selection of dofs.)
  do j = 1, face_selector_ew(map_w3_2d(1)) + face_selector_ns(map_w3_2d(1))

    df = j
    if (j == 3 .and. face_selector_ns(map_w3_2d(1)) == 2 .and. face_selector_ew(map_w3_2d(1)) == 1) df = N

    a0 = 0.0_r_def
    a1 = 0.0_r_def
    a2 = 0.0_r_def
    u_rhs = 0.0_r_def
    u_out = 0.0_r_def
    factor_u = 0.0_r_def

    ! Set up coeffs a0, a1, a2, u_rhs
    do k = 1, BLevs_m
      if (k == 1) then
        a0(1)    = 1.0_r_def + (Auv(map_w2(df) + 1) + Auv(map_w2(df) + 0)) / Buv_inv(map_w2(df) + 1)
        a1(k)    = -Auv(map_w2(df) + 1) / Buv_inv(map_w2(df) + 1)
        u_rhs(1) = (Auv(map_w2(df) + 1) &
          * (u(map_w2(df) + 1) - u(map_w2(df) + 0)) - Auv(map_w2(df) + 0) * u(map_w2(df) + 0)) / Buv_inv(map_w2(df) + 1)
      else if (k >  1 .and. k <  BLevs_m) then
        a0(k)    = 1.0_r_def + (Auv(map_w2(df) + k) + Auv(map_w2(df) + k - 1)) / Buv_inv(map_w2(df) + k)
        a2(k)    = -Auv(map_w2(df) + k - 1) / Buv_inv(map_w2(df) + k)
        a1(k)    = -Auv(map_w2(df) + k) / Buv_inv(map_w2(df) + k)
        u_rhs(k) = (Auv(map_w2(df) + k)                 &
          * (u(map_w2(df) + k) - u(map_w2(df) + k - 1)) &
          - Auv(map_w2(df) + k - 1) * (u(map_w2(df) + k - 1) - u(map_w2(df) + k - 2))) / Buv_inv(map_w2(df) + k)
      else
        a0(k)    = 1.0_r_def + Auv(map_w2(df) + k - 1) / Buv_inv(map_w2(df) + k)
        a2(k)    = -Auv(map_w2(df) + k - 1) / Buv_inv(map_w2(df) + k)
        u_rhs(k) = -(Auv(map_w2(df) + k - 1) * (u(map_w2(df) + k - 1) - u(map_w2(df) + k - 2))) / Buv_inv(map_w2(df) + k)
      end if
    end do

    ! Transform to upper triangular form
    do k = 1, BLevs_m
      if (k == 1) then
        a0(1) = 1.0_r_def / a0(1)
      else
        factor_u(k) = a2(k) * a0(k - 1)
        a0(k)       = 1.0_r_def / (a0(k) - factor_u(k) * a1(k - 1))
        u_rhs(k)    = u_rhs(k) - factor_u(k) * u_rhs(k - 1)
      end if
    end do

    ! Solve for u_inc
    u_out(BLevs_m) = a0(BLevs_m) * u_rhs(BLevs_m)
    u_inc(map_w2(df) + BLevs_m - 1) = u_out(BLevs_m)
    do k = BLevs_m - 1, 1, -1
      u_out(k) = a0(k) * (u_rhs(k) - a1(k) * u_out(k + 1))
      u_inc(map_w2(df) + k - 1) = u_out(k)
    end do

  end do ! Loop over horizontal W2 DoFs

end subroutine tl_bl_inc_code

end module tl_bl_inc_kernel_mod
