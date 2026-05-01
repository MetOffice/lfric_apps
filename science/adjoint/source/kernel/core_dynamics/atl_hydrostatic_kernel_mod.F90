!-----------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be brief.
!-----------------------------------------------------------------------------

!> @brief Computes the adjoint for hydrostatic balance.
module atl_hydrostatic_kernel_mod

  use argument_mod,      only : arg_type, func_type,     &
                                GH_FIELD, GH_REAL,       &
                                GH_SCALAR,               &
                                GH_READ, GH_READWRITE,   &
                                GH_BASIS, GH_DIFF_BASIS, &
                                CELL_COLUMN,             &
                                GH_QUADRATURE_XYoZ,      &
                                ANY_W2
  use constants_mod,     only : r_def
  use fs_continuity_mod, only : W3, Wtheta
  use kernel_mod,        only : kernel_type

  implicit none
  private
  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: atl_hydrostatic_kernel_type
  private
    type(arg_type) :: meta_args(8) = (/                           &
        arg_type(GH_FIELD,   GH_REAL, GH_READ,      ANY_W2),      &
        arg_type(GH_FIELD,   GH_REAL, GH_READWRITE, W3),          &
        arg_type(GH_FIELD,   GH_REAL, GH_READWRITE, Wtheta),      &
        arg_type(GH_FIELD*3, GH_REAL, GH_READWRITE, Wtheta),      &
        arg_type(GH_FIELD,   GH_REAL, GH_READ,      W3),          &
        arg_type(GH_FIELD,   GH_REAL, GH_READ,      Wtheta),      &
        arg_type(GH_FIELD*3, GH_REAL, GH_READ,      Wtheta),      &
        arg_type(GH_SCALAR,  GH_REAL, GH_READ)                    &
        /)
    type(func_type) :: meta_funcs(3) = (/                &
        func_type(ANY_W2,      GH_BASIS, GH_DIFF_BASIS), &
        func_type(W3,          GH_BASIS),                &
        func_type(Wtheta,      GH_BASIS, GH_DIFF_BASIS)  &
        /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: atl_hydrostatic_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: atl_hydrostatic_code

  contains

!> @brief Compute the tangent linear of the hydrostatic term.
!! @param[in]      nlayers       Number of layers
!! @param[in]      r_u               ACTIVE Change in Momentum equation right hand side
!! @param[in,out]  exner             ACTIVE Change in Exner pressure field
!! @param[in,out]  theta             ACTIVE Change in Potential temperature field
!! @param[in,out]  moist_dyn_gas     ACTIVE Change in Moist dynamics factor in gas law
!! @param[in,out]  moist_dyn_tot     ACTIVE Change in Moist dynamics total mass factor
!! @param[in]      moist_dyn_fac     ACTIVE Change in Moist dynamics water factor
!! @param[in]      ls_exner          Lin state Exner pressure field
!! @param[in]      ls_theta          Lin state Potential temperature field
!! @param[in]      ls_moist_dyn_gas  Lin state Moist dynamics factor in gas law
!! @param[in]      ls_moist_dyn_tot  Lin state Moist dynamics total mass factor
!! @param[in]      ls_moist_dyn_fac  Lin state Moist dynamics water factor
!! @param[in]      cp                Specific heat of dry air at constant pressure
!! @param[in]      ndf_w2            Number of degrees of freedom per cell for w2
!! @param[in]      undf_w2           Number of unique degrees of freedom  for w2
!! @param[in]      map_w2            Dofmap for the cell at the base of the column for w2
!! @param[in]      w2_basis          Basis functions evaluated at quadrature points
!! @param[in]      w2_diff_basis     Differential of the basis functions evaluated
!!                                   at quadrature points
!! @param[in]      ndf_w3            Number of degrees of freedom per cell for w3
!! @param[in]      undf_w3           Number of unique degrees of freedom  for w3
!! @param[in]      map_w3            Dofmap for the cell at the base of the column for w3
!! @param[in]      w3_basis          Basis functions evaluated at gaussian quadrature
!!                                   points
!! @param[in]      ndf_wt            Number of degrees of freedom per cell for wt
!! @param[in]      undf_wt           Number of unique degrees of freedom  for wt
!! @param[in]      map_wt            Dofmap for the cell at the base of the column for wt
!! @param[in]      wt_basis          Basis functions evaluated at gaussian quadrature
!!                                   points
!! @param[in]      wt_diff_basis     Differential of the basis functions evaluated at
!!                                   quadrature points
!! @param[in]      nqp_h             Number of quadrature points in the horizontal
!! @param[in]      nqp_v             Number of quadrature points in the vertical
!! @param[in]      wqp_h             Horizontal quadrature weights
!! @param[in]      wqp_v             Vertical quadrature weights
  subroutine atl_hydrostatic_code(nlayers,          &
                                  r_u,              &
                                  exner,            &
                                  theta,            &
                                  moist_dyn_gas,    &
                                  moist_dyn_tot,    &
                                  moist_dyn_fac,    &
                                  ls_exner,         &
                                  ls_theta,         &
                                  ls_moist_dyn_gas, &
                                  ls_moist_dyn_tot, &
                                  ls_moist_dyn_fac, &
                                  cp,               &
                                  ndf_w2,           &
                                  undf_w2,          &
                                  map_w2,           &
                                  w2_basis,         &
                                  w2_diff_basis,    &
                                  ndf_w3,           &
                                  undf_w3,          &
                                  map_w3,           &
                                  w3_basis,         &
                                  ndf_wt,           &
                                  undf_wt,          &
                                  map_wt,           &
                                  wt_basis,         &
                                  wt_diff_basis,    &
                                  nqp_h,            &
                                  nqp_v,            &
                                  wqp_h,            &
                                  wqp_v             &
                                  )

    implicit none

  ! Arguments
  integer, intent(in) :: nlayers,nqp_h, nqp_v
  integer, intent(in) :: ndf_wt, ndf_w2, ndf_w3
  integer, intent(in) :: undf_wt, undf_w2, undf_w3
  integer, dimension(ndf_wt), intent(in) :: map_wt
  integer, dimension(ndf_w2), intent(in) :: map_w2
  integer, dimension(ndf_w3), intent(in) :: map_w3

  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v), intent(in) :: w3_basis
  real(kind=r_def), dimension(3,ndf_w2,nqp_h,nqp_v), intent(in) :: w2_basis
  real(kind=r_def), dimension(1,ndf_wt,nqp_h,nqp_v), intent(in) :: wt_basis
  real(kind=r_def), dimension(1,ndf_w2,nqp_h,nqp_v), intent(in) :: w2_diff_basis
  real(kind=r_def), dimension(3,ndf_wt,nqp_h,nqp_v), intent(in) :: wt_diff_basis

  real(kind=r_def), dimension(undf_w2), intent(in)    :: r_u
  real(kind=r_def), dimension(undf_w3), intent(inout) :: exner
  real(kind=r_def), dimension(undf_wt), intent(inout) :: theta
  real(kind=r_def), dimension(undf_wt), intent(inout) :: moist_dyn_gas, &
                                                         moist_dyn_tot
  real(kind=r_def), dimension(undf_wt), intent(in)    :: moist_dyn_fac
  real(kind=r_def), dimension(undf_w3), intent(in)    :: ls_exner
  real(kind=r_def), dimension(undf_wt), intent(in)    :: ls_theta
  real(kind=r_def), dimension(undf_wt), intent(in)    :: ls_moist_dyn_gas, &
                                                         ls_moist_dyn_tot, &
                                                         ls_moist_dyn_fac
  real(kind=r_def),                     intent(in)    :: cp

  real(kind=r_def), dimension(nqp_h), intent(in)      ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)      ::  wqp_v

  ! Internal variables
  integer               :: df, k
  integer               :: qp1, qp2

  real(kind=r_def), dimension(ndf_w3)          :: exner_e
  real(kind=r_def), dimension(ndf_wt)          :: theta_v_e
  real(kind=r_def), dimension(ndf_w3)          :: ls_exner_e
  real(kind=r_def), dimension(ndf_wt)          :: ls_theta_v_e

  real(kind=r_def) :: grad_theta_v_at_quad(3), ls_grad_theta_v_at_quad(3), v(3)
  real(kind=r_def) :: exner_at_quad, theta_v_at_quad,       &
                      ls_exner_at_quad, ls_theta_v_at_quad, &
                      grad_term, dv

  theta_v_e = 0.0_r_def
  grad_theta_v_at_quad = 0.0_r_def
  exner_at_quad = 0.0_r_def
  theta_v_at_quad = 0.0_r_def
  exner_e = 0.0_r_def
  do k = nlayers - 1, 0, -1

    ! Linearisation state
    do df = 1, ndf_w3, 1
      ls_exner_e(df) = ls_exner( map_w3(df) + k )
    end do
    do df = 1, ndf_wt
      ls_theta_v_e(df) = ls_theta( map_wt(df) + k ) * &
                         ls_moist_dyn_gas( map_wt(df) + k ) / &
                         ls_moist_dyn_tot( map_wt(df) + k )
    end do

    do qp2 = nqp_v, 1, -1
      do qp1 = nqp_h, 1, -1

        ! Linearisation state
        ls_exner_at_quad = 0.0_r_def
        do df = 1, ndf_w3, 1
          ls_exner_at_quad = ls_exner_at_quad + ls_exner_e(df) * w3_basis(1,df,qp1,qp2)
        enddo
        ls_theta_v_at_quad = 0.0_r_def
        ls_grad_theta_v_at_quad(:) = 0.0_r_def
        do df = 1, ndf_wt
          ls_theta_v_at_quad   = ls_theta_v_at_quad                                 &
                            + ls_theta_v_e(df)*wt_basis(1,df,qp1,qp2)
          ls_grad_theta_v_at_quad(:) = ls_grad_theta_v_at_quad(:)                   &
                                  + ls_theta_v_e(df)*wt_diff_basis(:,df,qp1,qp2)
        end do

        ! Perturbation
        do df = ndf_w2, 1, -1
          v = w2_basis(:,df,qp1,qp2)
          dv = w2_diff_basis(1,df,qp1,qp2)
          grad_term = r_u(map_w2(df) + k) * wqp_h(qp1) * wqp_v(qp2)
          exner_at_quad = exner_at_quad + grad_term * cp * &
                          (ls_theta_v_at_quad * dv + dot_product(ls_grad_theta_v_at_quad, v))
          theta_v_at_quad = theta_v_at_quad + cp * dv * ls_exner_at_quad * grad_term
          grad_theta_v_at_quad = grad_theta_v_at_quad + cp * ls_exner_at_quad * grad_term * v
        end do

        do df = ndf_wt, 1, -1
          theta_v_e(df) = theta_v_e(df) + dot_product(grad_theta_v_at_quad(:),     &
                                                      wt_diff_basis(:,df,qp1,qp2)) &
                                        + theta_v_at_quad * wt_basis(1,df,qp1,qp2)
        end do
        grad_theta_v_at_quad = 0.0_r_def
        theta_v_at_quad = 0.0_r_def
        exner_e = exner_e + exner_at_quad * w3_basis(1,:,qp1,qp2)
        exner_at_quad = 0.0_r_def
      end do ! qp1
    end do   ! qp2

    do df = ndf_wt, 1, -1
      theta(k + map_wt(df)) = theta(k + map_wt(df)) +                     &
                              ls_theta_v_e(df) * theta_v_e(df) /          &
                              ls_theta(k + map_wt(df))
      moist_dyn_tot(k + map_wt(df)) = moist_dyn_tot(k + map_wt(df)) -     &
                                      ls_theta_v_e(df) * theta_v_e(df) /  &
                                      ls_moist_dyn_tot(k + map_wt(df))
      moist_dyn_gas(k + map_wt(df)) = moist_dyn_gas(k + map_wt(df)) +     &
                                      ls_theta_v_e(df) * theta_v_e(df) /  &
                                      ls_moist_dyn_gas(k + map_wt(df))
    end do
    theta_v_e = 0.0_r_def

    do df = ndf_w3, 1, -1
      exner(map_w3(df) + k) = exner(map_w3(df) + k) + exner_e(df)
      exner_e(df) = 0.0_r_def
    end do
  end do ! k

  end subroutine atl_hydrostatic_code

end module atl_hydrostatic_kernel_mod
