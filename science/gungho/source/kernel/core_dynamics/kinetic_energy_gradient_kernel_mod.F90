!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Computes the kinetic gradient component of the rhs of the momentum
!>        equation for the nonlinear equations.
!>
!> The kernel computes the kinetic gradient component of the rhs of the
!> momentum equation for the nonlinear equations, written in the vector
!> invariant form.
!>
!> This consists of four terms:
!> Pressure gradient: \f[ cp*\theta*\nabla(\Pi)\f]
!> geopotential gradient: \f[ \nabla(\Phi) ( \equiv g for some domains)\f]
!> gradient of kinetic energy: \f[ \nabla(1/2*u.u) \f]
!> vorticity advection: \f[ \xi/\rho \times F (with vorticity \xi and mass flux F) \f]
!>
!> This results in:
!> \f[ r_u = -\xi/\rho \times F - \nabla(\Phi + 1/2*u.u) - cp*\theta*\nabla(\Pi) \f]
!>
module kinetic_energy_gradient_kernel_mod

  use argument_mod,      only : arg_type, func_type,       &
                                GH_FIELD, GH_REAL,         &
                                GH_READ, GH_INC,           &
                                ANY_SPACE_9,               &
                                ANY_DISCONTINUOUS_SPACE_3, &
                                GH_BASIS, GH_DIFF_BASIS,   &
                                CELL_COLUMN, GH_QUADRATURE_XYoZ
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W2
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  type, public, extends(kernel_type) :: kinetic_energy_gradient_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                    &
         arg_type(GH_FIELD,   GH_REAL, GH_INC,  W2),                       &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, W2),                       &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ, ANY_SPACE_9),              &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_3) &
         /)
    type(func_type) :: meta_funcs(2) = (/                                  &
         func_type(W2,          GH_BASIS, GH_DIFF_BASIS),                  &
         func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                   &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: kinetic_energy_gradient_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: kinetic_energy_gradient_code

contains

!> @brief Computes the kinetic gradient component of the rhs of the momentum equation
!! @param[in] nlayers Number of layers
!! @param[in,out] r_u Right hand side of momentum equation
!! @param[in] u Velocity
!! @param[in] chi_1 1st coordinate field in Wchi
!! @param[in] chi_2 2nd coordinate field in Wchi
!! @param[in] chi_3 3rd coordinate field in Wchi
!! @param[in] panel_id Field giving the ID for mesh panels
!! @param[in] ndf_w2 Number of degrees of freedom per cell for w2
!! @param[in] undf_w2 Number of unique degrees of freedom  for w2
!! @param[in] map_w2 Dofmap for the cell at the base of the column for w2
!! @param[in] w2_basis Basis functions evaluated at quadrature points
!! @param[in] w2_diff_basis Differential of the W2 basis functions evaluated
!!                          at gaussian quadrature points
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi
!! @param[in] undf_chi Number of unique degrees of freedom  for chi
!! @param[in] map_chi Dofmap for the cell at the base of the column for chi
!! @param[in] chi_basis Wchi basis functions evaluated
!!                      at gaussian quadrature points
!! @param[in] chi_diff_basis Differential of the Wchi basis functions evaluated
!!                           at gaussian quadrature points
!! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
!! @param[in] undf_pid Number of unique degrees of freedom for panel_id
!! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
!! @param[in] nqp_h Number of quadrature points in the horizontal
!! @param[in] nqp_v Number of quadrature points in the vertical
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine kinetic_energy_gradient_code(nlayers,                       &
                                        r_u, u,                        &
                                        chi_1, chi_2, chi_3, panel_id, &
                                        ndf_w2, undf_w2, map_w2,       &
                                        w2_basis, w2_diff_basis,       &
                                        ndf_chi, undf_chi, map_chi,    &
                                        chi_basis, chi_diff_basis,     &
                                        ndf_pid, undf_pid, map_pid,    &
                                        nqp_h, nqp_v, wqp_h, wqp_v     &
                                        )

  use sci_coordinate_jacobian_mod, only: coordinate_jacobian

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers,nqp_h, nqp_v
  integer(kind=i_def), intent(in) :: ndf_chi, ndf_w2, ndf_pid
  integer(kind=i_def), intent(in) :: undf_chi, undf_w2, undf_pid
  integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_w2),  intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_pid), intent(in) :: map_pid

  real(kind=r_def), dimension(3,ndf_w2,nqp_h,nqp_v),  intent(in) :: w2_basis
  real(kind=r_def), dimension(1,ndf_w2,nqp_h,nqp_v),  intent(in) :: w2_diff_basis
  real(kind=r_def), dimension(1,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_basis
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_diff_basis

  real(kind=r_def), dimension(undf_w2),  intent(inout) :: r_u
  real(kind=r_def), dimension(undf_w2),  intent(in)    :: u
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(undf_pid), intent(in)    :: panel_id

  real(kind=r_def), dimension(nqp_h), intent(in)      ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)      ::  wqp_v

  ! Internal variables
  integer(kind=i_def) :: df, k, loc, ipanel
  integer(kind=i_def) :: qp1, qp2

  real(kind=r_def), dimension(ndf_chi)         :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac
  real(kind=r_def), dimension(ndf_w2)          :: ru_e, u_e

  real(kind=r_def) :: u_at_quad(3)
  real(kind=r_def) :: ke_at_quad, dv

  ipanel = int(panel_id(map_pid(1)), i_def)

  do k = 0, nlayers-1
  ! Extract element arrays of chi
    do df = 1, ndf_chi
      loc = map_chi(df) + k
      chi_1_e(df) = chi_1( loc )
      chi_2_e(df) = chi_2( loc )
      chi_3_e(df) = chi_3( loc )
    end do
    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi_1_e, chi_2_e, chi_3_e,  &
                             ipanel, chi_basis, chi_diff_basis, jac, dj)

    do df = 1, ndf_w2
      u_e(df) = u( map_w2(df) + k )
      ru_e(df) = 0.0_r_def
    end do
  ! compute the RHS integrated over one cell
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
  ! k.e term
        u_at_quad(:) = 0.0_r_def
        do df = 1, ndf_w2
          u_at_quad(:) = u_at_quad(:) &
                       + u_e(df)*w2_basis(:,df,qp1,qp2)
        end do
        ke_at_quad = 0.5_r_def*dot_product(matmul(jac(:,:,qp1,qp2),u_at_quad), &
                                           matmul(jac(:,:,qp1,qp2),u_at_quad))/(dj(qp1,qp2)**2)

        do df = 1, ndf_w2
          dv = w2_diff_basis(1,df,qp1,qp2)
          ru_e(df) = ru_e(df) +  wqp_h(qp1)*wqp_v(qp2)*dv*ke_at_quad

        end do
      end do
    end do
    do df = 1, ndf_w2
      r_u( map_w2(df) + k ) =  r_u( map_w2(df) + k ) + ru_e(df)
    end do
  end do

end subroutine kinetic_energy_gradient_code

end module kinetic_energy_gradient_kernel_mod
