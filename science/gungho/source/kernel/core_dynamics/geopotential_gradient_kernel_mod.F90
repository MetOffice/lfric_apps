!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Computes rhs of the momentum equation for the nonlinear equations.
!>
!> The kernel computes the rhs of the momentum equation for the nonlinear
!> equations, written in the vector invariant form.
!>
!> This consists of four terms:
!>
!> Pressure gradient: \f$ cp\theta\nabla(\Pi) \f$
!>
!> geopotential gradient: \f$ \nabla(\Phi) \f$ ( \f$\equiv g\f$ for some domains)
!>
!> gradient of kinetic energy: \f$ \nabla(1/2u.u) \f$
!>
!> vorticity advection: \f$ \xi/\rho \times F\f$ (with vorticity \f$\xi\f$ and mass flux F)
!>
!> This results in:
!> \f[ r_u = -\xi/\rho \times F - \nabla(\Phi + 1/2u.u) - cp*\theta*\nabla(\Pi) \f]
!>
module geopotential_gradient_kernel_mod

  use argument_mod,      only : arg_type, func_type,     &
                                GH_FIELD, GH_REAL,       &
                                GH_READ, GH_INC,         &
                                GH_BASIS, GH_DIFF_BASIS, &
                                CELL_COLUMN, GH_QUADRATURE_XYoZ
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W0, W2
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: geopotential_gradient_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/            &
         arg_type(GH_FIELD, GH_REAL, GH_INC,  W2), &
         arg_type(GH_FIELD, GH_REAL, GH_READ, W0)  &
         /)
    type(func_type) :: meta_funcs(2) = (/          &
         func_type(W2, GH_BASIS),                  &
         func_type(W0, GH_DIFF_BASIS)              &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: geopotential_gradient_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: geopotential_gradient_code

contains

!> @brief Kernel which computes rhs of the momentum equation for the nonlinear
!>        equations, written in the vector invariant form
!! @param[in] nlayers Number of layers
!! @param[in] ndf_w2 Number of degrees of freedom per cell for w2
!! @param[in] undf_w2 Number of unique degrees of freedom  for w2
!! @param[in] map_w2 Dofmap for the cell at the base of the column for w2
!! @param[in] w2_basis Basis functions evaluated at quadrature points
!! @param[in,out] r_u Right hand side of the momentum equation
!! @param[in] ndf_w0 Number of degrees of freedom per cell for w0
!! @param[in] undf_w0 Number of unique degrees of freedom  for w0
!! @param[in] map_w0 Dofmap for the cell at the base of the column for w0
!! @param[in] w0_diff_basis Differntial of the basis functions evaluated at
!!                          gaussian quadrature point
!! @param[in] phi Geopotential
!! @param[in] nqp_h Number of quadrature points in the horizontal
!! @param[in] nqp_v Number of quadrature points in the vertical
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine geopotential_gradient_code(nlayers,                                 &
                                      r_u, phi,                                &
                                      ndf_w2, undf_w2, map_w2, w2_basis,       &
                                      ndf_w0, undf_w0, map_w0, w0_diff_basis,  &
                                      nqp_h, nqp_v, wqp_h, wqp_v               &
                                      )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, nqp_h, nqp_v
  integer(kind=i_def), intent(in) :: ndf_w0, ndf_w2
  integer(kind=i_def), intent(in) :: undf_w0, undf_w2
  integer(kind=i_def), dimension(ndf_w0), intent(in) :: map_w0
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2

  real(kind=r_def), dimension(3,ndf_w2,nqp_h,nqp_v), intent(in) :: w2_basis
  real(kind=r_def), dimension(3,ndf_w0,nqp_h,nqp_v), intent(in) :: w0_diff_basis

  real(kind=r_def), dimension(undf_w2), intent(inout) :: r_u
  real(kind=r_def), dimension(undf_w0), intent(in)    :: phi

  real(kind=r_def), dimension(nqp_h), intent(in)      ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)      ::  wqp_v

  ! Internal variables
  integer(kind=i_def) :: df, k
  integer(kind=i_def) :: qp1, qp2

  real(kind=r_def), dimension(ndf_w2)          :: ru_e
  real(kind=r_def), dimension(ndf_w0)          :: phi_e

  real(kind=r_def) :: grad_phi_at_quad(3)
  real(kind=r_def) :: geo_term

  do k = 0, nlayers-1
    do df = 1, ndf_w0
      phi_e(df)   = phi(map_w0(df) + k)
    end do
    do df = 1, ndf_w2
      ru_e(df) = 0.0_r_def
    end do
  ! compute the RHS integrated over one cell
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        grad_phi_at_quad(:) = 0.0_r_def
        do df = 1, ndf_w0
          grad_phi_at_quad(:)   = grad_phi_at_quad(:) &
                                + phi_e(df)*w0_diff_basis(:,df,qp1,qp2)
        end do

        do df = 1, ndf_w2
! geopotential term
          geo_term = dot_product( grad_phi_at_quad(:), w2_basis(:,df,qp1,qp2))

          ru_e(df) = ru_e(df) -  wqp_h(qp1)*wqp_v(qp2)*geo_term

        end do
      end do
    end do
    do df = 1, ndf_w2
      r_u( map_w2(df) + k ) =  r_u( map_w2(df) + k ) + ru_e(df)
    end do
  end do

end subroutine geopotential_gradient_code

end module geopotential_gradient_kernel_mod
