!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes the vorticity component of the rhs of the momentum
!>        equation for the nonlinear equations, written in the vector invariant form


!> @details The kernel computes the  vorticity component of the rhs of the momentum equation
!>         for the nonlinear equations, written in the vector invariant form
!>         This consists of four terms:
!>         Pressure gradient: \f[ cp*\theta*\nabla(\Pi) \f]
!>         geopotential gradient: \f[ \nabla(\Phi) ( \equiv g for some domains) \f]
!>         gradient of kinetic energy: \f[ \nabla(1/2*u.u) \f]
!>         vorticity advection: \f[ \xi/\rho \times F (with vorticity \xi and mass flux F) \f]
!>         This results in:
!>         \f[ r_u = -\xi/\rho \times F - \nabla(\Phi + 1/2*u.u) - cp*\theta*\nabla(\Pi) \f]
!>         Where \f[ \xi \f] is the Vorticity
module vorticity_advection_kernel_mod

use kernel_mod,              only: kernel_type
use argument_mod,            only: arg_type, func_type,       &
                                   GH_FIELD, GH_REAL,         &
                                   GH_READ, GH_INC,           &
                                   ANY_SPACE_9,               &
                                   ANY_DISCONTINUOUS_SPACE_3, &
                                   GH_BASIS, GH_DIFF_BASIS,   &
                                   CELL_COLUMN, GH_QUADRATURE_XYoZ
use constants_mod,           only: r_def, i_def
use fs_continuity_mod,       only: W1, W2
use cross_product_mod,       only: cross_product

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: vorticity_advection_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                                    &
       arg_type(GH_FIELD,   GH_REAL, GH_INC,  W2),                       &
       arg_type(GH_FIELD,   GH_REAL, GH_READ, W2),                       &
       arg_type(GH_FIELD,   GH_REAL, GH_READ, W1),                       &
       arg_type(GH_FIELD*3, GH_REAL, GH_READ, ANY_SPACE_9),              &
       arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_3) &
       /)
  type(func_type) :: meta_funcs(3) = (/                                  &
       func_type(W2,          GH_BASIS),                                 &
       func_type(W1,          GH_BASIS),                                 &
       func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                   &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_QUADRATURE_XYoZ
contains
  procedure, nopass :: vorticity_advection_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: vorticity_advection_code

contains

!> @brief Compute the advection of the wind field by the vorticity
!! @param[in] nlayers Number of layers
!! @param[in,out] r_u Right hand side of the momentum equation
!! @param[in] wind Advecting wind field
!! @param[in] vorticity Vorticity field = curl(u)
!! @param[in] chi_1 1st coordinate field in Wchi
!! @param[in] chi_2 2nd coordinate field in Wchi
!! @param[in] chi_3 3rd coordinate field in Wchi
!! @param[in] panel_id Field giving the ID for mesh panels
!! @param[in] ndf_w2 Number of degrees of freedom per cell for W2
!! @param[in] undf_w2 Number of unique degrees of freedom  for W2
!! @param[in] map_w2 Dofmap for the cell at the base of the column for W2
!! @param[in] w2_basis Basis functions evaluated at quadrature points
!! @param[in] ndf_w1 Number of degrees of freedom per cell for W1
!! @param[in] undf_w1 Number of unique degrees of freedom  for W1
!! @param[in] map_w1 Dofmap for the cell at the base of the column for W1
!! @param[in] w1_basis Basis functions evaluated at Gaussian quadrature points
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi
!! @param[in] undf_chi Number of unique degrees of freedom  for chi
!! @param[in] map_chi Dofmap for the cell at the base of the column for chi
!! @param[in] chi_basis Wchi basis functions evaluated at Gaussian
!!                      quadrature points
!! @param[in] chi_diff_basis Derivatives of Wchi basis functions
!!                           evaluated at Gaussian quadrature points
!! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
!! @param[in] undf_pid Number of unique degrees of freedom for panel_id
!! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
!! @param[in] nqp_h Number of quadrature points in the horizontal
!! @param[in] nqp_v Number of quadrature points in the vertical
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine vorticity_advection_code(nlayers,                                   &
                                    r_u, wind, vorticity,                      &
                                    chi_1, chi_2, chi_3, panel_id,             &
                                    ndf_w2, undf_w2, map_w2, w2_basis,         &
                                    ndf_w1, undf_w1, map_w1, w1_basis,         &
                                    ndf_chi, undf_chi, map_chi,                &
                                    chi_basis, chi_diff_basis,                 &
                                    ndf_pid, undf_pid, map_pid,                &
                                    nqp_h, nqp_v, wqp_h, wqp_v                 &
                                    )

  use sci_coordinate_jacobian_mod, only: pointwise_coordinate_jacobian, &
                                         pointwise_coordinate_jacobian_inverse

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers,nqp_h, nqp_v
  integer(kind=i_def), intent(in) :: ndf_chi, ndf_w1, ndf_w2, ndf_pid
  integer(kind=i_def), intent(in) :: undf_chi, undf_w1, undf_w2, undf_pid
  integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_w1),  intent(in) :: map_w1
  integer(kind=i_def), dimension(ndf_w2),  intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_pid), intent(in) :: map_pid

  real(kind=r_def), dimension(3,ndf_w2,nqp_h,nqp_v),  intent(in) :: w2_basis
  real(kind=r_def), dimension(3,ndf_w1,nqp_h,nqp_v),  intent(in) :: w1_basis
  real(kind=r_def), dimension(1,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_basis
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_diff_basis

  real(kind=r_def), dimension(undf_w2),  intent(inout) :: r_u
  real(kind=r_def), dimension(undf_w2),  intent(in)    :: wind
  real(kind=r_def), dimension(undf_w1),  intent(in)    :: vorticity
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(undf_pid), intent(in)    :: panel_id

  real(kind=r_def), dimension(nqp_h),    intent(in)    ::  wqp_h
  real(kind=r_def), dimension(nqp_v),    intent(in)    ::  wqp_v

  ! Internal variables
  integer(kind=i_def) :: df, k, loc, ipanel
  integer(kind=r_def) :: qp1, qp2

  real(kind=r_def), dimension(ndf_chi)  :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def)                      :: dj
  real(kind=r_def), dimension(3,3)      :: jac, jac_inv
  real(kind=r_def), dimension(3)        :: vorticity_at_quad, u_at_quad, &
                                           vorticity_term, j_vorticity

  ipanel = int(panel_id(map_pid(1)), i_def)

  do k = 0, nlayers-1
  ! Extract element arrays of chi
    do df = 1, ndf_chi
      loc = map_chi(df) + k
      chi_1_e(df) = chi_1( loc )
      chi_2_e(df) = chi_2( loc )
      chi_3_e(df) = chi_3( loc )
    end do

    ! The term to be computed is:
    ! J*v . ( J^-T*vorticity cross J*u ) /det(J)
    ! This can be simplified through the vector triple product to
    ! J^-T*vorticity . ( J*u cross J*v ) /det(J)
    ! Which can again be simplified by pulling out the J factor to
    ! J^-T*vorticity . J^-T*( u cross v )
    ! Or equivalently
    ! [(J^-1 * J^-T)*vorticity].( u cross v ) = v.([(J^-1 * J^-T)*vorticity] cross u)
    ! This aviods having any matrix multiplies of cross products inside a dof loop
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h

        vorticity_at_quad(:) = 0.0_r_def
        do df = 1, ndf_w1
          vorticity_at_quad(:) = vorticity_at_quad(:) &
                               + vorticity( map_w1(df) + k )*w1_basis(:,df,qp1,qp2)
        end do
        call pointwise_coordinate_jacobian(ndf_chi, chi_1_e, chi_2_e, chi_3_e,  &
                                           ipanel, chi_basis(:,:,qp1,qp2),      &
                                           chi_diff_basis(:,:,qp1,qp2), jac, dj)
        jac_inv =  pointwise_coordinate_jacobian_inverse(jac, dj)
        jac = matmul(jac_inv,transpose(jac_inv))
        j_vorticity = wqp_h(qp1)*wqp_v(qp2)*matmul(jac,vorticity_at_quad)

        u_at_quad(:) = 0.0_r_def
        do df = 1, ndf_w2
          u_at_quad(:) = u_at_quad(:) + wind( map_w2(df) + k )*w2_basis(:,df,qp1,qp2)
        end do
        vorticity_term = cross_product(j_vorticity,u_at_quad)

        do df = 1, ndf_w2
          r_u( map_w2(df) + k ) = r_u( map_w2(df) + k ) &
                                - dot_product(w2_basis(:,df,qp1,qp2), &
                                              vorticity_term)
        end do
      end do
    end do
  end do

end subroutine vorticity_advection_code

end module vorticity_advection_kernel_mod
