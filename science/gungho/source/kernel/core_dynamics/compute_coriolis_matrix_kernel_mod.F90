!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Compute the Coriolis operator to apply the rotation vector Omega to
!>        the wind fields.
!> @details The form of the Coriolis operator is:
!!          \f[ < v, 2\Omega \times u > \f] where v is the test function for
!!          the velocity space W2, u corresponds to the wind field in W2 to
!!          which this operator will be applied, and Omega is the rotation
!!          vector of the domain.
!!          Multiplying this operator by a wind field will assemble the weak
!!          form of the Coriolis term.

module compute_coriolis_matrix_kernel_mod

use constants_mod,           only: i_def, r_def
use kernel_mod,              only: kernel_type
use argument_mod,            only: arg_type, func_type,       &
                                   GH_OPERATOR, GH_FIELD,     &
                                   GH_READ, GH_WRITE,         &
                                   GH_REAL, GH_SCALAR,        &
                                   ANY_SPACE_9,               &
                                   ANY_DISCONTINUOUS_SPACE_3, &
                                   GH_BASIS, GH_DIFF_BASIS,   &
                                   CELL_COLUMN, GH_QUADRATURE_XYoZ
use fs_continuity_mod,       only: W2

use sci_coordinate_jacobian_mod, only: coordinate_jacobian
use base_mesh_config_mod,    only: geometry,           &
                                   geometry_spherical
use rotation_vector_mod,     only: rotation_vector_fplane,  &
                                   rotation_vector_sphere
use cross_product_mod,       only: cross_product

implicit none
private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: compute_coriolis_matrix_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                                       &
       arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W2, W2),                    &
       arg_type(GH_FIELD*3,  GH_REAL, GH_READ,  ANY_SPACE_9),               &
       arg_type(GH_FIELD,    GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), &
       arg_type(GH_SCALAR,   GH_REAL, GH_READ),                             &
       arg_type(GH_SCALAR,   GH_REAL, GH_READ)                              &
       /)
  type(func_type) :: meta_funcs(2) = (/                                    &
       func_type(W2,          GH_BASIS),                                   &
       func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                     &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_QUADRATURE_XYoZ
contains
  procedure, nopass :: compute_coriolis_matrix_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: compute_coriolis_matrix_code
contains

!> @brief Compute the Coriolis operator to apply the rotation vector Omega to
!>        the wind fields.
!!
!! @param[in] cell     Identifying number of cell.
!! @param[in] nlayers  Number of layers.
!! @param[in] ncell_3d ncell*ndf
!! @param[in,out] matrix   Local stencil or Coriolis operator.
!! @param[in] chi_1 1st coordinate field
!! @param[in] chi_2 2nd coordinate field
!! @param[in] chi_3 3rd coordinate field
!! @param[in] panel_id A field giving the ID for mesh panels.
!! @param[in] omega    Planet angular velocity
!! @param[in] f_lat    F-plane latitude
!! @param[in] ndf      Degrees of freedom per cell.
!! @param[in] basis    Vector basis functions evaluated at quadrature points.
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi
!! @param[in] undf_chi Number of unique degrees of freedom  for chi
!! @param[in] map_chi Dofmap for the cell at the base of the column for chi
!! @param[in] chi_basis Basis functions for Wchi evaluated at
!!                      gaussian quadrature points
!! @param[in] chi_diff_basis Differential of the Wchi basis functions
!!                           evaluated at gaussian quadrature point
!! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
!! @param[in] undf_pid Number of unique degrees of freedom for panel_id
!! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
!! @param[in] nqp_h    Number of horizontal quadrature points.
!! @param[in] nqp_v    Number of vertical quadrature points.
!! @param[in] wqp_h    Horizontal quadrature weights.
!! @param[in] wqp_v    Vertical quadrature weights.
subroutine compute_coriolis_matrix_code(cell, nlayers, ncell_3d,           &
                                        matrix,                            &
                                        chi_1, chi_2, chi_3,               &
                                        panel_id,                          &
                                        omega, f_lat,                      &
                                        ndf, basis,                        &
                                        ndf_chi, undf_chi,                 &
                                        map_chi,                           &
                                        basis_chi, diff_basis_chi,         &
                                        ndf_pid, undf_pid, map_pid,        &
                                        nqp_h, nqp_v, wqp_h, wqp_v)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: cell, ndf, ndf_pid, ndf_chi
  integer(kind=i_def), intent(in)    :: undf_pid, undf_chi
  integer(kind=i_def), intent(in)    :: nqp_h, nqp_v
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ncell_3d

  real(kind=r_def), dimension(3,ndf,nqp_h,nqp_v), intent(in) :: basis

  integer(kind=i_def), intent(in)    :: map_chi(ndf_chi)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
  real(kind=r_def),    intent(inout) :: matrix(ncell_3d,ndf,ndf)
  real(kind=r_def),    intent(in)    :: basis_chi(1,ndf_chi,nqp_h,nqp_v)
  real(kind=r_def),    intent(in)    :: diff_basis_chi(3,ndf_chi,nqp_h,nqp_v)
  real(kind=r_def),    intent(in)    :: chi_1(undf_chi)
  real(kind=r_def),    intent(in)    :: chi_2(undf_chi)
  real(kind=r_def),    intent(in)    :: chi_3(undf_chi)
  real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
  real(kind=r_def),    intent(in)    :: wqp_h(nqp_h)
  real(kind=r_def),    intent(in)    :: wqp_v(nqp_v)
  real(kind=r_def),    intent(in)    :: omega
  real(kind=r_def),    intent(in)    :: f_lat

  ! Internal variables
  integer(kind=i_def)                          :: df, df2, k, ik
  integer(kind=i_def)                          :: qp1, qp2

  real(kind=r_def), dimension(ndf_chi)         :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac
  real(kind=r_def), dimension(3,nqp_h,nqp_v)   :: rotation_vector
  real(kind=r_def), dimension(3)               :: omega_cross_u
  real(kind=r_def), dimension(3)               :: jac_u, jac_v

  integer(kind=i_def) :: ipanel

  ipanel = int(panel_id(map_pid(1)), i_def)

  ! Loop over layers: Start from 1 as in this loop k is not an offset
  do k = 1, nlayers

     ! Indirect the chi coord field here
     do df = 1, ndf_chi
        chi_1_e(df) = chi_1(map_chi(df) + k - 1)
        chi_2_e(df) = chi_2(map_chi(df) + k - 1)
        chi_3_e(df) = chi_3(map_chi(df) + k - 1)
     end do

    ! Calculate planet's rotation vector
    if ( geometry == geometry_spherical ) then
      call rotation_vector_sphere(ndf_chi, nqp_h, nqp_v, chi_1_e, chi_2_e,     &
                                  chi_3_e, ipanel, basis_chi, rotation_vector)
    else
      call rotation_vector_fplane(nqp_h, nqp_v, omega, f_lat, rotation_vector)
    end if

    ! Calculate the Jacobian and its determinant
    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v,                            &
                             chi_1_e, chi_2_e, chi_3_e, ipanel,                &
                             basis_chi, diff_basis_chi, jac, dj)


    ! To convert from reference space to physical space:
    ! Let the volume element be dV and the unit volume be dV_hat, then:
    ! u = J u_hat/detJ, v = J v_hat/detJ, and dV = detJ dV_hat,
    ! where hats denote the reference element space, J is the Jacobian, and
    ! detJ its determinant. This gives
    ! \int (J w_hat) . (2 Omega x J u_hat) / detJ dV_hat
    ! Note that some of the detJ factors cancel
    ik = k + (cell-1)*nlayers
    matrix(ik,:,:) = 0.0_r_def
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        do df2 = 1, ndf
          jac_u = matmul(jac(:,:,qp1,qp2),basis(:,df2,qp1,qp2))
          omega_cross_u = wqp_h(qp1) * wqp_v(qp2)                              &
              * cross_product(rotation_vector(:,qp1,qp2), jac_u) / dj(qp1,qp2)
          do df = 1, ndf
            jac_v = matmul(jac(:,:,qp1,qp2), basis(:,df,qp1,qp2))
            matrix(ik,df,df2) = matrix(ik,df,df2) - dot_product(jac_v, omega_cross_u)
          end do
        end do
      end do
    end do
  end do ! end of k loop

end subroutine compute_coriolis_matrix_code

end module compute_coriolis_matrix_kernel_mod
