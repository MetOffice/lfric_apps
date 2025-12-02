!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Compute the vertical part of the Coriolis operator to apply the
!!        rotation vector Omega to the wind fields, leaving only the vertical
!!        component in Wtheta.
!> @details The form of vertical component of the Coriolis operator is:
!!          \f[ < w, k . (2\Omega \times u) > \f] where w is the test function
!!          in Wtheta, u corresponds to the wind field in W2 to which this
!!          operator will be applied, k is the vertical unit vector and
!!          Omega is the rotation vector of the domain.
!!          Multiplying this operator by a wind field will assemble the weak
!!          form of the vertical part of the Coriolis term.

module compute_vert_coriolis_matrix_kernel_mod

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
use fs_continuity_mod,       only: W2, Wtheta

use sci_coordinate_jacobian_mod, only: coordinate_jacobian
use base_mesh_config_mod,    only: geometry,           &
                                   geometry_spherical
use rotation_vector_mod,     only: rotation_vector_fplane,  &
                                   rotation_vector_sphere,  &
                                   vert_vector_sphere
use cross_product_mod,       only: cross_product

implicit none
private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: compute_vert_coriolis_matrix_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                                          &
      arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, Wtheta, W2),                    &
      arg_type(GH_FIELD*3,  GH_REAL, GH_READ,  ANY_SPACE_9),                   &
      arg_type(GH_FIELD,    GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),     &
      arg_type(GH_SCALAR,   GH_REAL, GH_READ),                                 &
      arg_type(GH_SCALAR,   GH_REAL, GH_READ)                                  &
  /)
  type(func_type) :: meta_funcs(3) = (/                                        &
      func_type(Wtheta,      GH_BASIS),                                        &
      func_type(W2,          GH_BASIS),                                        &
      func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                          &
  /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_QUADRATURE_XYoZ
contains
  procedure, nopass :: compute_vert_coriolis_matrix_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: compute_vert_coriolis_matrix_code
contains

!> @brief Compute the vertical part of the Coriolis operator to apply the
!>        rotation vector Omega to the wind fields, leaving only the vertical
!>        component in Wtheta.
!> @param[in]     col_idx         Index of column
!> @param[in]     nlayers         Number of layers in mesh
!> @param[in]     ncell_3d        Total number of 3D cells in mesh
!> @param[in,out] matrix          Vertical Coriolis operator
!> @param[in]     chi_1           1st coordinate field
!> @param[in]     chi_2           2nd coordinate field
!> @param[in]     chi_3           3rd coordinate field
!> @param[in]     panel_id        A field giving the ID for mesh panels
!> @param[in]     omega           Planet angular velocity
!> @param[in]     f_lat           F-plane latitude
!> @param[in]     ndf_wt          Num of DoFs per cell for Wtheta
!> @param[in]     basis_wt        Wtheta basis funcs evaluated at quad points
!> @param[in]     ndf_w2          Num of DoFs per cell for W2
!> @param[in]     basis_w2        W2 basis functions evaluated at quad points
!> @param[in]     ndf_chi         Num of DoFs per cell for Wchi
!> @param[in]     undf_chi        Num of DoFs in this partition for Wchi
!> @param[in]     map_chi         Dofmap for Wchi
!> @param[in]     basis_chi       Wchi basis functions evaluated at quad points
!> @param[in]     diff_basis_chi  Wchi differential basis funcs at quad points
!> @param[in]     ndf_chi         Num of DoFs per cell for Wchi
!> @param[in]     undf_chi        Num of DoFs in this partition for Wchi
!> @param[in]     map_chi         Dofmap for Wchi
!> @param[in]     nqp_h           Number of horizontal quadrature points
!> @param[in]     nqp_v           Number of vertical quadrature points
!> @param[in]     wqp_h           Horizontal quadrature weights
!> @param[in]     wqp_v           Vertical quadrature weights
subroutine compute_vert_coriolis_matrix_code(col_idx, nlayers, ncell_3d,       &
                                             matrix,                           &
                                             chi_1, chi_2, chi_3,              &
                                             panel_id,                         &
                                             omega, f_lat,                     &
                                             ndf_wt, basis_wt,                 &
                                             ndf_w2, basis_w2,                 &
                                             ndf_chi, undf_chi,                &
                                             map_chi,                          &
                                             basis_chi, diff_basis_chi,        &
                                             ndf_pid, undf_pid, map_pid,       &
                                             nqp_h, nqp_v, wqp_h, wqp_v)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: col_idx, ndf_pid, ndf_chi
  integer(kind=i_def), intent(in)    :: undf_pid, undf_chi
  integer(kind=i_def), intent(in)    :: ndf_wt, ndf_w2
  integer(kind=i_def), intent(in)    :: nqp_h, nqp_v
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ncell_3d

  real(kind=r_def),    intent(in)    :: basis_w2(3,ndf_w2,nqp_h,nqp_v)
  real(kind=r_def),    intent(in)    :: basis_wt(1,ndf_wt,nqp_h,nqp_v)
  real(kind=r_def),    intent(in)    :: basis_chi(1,ndf_chi,nqp_h,nqp_v)
  real(kind=r_def),    intent(in)    :: diff_basis_chi(3,ndf_chi,nqp_h,nqp_v)

  integer(kind=i_def), intent(in)    :: map_chi(ndf_chi)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
  real(kind=r_def),    intent(inout) :: matrix(ncell_3d,ndf_wt,ndf_w2)
  real(kind=r_def),    intent(in)    :: chi_1(undf_chi)
  real(kind=r_def),    intent(in)    :: chi_2(undf_chi)
  real(kind=r_def),    intent(in)    :: chi_3(undf_chi)
  real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
  real(kind=r_def),    intent(in)    :: wqp_h(nqp_h)
  real(kind=r_def),    intent(in)    :: wqp_v(nqp_v)
  real(kind=r_def),    intent(in)    :: omega
  real(kind=r_def),    intent(in)    :: f_lat

  ! Internal variables
  integer(kind=i_def) :: df_wt, df_w2, df_chi, k, ik
  integer(kind=i_def) :: qp1, qp2

  real(kind=r_def) :: chi_1_e(ndf_chi), chi_2_e(ndf_chi), chi_3_e(ndf_chi)
  real(kind=r_def) :: dj(nqp_h,nqp_v)
  real(kind=r_def) :: jac(3,3,nqp_h,nqp_v)
  real(kind=r_def) :: rotation_vector(3,nqp_h,nqp_v)
  real(kind=r_def) :: vert_vec(3,nqp_h,nqp_v)
  real(kind=r_def) :: k_dot_omega_cross_u
  real(kind=r_def) :: jac_u(3)

  integer(kind=i_def) :: ipanel

  ipanel = int(panel_id(map_pid(1)), i_def)

  ! Loop over layers: Start from 1 as in this loop k is not an offset
  do k = 1, nlayers
    ! Extract the coordinates for this cell
    do df_chi = 1, ndf_chi
      chi_1_e(df_chi) = chi_1(map_chi(df_chi) + k - 1)
      chi_2_e(df_chi) = chi_2(map_chi(df_chi) + k - 1)
      chi_3_e(df_chi) = chi_3(map_chi(df_chi) + k - 1)
    end do

    ! Calculate planet's rotation vector and the vertical vector
    if ( geometry == geometry_spherical ) then
      call rotation_vector_sphere(ndf_chi, nqp_h, nqp_v, chi_1_e, chi_2_e,     &
                                  chi_3_e, ipanel, basis_chi, rotation_vector)
      call vert_vector_sphere(ndf_chi, nqp_h, nqp_v, chi_1_e, chi_2_e,         &
                              chi_3_e, ipanel, basis_chi, vert_vec)
    else
      call rotation_vector_fplane(nqp_h, nqp_v, omega, f_lat, rotation_vector)
      do qp2 = 1, nqp_v
        do qp1 = 1, nqp_h
          vert_vec(:,qp1,qp2) = (/ 0.0_r_def, 0.0_r_def, 1.0_r_def /)
        end do
      end do
    end if

    ! Calculate the Jacobian
    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v,                            &
                             chi_1_e, chi_2_e, chi_3_e, ipanel,                &
                             basis_chi, diff_basis_chi, jac, dj)

    ! To convert from reference space to physical space:
    ! Let the volume element be dV and the unit volume be dV_hat, then:
    ! w = w_hat, v = J v_hat/detJ, and dV = detJ dV_hat,
    ! where hats denote the reference element space, J is the Jacobian, and
    ! detJ its determinant. This gives:
    ! \int w_hat * k . (2 Omega x J u_hat) dV_hat
    ! Note that some of the detJ factors cancel
    ik = k + (col_idx-1)*nlayers
    matrix(ik,:,:) = 0.0_r_def
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        do df_w2 = 1, ndf_w2
          jac_u = matmul(jac(:,:,qp1,qp2), basis_w2(:,df_w2,qp1,qp2))
          k_dot_omega_cross_u = wqp_h(qp1) * wqp_v(qp2) * dot_product(         &
            vert_vec(:,qp1,qp2),                                               &
            cross_product(rotation_vector(:,qp1,qp2), jac_u)                   &
          )
          do df_wt = 1, ndf_wt
            matrix(ik,df_wt,df_w2) = matrix(ik,df_wt,df_w2)                    &
              - basis_wt(1,df_wt,qp1,qp2) * k_dot_omega_cross_u
          end do
        end do
      end do
    end do
  end do ! end of k loop

end subroutine compute_vert_coriolis_matrix_code

end module compute_vert_coriolis_matrix_kernel_mod
