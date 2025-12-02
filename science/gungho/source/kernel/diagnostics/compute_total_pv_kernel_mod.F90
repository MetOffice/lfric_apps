!-----------------------------------------------------------------------------
! Copyright (c) 2022,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Computes the cell integrated potential vorticity.
!>
!> \f$ \int( (\xi + Omega) \cdot \nabla(\theta)/\rho dV ) \f$
!>
module compute_total_pv_kernel_mod

  use argument_mod,      only : arg_type, func_type,         &
                                GH_FIELD, GH_WRITE, GH_READ, &
                                GH_REAL, ANY_SPACE_9,        &
                                ANY_DISCONTINUOUS_SPACE_3,   &
                                GH_BASIS, GH_DIFF_BASIS,     &
                                CELL_COLUMN, GH_SCALAR,       &
                                GH_QUADRATURE_XYoZ
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W0, W1, W3
  use kernel_mod,        only : kernel_type
  use base_mesh_config_mod, &
                         only: geometry, &
                               geometry_spherical
  use rotation_vector_mod, &
                         only: rotation_vector_fplane, &
                               rotation_vector_sphere

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: compute_total_pv_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                                     &
        arg_type(GH_FIELD,   GH_REAL, GH_WRITE, W3),                        &
        arg_type(GH_FIELD,   GH_REAL, GH_READ,  W1),                        &
        arg_type(GH_FIELD,   GH_REAL, GH_READ,  W0),                        &
        arg_type(GH_FIELD,   GH_REAL, GH_READ, W3),                         &
        arg_type(GH_FIELD*3, GH_REAL, GH_READ,  ANY_SPACE_9),               &
        arg_type(GH_FIELD,   GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), &
        arg_type(GH_SCALAR,   GH_REAL, GH_READ),                            &
        arg_type(GH_SCALAR,   GH_REAL, GH_READ)                             &
        /)
    type(func_type) :: meta_funcs(4) = (/                                  &
        func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS),                   &
        func_type(W0,          GH_DIFF_BASIS),                             &
        func_type(W1,          GH_BASIS),                                  &
        func_type(W3,          GH_BASIS)                                   &
        /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: compute_total_pv_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: compute_total_pv_code

contains

!> @brief The kernel computes the cell integrated potential vorticity
!! @param[in] nlayers   Number of layers
!! @param[in,out] pv Cell integrated potential vorticity
!! @param[in] xi        Absolute vorticity
!! @param[in] theta     Potential temperature
!! @param[in] rho       Density
!! @param[in] chi1      1st coordinate field in Wchi
!! @param[in] chi2      2nd coordinate field in Wchi
!! @param[in] chi3      3rd coordinate field in Wchi
!! @param[in] panel_id  Field giving the ID for mesh panels.
!! @param[in] omega     Planet angular velocity
!! @param[in] f_lat     F-plane latitude
!! @param[in] ndf_w3    Number of degrees of freedom per cell for w3
!! @param[in] undf_w3   Number of unique degrees of freedom  for w3
!! @param[in] map_w3    Dofmap for the cell at the base of the column for w3
!! @param[in] w3_basis  Basis functions evaluated at gaussian quadrature points
!! @param[in] ndf_w1    Number of degrees of freedom per cell for w1
!! @param[in] undf_w1   Number of unique degrees of freedom  for w1
!! @param[in] map_w1    Dofmap for the cell at the base of the column for w1
!! @param[in] w1_basis  Basis functions evaluated at gaussian quadrature points
!! @param[in] ndf_w0    Number of degrees of freedom per cell for w0
!! @param[in] undf_w0   Number of unique degrees of freedom  for w0
!! @param[in] map_w0    Dofmap for the cell at the base of the column for w0
!! @param[in] w0_diff_basis Differential basis functions evaluated at gaussian quadrature points
!! @param[in] ndf_chi   Number of degrees of freedom per cell for chi
!! @param[in] undf_chi  Number of unique degrees of freedom  for chi
!! @param[in] map_chi   Dofmap for the cell at the base of the column for chi
!! @param[in] chi_basis Wchi basis functions evaluated at gaussian quadrature points
!! @param[in] chi_diff_basis Derivatives of Wchi basis functions
!!                           evaluated at gaussian quadrature points
!! @param[in] ndf_pid   Number of degrees of freedom per cell for panel_id
!! @param[in] undf_pid  Number of unique degrees of freedom for panel_id
!! @param[in] map_pid   Dofmap for the cell at the base of the column for panel_id
!! @param[in] nqp_h     Number of horizontal quadrature points
!! @param[in] nqp_v     Number of vertical quadrature points
!! @param[in] wqp_h     Weights of the horizontal quadrature points
!! @param[in] wqp_v     Weights of the vertical quadrature points
subroutine compute_total_pv_code(                                                        &
                                 nlayers,                                                &
                                 pv,                                                     &
                                 xi,                                                     &
                                 theta,                                                  &
                                 rho,                                                    &
                                 chi1, chi2, chi3, panel_id,                             &
                                 omega, f_lat,                                           &
                                 ndf_w3, undf_w3, map_w3, w3_basis,                      &
                                 ndf_w1, undf_w1, map_w1, w1_basis,                      &
                                 ndf_w0, undf_w0, map_w0, w0_diff_basis,                 &
                                 ndf_chi, undf_chi, map_chi, chi_basis, chi_diff_basis,  &
                                 ndf_pid, undf_pid, map_pid,                             &
                                 nqp_h, nqp_v, wqp_h, wqp_v )

  use sci_coordinate_jacobian_mod, only: coordinate_jacobian, &
                                         coordinate_jacobian_inverse

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, nqp_h, nqp_v
  integer(kind=i_def), intent(in) :: ndf_w0, ndf_w1, ndf_w3, ndf_chi, ndf_pid
  integer(kind=i_def), intent(in) :: undf_w0, undf_w1, undf_w3, undf_chi, undf_pid

  integer(kind=i_def), dimension(ndf_w3),  intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_w0),  intent(in) :: map_w0
  integer(kind=i_def), dimension(ndf_w1),  intent(in) :: map_w1
  integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_pid), intent(in) :: map_pid

  real(kind=r_def), dimension(3,ndf_w0,nqp_h,nqp_v),  intent(in) :: w0_diff_basis
  real(kind=r_def), dimension(1,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_basis
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_diff_basis
  real(kind=r_def), dimension(3,ndf_w1,nqp_h,nqp_v),  intent(in) :: w1_basis
  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v),  intent(in) :: w3_basis

  real(kind=r_def), dimension(undf_w3),  intent(inout) :: pv
  real(kind=r_def), dimension(undf_w0),  intent(in) :: theta
  real(kind=r_def), dimension(undf_w3),  intent(in) :: rho
  real(kind=r_def), dimension(undf_chi), intent(in) :: chi1, chi2, chi3
  real(kind=r_def), dimension(undf_pid), intent(in) :: panel_id
  real(kind=r_def), dimension(undf_w1),  intent(in) :: xi

  real(kind=r_def), dimension(nqp_h),    intent(in) ::  wqp_h
  real(kind=r_def), dimension(nqp_v),    intent(in) ::  wqp_v

  real(kind=r_def),    intent(in)    :: omega
  real(kind=r_def),    intent(in)    :: f_lat

  ! Internal variables
  integer(kind=i_def) :: df, k, ipanel
  integer(kind=i_def) :: qp1, qp2

  real(kind=r_def), dimension(ndf_chi)         :: chi1_e, chi2_e, chi3_e
  real(kind=r_def), dimension(ndf_w0)          :: theta_e
  real(kind=r_def), dimension(ndf_w3)          :: rho_e
  real(kind=r_def), dimension(ndf_w1)          :: xi_e
  real(kind=r_def), dimension(ndf_w3)          :: pv_e
  real(kind=r_def), dimension(3)               :: xi_at_quad
  real(kind=r_def), dimension(3)               :: grad_theta_at_quad
  real(kind=r_def)                             :: rho_at_quad
  real(kind=r_def)                             :: pv_at_quad
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac, jac_inv
  real(kind=r_def), dimension(3,nqp_h,nqp_v)   :: rotation_vector

  ipanel = int(panel_id(map_pid(1)), i_def)

  do k = 0, nlayers-1
  ! Extract element arrays of chi and theta
    do df = 1, ndf_chi
      chi1_e(df) = chi1( map_chi(df) + k )
      chi2_e(df) = chi2( map_chi(df) + k )
      chi3_e(df) = chi3( map_chi(df) + k )
    end do
    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi1_e, chi2_e, chi3_e,  &
                             ipanel, chi_basis, chi_diff_basis, jac, dj)
    call coordinate_jacobian_inverse(nqp_h, nqp_v, jac, dj, jac_inv)

    ! Calculate rotation vector Omega = (0, 2*cos(lat), 2*sin(lat))
    if ( geometry == geometry_spherical ) then
      call rotation_vector_sphere(ndf_chi, nqp_h, nqp_v, chi1_e, chi2_e,    &
                                  chi3_e, ipanel, chi_basis, rotation_vector)
    else
      call rotation_vector_fplane(nqp_h, nqp_v, omega, f_lat, &
                                  rotation_vector)
    end if

    do df = 1, ndf_w0
      theta_e(df)  = theta( map_w0(df) + k )
    end do

    do df = 1, ndf_w1
      xi_e(df) = xi( map_w1(df) + k )
    end do
    do df = 1, ndf_w3
      rho_e(df)  = rho( map_w3(df) + k )
    end do
    pv_e(:) = 0.0_r_def
  ! compute the pv integrated over one cell
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        xi_at_quad(:) = 0.0_r_def
        do df = 1, ndf_w1
          xi_at_quad(:)  = xi_at_quad(:)  + xi_e(df)*w1_basis(:,df,qp1,qp2)
        end do
        grad_theta_at_quad(:) = 0.0_r_def
        do df = 1, ndf_w0
          grad_theta_at_quad(:) = grad_theta_at_quad(:) &
                                + theta_e(df)*w0_diff_basis(:,df,qp1,qp2)
        end do
        rho_at_quad = 0.0_r_def
        do df = 1, ndf_w3
          rho_at_quad= rho_at_quad &
                                + rho_e(df)*w3_basis(1,df,qp1,qp2)
        end do
        pv_at_quad = wqp_h(qp1)*wqp_v(qp2)*dj(qp1,qp2) &
        * dot_product(matmul(transpose(jac_inv(:,:,qp1,qp2)),xi_at_quad) + rotation_vector(:,qp1,qp2), &
                      matmul(transpose(jac_inv(:,:,qp1,qp2)),grad_theta_at_quad))/rho_at_quad
        do df = 1,ndf_w3
          pv_e(df)= pv_e(df)+ pv_at_quad*w3_basis(1,df,qp1,qp2)
        end do
      end do
    end do
    do df = 1,ndf_w3
      pv(map_w3(df)+k)= pv_e(df)
    end do
  end do
end subroutine compute_total_pv_code

end module compute_total_pv_kernel_mod
