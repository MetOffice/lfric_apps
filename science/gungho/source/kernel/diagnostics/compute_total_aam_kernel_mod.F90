!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Computes the cell integrated axial angular momentum.
!>
!> @details The kernel computes the  cell integrated axial angular momentum:
!> \f[ \int ( \rho  [u + \Omega r \cos(\phi) ] r \cos(\phi) ) dV \f]
!> where \phi is the latitude and u is the zonal velocity.
!> From vector components this can be calculated from
!> \f[ \int ( \rho * \hat{z} . [\textbf{r} \times { \textbf{v} + \Omega \times \textbf{r} } ] ) dV \f]
!> where \textbf{v} is now the velocity vector and \hat{z} is aligned with the pole.
!>
module compute_total_aam_kernel_mod

  use argument_mod,      only : arg_type, func_type,         &
                                GH_FIELD, GH_READ, GH_WRITE, &
                                GH_REAL, GH_SCALAR,          &
                                ANY_SPACE_9,                 &
                                ANY_DISCONTINUOUS_SPACE_3,   &
                                GH_BASIS, GH_DIFF_BASIS,     &
                                CELL_COLUMN, GH_QUADRATURE_XYoZ
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W2, W3
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: compute_total_aam_kernel_type
    private
    type(arg_type) :: meta_args(6) = (/                                      &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, W3),                        &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  W2),                        &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),                        &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ,  ANY_SPACE_9),               &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ)                              &
         /)
    type(func_type) :: meta_funcs(3) = (/                                   &
         func_type(W2,          GH_BASIS),                                  &
         func_type(W3,          GH_BASIS),                                  &
         func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                    &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: compute_total_aam_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: compute_total_aam_code

contains

!> @brief The subroutine to compute the total axial angular momentum
!! @param[in] nlayers Number of layers
!! @param[in,out] aam Cell integrated axial angular momentum
!! @param[in] u Velocity array
!! @param[in] rho density
!! @param[in] chi_1 1st coordinate field
!! @param[in] chi_2 2nd coordinate field
!! @param[in] chi_3 3rd coordinate field
!! @param[in] panel_id Field giving the ID for mesh panels
!! @param[in] omega Planet angular velocity
!! @param[in] ndf_w2 Number of degrees of freedom per cell for w2
!! @param[in] undf_w2 Number of unique degrees of freedom  for w2
!! @param[in] map_w2 Dofmap for the cell at the base of the column for w2
!! @param[in] w2_basis Basis functions evaluated at quadrature points
!! @param[in] ndf_w3 Number of degrees of freedom per cell for w3
!! @param[in] undf_w3 Number of unique degrees of freedom  for w3
!! @param[in] map_w3 Dofmap for the cell at the base of the column for w3
!! @param[in] w3_basis Basis functions evaluated at gaussian quadrature points
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi
!! @param[in] undf_chi Number of unique degrees of freedom for chi
!! @param[in] map_chi Dofmap for the cell at the base of the column for chi
!! @param[in] chi_basis Basis functions for Wchi evaluated at
!!                      gaussian quadrature points
!! @param[in] chi_diff_basis Differential of the Wchi basis functions
!!                           evaluated at gaussian quadrature points
!! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
!! @param[in] undf_pid Number of unique degrees of freedom for panel_id
!! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
!! @param[in] nqp_h Number of quadrature points in the horizontal
!! @param[in] nqp_v Number of quadrature points in the vertical
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine compute_total_aam_code(                                           &
                                  nlayers,                                   &
                                  aam, u, rho,                               &
                                  chi_1, chi_2, chi_3, panel_id,             &
                                  omega,                                     &
                                  ndf_w3, undf_w3, map_w3, w3_basis,         &
                                  ndf_w2, undf_w2, map_w2, w2_basis,         &
                                  ndf_chi, undf_chi, map_chi,                &
                                  chi_basis, chi_diff_basis,                 &
                                  ndf_pid, undf_pid, map_pid,                &
                                  nqp_h, nqp_v, wqp_h, wqp_v                 &
                                 )

  use sci_coordinate_jacobian_mod, only: coordinate_jacobian
  use sci_chi_transform_mod,       only: chi2llr, chi2xyz
  use coord_transform_mod,         only: cart2sphere_vector
  use cross_product_mod,           only: cross_product

  implicit none

  ! Arguments
  integer(kind=i_def),             intent(in) :: nlayers, nqp_h, nqp_v
  integer(kind=i_def),             intent(in) :: ndf_w2, ndf_w3, ndf_pid
  integer(kind=i_def),             intent(in) :: ndf_chi
  integer(kind=i_def),             intent(in) :: undf_chi
  integer(kind=i_def),             intent(in) :: undf_w2, undf_w3, undf_pid
  integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_w2),  intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_w3),  intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_pid), intent(in) :: map_pid

  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v),  intent(in) :: w3_basis
  real(kind=r_def), dimension(3,ndf_w2,nqp_h,nqp_v),  intent(in) :: w2_basis
  real(kind=r_def), dimension(1,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_basis
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_diff_basis

  real(kind=r_def), dimension(undf_w3),      intent(inout) :: aam
  real(kind=r_def), dimension(undf_w2),      intent(in)    :: u
  real(kind=r_def), dimension(undf_w3),      intent(in)    :: rho
  real(kind=r_def), dimension(undf_chi),     intent(in)    :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(undf_pid),     intent(in)    :: panel_id
  real(kind=r_def),                          intent(in)    :: omega

  real(kind=r_def), dimension(nqp_h), intent(in) ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in) ::  wqp_v

  ! Internal variables
  integer(kind=i_def) :: df, k, loc
  integer(kind=i_def) :: qp1, qp2
  integer(kind=i_def) :: ipanel

  real(kind=r_def), dimension(ndf_chi)         :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac
  real(kind=r_def), dimension(ndf_w3)          :: rho_e, aam_e
  real(kind=r_def), dimension(ndf_w2)          :: u_e

  real(kind=r_def) :: u_at_quad(3), omega_vec(3), llr_vec(3), coords(3)
  real(kind=r_def) :: am(3), x_vec(3), u_vec(3), j_u(3), r_vec(3), spherical_z_hat(3)
  real(kind=r_def) :: rho_at_quad
  real(kind=r_def), parameter :: z_hat(3) = (/ 0.0_r_def, 0.0_r_def, 1.0_r_def /)

  ipanel = int(panel_id(map_pid(1)), i_def)

  do k = 0, nlayers-1
  ! Extract element arrays of chi
    do df = 1, ndf_chi
      loc = map_chi(df) + k
      chi_1_e(df) = chi_1( loc )
      chi_2_e(df) = chi_2( loc )
      chi_3_e(df) = chi_3( loc )
    end do

    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v,             &
                             chi_1_e, chi_2_e, chi_3_e, ipanel, &
                             chi_basis, chi_diff_basis, jac, dj)
    do df = 1, ndf_w3
      rho_e(df) = rho( map_w3(df) + k )
      aam_e(df) = 0.0_r_def
    end do
    do df = 1, ndf_w2
      u_e(df) = u( map_w2(df) + k )
    end do
  ! compute the aam integrated over one cell
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        coords(:) = 0.0_r_def
        do df = 1, ndf_chi
          coords(1) = coords(1) + chi_1_e(df)*chi_basis(1,df,qp1,qp2)
          coords(2) = coords(2) + chi_2_e(df)*chi_basis(1,df,qp1,qp2)
          coords(3) = coords(3) + chi_3_e(df)*chi_basis(1,df,qp1,qp2)
        end do

        ! Obtain (X,Y,Z) and (lon,lat,r) coords
        call chi2xyz(coords(1), coords(2), coords(3),  &
                     ipanel, x_vec(1), x_vec(2), x_vec(3))
        call chi2llr(coords(1), coords(2), coords(3),  &
                     ipanel, llr_vec(1), llr_vec(2), llr_vec(3))

        ! get position vector with spherical components
        r_vec(:) = cart2sphere_vector(x_vec, x_vec)
        spherical_z_hat(:) = cart2sphere_vector(x_vec, z_hat)

        ! get Omega vector expressed in (long,lat,r) components
        omega_vec(1) = 0.0_r_def
        omega_vec(2) = omega*cos(llr_vec(2))
        omega_vec(3) = omega*sin(llr_vec(2))

        rho_at_quad = 0.0_r_def
        do df = 1, ndf_w3
          rho_at_quad  = rho_at_quad + rho_e(df)*w3_basis(1,df,qp1,qp2)
        end do

        u_at_quad(:) = 0.0_r_def
        do df = 1, ndf_w2
          u_at_quad(:) = u_at_quad(:) &
                       + u_e(df)*w2_basis(:,df,qp1,qp2)
        end do

        ! integral is
        ! rho * cross(r, u + cross(Omega, r)) dV
        ! transforming to reference element this gives
        ! rho' * cross(r, J*u/dj + cross(Omega, r)) * dj * dV
        j_u = matmul(jac(:,:,qp1,qp2),u_at_quad)
        u_vec(:) = cart2sphere_vector(x_vec,j_u) &
                 + cross_product(omega_vec,r_vec)*dj(qp1,qp2)

        am(:) = cross_product(r_vec,u_vec)

        do df = 1, ndf_w3
          aam_e(df) = aam_e(df) &
                      + wqp_h(qp1)*wqp_v(qp2)*rho_at_quad*dot_product(spherical_z_hat,am)
        end do
      end do
    end do
    do df = 1,ndf_w3
      aam(map_w3(df)+k) = aam_e(df)
    end do
  end do

end subroutine compute_total_aam_code

end module compute_total_aam_kernel_mod
