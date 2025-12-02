!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the potential vorticity rhs for vorticity as a field
!!        in V2 (stored in W3).
!!
!> @details The potential vorticity is computed using a projection:
!!          \f[ q = P_{W_2}((\nabla \times P_{W_1}(u) + Omega)/geopot) \f]
!!
!!
module initial_vorticity_v2_kernel_mod

  use argument_mod,            only: arg_type, func_type,         &
                                     GH_FIELD, GH_READ, GH_WRITE, &
                                     GH_BASIS, GH_DIFF_BASIS,     &
                                     GH_INTEGER, CELL_COLUMN,     &
                                     GH_QUADRATURE_XYoZ,          &
                                     ANY_SPACE_9,                 &
                                     ANY_DISCONTINUOUS_SPACE_3,   &
                                     GH_REAL
  use constants_mod,           only: r_def, i_def
  use fs_continuity_mod,       only: W2, W3
  use kernel_mod,              only: kernel_type
  use sci_coordinate_jacobian_mod, only: coordinate_jacobian, &
                                         coordinate_jacobian_inverse
  use base_mesh_config_mod,    only: geometry,           &
                                     geometry_spherical, &
                                     f_lat
  use rotation_vector_mod,     only: rotation_vector_fplane,  &
                                     rotation_vector_sphere
  use planet_config_mod,       only: scaled_omega

  implicit none

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: initial_vorticity_v2_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                    &
        arg_type(GH_FIELD,   GH_REAL, GH_WRITE, W3),                       &
        arg_type(GH_FIELD,   GH_REAL, GH_READ,  W2),                       &
        arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),                       &
        arg_type(GH_FIELD*3, GH_REAL, GH_READ,  ANY_SPACE_9),              &
        arg_type(GH_FIELD,   GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3) &
        /)
    type(func_type) :: meta_funcs(3) = (/               &
        func_type(W3, GH_BASIS),                        &
        func_type(W2, GH_BASIS),                        &
        func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS) &
        /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, public, nopass :: initial_vorticity_v2_code
  end type

!---------------------------------------------------------------------------
! Contained functions/subroutines
!---------------------------------------------------------------------------
public initial_vorticity_v2_code

contains

!> @brief Compute the right-hand side of the potential vorticity diagnostic equation
!> @param[in]     nlayers        Number of layers
!> @param[in,out] r_q            Potential vorticity right hand side
!> @param[in]     curl_u         Strong curl of wind field projected into W1 space
!> @param[in]     geopot         Geopotential field
!> @param[in]     chi_1          1st coordinate field
!> @param[in]     chi_2          2nd coordinate field
!> @param[in]     chi_3          3rd coordinate field
!> @param[in]     panel_id       Field containing the ID of the mesh panel
!> @param[in]     ndf_w2         Number of degrees of freedom per cell for w2
!> @param[in]     undf_w2        Number unique of degrees of freedom  for w2
!> @param[in]     map_w2         Dofmap for the cell at the base of the column for w2
!> @param[in]     w2_basis       Basis functions evaluated at quadrature points
!> @param[in]     ndf_chi        Number of degrees of freedom per cell for chi
!> @param[in]     undf_chi       Number unique of degrees of freedom  for chi
!> @param[in]     map_chi        Dofmap for the cell at the base of the column for chi
!> @param[in]     chi_basis      Basis functions evaluated at quadrature points.
!> @param[in]     chi_diff_basis Differential of the basis functions evaluated at gaussian quadrature point
!> @param[in]     ndf_pid        The number of DoFs per cell for the panel ID
!> @param[in]     undf_pid       The number of DoFs for this partition for the panel ID
!> @param[in]     map_pid        DoF-map for the panel ID
!> @param[in]     nqp_h          Number of quadrature points in the horizontal
!> @param[in]     nqp_v          Number of quadrature points in the vertical
!> @param[in]     wqp_h          Horizontal quadrature weights
!> @param[in]     wqp_v          Vertical quadrature weights
subroutine initial_vorticity_v2_code(nlayers, r_q, curl_u, geopot,      &
                                     chi_1, chi_2, chi_3, panel_id,     &
                                     ndf_w3, undf_w3, map_w3, w3_basis, &
                                     ndf_w2, undf_w2, map_w2, w2_basis, &
                                     ndf_chi, undf_chi, map_chi,        &
                                     chi_basis, chi_diff_basis,         &
                                     ndf_pid, undf_pid, map_pid,        &
                                     nqp_h, nqp_v, wqp_h, wqp_v )

  implicit none

  ! Arguments
  integer, intent(in) :: nlayers, nqp_h, nqp_v
  integer, intent(in) :: ndf_chi, ndf_w2, ndf_w3, ndf_pid
  integer, intent(in) :: undf_chi, undf_w2, undf_w3, undf_pid
  integer, dimension(ndf_chi), intent(in)  :: map_chi
  integer, dimension(ndf_w2),  intent(in)  :: map_w2
  integer, dimension(ndf_w3),  intent(in)  :: map_w3
  integer(kind=i_def), dimension(ndf_pid),     intent(in) :: map_pid

  real(kind=r_def), dimension(1,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_basis
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_diff_basis
  real(kind=r_def), dimension(3,ndf_w2,nqp_h,nqp_v),  intent(in) :: w2_basis
  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v),  intent(in) :: w3_basis

  real(kind=r_def), dimension(undf_w3),  intent(inout) :: r_q
  real(kind=r_def), dimension(undf_w2),  intent(in)    :: curl_u
  real(kind=r_def), dimension(undf_w3),  intent(in)    :: geopot
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi_1, chi_2, chi_3

  real(kind=r_def), dimension(nqp_h), intent(in)      :: wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)      :: wqp_v

  ! Internal variables
  integer             :: df, loc
  integer             :: qp1, qp2
  real(kind=r_def), dimension(ndf_chi)         :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(undf_pid),  intent(in)    :: panel_id
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac, jac_inv

  real(kind=r_def), dimension(ndf_w3)          :: r_q_e
  real(kind=r_def), dimension(ndf_w2)          :: curl_u_e
  real(kind=r_def), dimension(ndf_w3)          :: geopot_e

  real(kind=r_def)                           :: geopot_at_quad
  real(kind=r_def), dimension(3)             :: curl_u_at_quad, jac_rot
  real(kind=r_def), dimension(3,nqp_h,nqp_v) :: rotation_vector

  integer(kind=i_def) :: ipanel

  ipanel = int(panel_id(map_pid(1)), i_def)

  ! Avoid compile warning for unused variables
  df = nlayers

  ! Extract element arrays of chi
  do df = 1, ndf_chi
    loc = map_chi(df)
    chi_1_e(df) = chi_1( loc )
    chi_2_e(df) = chi_2( loc )
    chi_3_e(df) = chi_3( loc )
  end do

  ! Calculate rotation vector Omega = (0, 2*cos(lat), 2*sin(lat))
  if ( geometry == geometry_spherical ) then
    call rotation_vector_sphere(ndf_chi, nqp_h, nqp_v, chi_1_e, chi_2_e, &
                                chi_3_e, ipanel, chi_basis, rotation_vector)
  else
    call rotation_vector_fplane(nqp_h, nqp_v, scaled_omega, f_lat, &
                                rotation_vector)
  end if

  call coordinate_jacobian(ndf_chi,        &
                           nqp_h,          &
                           nqp_v,          &
                           chi_1_e,        &
                           chi_2_e,        &
                           chi_3_e,        &
                           ipanel,         &
                           chi_basis,      &
                           chi_diff_basis, &
                           jac,            &
                           dj)

  call coordinate_jacobian_inverse(nqp_h, nqp_v, jac, dj, jac_inv)

  do df = 1, ndf_w2
    curl_u_e(df) = curl_u( map_w2(df) )
  end do
  do df = 1, ndf_w3
    r_q_e(df) = 0.0_r_def
    geopot_e(df) = geopot( map_w3(df) )
  end do

  ! Compute the RHS integrated over one cell
  do qp2 = 1, nqp_v
    do qp1 = 1, nqp_h
      curl_u_at_quad(:) = 0.0_r_def
      do df = 1, ndf_w2
        curl_u_at_quad(3) = curl_u_at_quad(3) + curl_u_e(df)*w2_basis(3,df,qp1,qp2)
      end do

      ! For test function eta in W1, the Coriolis term is J^-T eta dot rot_vec.
      ! For test function sigma in W3, move Jacobian over to rot_vec to get
      ! (sigma z_hat) dot J^-1 rot_vec.
      jac_rot(:) = matmul(jac_inv(:,:,qp1,qp2), rotation_vector(:,qp1,qp2))

      geopot_at_quad = 0.0_r_def
      do df = 1, ndf_w3
        geopot_at_quad = geopot_at_quad + geopot_e(df)*w3_basis(1,df,qp1,qp2)
      end do

      do df = 1, ndf_w3
        ! Strong curl and Coriolis, divided by geopotential
        ! only implemented with rehabilitation
        r_q_e(df) = r_q_e(df)                                         &
                     + ( curl_u_at_quad(3) + jac_rot(3)*dj(qp1,qp2) ) &
                        /geopot_at_quad*w3_basis(1,df,qp1,qp2)*wqp_h(qp1)*wqp_v(qp2)
      end do
    end do
  end do
  do df = 1, ndf_w3
    r_q( map_w3(df) ) =  r_q( map_w3(df) ) + r_q_e(df)
  end do

end subroutine initial_vorticity_v2_code

end module initial_vorticity_v2_kernel_mod
