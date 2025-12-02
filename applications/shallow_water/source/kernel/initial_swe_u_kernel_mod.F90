!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the rhs for the initialisation of the wind field.
!!
!> @details The kernel computes the rhs of the equation u = u0 where u0 is the
!!          analytically defined wind field. The analytic wind field is projected onto
!!          using Galerkin projection.
!!
module initial_swe_u_kernel_mod

  use argument_mod,            only : arg_type, func_type,             &
                                      GH_FIELD, GH_INC, GH_READ,       &
                                      ANY_SPACE_9, GH_REAL, GH_SCALAR, &
                                      GH_BASIS, GH_DIFF_BASIS,         &
                                      CELL_COLUMN, GH_QUADRATURE_XYoZ, &
                                      ANY_DISCONTINUOUS_SPACE_3
  use constants_mod,           only : r_def, PI, i_def
  use fs_continuity_mod,       only : W2
  use initial_wind_config_mod, only : profile_sin_uv,                        &
                                      profile, sbr_angle_lat, sbr_angle_lon, &
                                      u0, v0, shear, wavelength

  use kernel_mod,              only : kernel_type
  use shallow_water_settings_config_mod, &
                               only : swe_test

  implicit none

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: initial_swe_u_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                    &
        arg_type(GH_FIELD,   GH_REAL, GH_INC,  W2),                        &
        ARG_TYPE(GH_FIELD*3, GH_REAL, GH_READ, ANY_SPACE_9),               &
        arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_3), &
        arg_type(GH_SCALAR,  GH_REAL, GH_READ)                             &
        /)
    type(func_type) :: meta_funcs(2) = (/               &
        func_type(W2, GH_BASIS),                        &
        func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS) &
        /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, public, nopass :: initial_swe_u_code
  end type

!-----------------------------------------------------------------------------
! Contained functions/subroutines
!-----------------------------------------------------------------------------
public initial_swe_u_code

contains

!> @brief Compute the right hand side to initialise the wind field.
!> @param[in]     nlayers        Number of layers
!> @param[in,out] rhs            Right hand side field to compute
!> @param[in]     chi_1          1st component of the coordinate field
!> @param[in]     chi_2          2nd component of the coordinate field
!> @param[in]     chi_3          3rd component of the coordinate field
!> @param[in]     panel_id       Field containing the ID of the mesh panel
!> @param[in]     domain_x       domain_size in x direction
!> @param[in]     ndf            Number of degrees of freedom per cell
!> @param[in]     undf           Total number of degrees of freedom
!> @param[in]     map            Dofmap for the cell at the base of the column
!> @param[in]     basis          Basis functions evaluated at gaussian quadrature points
!> @param[in]     ndf_chi        Number of dofs per cell for the coordinate field
!> @param[in]     undf_chi       Total number of degrees of freedom
!> @param[in]     map_chi        Dofmap for the coordinate field
!> @param[in]     chi_basis      Basis functions evaluated at gaussian quadrature points
!> @param[in]     chi_diff_basis Basis functions evaluated at gaussian quadrature points
!> @param[in]     ndf_pid        The number of DoFs per cell for the panel ID
!> @param[in]     undf_pid       The number of DoFs for this partition for the panel ID
!> @param[in]     map_pid        DoF-map for the panel ID
!> @param[in]     nqp_h          Number of quadrature points in the horizontal
!> @param[in]     nqp_v          Number of quadrature points in the vertical
!> @param[in]     wqp_h          Horizontal quadrature weights
!> @param[in]     wqp_v          Vertical quadrature weights
subroutine initial_swe_u_code( nlayers, rhs,                       &
                               chi_1, chi_2, chi_3, panel_id,      &
                               domain_x,                           &
                               ndf, undf, map, basis,              &
                               ndf_chi, undf_chi,                  &
                               map_chi, chi_basis, chi_diff_basis, &
                               ndf_pid, undf_pid, map_pid,         &
                               nqp_h, nqp_v, wqp_h, wqp_v )

  use analytic_swe_wind_profiles_mod, only : analytic_swe_wind
  use base_mesh_config_mod,           only : geometry, &
                                             geometry_spherical
  use sci_coordinate_jacobian_mod,    only : coordinate_jacobian
  use coord_transform_mod,            only : sphere2cart_vector
  use sci_chi_transform_mod,          only : chi2llr, chi2xyz

  implicit none

  ! Arguments
  integer, intent(in) :: nlayers, ndf, ndf_chi, ndf_pid
  integer, intent(in) :: undf, undf_chi, undf_pid
  integer, intent(in) :: nqp_h, nqp_v

  integer, dimension(ndf),     intent(in) :: map
  integer, dimension(ndf_chi), intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_pid),     intent(in) :: map_pid

  real(kind=r_def), intent(in), dimension(3,ndf,    nqp_h,nqp_v) :: basis
  real(kind=r_def), intent(in), dimension(3,ndf_chi,nqp_h,nqp_v) :: chi_diff_basis
  real(kind=r_def), intent(in), dimension(1,ndf_chi,nqp_h,nqp_v) :: chi_basis

  real(kind=r_def), dimension(undf),     intent(inout) :: rhs
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(undf_pid), intent(in)    :: panel_id
  real(kind=r_def),                      intent(in)    :: domain_x

  real(kind=r_def), dimension(nqp_h), intent(in)      ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)      ::  wqp_v

  ! Internal variables
  integer                                      :: df, qp1, qp2
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jacobian
  real(kind=r_def), dimension(ndf_chi)         :: chi_1_cell, chi_2_cell, chi_3_cell
  real(kind=r_def), dimension(3)               :: u_physical, u_spherical
  real(kind=r_def), dimension(3)               :: coord, xyz, llr
  real(kind=r_def)                             :: integrand

  integer(kind=i_def) :: ipanel

  ipanel = int(panel_id(map_pid(1)), i_def)

  do df = 1, ndf_chi
    chi_1_cell(df) = chi_1( map_chi(df) )
    chi_2_cell(df) = chi_2( map_chi(df) )
    chi_3_cell(df) = chi_3( map_chi(df) )
  end do

  call coordinate_jacobian(ndf_chi,         &
                            nqp_h,          &
                            nqp_v,          &
                            chi_1_cell,     &
                            chi_2_cell,     &
                            chi_3_cell,     &
                            ipanel,         &
                            chi_basis,      &
                            chi_diff_basis, &
                            jacobian,       &
                            dj)
  do qp2 = 1, nqp_v
    do qp1 = 1, nqp_h
      ! Compute analytical vector wind in physical space
      coord(:) = 0.0_r_def
      do df = 1, ndf_chi
        coord(1) = coord(1) + chi_1_cell(df)*chi_basis(1,df,qp1,qp2)
        coord(2) = coord(2) + chi_2_cell(df)*chi_basis(1,df,qp1,qp2)
        coord(3) = coord(3) + chi_3_cell(df)*chi_basis(1,df,qp1,qp2)
      end do
      if ( geometry == geometry_spherical ) then
        call chi2llr(coord(1), coord(2), coord(3), ipanel, llr(1), llr(2), llr(3))
        u_spherical = analytic_swe_wind(llr, swe_test, domain_x)
        u_physical = sphere2cart_vector(u_spherical,llr)
      else
        call chi2xyz(coord(1), coord(2), coord(3), ipanel, xyz(1), xyz(2), xyz(3))
        u_physical = analytic_swe_wind(xyz, swe_test, domain_x)
      end if
      do df = 1, ndf
        integrand = dot_product(matmul(jacobian(:,:,qp1,qp2),&
                                       basis(:,df,qp1,qp2)),u_physical)
        rhs(map(df)) = rhs(map(df)) &
                         + wqp_h(qp1)*wqp_v(qp2)*integrand
      end do
    end do
  end do

end subroutine initial_swe_u_code

end module initial_swe_u_kernel_mod
