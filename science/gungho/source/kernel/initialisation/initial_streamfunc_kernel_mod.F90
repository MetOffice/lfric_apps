!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Compute a stream function psi such that u = curl(psi)

module initial_streamfunc_kernel_mod

use argument_mod,              only : arg_type, func_type,       &
                                      GH_FIELD, GH_SCALAR,       &
                                      GH_INC, GH_READ,           &
                                      GH_REAL, ANY_SPACE_9,      &
                                      ANY_DISCONTINUOUS_SPACE_3, &
                                      GH_BASIS, GH_DIFF_BASIS,   &
                                      CELL_COLUMN, GH_QUADRATURE_XYoZ
use constants_mod,             only : r_def, i_def, PI
use fs_continuity_mod,         only : W1
use kernel_mod,                only : kernel_type
use log_mod,                   only : log_event, LOG_LEVEL_ERROR
use initial_wind_config_mod,   only : profile

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: initial_streamfunc_kernel_type
  private
  type(arg_type) :: meta_args(9) = (/                                     &
       arg_type(GH_FIELD,   GH_REAL, GH_INC,  W1),                        &
       arg_type(GH_FIELD*3, GH_REAL, GH_READ, ANY_SPACE_9),               &
       arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_3), &
       arg_type(GH_SCALAR,  GH_REAL, GH_READ),                            &
       arg_type(GH_SCALAR,  GH_REAL, GH_READ),                            &
       arg_type(GH_SCALAR,  GH_REAL, GH_READ),                            &
       arg_type(GH_SCALAR,  GH_REAL, GH_READ),                            &
       arg_type(GH_SCALAR,  GH_REAL, GH_READ),                            &
       arg_type(GH_SCALAR,  GH_REAL, GH_READ)                             &
       /)
  type(func_type) :: meta_funcs(2) = (/                                   &
       func_type(W1,          GH_BASIS),                                  &
       func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                    &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_QUADRATURE_XYoZ
contains
  procedure, nopass :: initial_streamfunc_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: initial_streamfunc_code
contains

!> @brief Computes the righthand side of the galerkin projection of a
!>        stream function given by an analytical expression
!> @details Computes rhs = int (c . psi dV) for a vector field psi whose
!!          vertical components contain the values of a stream function given by a
!!          chosen analytic expression
!! @param[in] nlayers Number of layers
!! @param[in,out] rhs Right hand side field to compute
!! @param[in] chi_1 1st coordinate field
!! @param[in] chi_2 2nd coordinate field
!! @param[in] chi_3 3rd coordinate field
!! @param[in] panel_id A field giving the ID for mesh panels
!! @param[in] time  The time to be passed to the analytic stream function
!! @param[in] domain_max_x  Domain maximum x-coordinate.
!! @param[in] sbr_angle_lat  SBR angle latitude
!! @param[in] sbr_angle_lon  SBR angle longitude
!! @param[in] u0  Initial u-wind
!! @param[in] v0  Initial v-wind
!! @param[in] ndf Number of degrees of freedom per cell
!! @param[in] undf Total number of degrees of freedom
!! @param[in] map Dofmap for the cell at the base of the column for W1
!! @param[in] basis Basis functions evaluated at gaussian quadrature points
!! @param[in] ndf_chi Number of degrees of freedom per cell for spherical chi
!! @param[in] undf_chi Number of unique degrees of freedom for spherical chi
!! @param[in] map_chi Dofmap for the cell at the base of the column for spherical chi
!! @param[in] chi_basis Basis functions for spherical Wchi evaluated at
!!                      gaussian quadrature points
!! @param[in] chi_diff_basis Differential of the spherical Wchi basis functions
!!                           evaluated at gaussian quadrature point
!! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
!! @param[in] undf_pid Number of unique degrees of freedom for panel_id
!! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
!! @param[in] nqp_h Number of quadrature points in the horizontal
!! @param[in] nqp_v Number of quadrature points in the vertical
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine initial_streamfunc_code(nlayers,                         &
                                   rhs,                             &
                                   chi_1, chi_2, chi_3,             &
                                   panel_id,                        &
                                   time, domain_max_x,              &
                                   sbr_angle_lat, sbr_angle_lon,    &
                                   u0, v0,                          &
                                   ndf, undf, map, basis,           &
                                   ndf_chi, undf_chi,               &
                                   map_chi, chi_basis,              &
                                   chi_diff_basis,                  &
                                   ndf_pid, undf_pid, map_pid,      &
                                   nqp_h, nqp_v, wqp_h, wqp_v       &
                                   )

  use analytic_streamfunction_profiles_mod, only: analytic_streamfunction
  use base_mesh_config_mod,                 only: geometry,           &
                                                  geometry_planar,    &
                                                  geometry_spherical
  use chi_transform_mod,                    only: chi2llr
  use coordinate_jacobian_mod,              only: coordinate_jacobian, &
                                                  coordinate_jacobian_inverse
  use coord_transform_mod,                  only: sphere2cart_vector

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, ndf, ndf_chi, ndf_pid
  integer(kind=i_def), intent(in) :: undf, undf_chi, undf_pid
  integer(kind=i_def), intent(in) :: nqp_h, nqp_v

  integer(kind=i_def), dimension(ndf),     intent(in) :: map
  integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_pid), intent(in) :: map_pid


  real(kind=r_def), intent(in), dimension(3,ndf,nqp_h,nqp_v)     :: basis
  real(kind=r_def), intent(in), dimension(3,ndf_chi,nqp_h,nqp_v) :: chi_diff_basis
  real(kind=r_def), intent(in), dimension(1,ndf_chi,nqp_h,nqp_v) :: chi_basis

  real(kind=r_def), dimension(undf),     intent(inout) :: rhs
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(undf_pid), intent(inout) :: panel_id

  real(kind=r_def), dimension(nqp_h), intent(in) ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in) ::  wqp_v
  real(kind=r_def),                   intent(in) ::  time
  real(kind=r_def),                   intent(in) ::  domain_max_x
  real(kind=r_def),                   intent(in) ::  sbr_angle_lat
  real(kind=r_def),                   intent(in) ::  sbr_angle_lon
  real(kind=r_def),                   intent(in) ::  u0
  real(kind=r_def),                   intent(in) ::  v0

  ! Internal variables
  integer(kind=i_def)                          :: df, k, qp1, qp2, ipanel
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jacobian, jac_inv
  real(kind=r_def), dimension(ndf_chi)         :: chi_1_cell, chi_2_cell, chi_3_cell
  real(kind=r_def), dimension(3)               :: psi_physical, psi_spherical, coords, llr
  real(kind=r_def)                             :: integrand
  real(kind=r_def), dimension(2)               :: option2
  real(kind=r_def), dimension(3)               :: option3

  ipanel = int(panel_id(map_pid(1)), i_def)

  option3 = (/ U0, sbr_angle_lat, sbr_angle_lon /)
  option2 = (/ U0, V0 /)

  do k = 0, nlayers-1

    do df = 1, ndf_chi
      chi_1_cell(df) = chi_1( map_chi(df) + k)
      chi_2_cell(df) = chi_2( map_chi(df) + k)
      chi_3_cell(df) = chi_3( map_chi(df) + k)
    end do


    call coordinate_jacobian(ndf_chi,        &
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

    call coordinate_jacobian_inverse(nqp_h, nqp_v, jacobian, dj, jac_inv)
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        ! Compute analytical vector streamfunction in physical space
        coords(:) = 0.0_r_def
        do df = 1, ndf_chi
          coords(1) = coords(1) + chi_1_cell(df)*chi_basis(1,df,qp1,qp2)
          coords(2) = coords(2) + chi_2_cell(df)*chi_basis(1,df,qp1,qp2)
          coords(3) = coords(3) + chi_3_cell(df)*chi_basis(1,df,qp1,qp2)
        end do

        if ( geometry == geometry_spherical ) then
          ! Need (lon,lat,r) coordinates
          call chi2llr(coords(1), coords(2), coords(3), &
                       ipanel, llr(1), llr(2), llr(3))

          psi_spherical = analytic_streamfunction( llr, profile, 3, option3, &
                                                   time, domain_max_x )
          psi_physical = sphere2cart_vector( psi_spherical,llr )
        else if ( geometry == geometry_planar ) then
          ! Coords must be (X,Y,Z)
          psi_physical = analytic_streamfunction( coords, profile, 2, option2, &
                                                  time, domain_max_x )
        else
          call log_event('initial_streamfunc_kernel is not implemented ' // &
                         'with your geometry and topology',                 &
                         LOG_LEVEL_ERROR)
        end if

        do df = 1, ndf
          integrand = dot_product(matmul(transpose(jac_inv(:,:,qp1,qp2)),&
                                         basis(:,df,qp1,qp2)),psi_physical)*dj(qp1,qp2)
          rhs(map(df) + k) = rhs(map(df) + k) &
                           + wqp_h(qp1)*wqp_v(qp2)*integrand
        end do
      end do
    end do
  end do

end subroutine initial_streamfunc_code

end module initial_streamfunc_kernel_mod
