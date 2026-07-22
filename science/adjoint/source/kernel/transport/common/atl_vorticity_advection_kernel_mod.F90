!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Adjoint to the vorticity advection component
!>        of the rhs of the momentum equation.
module atl_vorticity_advection_kernel_mod
use kernel_mod,                only: kernel_type
use argument_mod,              only: arg_type, func_type,                 &
                                     GH_FIELD, GH_REAL,                   &
                                     GH_READ, GH_INC,                     &
                                     ANY_SPACE_9,                         &
                                     ANY_DISCONTINUOUS_SPACE_3,           &
                                     GH_BASIS, GH_DIFF_BASIS,             &
                                     CELL_COLUMN, GH_QUADRATURE_XYoZ
use constants_mod,             only: r_def, i_def
use fs_continuity_mod,         only: W1, W2, Wchi

use base_mesh_config_mod,      only: geometry, topology
use finite_element_config_mod, only: coord_system
use planet_config_mod,         only: scaled_radius

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: atl_vorticity_advection_kernel_type
  private
  type(arg_type) :: meta_args(7) = (/                                    &
       arg_type(GH_FIELD,   GH_REAL, GH_READ, W2),                       &
       arg_type(GH_FIELD,   GH_REAL, GH_INC,  W2),                       &
       arg_type(GH_FIELD,   GH_REAL, GH_INC,  W1),                       &
       arg_type(GH_FIELD,   GH_REAL, GH_READ, W2),                       &
       arg_type(GH_FIELD,   GH_REAL, GH_READ, W1),                       &
       arg_type(GH_FIELD*3, GH_REAL, GH_READ, WChi),                     &
       arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_3) &
       /)
  type(func_type) :: meta_funcs(3) = (/                                &
       func_type(W2, GH_BASIS),                                        &
       func_type(W1, GH_BASIS),                                        &
       func_type(WChi, GH_BASIS, GH_DIFF_BASIS)                        &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_QUADRATURE_XYoZ
contains
  procedure, nopass :: atl_vorticity_advection_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public atl_vorticity_advection_code
contains

!> @brief Compute the adjoint for vorticity advection.
!! @param[in] nlayers      Number of layers
!! @param[in] r_u          ACTIVE Change in right hand side of the mom equation
!! @param[inout] wind      ACTIVE Change in advecting wind field
!! @param[inout] vorticity ACTIVE Change in vorticity field = curl(u)
!! @param[in] ls_wind      Lin state advecting wind field
!! @param[in] ls_vorticity Lin state vorticity field = curl(ls_u)
!! @param[in] chi_1        1st (spherical) coordinate field in Wchi
!! @param[in] chi_2        2nd (spherical) coordinate field in Wchi
!! @param[in] chi_3        3rd (spherical) coordinate field in Wchi
!! @param[in] panel_id     Field giving the ID for mesh panels.
!! @param[in] ndf_w2       Number of degrees of freedom per cell for w2
!! @param[in] undf_w2      Number of unique degrees of freedom  for w2
!! @param[in] map_w2       Dofmap for the cell at the base of the column for w2
!! @param[in] w2_basis     Basis functions evaluated at quadrature points
!! @param[in] ndf_w1       Number of degrees of freedom per cell for w1
!! @param[in] undf_w1      Number of unique degrees of freedom  for w1
!! @param[in] map_w1       Dofmap for the cell at the base of the column for w1
!! @param[in] w1_basis     Basis functions evaluated at gaussian quadrature points
!! @param[in] ndf_chi      Number of degrees of freedom per cell for chi
!! @param[in] undf_chi     Number of unique degrees of freedom  for chi
!! @param[in] map_chi      Dofmap for the cell at the base of the column for chi
!! @param[in] chi_basis    Wchi basis functions evaluated at gaussian quadrature points.
!! @param[in] chi_diff_basis Derivatives of Wchi basis functions
!!                           evaluated at gaussian quadrature points
!! @param[in] ndf_pid      Number of degrees of freedom per cell for panel_id
!! @param[in] undf_pid     Number of unique degrees of freedom for panel_id
!! @param[in] map_pid      Dofmap for the cell at the base of the column for panel_id
!! @param[in] nqp_h        Number of quadrature points in the horizontal
!! @param[in] nqp_v        Number of quadrature points in the vertical
!! @param[in] wqp_h        Horizontal quadrature weights
!! @param[in] wqp_v        Vertical quadrature weights
subroutine atl_vorticity_advection_code(nlayers,         &
                                        r_u,             &
                                        wind,            &
                                        vorticity,       &
                                        ls_wind,         &
                                        ls_vorticity,    &
                                        chi_1,           &
                                        chi_2,           &
                                        chi_3,           &
                                        panel_id,        &
                                        ndf_w2,          &
                                        undf_w2,         &
                                        map_w2,          &
                                        w2_basis,        &
                                        ndf_w1,          &
                                        undf_w1,         &
                                        map_w1,          &
                                        w1_basis,        &
                                        ndf_chi,         &
                                        undf_chi,        &
                                        map_chi,         &
                                        chi_basis,       &
                                        chi_diff_basis,  &
                                        ndf_pid,         &
                                        undf_pid,        &
                                        map_pid,         &
                                        nqp_h,           &
                                        nqp_v,           &
                                        wqp_h,           &
                                        wqp_v            &
                                        )

  use sci_coordinate_jacobian_mod, only: pointwise_coordinate_jacobian, &
                                         pointwise_coordinate_jacobian_inverse
  implicit none

  ! Arguments
  integer, intent(in) :: nlayers,nqp_h, nqp_v
  integer, intent(in) :: ndf_chi, ndf_w1, ndf_w2, ndf_pid
  integer, intent(in) :: undf_chi, undf_w1, undf_w2, undf_pid
  integer, dimension(ndf_chi), intent(in) :: map_chi
  integer, dimension(ndf_w1),  intent(in) :: map_w1
  integer, dimension(ndf_w2),  intent(in) :: map_w2
  integer, dimension(ndf_pid), intent(in) :: map_pid

  real(kind=r_def), dimension(3,ndf_w2,nqp_h,nqp_v),  intent(in) :: w2_basis
  real(kind=r_def), dimension(3,ndf_w1,nqp_h,nqp_v),  intent(in) :: w1_basis
  real(kind=r_def), dimension(1,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_basis
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_diff_basis

  real(kind=r_def), dimension(undf_w2),  intent(in)    :: r_u
  real(kind=r_def), dimension(undf_w2),  intent(inout) :: wind
  real(kind=r_def), dimension(undf_w1),  intent(inout) :: vorticity
  real(kind=r_def), dimension(undf_w2),  intent(in)    :: ls_wind
  real(kind=r_def), dimension(undf_w1),  intent(in)    :: ls_vorticity
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
  real(kind=r_def), dimension(3)        :: vorticity_at_quad, u_at_quad,       &
                                           ls_vorticity_at_quad, u_ls_at_quad, &
                                           j_vorticity, j_ls_vorticity,        &
                                           vorticity_term

  ipanel = int(panel_id(map_pid(1)), i_def)

  do k = 0, nlayers-1
  ! Extract element arrays of chi
    do df = 1, ndf_chi
      loc = map_chi(df) + k
      chi_1_e(df) = chi_1( loc )
      chi_2_e(df) = chi_2( loc )
      chi_3_e(df) = chi_3( loc )
    end do

    ! The nonlinear term is:
    ! v.( [(J^-1 * J^-T) * vorticity] cross u )
    ! The corresponding linear term is:
    ! v.( ( [(J^-1 * J^-T) * ls_vorticity] cross u )
    ! + ( [(J^-1 * J^-T) * vorticity] cross u_ls ) )

    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h

        ! Constants
        call pointwise_coordinate_jacobian(coord_system, geometry,             &
                                           topology, scaled_radius,            &
                                           ndf_chi, chi_1_e, chi_2_e, chi_3_e, &
                                           ipanel, chi_basis(:,:,qp1,qp2),     &
                                           chi_diff_basis(:,:,qp1,qp2), jac, dj)
        jac_inv =  pointwise_coordinate_jacobian_inverse(jac, dj)
        jac = matmul(jac_inv,transpose(jac_inv))

        ! Linearisation state
        ls_vorticity_at_quad(:) = 0.0_r_def
        do df = 1, ndf_w1
          ls_vorticity_at_quad(:) = ls_vorticity_at_quad(:)        &
                                  + ls_vorticity( map_w1(df) + k ) &
                                  * w1_basis(:,df,qp1,qp2)
        end do
        u_ls_at_quad(:) = 0.0_r_def
        do df = 1, ndf_w2
          u_ls_at_quad(:) = u_ls_at_quad(:)           &
                          + ls_wind( map_w2(df) + k ) &
                          * w2_basis(:,df,qp1,qp2)
        end do

        j_ls_vorticity(:) = wqp_h(qp1) * wqp_v(qp2) * matmul( jac, ls_vorticity_at_quad )

        ! Perturbation
        vorticity_term(:) = 0.0_r_def
        do df = ndf_w2, 1, -1
           vorticity_term(:) = vorticity_term(:) - w2_basis(:,df,qp1,qp2) * r_u( map_w2(df) + k )
        end do

        j_vorticity(:) = 0.0_r_def
        u_at_quad(:) = 0.0_r_def
        j_vorticity(1) = j_vorticity(1) + vorticity_term(3) * u_ls_at_quad(2)
        j_vorticity(2) = j_vorticity(2) - vorticity_term(3) * u_ls_at_quad(1)
        u_at_quad(2) = u_at_quad(2) + vorticity_term(3) * j_ls_vorticity(1)
        u_at_quad(1) = u_at_quad(1) - vorticity_term(3) * j_ls_vorticity(2)

        j_vorticity(3) = j_vorticity(3) + vorticity_term(2) * u_ls_at_quad(1)
        j_vorticity(1) = j_vorticity(1) - vorticity_term(2) * u_ls_at_quad(3)
        u_at_quad(1) = u_at_quad(1) + vorticity_term(2) * j_ls_vorticity(3)
        u_at_quad(3) = u_at_quad(3) - vorticity_term(2) * j_ls_vorticity(1)

        j_vorticity(2) = j_vorticity(2) + vorticity_term(1) * u_ls_at_quad(3)
        j_vorticity(3) = j_vorticity(3) - vorticity_term(1) * u_ls_at_quad(2)
        u_at_quad(3) = u_at_quad(3) + vorticity_term(1) * j_ls_vorticity(2)
        u_at_quad(2) = u_at_quad(2) - vorticity_term(1) * j_ls_vorticity(3)

        vorticity_at_quad = wqp_h(qp1) * wqp_v(qp2) * matmul( transpose(jac), j_vorticity )

        do df = ndf_w2, 1, -1
          wind( map_w2(df) + k ) = wind( map_w2(df) + k ) &
                                 + dot_product( u_at_quad, w2_basis(:,df,qp1,qp2) )
        end do

        do df = ndf_w1, 1, -1
          vorticity( map_w1(df) + k ) = vorticity( map_w1(df) + k ) &
                                      + dot_product( vorticity_at_quad, w1_basis(:,df,qp1,qp2) )
        end do
      end do
    end do
  end do

end subroutine atl_vorticity_advection_code

end module atl_vorticity_advection_kernel_mod
