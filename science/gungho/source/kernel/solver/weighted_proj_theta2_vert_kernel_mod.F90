!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Compute the projection operator from the vertical part of the
!>        velocity space to the potential temperature space weighted by the
!>        potential temperature gradient.
!>
!> @details Compute the projection operator \f[<\gamma,\nabla(\theta*)v>\f]
!!          where v is in W2 and gamma is in the potential temperature space
!> @todo Create unit test for this kernel, see #2935
module weighted_proj_theta2_vert_kernel_mod

  use argument_mod,      only : arg_type, func_type,     &
                                GH_OPERATOR, GH_FIELD,   &
                                GH_SCALAR, GH_REAL,      &
                                GH_INTEGER,              &
                                GH_READ, GH_WRITE,       &
                                GH_BASIS, GH_DIFF_BASIS, &
                                CELL_COLUMN, GH_QUADRATURE_XYoZ
  use constants_mod,     only : r_def, i_def, r_solver
  use fs_continuity_mod, only : W2, Wtheta, W2H
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: weighted_proj_theta2_vert_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                        &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, Wtheta, W2), &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  Wtheta),     &
         arg_type(GH_SCALAR,   GH_REAL, GH_READ),              &
         arg_type(GH_SCALAR,   GH_INTEGER, GH_READ),           &
         arg_type(GH_SCALAR,   GH_INTEGER, GH_READ)            &
         /)
    type(func_type) :: meta_funcs(2) = (/                      &
         func_type(Wtheta, GH_BASIS, GH_DIFF_BASIS),           &
         func_type(W2,     GH_BASIS)                           &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: weighted_proj_theta2_vert_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: weighted_proj_theta2_vert_code

contains

!> @brief Compute the weighted projection operator that maps from W2 to Wtheta
!! @param[in] cell Current cell index
!! @param[in] nlayers Number of layers
!! @param[in] ncell_3d Total number of cells in the 3d mesh
!! @param[in,out] projection Locally assembled projection operator
!! @param[in] theta Potential temperature array
!! @param[in] scalar Real to scale matrix by
!! @param[in] element_order_h Horizontal element order of the function space
!! @param[in] element_order_v Vertical element order of the function space
!! @param[in] ndf_wtheta Number of degrees of freedom per cell for Wtheta
!! @param[in] undf_wtheta Number of unique degrees of freedom for Wtheta
!! @param[in] map_wtheta Dofmap for the cell at the base of the column for Wtheta
!! @param[in] wtheta_basis Basis functions evaluated at Gaussian quadrature points
!! @param[in] wtheta_diff_basis Differential basis functions evaluated at
!!                              Gaussian quadrature points
!! @param[in] ndf_w2 Number of degrees of freedom per cell for W2
!! @param[in] w2_basis Basis functions evaluated at Gaussian quadrature points
!! @param[in] nqp_h Number of horizontal quadrature points
!! @param[in] nqp_v Number of vertical quadrature points
!! @param[in] wqp_h Weights of the horizontal quadrature points
!! @param[in] wqp_v Weights of the vertical quadrature points
subroutine weighted_proj_theta2_vert_code(cell, nlayers, ncell_3d,              &
                                          projection,                           &
                                          theta,                                &
                                          scalar,                               &
                                          element_order_h,                      &
                                          element_order_v,                      &
                                          ndf_wtheta, undf_wtheta, map_wtheta,  &
                                          wtheta_basis, wtheta_diff_basis,      &
                                          ndf_w2, w2_basis,                     &
                                          nqp_h, nqp_v, wqp_h, wqp_v)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: cell, nlayers, ncell_3d, nqp_h, nqp_v
  integer(kind=i_def), intent(in) :: ndf_wtheta, ndf_w2, undf_wtheta
  integer(kind=i_def), intent(in) :: element_order_h, element_order_v

  integer(kind=i_def), dimension(ndf_wtheta), intent(in) :: map_wtheta

  real(kind=r_def), dimension(1,ndf_wtheta,nqp_h,nqp_v), intent(in) :: wtheta_basis
  real(kind=r_def), dimension(3,ndf_wtheta,nqp_h,nqp_v), intent(in) :: wtheta_diff_basis
  real(kind=r_def), dimension(3,ndf_w2,nqp_h,nqp_v),     intent(in) :: w2_basis

  real(kind=r_solver), dimension(ncell_3d,ndf_wtheta,ndf_w2), intent(inout) :: projection
  real(kind=r_solver), dimension(undf_wtheta),                intent(in)    :: theta
  real(kind=r_solver),                                        intent(in)    :: scalar

  real(kind=r_def), dimension(nqp_h), intent(in) ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in) ::  wqp_v

  ! Internal variables
  integer(kind=i_def) :: df, k, ik, dft, df2
  integer(kind=i_def) :: qp1, qp2
  integer(kind=i_def) :: k_h, k_v
  integer(kind=i_def) :: ndf_w2h, ndf_w2h_vol, ndf_w2v_vol, ndf_w2_vol

  real(kind=r_solver), dimension(ndf_wtheta) :: theta_e
  real(kind=r_solver) :: grad_theta_at_quad(3)
  real(kind=r_solver) :: integrand, i1(3), i2, wt

  real(kind=r_solver), dimension(1,ndf_wtheta,nqp_h,nqp_v) :: rsol_wtheta_basis
  real(kind=r_solver), dimension(3,ndf_wtheta,nqp_h,nqp_v) :: rsol_wtheta_diff_basis
  real(kind=r_solver), dimension(3,ndf_w2,nqp_h,nqp_v)     :: rsol_w2_basis

  rsol_wtheta_basis      = real(wtheta_basis, r_solver)
  rsol_wtheta_diff_basis = real(wtheta_diff_basis, r_solver)
  rsol_w2_basis          = real(w2_basis, r_solver)

  k_h = element_order_h
  k_v = element_order_v

  ! Calculate number of volume dofs of W2, W2H and W2V
  ndf_w2h     = 2*(k_h + 1)*(k_h + 2)*(k_v + 1)
  ndf_w2h_vol = 2*k_h*(k_h + 1)*(k_v + 1)
  ndf_w2v_vol = (k_h + 1)*(k_h + 1)*k_v
  ndf_w2_vol  = ndf_w2h_vol + ndf_w2v_vol

  do k = 0, nlayers-1
    ik = k + 1 + (cell-1)*nlayers
    do df = 1, ndf_wtheta
      theta_e(df)  = theta( map_wtheta(df) + k )
    end do
    projection(ik,:,:) = 0.0_r_solver
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        grad_theta_at_quad = 0.0_r_solver
        do df = 1, ndf_wtheta
          grad_theta_at_quad(:) = grad_theta_at_quad(:) &
                                + theta_e(df)*rsol_wtheta_diff_basis(:,df,qp1,qp2)
        end do
        ! Avoid issues when theta gradient is zero or negative:
        ! The helmholtz equation is (vastly simplified)
        ! (1 + div[g/theta * dtheta/dz * grad])Pi = rhs
        ! and so if dtheta/dz =< 0 then the characteristic of this
        ! equation changes and this can cause issues for iterative solvers
        ! so the pragmatic solution is to enfoce dtheta/dz > 0 here
        grad_theta_at_quad(3) = max(1.0_r_solver, grad_theta_at_quad(3))
        wt = real(wqp_h(qp1)*wqp_v(qp2), r_solver)
        i1 = scalar*grad_theta_at_quad*wt
        do df2 = ndf_w2h_vol+1,ndf_w2
          ! W2 dofs are ordered:
          !   a) Horizontal volume dofs
          !   b) Vertical volume dofs
          !   c) Horizontal face dofs
          !   d) Vertical face dofs
          ! Only use the vertical dofs of W2, which satisfy one of two
          ! conditions:
          !   b) df2 is a vertical volume dof of W2 (higher order spaces
          !      only), so
          !        ndf_w2h_vol < df2 <= ndf_w2_vol
          !   d) df2 is a vertical face dof of W2 (lives on the top or
          !      bottom face of a cell), so
          !        ndf_w2h + ndf_w2v_vol < df2
          if ( df2 <= ndf_w2_vol .or. ndf_w2h + ndf_w2v_vol < df2 ) then
            i2 = dot_product(i1, rsol_w2_basis(:,df2,qp1,qp2))
            do dft = 1,ndf_wtheta
              integrand = rsol_wtheta_basis(1,dft,qp1,qp2)*i2
              projection(ik,dft,df2) = projection(ik,dft,df2) + integrand
            end do
          end if
        end do
      end do
    end do
  end do

end subroutine weighted_proj_theta2_vert_code

end module weighted_proj_theta2_vert_kernel_mod
