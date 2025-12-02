!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Computes the absolute vorticity.
!>
!> Compute the projection of the absolute vorticity into the w1 space.
!> That is: curl(u) + 2*Omega -> curl(c).u + 2*c.Omega
!>
module vorticity_rhs_kernel_mod

  use argument_mod,      only : arg_type, func_type,       &
                                GH_FIELD, GH_REAL,         &
                                GH_READ, GH_INC,           &
                                ANY_SPACE_9,               &
                                ANY_DISCONTINUOUS_SPACE_3, &
                                GH_BASIS, GH_DIFF_BASIS,   &
                                CELL_COLUMN, GH_QUADRATURE_XYoZ
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W1, W2
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: vorticity_rhs_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                    &
         arg_type(GH_FIELD,   GH_REAL, GH_INC,  W1),                       &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, W2),                       &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ, ANY_SPACE_9),              &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_3) &
         /)
    type(func_type) :: meta_funcs(3) = (/                                  &
         func_type(W1,          GH_DIFF_BASIS),                            &
         func_type(W2,          GH_BASIS),                                 &
         func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                   &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: vorticity_rhs_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: vorticity_rhs_code

contains

!> @brief Compute the projection of curl(u) into the vorticity function space
!! @param[in] nlayers Number of layers
!! @param[in,out] rhs Right hand side to be computed
!! @param[in] u Velocity field
!! @param[in] chi_1 1st coordinate field in Wchi
!! @param[in] chi_2 2nd coordinate field in Wchi
!! @param[in] chi_3 3rd coordinate field in Wchi
!! @param[in] panel_id Field giving the ID for mesh panels
!! @param[in] ndf_xi Number of degrees of freedom per cell for W1
!! @param[in] undf_xi Unique number of degrees of freedom for W1
!! @param[in] map_xi Dofmap for the cell at the base of the column for W1
!! @param[in] diff_basis_xi Differential of the basis functions evaluated
!!                          at Gaussian quadrature point
!! @param[in] ndf_u Number of degrees of freedom per cell for the velocity field
!! @param[in] undf_u Unique number of degrees of freedom for the velocity field
!! @param[in] map_u Dofmap for the cell at the base of the column for the velocity field
!! @param[in] basis_u Basis functions evaluated at Gaussian quadrature points
!! @param[in] ndf_chi Number of degrees of freedom per cell for the function
!!                    space containing chi
!! @param[in] undf_chi Unique number of degrees of freedom for the chi arrays
!! @param[in] map_chi Dofmap for the cell at the base of the column for the
!!                    function space containing chi
!! @param[in] basis_chi Wchi basis functions evaluated at Gaussian quadrature points
!! @param[in] diff_basis_chi Derivatives of Wchi basis functions
!!                           evaluated at Gaussian quadrature points
!! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
!! @param[in] undf_pid Number of unique degrees of freedom for panel_id
!! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
!! @param[in] nqp_h Number of horizontal quadrature points
!! @param[in] nqp_v Number of vertical quadrature points
!! @param[in] wqp_h Weights of the horizontal quadrature points
!! @param[in] wqp_v Weights of the vertical quadrature points
subroutine vorticity_rhs_code(nlayers,                                &
                              rhs, u, chi_1, chi_2, chi_3, panel_id,  &
                              ndf_xi, undf_xi, map_xi, diff_basis_xi, &
                              ndf_u, undf_u, map_u, basis_u,          &
                              ndf_chi, undf_chi, map_chi,             &
                              basis_chi, diff_basis_chi,              &
                              ndf_pid, undf_pid, map_pid,             &
                              nqp_h, nqp_v, wqp_h, wqp_v              &
                              )

  use sci_coordinate_jacobian_mod, only: coordinate_jacobian

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_chi, ndf_u, ndf_xi, ndf_pid
  integer(kind=i_def), intent(in) :: undf_chi, undf_u, undf_xi, undf_pid
  integer(kind=i_def), intent(in) :: nqp_h, nqp_v
  integer(kind=i_def), dimension(ndf_xi),  intent(in) :: map_xi
  integer(kind=i_def), dimension(ndf_u),   intent(in) :: map_u
  integer(kind=i_def), dimension(ndf_pid), intent(in) :: map_pid
  integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi

  real(kind=r_def), dimension(3,ndf_u,  nqp_h,nqp_v), intent(in) :: basis_u
  real(kind=r_def), dimension(3,ndf_xi, nqp_h,nqp_v), intent(in) :: diff_basis_xi
  real(kind=r_def), dimension(1,ndf_chi,nqp_h,nqp_v), intent(in) :: basis_chi
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in) :: diff_basis_chi

  real(kind=r_def), dimension(undf_xi),               intent(inout) :: rhs
  real(kind=r_def), dimension(undf_u),                intent(in)    :: u
  real(kind=r_def), dimension(undf_chi),              intent(in)    :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(undf_pid),              intent(in)    :: panel_id

  real(kind=r_def), dimension(nqp_h), intent(in) :: wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in) :: wqp_v

  ! Internal variables
  integer(kind=i_def) :: df, k, loc, ipanel
  integer(kind=i_def) :: qp1, qp2

  real(kind=r_def), dimension(ndf_chi) :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(nqp_h,nqp_v)        :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v)    :: jac
  real(kind=r_def), dimension(ndf_u)  :: u_cell(ndf_u)
  real(kind=r_def), dimension(ndf_xi) :: rhs_cell(ndf_xi)
  real(kind=r_def) :: u_at_quad(3), &
                      dc(3)

  ipanel = int(panel_id(map_pid(1)), i_def)

  do k = 0, nlayers-1
  ! Extract element arrays of chi
    do df = 1, ndf_chi
      loc = map_chi(df) + k
      chi_1_e(df) = chi_1( loc )
      chi_2_e(df) = chi_2( loc )
      chi_3_e(df) = chi_3( loc )
    end do
    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi_1_e, chi_2_e, chi_3_e,  &
                             ipanel, basis_chi, diff_basis_chi, jac, dj)
    do df = 1, ndf_u
      u_cell(df) = u( map_u(df) + k )
    end do
  ! compute the RHS integrated over one cell
    rhs_cell(:) = 0.0_r_def
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        u_at_quad(:) = 0.0_r_def
        do df = 1, ndf_u
          u_at_quad(:) = u_at_quad(:) + u_cell(df)*basis_u(:,df,qp1,qp2)
        end do
        u_at_quad(:) = wqp_h(qp1)*wqp_v(qp2)*matmul(jac(:,:,qp1,qp2),u_at_quad(:))/dj(qp1,qp2)
        do df = 1, ndf_xi
          dc = matmul(jac(:,:,qp1,qp2), diff_basis_xi(:,df,qp1,qp2))
          rhs_cell(df) = rhs_cell(df) + dot_product(dc,u_at_quad(:))
        end do
      end do
    end do
    do df = 1, ndf_xi
      rhs( map_xi(df) + k ) =  rhs( map_xi(df) + k ) + rhs_cell(df)
    end do
  end do
end subroutine vorticity_rhs_code

end module vorticity_rhs_kernel_mod
