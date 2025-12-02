!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Computes LHS of Galerkin projection and solves equation in W3 space.
!>
module set_exner_kernel_mod

  use argument_mod,         only : arg_type, func_type,       &
                                   GH_FIELD, GH_SCALAR,       &
                                   GH_READ, GH_WRITE,         &
                                   GH_REAL, ANY_SPACE_9,      &
                                   ANY_DISCONTINUOUS_SPACE_3, &
                                   GH_BASIS, GH_DIFF_BASIS,   &
                                   CELL_COLUMN, GH_QUADRATURE_XYoZ
  use constants_mod,        only : r_def, i_def
  use fs_continuity_mod,    only : W3
  use idealised_config_mod, only : test
  use kernel_mod,           only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: set_exner_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                      &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, W3),                        &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ,  ANY_SPACE_9),               &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ)                              &
         /)
    type(func_type) :: meta_funcs(2) = (/                                    &
         func_type(W3,          GH_BASIS),                                   &
         func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                     &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: set_exner_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: set_exner_code

contains

!> @brief Computes LHS of Galerkin projection and solves equation in W3 space
!! @param[in] nlayers Number of layers
!! @param[in,out] exner Pressure field
!! @param[in] chi_1 1st coordinate field
!! @param[in] chi_2 2nd coordinate field
!! @param[in] chi_3 3rd coordinate field
!! @param[in] panel_id Field giving the ID for mesh panels
!! @param[in] time Time evaluated as a real value
!! @param[in] ndf_w3 Number of degrees of freedom per cell
!! @param[in] undf_w3 Total number of degrees of freedom
!! @param[in] map_w3 Dofmap for the cell at the base of the column
!! @param[in] w3_basis Basis functions evaluated at Gaussian quadrature points
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi
!! @param[in] undf_chi Number of unique degrees of freedom for chi
!! @param[in] map_chi Dofmap for the cell at the base of the column for chi
!! @param[in] chi_basis Basis functions for Wchi evaluated at
!!                      Gaussian quadrature points
!! @param[in] chi_diff_basis Differential of the Wchi basis functions
!!                           evaluated at Gaussian quadrature points
!! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
!! @param[in] undf_pid Number of unique degrees of freedom for panel_id
!! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
!! @param[in] nqp_h Number of horizontal quadrature points
!! @param[in] nqp_v Number of vertical quadrature points
!! @param[in] wqp_h Weights of horizontal quadrature points
!! @param[in] wqp_v Weights of vertical quadrature points
subroutine set_exner_code(nlayers,                                    &
                          exner,                                      &
                          chi_1, chi_2, chi_3, panel_id,              &
                          time,                                       &
                          ndf_w3, undf_w3, map_w3,                    &
                          w3_basis,                                   &
                          ndf_chi, undf_chi, map_chi,                 &
                          chi_basis, chi_diff_basis,                  &
                          ndf_pid, undf_pid, map_pid,                 &
                          nqp_h, nqp_v, wqp_h, wqp_v)

  use matrix_invert_mod,              only: matrix_invert_lu
  use sci_coordinate_jacobian_mod,    only: coordinate_jacobian
  use sci_chi_transform_mod,          only: chi2xyz
  use analytic_pressure_profiles_mod, only: analytic_pressure

  ! Needs to compute the integral of exner_df * P
  ! P_analytic over a single column

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w3, ndf_chi, ndf_pid
  integer(kind=i_def), intent(in) :: undf_w3, undf_chi, undf_pid
  integer(kind=i_def), intent(in) :: nqp_h, nqp_v

  integer(kind=i_def), dimension(ndf_w3),  intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_pid), intent(in) :: map_pid
  integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi

  real(kind=r_def), dimension(undf_w3),             intent(inout) :: exner
  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v),   intent(in) :: w3_basis
  real(kind=r_def),                                    intent(in) :: time
  real(kind=r_def), dimension(undf_chi),               intent(in) :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v),  intent(in) :: chi_diff_basis
  real(kind=r_def), dimension(1,ndf_chi,nqp_h,nqp_v),  intent(in) :: chi_basis
  real(kind=r_def), dimension(undf_pid),               intent(in) :: panel_id
  real(kind=r_def), dimension(nqp_h),                  intent(in) :: wqp_h
  real(kind=r_def), dimension(nqp_v),                  intent(in) :: wqp_v

  ! Internal variables
  integer(kind=i_def) :: df1, df2, k, ipanel
  integer(kind=i_def) :: qp1, qp2

  real(kind=r_def), dimension(ndf_w3)          :: exner_e, rhs_e
  real(kind=r_def), dimension(ndf_w3,ndf_w3)   :: mass_matrix_w3, inv_mass_matrix_w3
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac
  real(kind=r_def), dimension(ndf_chi)         :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def)                             :: exner_ref, integrand
  real(kind=r_def)                             :: xyz(3), coords(3)

  ipanel = int(panel_id(map_pid(1)), i_def)

  ! Compute the RHS & LHS integrated over one cell and solve
  do k = 0, nlayers-1
    do df1 = 1, ndf_chi
      chi_1_e(df1) = chi_1(map_chi(df1) + k)
      chi_2_e(df1) = chi_2(map_chi(df1) + k)
      chi_3_e(df1) = chi_3(map_chi(df1) + k)
    end do
    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v,             &
                             chi_1_e, chi_2_e, chi_3_e,         &
                             ipanel, chi_basis, chi_diff_basis, &
                             jac, dj)
    ! Compute RHS
    do df1 = 1, ndf_w3
      rhs_e(df1) = 0.0_r_def
      do qp2 = 1, nqp_v
        do qp1 = 1, nqp_h
          coords(:) = 0.0_r_def
          do df2 = 1, ndf_chi
            coords(1) = coords(1) + chi_1_e(df2)*chi_basis(1,df2,qp1,qp2)
            coords(2) = coords(2) + chi_2_e(df2)*chi_basis(1,df2,qp1,qp2)
            coords(3) = coords(3) + chi_3_e(df2)*chi_basis(1,df2,qp1,qp2)
          end do

          ! Get (X,Y,Z) coordinates
          call chi2xyz(coords(1), coords(2), coords(3), &
                       ipanel, xyz(1), xyz(2), xyz(3))

          exner_ref = analytic_pressure(xyz, test, time)

          integrand =  w3_basis(1,df1,qp1,qp2) * exner_ref * dj(qp1,qp2)
          rhs_e(df1) = rhs_e(df1) + wqp_h(qp1)*wqp_v(qp2)*integrand
        end do
      end do
    end do
    ! Compute LHS
    do df1 = 1, ndf_w3
       do df2 = 1, ndf_w3
          mass_matrix_w3(df1,df2) = 0.0_r_def
          do qp2 = 1, nqp_v
             do qp1 = 1, nqp_h
                integrand =  w3_basis(1,df1,qp1,qp2) * &
                             w3_basis(1,df2,qp1,qp2) * dj(qp1,qp2)
                 mass_matrix_w3(df1,df2) = mass_matrix_w3(df1,df2) &
                                         + wqp_h(qp1)*wqp_v(qp2)*integrand
             end do
          end do
       end do
    end do
    ! Solve
    call matrix_invert_lu(mass_matrix_w3,inv_mass_matrix_w3,ndf_w3)
    exner_e = matmul(inv_mass_matrix_w3,rhs_e)
    do df1 = 1,ndf_w3
      exner(map_w3(df1)+k) = exner_e(df1)
    end do
  end do

end subroutine set_exner_code

end module set_exner_kernel_mod
