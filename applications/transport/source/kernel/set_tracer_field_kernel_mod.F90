!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Initialise a tracer field at either rho or theta spaces.
!>
module set_tracer_field_kernel_mod

  use argument_mod,         only : arg_type, func_type,       &
                                   GH_FIELD, GH_SCALAR,       &
                                   GH_READ, GH_WRITE,         &
                                   GH_REAL, ANY_SPACE_9,      &
                                   ANY_DISCONTINUOUS_SPACE_1, &
                                   ANY_DISCONTINUOUS_SPACE_3, &
                                   GH_BASIS, GH_DIFF_BASIS,   &
                                   GH_INTEGER, CELL_COLUMN,   &
                                   GH_QUADRATURE_XYoZ
  use fs_continuity_mod,    only : Wchi
  use constants_mod,        only : r_def, i_def
  use idealised_config_mod, only : test
  use kernel_mod,           only : kernel_type
  use log_mod,              only : log_event, LOG_LEVEL_ERROR

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: set_tracer_field_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                         &
         arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  Wchi),                      &
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3), &
         arg_type(GH_SCALAR,  GH_REAL,    GH_READ),                             &
         arg_type(GH_SCALAR,  GH_INTEGER, GH_READ)                              &
         /)
    type(func_type) :: meta_funcs(2) = (/                                       &
         func_type(ANY_DISCONTINUOUS_SPACE_1, GH_BASIS),                        &
         func_type(Wchi,                      GH_BASIS, GH_DIFF_BASIS)          &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: set_tracer_field_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: set_tracer_field_code

contains

!> @brief Initialise a tracer field at either rho or theta spaces.
!> @param[in] nlayers Number of layers
!> @param[in,out] tracer     Tracer field
!> @param[in] chi_1 1st      Coordinate field
!> @param[in] chi_2 2nd      Coordinate field
!> @param[in] chi_3 3rd      Coordinate field
!> @param[in] panel_id       Field giving the ID for mesh panels
!> @param[in] domain_max_x   Domain maximum extent in x direction.
!> @param[in] const_flag     Flag to set field to constant value (1) or use namelist (0)
!> @param[in] ndf_trc        Number of degrees of freedom per cell
!> @param[in] undf_trc       Total number of degrees of freedom
!> @param[in] map_trc        Dofmap for the cell at the base of the column
!> @param[in] trc_basis      Basis functions evaluated at Gaussian quadrature points
!> @param[in] ndf_chi        Number of degrees of freedom per cell for chi
!> @param[in] undf_chi       Number of unique degrees of freedom for chi
!> @param[in] map_chi        Dofmap for the cell at the base of the column for chi
!> @param[in] chi_basis      Basis functions for Wchi evaluated at
!!                           Gaussian quadrature points
!> @param[in] chi_diff_basis Differential of the Wchi basis functions
!!                           evaluated at Gaussian quadrature points
!> @param[in] ndf_pid        Number of degrees of freedom per cell for panel_id
!> @param[in] undf_pid       Number of unique degrees of freedom for panel_id
!> @param[in] map_pid        Dofmap for the panel_id function space
!> @param[in] nqp_h          Number of horizontal quadrature points
!> @param[in] nqp_v          Number of vertical quadrature points
!> @param[in] wqp_h          Weights of horizontal quadrature points
!> @param[in] wqp_v          Weights of vertical quadrature points
subroutine set_tracer_field_code(nlayers, tracer,                        &
                                 chi_1, chi_2, chi_3,                    &
                                 panel_id,                               &
                                 domain_max_x,                           &
                                 const_flag,                             &
                                 ndf_trc, undf_trc, map_trc, trc_basis,  &
                                 ndf_chi, undf_chi, map_chi,             &
                                 chi_basis, chi_diff_basis,              &
                                 ndf_pid, undf_pid, map_pid,             &
                                 nqp_h, nqp_v, wqp_h, wqp_v )

  use matrix_invert_mod,                  only : matrix_invert
  use coordinate_jacobian_mod,            only : coordinate_jacobian
  use chi_transform_mod,                  only : chi2xyz
  use analytic_tracer_field_profiles_mod, only : analytic_tracer_field

  ! needs to compute the integral of trc_df * P
  ! P_analytic over a single column

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, nqp_h, nqp_v
  integer(kind=i_def), intent(in) :: ndf_trc, ndf_chi, ndf_pid
  integer(kind=i_def), intent(in) :: undf_trc, undf_chi, undf_pid
  integer(kind=i_def), dimension(ndf_trc), intent(in) :: map_trc
  integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_pid), intent(in) :: map_pid

  real(kind=r_def), dimension(1,ndf_trc,nqp_h,nqp_v), intent(in)    :: trc_basis
  real(kind=r_def), dimension(undf_trc),              intent(inout) :: tracer
  real(kind=r_def),                                   intent(in)    :: domain_max_x
  integer(kind=i_def),                                intent(in)    :: const_flag

  real(kind=r_def), dimension(undf_chi), intent(in) :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(undf_pid), intent(in) :: panel_id
  real(kind=r_def), dimension(nqp_h),    intent(in) :: wqp_h
  real(kind=r_def), dimension(nqp_v),    intent(in) :: wqp_v

  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_diff_basis
  real(kind=r_def), dimension(1,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_basis

  ! Internal variables
  integer(kind=i_def) :: df1, df2, k, ipanel
  integer(kind=i_def) :: qp1, qp2

  real(kind=r_def), dimension(ndf_trc)         :: trc_e, rhs_e
  real(kind=r_def), dimension(ndf_trc,ndf_trc) :: mass_matrix_trc, inv_mass_matrix_trc
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac
  real(kind=r_def), dimension(ndf_chi)         :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def)                             :: tracer_ref, integrand
  real(kind=r_def)                             :: xyz(3), coords(3)

  ipanel = int(panel_id(map_pid(1)), i_def)

  ! compute the RHS & LHS integrated over one cell and solve
  do k = 0, nlayers-1
    do df1 = 1, ndf_chi
      chi_1_e(df1) = chi_1( map_chi(df1) + k )
      chi_2_e(df1) = chi_2( map_chi(df1) + k )
      chi_3_e(df1) = chi_3( map_chi(df1) + k )
    end do
    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v,             &
                             chi_1_e, chi_2_e, chi_3_e,         &
                             ipanel, chi_basis, chi_diff_basis, &
                             jac, dj)

! Compute RHS
    do df1 = 1, ndf_trc
      rhs_e(df1) = 0.0_r_def
      do qp2 = 1, nqp_v
        do qp1 = 1, nqp_h
          coords(:) = 0.0_r_def
          do df2 = 1, ndf_chi
            coords(1) = coords(1) + chi_1_e(df2)*chi_basis(1,df2,qp1,qp2)
            coords(2) = coords(2) + chi_2_e(df2)*chi_basis(1,df2,qp1,qp2)
            coords(3) = coords(3) + chi_3_e(df2)*chi_basis(1,df2,qp1,qp2)
          end do

          ! Need (X,Y,Z) coordinate
          call chi2xyz(coords(1), coords(2), coords(3), &
                       ipanel, xyz(1), xyz(2), xyz(3))

          if (const_flag == 1_i_def) then
            ! Set tracer field to be 1 everywhere
            tracer_ref = 1.0_r_def
          else
            ! Set tracer field based on namelist
            tracer_ref = analytic_tracer_field(xyz, test, domain_max_x)
          end if

          integrand =  trc_basis(1,df1,qp1,qp2) * tracer_ref * dj(qp1,qp2)
          rhs_e(df1) = rhs_e(df1) + wqp_h(qp1)*wqp_v(qp2)*integrand
        end do
      end do
    end do
! Compute LHS
    do df1 = 1, ndf_trc
       do df2 = 1, ndf_trc
          mass_matrix_trc(df1,df2) = 0.0_r_def
          do qp2 = 1, nqp_v
             do qp1 = 1, nqp_h
                integrand =  trc_basis(1,df1,qp1,qp2) * &
                             trc_basis(1,df2,qp1,qp2) * dj(qp1,qp2)
                 mass_matrix_trc(df1,df2) = mass_matrix_trc(df1,df2) &
                                          + wqp_h(qp1)*wqp_v(qp2)*integrand
             end do
          end do
       end do
    end do
! Solve
    call matrix_invert(mass_matrix_trc,inv_mass_matrix_trc,ndf_trc)
    trc_e = matmul(inv_mass_matrix_trc,rhs_e)
    do df1 = 1,ndf_trc
      tracer(map_trc(df1)+k) = trc_e(df1)
    end do
  end do

end subroutine set_tracer_field_code

end module set_tracer_field_kernel_mod
