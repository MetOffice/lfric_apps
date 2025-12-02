!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the initial (unperturbed) reference theta field.
!> @details The kernel computes initial theta reference field for theta in the
!!          space of horizontally discontinuous, vertically continuous
!!          polynomials. This is based on an analytical unperturbed reference
!!          vertical profile.
module initial_theta_ref_kernel_mod

  use argument_mod,      only: arg_type, func_type,       &
                               GH_FIELD, GH_REAL,         &
                               GH_INTEGER, GH_SCALAR,     &
                               GH_WRITE, GH_READ,         &
                               GH_BASIS,                  &
                               ANY_DISCONTINUOUS_SPACE_3, &
                               CELL_COLUMN, GH_EVALUATOR
  use constants_mod,     only: r_def, i_def
  use fs_continuity_mod, only: Wtheta, Wchi
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !----------------------------------------------------------------------------
  ! Public types
  !----------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: initial_theta_ref_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                      &
      arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, Wtheta),                    &
      arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  Wchi),                      &
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3), &
      arg_type(GH_SCALAR,  GH_INTEGER, GH_READ )                             &
      /)
    type(func_type) :: meta_funcs(1) = (/                                    &
      func_type(Wchi, GH_BASIS)                                              &
      /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: initial_theta_ref_code
  end type

  !----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !----------------------------------------------------------------------------
  public :: initial_theta_ref_code

contains

!> @brief Computes the initial theta field.
!> @param[in] nlayers     Number of layers
!> @param[in,out] theta   Potential temperature
!> @param[in] chi_1       First component of the chi coordinate field
!> @param[in] chi_2       Second component of the chi coordinate field
!> @param[in] chi_3       Third component of the chi coordinate field
!> @param[in] panel_id    A field giving the ID for mesh panels
!> @param[in] test        Analytical profile configuration
!> @param[in] ndf_wtheta  Number of degrees of freedom per cell for Wtheta
!> @param[in] undf_wtheta Number of total degrees of freedom for Wtheta
!> @param[in] map_wtheta  Dofmap for the cell at the base of the column for Wtheta
!> @param[in] ndf_chi     Number of degrees of freedom per cell for Wchi
!> @param[in] undf_chi    Number of total degrees of freedom for Wchi
!> @param[in] map_chi     Dofmap for the cell at the base of the column for Wchi
!> @param[in] chi_basis   Basis functions evaluated at Wchi points
!> @param[in] ndf_pid     Number of degrees of freedom per cell for panel_id
!> @param[in] undf_pid    Number of unique degrees of freedom for panel_id
!> @param[in] map_pid     Dofmap for the cell at the base of the column for panel_id
subroutine initial_theta_ref_code(nlayers,                               &
                                  theta,                                 &
                                  chi_1, chi_2, chi_3,                   &
                                  panel_id,                              &
                                  test,                                  &
                                  ndf_wtheta, undf_wtheta, map_wtheta,   &
                                  ndf_chi, undf_chi, map_chi, chi_basis, &
                                  ndf_pid, undf_pid, map_pid             )

  use analytic_ref_temperature_profiles_mod, only: analytic_ref_temperature
  use sci_chi_transform_mod,                 only: chi2xyz

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, ndf_wtheta, ndf_chi, ndf_pid
  integer(kind=i_def), intent(in) :: undf_wtheta, undf_chi, undf_pid
  integer(kind=i_def), dimension(ndf_wtheta), intent(in) :: map_wtheta
  integer(kind=i_def), dimension(ndf_chi),    intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_pid),    intent(in) :: map_pid

  real(kind=r_def),    dimension(undf_wtheta),          intent(inout) :: theta
  real(kind=r_def),    dimension(undf_chi),             intent(in)    :: chi_1, chi_2, chi_3
  real(kind=r_def),    dimension(undf_pid),             intent(in)    :: panel_id
  real(kind=r_def),    dimension(1,ndf_chi,ndf_wtheta), intent(in)    :: chi_basis

  integer(kind=i_def), intent(in) :: test

  ! Internal variables
  real(kind=r_def),   dimension(ndf_chi) :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def)                       :: coords(3), xyz(3)
  integer(kind=i_def)                    :: df, dfc, k, ipanel

  ipanel = int(panel_id(map_pid(1)), i_def)

  ! Compute the pointwise theta profile

  do k = 0, nlayers-1
    do dfc = 1, ndf_chi
      chi_1_e(dfc) = chi_1( map_chi(dfc) + k)
      chi_2_e(dfc) = chi_2( map_chi(dfc) + k)
      chi_3_e(dfc) = chi_3( map_chi(dfc) + k)
    end do

    do df = 1, ndf_wtheta
      coords(:) = 0.0_r_def
      do dfc = 1, ndf_chi
        coords(1) = coords(1) + chi_1_e(dfc)*chi_basis(1,dfc,df)
        coords(2) = coords(2) + chi_2_e(dfc)*chi_basis(1,dfc,df)
        coords(3) = coords(3) + chi_3_e(dfc)*chi_basis(1,dfc,df)
      end do

      call chi2xyz(coords(1), coords(2), coords(3), &
                   ipanel, xyz(1), xyz(2), xyz(3))
      theta(map_wtheta(df) + k) = analytic_ref_temperature(xyz, test)

    end do
  end do

end subroutine initial_theta_ref_code

end module initial_theta_ref_kernel_mod
