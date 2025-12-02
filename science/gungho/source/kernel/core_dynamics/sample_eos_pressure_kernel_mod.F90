!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! For further details please refer to the file LICENCE which you should have
! received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Computes the pressure field by sampling of the density field at W3
!>        DoFs
!>
module sample_eos_pressure_kernel_mod

  use argument_mod,      only : arg_type, func_type,       &
                                GH_FIELD,                  &
                                GH_READ, GH_WRITE,         &
                                GH_REAL, GH_BASIS,         &
                                CELL_COLUMN, GH_EVALUATOR
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W3, Wtheta
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: sample_eos_pressure_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                       &
         arg_type(GH_FIELD,    GH_REAL, GH_WRITE, W3),                        &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  W3),                        &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  Wtheta),                    &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  Wtheta)                     &
         /)
    type(func_type) :: meta_funcs(2) = (/                                     &
         func_type(W3,          GH_BASIS),                                    &
         func_type(Wtheta,      GH_BASIS)                                     &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: sample_eos_pressure_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: sample_eos_pressure_code

contains

!> @brief Sample pressure pointwise using the equation of state
!! @param[in] nlayers Number of layers
!! @param[in,out] exner Pressure field
!! @param[in] rho Density
!! @param[in] theta Potential temperature
!! @param[in] moist_dyn_gas Moist dynamics factor in gas law (1+mv/epsilon)
!! @param[in] ndf_w3 Number of degrees of freedom per cell for w3
!! @param[in] undf_w3 Number of unique degrees of freedom for w3
!! @param[in] map_w3 Dofmap for the cell at the base of the column for w3
!! @param[in] w3_basis Basis functions evaluated at the W3 DoFs
!! @param[in] ndf_wt Number of degrees of freedom per cell for theta space
!! @param[in] undf_wt Number of unique degrees of freedom for theta space
!! @param[in] map_wt Dofmap for the cell at the base of the column for theta space
!! @param[in] wt_basis Basis functions evaluated at the W3 DoFs
subroutine sample_eos_pressure_code(nlayers,                           &
                                 exner, rho, theta, moist_dyn_gas, &
                                 ndf_w3, undf_w3, map_w3, w3_basis,&
                                 ndf_wt, undf_wt, map_wt, wt_basis)

  use calc_exner_pointwise_mod, only: calc_exner_pointwise

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_wt, ndf_w3
  integer(kind=i_def), intent(in) :: undf_wt, undf_w3
  integer(kind=i_def), dimension(ndf_wt),  intent(in) :: map_wt
  integer(kind=i_def), dimension(ndf_w3),  intent(in) :: map_w3


  real(kind=r_def), dimension(1,ndf_w3,ndf_w3), intent(in) :: w3_basis
  real(kind=r_def), dimension(1,ndf_wt,ndf_w3), intent(in) :: wt_basis

  real(kind=r_def), dimension(undf_w3),  intent(inout) :: exner
  real(kind=r_def), dimension(undf_w3),  intent(in)    :: rho
  real(kind=r_def), dimension(undf_wt),  intent(in)    :: theta
  real(kind=r_def), dimension(undf_wt),  intent(in)    :: moist_dyn_gas

  ! Internal variables
  integer(kind=i_def) :: df, k, dft, df3
  real(kind=r_def), dimension(ndf_w3)          :: rho_e
  real(kind=r_def), dimension(ndf_wt)          :: theta_vd_e

  real(kind=r_def) :: rho_cell, theta_vd_cell

  ! Sample Exner pointwise from equation of state
  do k = 0, nlayers-1

    do df = 1, ndf_w3
      rho_e(df) = rho( map_w3(df) + k )
    end do

    do df = 1, ndf_wt
      theta_vd_e(df) = theta( map_wt(df) + k ) * moist_dyn_gas( map_wt(df) + k )
    end do

    do df = 1, ndf_w3

      ! Evaluate theta and rho pointwise
      rho_cell = 0.0_r_def
      do df3 = 1, ndf_w3
        rho_cell = rho_cell + rho_e(df3)*w3_basis(1,df3,df)
      end do

      theta_vd_cell = 0.0_r_def
      do dft = 1, ndf_wt
        theta_vd_cell = theta_vd_cell + theta_vd_e(dft)*wt_basis(1,dft,df)
      end do

      ! Calcualte exner
      exner(map_w3(df)+k) = calc_exner_pointwise(rho_cell, theta_vd_cell)

    end do
  end do

end subroutine sample_eos_pressure_code

end module sample_eos_pressure_kernel_mod
