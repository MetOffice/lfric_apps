!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module hydrostatic_coriolis_kernel_mod

!> @brief   Computes exner from equation of state and vertical balance including
!!          Coriolis term.
!> @details Calculate exner on the top level, using equation of state. Then
!!          use finite differencing to calculate exner on successive levels
!!          below. Unlike hydrostatic_eos_kernel which balances upwards, the
!!          Coriolis term is also added to the vertical balance equation.

use argument_mod,               only : arg_type, func_type,      &
                                       GH_FIELD, GH_REAL,        &
                                       GH_SCALAR, GH_INTEGER,    &
                                       GH_READ, GH_WRITE,        &
                                       ANY_SPACE_9, ANY_SPACE_1, &
                                       GH_BASIS, CELL_COLUMN, GH_EVALUATOR
use constants_mod,              only : r_def, i_def
use fs_continuity_mod,          only : Wtheta, W3
use kernel_mod,                 only : kernel_type
use reference_element_mod,      only : B

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
type, public, extends(kernel_type) :: hydrostatic_coriolis_kernel_type
  private
  type(arg_type) :: meta_args(13) = (/                      &
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W3),      &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3),      &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  Wtheta),  &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  Wtheta),  &
       arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  Wtheta),  &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3),      &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3),      &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3),      &
       arg_type(GH_SCALAR,  GH_REAL,    GH_READ),           &
       arg_type(GH_SCALAR,  GH_REAL,    GH_READ),           &
       arg_type(GH_SCALAR,  GH_REAL,    GH_READ),           &
       arg_type(GH_SCALAR,  GH_REAL,    GH_READ),           &
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ)            &
       /)
  type(func_type) :: meta_funcs(2) = (/                     &
       func_type(W3,     GH_BASIS),                         &
       func_type(Wtheta, GH_BASIS)                          &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass :: hydrostatic_coriolis_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: hydrostatic_coriolis_code

contains

!! @param[in]  nlayers       Number of layers
!! @param[in,out] exner      Exner pressure field
!! @param[in]  rho           Density field
!! @param[in]  theta         Potential temperature field
!! @param[in]  coriolis_term Vertical component of the coriolis term
!! @param[in]  moist_dyn_gas Gas factor 1+ m_v/epsilon
!! @param[in]  moist_dyn_tot Total mass factor 1 + sum m_x
!! @param[in]  moist_dyn_fac Water factor
!! @param[in]  phi           Geopotential field
!! @param[in]  height_w3     Height coordinate in w3
!! @param[in]  w3_mask       LBC mask or Dummy mask for w3 space
!! @param[in]  p_zero        Reference surface pressure
!! @param[in]  kappa         Ratio of Rd and cp
!! @param[in]  rd            Gas constant for dry air
!! @param[in]  cp            Specific heat of dry air at constant pressure
!! @param[in]  eos_index     Vertical level at which the equation of state
!!                           is satisfied
!! @param[in]  ndf_w3        Number of degrees of freedom per cell for w3
!! @param[in]  undf_w3       Total number of degrees of freedom for w3
!! @param[in]  map_w3        Dofmap for the cell at column base for w3
!! @param[in]  basis_w3      Basis functions evaluated at w3 nodes
!! @param[in]  ndf_wt        Number of degrees of freedom per cell for wtheta
!! @param[in]  undf_wt       Total number of degrees of freedom for wtheta
!! @param[in]  map_wt        Dofmap for the cell at column base for wt
!! @param[in]  basis_wt      Basis functions evaluated at wt nodes
subroutine hydrostatic_coriolis_code( nlayers,       &
                                      exner,         &
                                      rho,           &
                                      theta,         &
                                      coriolis_term, &
                                      moist_dyn_gas, &
                                      moist_dyn_tot, &
                                      moist_dyn_fac, &
                                      phi,           &
                                      height_w3,     &
                                      w3_mask,       &
                                      p_zero,        &
                                      kappa,         &
                                      rd,            &
                                      cp,            &
                                      eos_index,     &
                                      ndf_w3,        &
                                      undf_w3,       &
                                      map_w3,        &
                                      basis_w3,      &
                                      ndf_wt,        &
                                      undf_wt,       &
                                      map_wt,        &
                                      basis_wt )

  implicit none

  ! Arguments
  integer(kind=i_def),                          intent(in) :: nlayers, &
                                                              ndf_w3,  &
                                                              undf_w3, &
                                                              ndf_wt,  &
                                                              undf_wt
  integer(kind=i_def), dimension(ndf_w3),       intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_wt),       intent(in) :: map_wt
  integer(kind=i_def),                          intent(in) :: eos_index

  real(kind=r_def), dimension(undf_w3),      intent(inout) :: exner
  real(kind=r_def), dimension(undf_w3),         intent(in) :: rho,           &
                                                              height_w3,     &
                                                              w3_mask,       &
                                                              phi
  real(kind=r_def), dimension(undf_wt),         intent(in) :: moist_dyn_gas, &
                                                              moist_dyn_tot, &
                                                              moist_dyn_fac
  real(kind=r_def), dimension(undf_wt),         intent(in) :: theta
  real(kind=r_def), dimension(undf_wt),         intent(in) :: coriolis_term
  real(kind=r_def), dimension(1,ndf_w3,ndf_w3), intent(in) :: basis_w3
  real(kind=r_def), dimension(1,ndf_wt,ndf_w3), intent(in) :: basis_wt
  real(kind=r_def),                             intent(in) :: p_zero
  real(kind=r_def),                             intent(in) :: kappa
  real(kind=r_def),                             intent(in) :: rd
  real(kind=r_def),                             intent(in) :: cp

  ! Internal variables
  integer(kind=i_def)                  :: k, df, dft, df3, eos_index_m1
  real(kind=r_def), dimension(ndf_w3)  :: rho_e
  real(kind=r_def), dimension(ndf_wt)  :: theta_moist_e
  real(kind=r_def)                     :: rho_cell, theta_moist, dz

  ! Return if the mask is 0 (with tolerance of 0.5 as mask is real, 0 or 1)
  if ( w3_mask( map_w3(1) ) < 0.5_r_def ) then
    return
  end if

  ! Compute exner from eqn of state at vertical-level eos_index
  ! Levels in the kernel are numbered 0 to nlayers-1, rather than 1 to nlayers

  eos_index_m1 = eos_index - 1

  k = eos_index_m1

  do df3 = 1, ndf_w3
    rho_e( df3 ) = rho( map_w3(df3) + k)
  end do

  do dft = 1, ndf_wt
    theta_moist_e( dft ) = theta( map_wt(dft) + k) * moist_dyn_gas( map_wt(dft) + k )
  end do

  do df = 1, ndf_w3

    rho_cell = 0.0_r_def
    do df3 = 1, ndf_w3
      rho_cell = rho_cell + rho_e( df3 )*basis_w3( 1,df3,df )
    end do

    theta_moist = 0.0_r_def
    do dft = 1, ndf_wt
      theta_moist = theta_moist + theta_moist_e( dft )*basis_wt( 1,dft,df )
    end do

    exner( map_w3(df) + k ) = ( rd*rho_cell*theta_moist/p_zero ) &
                              **( kappa/( 1.0_r_def-kappa ) )
  end do

  ! Exner on other levels from hydrostatic balance
  ! Compute exner at level k-1 from exner at level k plus other terms.
  ! Add the coriolis term using the bottom face from the cell at level k-1
  do k = eos_index_m1, 1, -1

    dz = height_w3( map_w3(1) + k ) - height_w3( map_w3(1) + k - 1 )
    theta_moist = moist_dyn_gas( map_wt(1) + k ) * theta( map_wt(1) + k ) /   &
                  moist_dyn_tot( map_wt(1) + k )
    exner( map_w3(1) + k - 1 ) = exner( map_w3 (1) + k )         &
       -  ( coriolis_term( map_wt(1) + k - 1 ) * dz              &
       +  ( phi( map_w3(1) + k - 1 ) - phi( map_w3(1) + k  ) ) ) &
       / ( cp * theta_moist )

  end do

  do k = eos_index_m1, nlayers-2

    dz = height_w3( map_w3(1) + k + 1 ) - height_w3( map_w3(1) + k )
    theta_moist = moist_dyn_gas( map_wt(1) + k + 1 ) * theta( map_wt(1) + k + 1) /   &
                  moist_dyn_tot( map_wt(1) + k + 1)
    exner( map_w3(1) + k + 1 ) = exner( map_w3 (1) + k )        &
       +  ( coriolis_term( map_wt(1) + k ) * dz                 &
       +  ( phi( map_w3(1) + k ) - phi( map_w3(1) + k + 1 ) ) ) &
       / ( cp * theta_moist )

  end do

end subroutine hydrostatic_coriolis_code

end module hydrostatic_coriolis_kernel_mod
