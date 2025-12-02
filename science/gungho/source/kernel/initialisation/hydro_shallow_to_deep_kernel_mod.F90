!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief   Defines exner pressure from vertical balance with a deep gravity,
!!          but using input data based on a shallow gravity.
!> @details Using the input exner pressure at level eos_index, integrate down
!!          to the surface using hydrostatic balance, plus Coriolis terms,
!!          using a discretisation using absolute temperature and constant
!!          gravity g (shallow definition). Then integrate from the surface to
!!          the top but using the model geopotential (may be defined as shallow
!!          or deep depending on configuration).
module hydro_shallow_to_deep_kernel_mod

use argument_mod,               only : arg_type, func_type,      &
                                       GH_FIELD, GH_REAL,        &
                                       GH_SCALAR, GH_INTEGER,    &
                                       GH_READ, GH_READWRITE,    &
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
type, public, extends(kernel_type) :: hydro_shallow_to_deep_kernel_type
  private
  type(arg_type) :: meta_args(11) = (/                          &
       arg_type(GH_FIELD,   GH_REAL,    GH_READWRITE, W3),      &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,      Wtheta),  &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,      Wtheta),  &
       arg_type(GH_FIELD*3, GH_REAL,    GH_READ,      Wtheta),  &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,      W3),      &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,      W3),      &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,      Wtheta),  &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,      W3),      &
       arg_type(GH_SCALAR,  GH_REAL,    GH_READ),               &
       arg_type(GH_SCALAR,  GH_REAL,    GH_READ),               &
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ)                &
       /)
  type(func_type) :: meta_funcs(2) = (/                         &
       func_type(W3,     GH_BASIS),                             &
       func_type(Wtheta, GH_BASIS)                              &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass :: hydro_shallow_to_deep_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: hydro_shallow_to_deep_code

contains

!! @param[in]  nlayers       Number of layers
!! @param[in,out] exner      Exner pressure field
!! @param[in]  temperature   Absolute temperature field
!! @param[in]  coriolis_term Vertical component of the coriolis term
!! @param[in]  moist_dyn_gas Gas factor 1+ m_v/epsilon
!! @param[in]  moist_dyn_tot Total mass factor 1 + sum m_x
!! @param[in]  moist_dyn_fac Water factor
!! @param[in]  phi           Geopotential field
!! @param[in]  height_w3     Height coordinate in w3
!! @param[in]  height_wth    Height coordinate in wth
!! @param[in]  w3_mask       LBC mask or Dummy mask for w3 space
!! @param[in]  gravity       The planet gravity
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
subroutine hydro_shallow_to_deep_code( &
                                      nlayers,       &
                                      exner,         &
                                      temperature,   &
                                      coriolis_term, &
                                      moist_dyn_gas, &
                                      moist_dyn_tot, &
                                      moist_dyn_fac, &
                                      phi,           &
                                      height_w3,     &
                                      height_wth,    &
                                      w3_mask,       &
                                      gravity,       &
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
  real(kind=r_def), dimension(undf_w3),         intent(in) :: height_w3,     &
                                                              w3_mask,       &
                                                              phi
  real(kind=r_def), dimension(undf_wt),         intent(in) :: moist_dyn_gas, &
                                                              moist_dyn_tot, &
                                                              moist_dyn_fac
  real(kind=r_def), dimension(undf_wt),         intent(in) :: temperature,   &
                                                              height_wth
  real(kind=r_def), dimension(undf_wt),         intent(in) :: coriolis_term
  real(kind=r_def), dimension(1,ndf_w3,ndf_w3), intent(in) :: basis_w3
  real(kind=r_def), dimension(1,ndf_wt,ndf_w3), intent(in) :: basis_wt
  real(kind=r_def),                             intent(in) :: gravity
  real(kind=r_def),                             intent(in) :: cp

  ! Internal variables
  integer(kind=i_def)                  :: k, eos_index_m1
  real(kind=r_def)                     :: temp_moist, dz, &
                                          weight1, weight2

  ! Return if the mask is 0 (with tolerance of 0.5 as mask is real, 0 or 1)
  ! setting exner to 1
  if ( w3_mask( map_w3(1) ) < 0.5_r_def ) then
    do k = 0, nlayers-1
      exner(map_w3(1) + k ) = 1.0_r_def
    enddo
    return
  end if

  ! Levels in the kernel are numbered 0 to nlayers-1, rather than 1 to nlayers
  eos_index_m1 = eos_index - 1

  ! Integrate from eos_level to surface using shallow gravity
  do k = eos_index_m1, 1, -1

    ! Absolute temperature at which dry air would have same pressure and density
    ! T_moist = T (1 + mr/eps)/(1 + sum mr)
    temp_moist = moist_dyn_gas( map_wt(1) + k ) * &
                 temperature( map_wt(1) + k ) / &
                 moist_dyn_tot( map_wt(1) + k )

    ! Discretisation weights e.g. weight 1 = w * dz where w = (z_k - z_k-1/2 ) / dz
    ! and weight 2 = (1 - w) * dz, dz = z_k - z_k-1
    !
    ! z_k      W3 k
    ! z_k-1/2  Wtheta k
    ! z_k-1    W3 k-1
    weight1 = height_w3( map_w3(1) + k ) - height_wth( map_wt(1) + k )
    weight2 = height_wth( map_wt(1) + k ) -  height_w3( map_w3(1) + k - 1 )

    ! Vertical momentum equation in hydrostatic balance (plus coriolis) with T = theta Pi
    ! cp T d Pi /dz = Pi ( F - d Phi/ dz )
    ! Discretised, using d Phi /dz = g as
    ! cp T_{k-1/2} ( Pi_k - Pi_k-1 ) / dz = ( w Pi_k + (1-w) Pi_k-1 ) ( F_k-1/2 - g)
    ! Rearranged as
    ! Pi_k-1 = Pi_k * ( cp T_{k-1/2} + (g - F_k-1/2) * (1-w) * dz ) /
    !                 ( cp T_{k-1/2} - (g - F_k-1/2) * w * dz )
    exner( map_w3 (1) + k-1 ) = exner( map_w3(1) + k ) * &
      ( cp * temp_moist + ( gravity - coriolis_term( map_wt(1) + k ) ) * weight2 ) / &
      ( cp * temp_moist - ( gravity - coriolis_term( map_wt(1) + k ) ) * weight1 )

  end do

  ! Integrate from surface to top using model geopotential phi
  do k = 1, nlayers-1

    temp_moist = moist_dyn_gas( map_wt(1) + k ) * temperature( map_wt(1) + k ) / &
         moist_dyn_tot( map_wt(1) + k )

    dz = height_w3( map_w3(1) + k ) - height_w3( map_w3(1) + k - 1 )
    weight1 = ( height_w3( map_w3(1) + k ) - height_wth( map_wt(1) + k ) ) / dz
    weight2 = ( height_wth( map_wt(1) + k ) -  height_w3( map_w3(1) + k - 1 ) ) /dz

    ! Pi_k = Pi_k-1 * ( cp T_k-1/2 - ( (phi_k - phi_k-1) - F_k-1/2 * dz) * (1-w) ) /
    !                 ( cp T_k-1/2 + ( (phi_k - phi_k-1) - F_k-1/2 * dz) * w )
    exner( map_w3(1) + k ) = exner( map_w3 (1) + k-1 ) * &
      ( cp * temp_moist - ( phi( map_w3(1) + k ) - phi( map_w3(1) + k - 1 ) &
                            - coriolis_term( map_wt(1) + k ) * dz ) * weight1 ) / &
      ( cp * temp_moist + ( phi( map_w3(1) + k ) - phi( map_w3(1) + k - 1 )  &
                            - coriolis_term( map_wt(1) + k ) * dz ) * weight2 )

  end do

end subroutine hydro_shallow_to_deep_code

end module hydro_shallow_to_deep_kernel_mod
