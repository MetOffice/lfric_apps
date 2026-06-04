!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief   Defines exner pressure from vertical balance with a deep gravity,
!>          but using input data based on a shallow gravity.
!> @details Using the input exner pressure at level eos_index, integrate down
!>          to the surface using hydrostatic balance, plus Coriolis terms,
!>          using a discretisation using absolute temperature and constant
!>          gravity g (shallow definition). Then integrate from the surface to
!>          the top but using the model geopotential (may be defined as shallow
!>          or deep depending on configuration).
module regrav_preserve_t_kernel_mod

use argument_mod,               only : arg_type,                 &
                                       GH_FIELD, GH_REAL,        &
                                       GH_SCALAR, GH_INTEGER,    &
                                       GH_READ, GH_READWRITE,    &
                                       CELL_COLUMN
use constants_mod,              only : r_def, i_def
use fs_continuity_mod,          only : Wtheta, W3
use kernel_mod,                 only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
type, public, extends(kernel_type) :: regrav_preserve_t_kernel_type
  private
  type(arg_type) :: meta_args(10) = (/                      &
       arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, W3),      &
       arg_type(GH_FIELD,  GH_REAL, GH_READ,      Wtheta),  &
       arg_type(GH_FIELD,  GH_REAL, GH_READ,      Wtheta),  &
       arg_type(GH_FIELD,  GH_REAL, GH_READ,      Wtheta),  &
       arg_type(GH_FIELD,  GH_REAL, GH_READ,      Wtheta),  &
       arg_type(GH_FIELD,  GH_REAL, GH_READ,      W3),      &
       arg_type(GH_FIELD,  GH_REAL, GH_READ,      W3),      &
       arg_type(GH_FIELD,  GH_REAL, GH_READ,      Wtheta),  &
       arg_type(GH_FIELD,  GH_REAL, GH_READ,      W3),      &
       arg_type(GH_SCALAR, GH_REAL, GH_READ)
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: regrav_preserve_t_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: regrav_preserve_t_code

contains

!> @param[in]  nlayers       Number of layers
!> @param[in,out] exner      Exner pressure field
!> @param[in]  temperature   Absolute temperature field
!> @param[in]  coriolis_term Vertical component of the coriolis term
!> @param[in]  moist_dyn_gas Gas factor 1+ m_v/epsilon
!> @param[in]  moist_dyn_tot Total mass factor 1 + sum m_x
!> @param[in]  phi           Geopotential field
!> @param[in]  height_w3     Height coordinate in w3
!> @param[in]  height_wth    Height coordinate in wth
!> @param[in]  w3_mask       LBC mask or Dummy mask for w3 space
!> @param[in]  cp            Specific heat capacity at constant pressure
!> @param[in]  ndf_w3        Number of degrees of freedom per cell for w3
!> @param[in]  undf_w3       Total number of degrees of freedom for w3
!> @param[in]  map_w3        Dofmap for the cell at column base for w3
!> @param[in]  ndf_wt        Number of degrees of freedom per cell for wtheta
!> @param[in]  undf_wt       Total number of degrees of freedom for wtheta
!> @param[in]  map_wt        Dofmap for the cell at column base for wt
subroutine regrav_preserve_t_code( nlayers,       &
                                 exner,         &
                                 temperature,   &
                                 coriolis_term, &
                                 moist_dyn_gas, &
                                 moist_dyn_tot, &
                                 phi,           &
                                 height_w3,     &
                                 height_wth,    &
                                 w3_mask,       &
                                 cp,            &
                                 ndf_w3,        &
                                 undf_w3,       &
                                 map_w3,        &
                                 ndf_wt,        &
                                 undf_wt,       &
                                 map_wt )

  implicit none

  ! Arguments

  integer(kind=i_def),                    intent(in) :: nlayers, &
                                                        ndf_w3,  &
                                                        undf_w3, &
                                                        ndf_wt,  &
                                                        undf_wt
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt

  real(kind=r_def), dimension(undf_w3), intent(inout) :: exner
  real(kind=r_def), dimension(undf_w3), intent(in)    :: height_w3, &
                                                         w3_mask,   &
                                                         phi
  real(kind=r_def), dimension(undf_wt), intent(in)    :: moist_dyn_gas, &
                                                         moist_dyn_tot
  real(kind=r_def), dimension(undf_wt), intent(in)    :: temperature, &
                                                         height_wth
  real(kind=r_def), dimension(undf_wt), intent(in)    :: coriolis_term
  real(kind=r_def),                     intent(in)    :: cp

  ! Internal variables
  integer(kind=i_def) :: k
  real(kind=r_def)    :: temp_virtual
  real(kind=r_def)    :: dz, weight1, weight2

  ! Return if the mask is 0 (with tolerance of 0.5 as mask is real, 0 or 1)
  ! setting exner to 1
  if ( w3_mask( map_w3(1) ) < 0.5_r_def ) then
    do k = 0, nlayers-1
      exner(map_w3(1) + k ) = 1.0_r_def
    enddo
    return
  end if

  ! Integrate from surface to top using model geopotential phi
  do k = 1, nlayers-1

    temp_virtual = moist_dyn_gas( map_wt(1) + k ) * temperature( map_wt(1) + k ) / &
                   moist_dyn_tot( map_wt(1) + k )

    dz = height_w3( map_w3(1) + k ) - height_w3( map_w3(1) + k - 1 )
    weight1 = ( height_w3( map_w3(1) + k ) - height_wth( map_wt(1) + k ) ) / dz
    weight2 = ( height_wth( map_wt(1) + k ) -  height_w3( map_w3(1) + k - 1 ) ) /dz

    ! Pi_k = Pi_k-1 * ( cp T_k-1/2 - ( (phi_k - phi_k-1) - F_k-1/2 * dz) * (1-w) ) /
    !                 ( cp T_k-1/2 + ( (phi_k - phi_k-1) - F_k-1/2 * dz) * w )
    exner( map_w3(1) + k ) = exner( map_w3 (1) + k-1 ) * &
      ( cp * temp_virtual - ( phi( map_w3(1) + k ) - phi( map_w3(1) + k - 1 )     &
                            - coriolis_term( map_wt(1) + k ) * dz ) * weight1 ) / &
      ( cp * temp_virtual + ( phi( map_w3(1) + k ) - phi( map_w3(1) + k - 1 )     &
                            - coriolis_term( map_wt(1) + k ) * dz ) * weight2 )

  end do

end subroutine regrav_preserve_t_code

end module regrav_preserve_t_kernel_mod
