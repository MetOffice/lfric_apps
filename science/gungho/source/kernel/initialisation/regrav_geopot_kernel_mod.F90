!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief   Need to say something here
!> @details More to be said here
module regrav_geopot_kernel_mod

use argument_mod,         only: arg_type,              &
                                GH_FIELD, GH_REAL,     &
                                GH_READ, GH_READWRITE, &
                                CELL_COLUMN
use constants_mod,        only: r_def, i_def
use extrusion_config_mod, only: planet_radius
use fs_continuity_mod,    only: Wtheta, W3
use kernel_mod,           only: kernel_type
use planet_config_mod,    only: gravity, cp

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
type, public, extends(kernel_type) :: regrav_geopot_kernel_type
  private
  type(arg_type) :: meta_args(9) = (/                     &
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta), &
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta), &
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, W3),     &
       arg_type(GH_FIELD, GH_REAL, GH_READ,      Wtheta), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,      Wtheta), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,      Wtheta), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,      W3),     &
       arg_type(GH_FIELD, GH_REAL, GH_READ,      Wtheta), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,      W3)      &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: regrav_geopot_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: regrav_geopot_code

contains

!> @param[in]  nlayers        Number of layers
!> @param[in,out] temperature Absolute temperature field
!> @param[in,out] theta       Potential temperature field
!> @param[in,out] exner       Exner pressure field
!> @param[in]  coriolis_term  Vertical component of the coriolis term
!> @param[in]  moist_dyn_gas  Gas factor 1+ m_v/epsilon
!> @param[in]  moist_dyn_tot  Total mass factor 1 + sum m_x
!> @param[in]  height_w3      Height coordinate in w3
!> @param[in]  height_wth     Height coordinate in wth
!> @param[in]  w3_mask        LBC mask or Dummy mask for w3 space
!> @param[in]  ndf_wt         Number of degrees of freedom per cell for wtheta
!> @param[in]  undf_wt        Total number of degrees of freedom for wtheta
!> @param[in]  map_wt         Dofmap for the cell at column base for wt
!> @param[in]  ndf_w3         Number of degrees of freedom per cell for w3
!> @param[in]  undf_w3        Total number of degrees of freedom for w3
!> @param[in]  map_w3         Dofmap for the cell at column base for w3
subroutine regrav_geopot_code( nlayers,       &
                               temperature,   &
                               theta,         &
                               exner,         &
                               coriolis_term, &
                               moist_dyn_gas, &
                               moist_dyn_tot, &
                               height_w3,     &
                               height_wth,    &
                               w3_mask,       &
                               ndf_wt,        &
                               undf_wt,       &
                               map_wt,        &
                               ndf_w3,        &
                               undf_w3,       &
                               map_w3 )

  implicit none

  ! Arguments
  integer(kind=i_def),                    intent(in) :: nlayers, &
                                                        ndf_w3,  &
                                                        undf_w3, &
                                                        ndf_wt,  &
                                                        undf_wt
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt

  real(kind=r_def), dimension(undf_wt),   intent(inout) :: theta, temperature
  real(kind=r_def), dimension(undf_w3),   intent(inout) :: exner


  real(kind=r_def), dimension(undf_w3),   intent(in) :: height_w3,     &
                                                        w3_mask
  real(kind=r_def), dimension(undf_wt),   intent(in) :: moist_dyn_gas, &
                                                        moist_dyn_tot
  real(kind=r_def), dimension(undf_wt),   intent(in) :: height_wth
  real(kind=r_def), dimension(undf_wt),   intent(in) :: coriolis_term

  ! Internal variables
  integer(kind=i_def) :: k
  real(kind=r_def)    :: temp_virt
  real(kind=r_def)    :: ht_wt(0:nlayers)
  real(kind=r_def)    :: th(0:nlayers)
  real(kind=r_def)    :: exner_surf
  real(kind=r_def)    :: weight1

  ! Return if the mask is 0 (with tolerance of 0.5 as mask is real, 0 or 1)
  ! setting exner to 1
  if ( w3_mask( map_w3(1) ) < 0.5_r_def ) then
    do k = 0, nlayers-1
      exner(map_w3(1) + k ) = 1.0_r_def
    enddo
    return
  end if

  ! Geopotential height
  do k = 0, nlayers
    ht_wt(k) = planet_radius * height_wth(map_wt(1)+k) / &
              ( planet_radius - height_wth(map_wt(1)+k) )
  end do

  temp_virt = moist_dyn_gas( map_wt(1) ) * temperature( map_wt(1) ) / &
              moist_dyn_tot( map_wt(1) )
  weight1 = height_w3( map_w3(1) ) - height_wth( map_wt(1) )
  exner_surf = exner( map_w3(1) ) * ( cp * temp_virt  ) / &
               ( cp * temp_virt - ( gravity - coriolis_term( map_wt(1) ) ) * weight1 )
  theta( map_wt(1) ) = temperature( map_wt(1) ) / exner_surf


! Reuse th as temporary storage for T
  do k = 0, nlayers
    th(k) = temperature( map_wt(1) + k )
  end do

  ! Map temperature from newly computed grid back onto current model grid
  call interp( nlayers+1, ht_wt, th, nlayers+1, height_wth(map_wt(1)), &
               temperature(map_wt(1)) )

end subroutine regrav_geopot_code

subroutine interp( nin, zin, fin, nout, zout, fout )

implicit none

integer(kind=i_def), intent(in)    :: nin, nout
real(kind=r_def),    intent(in)    :: zin(nin), fin(nin), zout(nout)
real(kind=r_def),    intent(inout) :: fout(nout)

integer(kind=i_def) :: kout, kin
real(kind=r_def)    :: xi, zmin, zmax, f_at_min, f_at_max
real(kind=r_def)    :: eta_out, eta_1, eta_2

if ( zin(1) > zin(nin) ) then
  zmax = zin(1)
  zmin = zin(nin)
  f_at_max = fin(1)
  f_at_min = fin(nin)
else
  zmax = zin(nin)
  zmin = zin(1)
  f_at_max = fin(nin)
  f_at_min = fin(1)
end if


do kout = 1, nout
  if ( zout(kout) <= zmin ) then
    fout(kout) = f_at_min
  else if ( zout(kout) >= zmax ) then
    fout(kout) = f_at_max
  else
    do kin = 1, nin - 1
      if ( ( zout(kout) - zin(kin) ) * &
           ( zin(kin+1) - zout(kout) ) >= 0.0_r_def ) then
        eta_out = sqrt(max(zout(kout), 0.0_r_def))
        eta_1 = sqrt(max(zin(kin), 0.0_r_def))
        eta_2 = sqrt(max(zin(kin+1), 0.0_r_def))
        eta_out = zout(kout)
        eta_1 = zin(kin)
        eta_2 = zin(kin+1)
        xi = ( eta_out - eta_1 ) / ( eta_2 - eta_1 )
        fout(kout) = ( 1.0_r_def - xi ) * log(fin(kin)) + xi * log(fin(kin+1))
        fout(kout) = exp(fout(kout))
      end if
    end do
  end if
end do

end subroutine interp

! Does the input UM data, T(z), pi(z), have the correct z information?
! If T and pi have been interpolated onto the LFRic grid by UM2LFRic, why
! is there a problem with the lower boundary? I'm assuming this problem is
! that the orographic height of UM and LFRic are different. but if they are,
! then the UM and LFRic grids are different and the T(z), pi(z) are not
! valid at LFRic grid heights. Need to look at precisely what data is
! coming from um2lfric.

end module regrav_geopot_kernel_mod
