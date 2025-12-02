! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_rho_dry_mod

implicit none

contains

! Subroutine to calculate the dry-mass density of air,
! based on the temperature and water-vapour mixing ratio
subroutine calc_rho_dry( n_points,                                             &
                         temperature, q_vap, pressure, rho_dry )

use comorph_constants_mod, only: R_dry, real_cvprec
use calc_virt_temp_dry_mod, only: calc_virt_temp_dry

implicit none

! Number of points
integer, intent(in) :: n_points

! Air temperature
real(kind=real_cvprec), intent(in) :: temperature(n_points)

! Water vapour mixing ratio
real(kind=real_cvprec), intent(in) :: q_vap(n_points)

! Ambient total pressure
real(kind=real_cvprec), intent(in) :: pressure(n_points)

! Dry-density in kg m-3
real(kind=real_cvprec), intent(out) :: rho_dry(n_points)

! Dry virtual temperature / K
real(kind=real_cvprec) :: virt_temp_dry(n_points)

! Loop counter
integer :: ic

! Calculate dry virtual temperature = T ( 1 + Rv/Rd q_vap )
!                                   = p / ( Rd rho_dry )
call calc_virt_temp_dry( n_points,                                             &
                         temperature, q_vap, virt_temp_dry )

! Rearrange to get actual dry-density
do ic = 1, n_points
  rho_dry(ic) = pressure(ic) / ( R_dry * virt_temp_dry(ic) )
end do

return
end subroutine calc_rho_dry

end module calc_rho_dry_mod
