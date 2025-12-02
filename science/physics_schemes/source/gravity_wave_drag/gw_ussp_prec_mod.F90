! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module gw_ussp_prec_mod

! Purpose: Holds variables and parameters used in the USSP scheme in
!          real_usprec precision.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: gravity_wave_drag

! parameters
use um_types, only: real_usprec
use conversions_mod, only: pi_orig => pi


implicit none

character(len=*), parameter, private :: ModuleName='GW_USSP_PREC_MOD'

! parameters
real (kind=real_usprec), parameter :: pi = real(pi_orig, kind=real_usprec)

! variables
real (kind=real_usprec) :: two_omega, ussp_launch_factor, wavelstar,           &
    cgw_scale_factor

contains

subroutine gw_ussp_prec_set_reals()

! Purpose: Set runtime constants that we need to be in real_usprec precision.

! variables from other modules. We alias them as *_orig.
use planet_constants_mod, only: two_omega_orig => two_omega
use g_wave_input_mod, only: ussp_launch_factor_orig => ussp_launch_factor,     &
    wavelstar_orig => wavelstar,                                               &
    cgw_scale_factor_orig => cgw_scale_factor

implicit none

two_omega = real(two_omega_orig, kind=real_usprec)
ussp_launch_factor = real(ussp_launch_factor_orig, kind=real_usprec)
wavelstar = real(wavelstar_orig, kind=real_usprec)
cgw_scale_factor = real(cgw_scale_factor_orig, kind=real_usprec)

end subroutine gw_ussp_prec_set_reals


end module gw_ussp_prec_mod
