!-------------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Map a prescribed (pseudo) profile to full photolysis rates array using
!>        a multiplicative factor.

module map_pseudophotol_kernel_mod

use argument_mod,      only: arg_type, GH_FIELD, GH_SCALAR, GH_INTEGER,        &
                             GH_REAL, GH_READWRITE,  GH_READ,                  &
                             ANY_DISCONTINUOUS_SPACE_1,                        &
                             CELL_COLUMN

use kernel_mod,        only: kernel_type
use fs_continuity_mod, only: WTHETA

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel.
!> Contains the metadata needed by the Psy layer

type, public, extends(kernel_type) :: map_pseudophotol_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/                                          &       
       arg_type( GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1 ),  & ! Photol_rates all species
       arg_type( GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                & ! Photol_rate single
       arg_type( GH_SCALAR, GH_INTEGER, GH_READ),                             & ! Index to process
       arg_type( GH_SCALAR, GH_REAL,    GH_READ)                              & ! Factor to apply
       
       /)
  integer :: operates_on = CELL_COLUMN

contains
  procedure, nopass :: map_pseudophotol_code
end type

public map_pseudophotol_code

contains

!> @brief UKCA Photolysis prescribed calculations.
!> @details Extend a single photolysis rate profile to given species in
!>          full array using a scaling factor derived from sample run

!> @param[in]      nlayers             Number of layers
!> @param[in, out] photol_rates        Photolysis rates all species (s-1)
!> @param[in]      photol_rate_single  Photolysis rate for single species (s-1)
!> @param[in]      jp2                 Index of species to process in full array
!> @param[in]      jrate_fac           Scaling factor to be applied to profile

subroutine map_pseudophotol_code( nlayers,                                     &
                        photol_rates,                                          &
                        photol_rate_single,                                    &
                        jp2,                                                   &
                        jrate_fac,                                             &
                        ndf_wth, undf_wth, map_wth,                            &
                        ndf_nphot, undf_nphot, map_nphot )

  use constants_mod,    only: r_def, i_def
  
  implicit none

  ! Arguments

  integer(kind=i_def), intent(in) :: nlayers  
  integer(kind=i_def), intent(in) :: ndf_wth
  integer(kind=i_def), intent(in) :: undf_wth
  integer(kind=i_def), dimension(ndf_wth), intent(in) :: map_wth
  integer(kind=i_def), intent(in) :: ndf_nphot
  integer(kind=i_def), intent(in) :: undf_nphot
  integer(kind=i_def), dimension(ndf_nphot), intent(in) :: map_nphot
    
  real(kind=r_def),    intent(in out), dimension(undf_nphot) :: photol_rates
  real(kind=r_def),    intent(in), dimension(undf_wth) :: photol_rate_single
  integer(kind=i_def), intent(in)                      :: jp2
  real(kind=r_def),    intent(in)                      :: jrate_fac
    
  ! Local variables for the kernel
  integer(kind=i_def) :: k

  ! Map rates from profile to full array
  do k = 1, nlayers
    photol_rates( map_nphot(1) + ( (jp2-1)*(nlayers+1) ) + k ) =              &
           photol_rate_single( map_wth(1) + k ) * jrate_fac
  end do  
  
return

end subroutine map_pseudophotol_code

end module map_pseudophotol_kernel_mod