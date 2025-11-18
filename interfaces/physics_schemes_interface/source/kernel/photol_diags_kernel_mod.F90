!-------------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Interface to Photolysis diagnostics calculations

module photol_diags_kernel_mod

use argument_mod,      only: arg_type, GH_FIELD, GH_SCALAR, GH_INTEGER,        &
                             GH_REAL, GH_READWRITE, GH_READ,                   &
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

type, public, extends(kernel_type) :: photol_diags_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                          &       
       arg_type( GH_FIELD,  GH_REAL,     GH_READWRITE, WTHETA),                & ! Photol_rate single
       arg_type( GH_SCALAR, GH_INTEGER,  GH_READ),                             & ! Index to extract
       arg_type( GH_FIELD,  GH_REAL,     GH_READ, ANY_DISCONTINUOUS_SPACE_1 )  & ! Photol_rates all species
       
       /)
  integer :: operates_on = CELL_COLUMN

contains
  procedure, nopass :: photol_diags_code
end type

public photol_diags_code

contains

!> @brief UKCA Photolysis diagnostic calculation.
!> @details Copy specified photolysis rates from full photol rates array to 
!           a 3-D diagnostic field.

!> @param[in]     nlayers             Number of layers
!> @param[in,out] photol_rate_single  Photolysis rate for single species (s-1)
!> @param[in]     photol_rates        Photolysis rates all species (s-1)
!> @param[in]     jp1                 Index of species to extract

subroutine photol_diags_code( nlayers,                                         &                        
                        photol_rate_single,                                    &
                        jp1,                                                   &
                        photol_rates,                                          &                        
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
  
  real(kind=r_def),    intent(in out), dimension(undf_wth)   :: photol_rate_single
  integer(kind=i_def), intent(in)                            :: jp1
  real(kind=r_def),    intent(in),     dimension(undf_nphot) :: photol_rates
    
  ! Local variables for the kernel
  integer(kind=i_def) :: k  

  ! Transfer rates from full array to diags array
  do k = 1, nlayers      
      photol_rate_single( map_wth(1) + k ) =                                   &
             photol_rates( map_nphot(1) + ( (jp1-1)*(nlayers+1) ) + k )
  end do  
  
return

end subroutine photol_diags_code

end module photol_diags_kernel_mod
