! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module set_dependent_constants_mod

implicit none

contains


!----------------------------------------------------------------
! Subroutine to set various constants and bits of infrastructure
! which depend on the values of other, user-defined, constants.
!----------------------------------------------------------------
! There are various things in the main constants module
! which need to be set during run-time:
subroutine set_dependent_constants()

use comorph_constants_mod, only: newline, l_init_constants,                    &
                     L_con_ref, L_sub_ref, L_fus_ref,                          &
                     L_con_0, L_sub_0, L_fus_0,                                &
                     melt_temp,                                                &
                     l_cv_rain, l_cv_cf, l_cv_snow, l_cv_graup,                &
                     n_cond_species,                                           &
                     n_cond_species_liq, n_cond_species_ice,                   &
                     i_cond_cl, i_cond_rain,                                   &
                     i_cond_cf, i_cond_snow, i_cond_graup,                     &
                     k_bot_conv, k_top_conv,                                   &
                     params_cl, params_rain,                                   &
                     params_cf, params_snow, params_graup,                     &
                     cond_params,                                              &
                     i_cl, i_rain, i_cf, i_snow, i_graup,                      &
                     cp_vap, cp_liq, cp_ice,                                   &
                     rho_liq, rho_ice, rho_rim,                                &
                     nx_full, ny_full, k_bot_conv, k_top_conv, k_top_init

use raise_error_mod, only: raise_fatal

implicit none

! Character string for error reporting
character(len=50) :: str

! Loop counter
integer :: i_cond

character(len=*), parameter :: routinename                                     &
                               = "SET_DEPENDENT_CONSTANTS"

! Set flag to indicate that we don't need to call this routine
! again!
l_init_constants = .true.

! Set newline character used in error messages etc
newline = new_line("a")

! Check array dimensions are consistent
if ( .not. nx_full > 0 ) then
  call raise_fatal( routinename,                                               &
                    "x-dimension array size nx_full is not a valid "         //&
                    "positive number." )
end if
if ( .not. ny_full > 0 ) then
  call raise_fatal( routinename,                                               &
                    "y-dimension array size ny_full is not a valid "         //&
                    "positive number." )
end if
if ( .not. k_top_conv > k_bot_conv ) then
  call raise_fatal( routinename,                                               &
                    "Highest convection model-level k_top_conv is not "      //&
                    "greater than the lowest convection model-level "        //&
                    "k_bot_conv." )
end if
if ( .not. ( k_top_init >= k_bot_conv .and. k_top_init <= k_top_conv ) ) then
  call raise_fatal( routinename,                                               &
                    "Highest convection initiation model-level k_top_init "  //&
                    "is not one of the allowed convection model-levels "     //&
                    "k_bot_conv : k_top_conv." )
end if

! Latent heat of sublimation is sum of condensation and fusion
L_sub_ref = L_con_ref + L_fus_ref

! Extrapolated latent heats at absolute zero
L_con_0 = L_con_ref + ( cp_liq - cp_vap ) * melt_temp
L_fus_0 = L_fus_ref + ( cp_ice - cp_liq ) * melt_temp
L_sub_0 = L_con_0 + L_fus_0


! Setup addresses of condensed water species in a condensed
! water mixing-ratio super-array, and count total number of
! condensed water fields actually in use...

! NOTE: If adding more species, all the liquid species MUST be
! assigned addresses first, with all the ice species addresses
! coming afterwards.  This is because various points in the code
! assume that the liquid and ice species each form contiguous
! blocks within the condensed water super-array.
n_cond_species = 0

! Liquid species:
! Liquid cloud is always used
n_cond_species = n_cond_species + 1
i_cond_cl = n_cond_species
! Other species are optional...
if ( l_cv_rain ) then
  n_cond_species = n_cond_species + 1
  i_cond_rain = n_cond_species
end if
! End of liquid species; count how many
n_cond_species_liq = n_cond_species

! Ice species:
if ( l_cv_cf ) then
  n_cond_species = n_cond_species + 1
  i_cond_cf = n_cond_species
end if
if ( l_cv_snow ) then
  n_cond_species = n_cond_species + 1
  i_cond_snow = n_cond_species
end if
if ( l_cv_graup ) then
  n_cond_species = n_cond_species + 1
  i_cond_graup = n_cond_species
end if
! End of ice species; count how many
n_cond_species_ice = n_cond_species - n_cond_species_liq


! Allocate list of pointers to properties of each condensed
! water species, and assign pointers
allocate( cond_params(n_cond_species) )
cond_params(i_cond_cl)%pt => params_cl
if ( l_cv_rain )  cond_params(i_cond_rain)%pt  => params_rain
if ( l_cv_cf )    cond_params(i_cond_cf)%pt    => params_cf
if ( l_cv_snow )  cond_params(i_cond_snow)%pt  => params_snow
if ( l_cv_graup ) cond_params(i_cond_graup)%pt => params_graup
! All the constants which are specific to each condensed
! water species can now be accessed from cond_params, with the
! appropriate condensate super-array address i_cond.


! Set condensed water species constants which depend on
! the values of other constants.
! Loop over active condensed water species:
do i_cond = 1, n_cond_species

  ! Store condensed water species name and i_frzmlt in a
  ! character string for error-reporting.
  write(str,"(A,I3)") "For condensed water species: "                        //&
                 trim(adjustl(cond_params(i_cond)%pt % cond_name))  //newline//&
               "i_frzmlt = ", cond_params(i_cond)%pt % i_frzmlt

  ! Set super-array address of the species which each one
  ! converts to when freezing or melting, using the set
  ! indicators i_frzmlt, converting to the corresponding
  ! super-array addresses.
  select case ( cond_params(i_cond)%pt % i_frzmlt )
  case ( i_cl )
    cond_params(i_cond)%pt % i_cond_frzmlt = i_cond_cl
  case ( i_rain )
    cond_params(i_cond)%pt % i_cond_frzmlt = i_cond_rain
  case ( i_cf )
    cond_params(i_cond)%pt % i_cond_frzmlt = i_cond_cf
  case ( i_snow )
    cond_params(i_cond)%pt % i_cond_frzmlt = i_cond_snow
  case ( i_graup )
    cond_params(i_cond)%pt % i_cond_frzmlt = i_cond_graup
  case DEFAULT
    call raise_fatal( routinename,                                             &
           "Freezing or melting condensed water species i_frzmlt "           //&
           "is set to an unrecognised species indicator."           //newline//&
           trim(adjustl(str)) )
  end select

  ! Check that the species that this is set to melt or freeze
  ! into is a valid species for such an action...
  if ( cond_params(i_cond)%pt % i_cond_frzmlt == 0 ) then
    call raise_fatal( routinename,                                             &
           "Freezing or melting condensed water species i_frzmlt "           //&
           "is set to a species which is turned off."               //newline//&
           trim(adjustl(str)) )
  else if ( cond_params(i_cond)%pt % l_ice .eqv.                               &
            cond_params(cond_params(i_cond)%pt % i_cond_frzmlt                 &
                               )%pt % l_ice ) then
    call raise_fatal( routinename,                                             &
           "Freezing or melting condensed water species i_frzmlt "           //&
           "is set to a species with the wrong phase "              //newline//&
           "(ice melting into ice, or liquid freezing into "                 //&
           "liquid)."                                               //newline//&
           trim(adjustl(str)) )
  end if

  ! Set specific heat capacity and density depending on whether
  ! the species is liquid or ice:
  if ( cond_params(i_cond)%pt % l_ice ) then
    ! Ice species
    cond_params(i_cond)%pt % cp = cp_ice
    cond_params(i_cond)%pt % rho = rho_ice
  else
    ! Liquid water species
    cond_params(i_cond)%pt % cp = cp_liq
    cond_params(i_cond)%pt % rho = rho_liq
  end if

end do  ! i_cond = 1, n_cond_species

! Set special reduced density of rimed ice for graupel
if ( i_cond_graup > 0 ) then
  cond_params(i_cond_graup)%pt % rho = rho_rim
end if


return
end subroutine set_dependent_constants


end module set_dependent_constants_mod
