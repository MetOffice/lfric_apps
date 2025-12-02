! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module turb_type_mod

use comorph_constants_mod, only: real_hmprec, name_length
use fields_type_mod, only: fields_list_type

implicit none


!----------------------------------------------------------------
! Type definition containing pointers to the turbulence fields
! which maybe used by CoMorph
!----------------------------------------------------------------
type :: turb_type

  ! These are all assumed to be defined on half-levels between
  ! the main thermodynamic fields, where e.g. w_var(:,:,k)
  ! is defined at height_half(:,:,k), where
  ! height_full(k-1) < height_half(k) < height_full(k)

  ! Turbulent Vertical Velocity Variance / m2 s-2
  real(kind=real_hmprec), pointer :: w_var(:,:,:) => null()

  ! Turbulent fluxes.
  ! These are all assumed to be defined in the upwards direction.
  ! Upwards momentum flux = minus the downwards turbulent stress.
  ! These are all asumed to be in units such that
  !   f_x = <w'x'> = -k dx/dz

  ! Flux of liquid water temperature / K m s-1
  real(kind=real_hmprec), pointer :: f_templ(:,:,:) => null()
  ! Flux of total-water mixing-ratio / kg kg-1 m s-1
  real(kind=real_hmprec), pointer :: f_q_tot(:,:,:) => null()
  ! Flux of u-wind momentum / m2 s-2
  real(kind=real_hmprec), pointer :: f_wind_u(:,:,:) => null()
  ! Flux of v-wind momentum / m2 s-2
  real(kind=real_hmprec), pointer :: f_wind_v(:,:,:) => null()

  ! List of pointers to all the above fields,
  ! for the purposes of looping over them to do
  ! the same thing to all the fields
  type(fields_list_type), allocatable :: list(:)

  ! Miscellaneous extra fields not included in the list...

  ! Turbulence length-scale, defined on full-levels
  real(kind=real_hmprec), pointer :: lengthscale(:,:,:) => null()

  ! Parcel radius amplification factor
  real(kind=real_hmprec), pointer :: par_radius_amp(:,:) => null()

  ! Boundary-layer top height from the boundary-layer scheme
  real(kind=real_hmprec), pointer :: z_bl_top(:,:) => null()

end type turb_type


!----------------------------------------------------------------
! Addresses of fields in the lists and in compressed super-arrays
!----------------------------------------------------------------

! Total number of fields in the list
integer, parameter :: n_turb = 5

! Address of each turb field
integer, parameter :: i_w_var = 1
integer, parameter :: i_f_templ = 2
integer, parameter :: i_f_q_tot = 3
integer, parameter :: i_f_wind_u = 4
integer, parameter :: i_f_wind_v = 5

! Array of strings storing name of each turbulence field,
! for error reporting
character(len=name_length), parameter :: turb_names(n_turb)                    &
              = [ "w_var   ", "f_templ ", "f_q_tot ", "f_wind_u", "f_wind_v" ]

! Flag for whether each turbulence field should only be positive
logical, parameter :: turb_positive(n_turb)                                    &
              = [ .true., .false., .false., .false., .false. ]
! TKE and diffusivities must be positive,
! but fluxes could go either way.

! Max allowed ratio of heat-flux over sqrt(TKE) / K
real(kind=real_hmprec), parameter :: max_templ = 20.0
! Max allowed ratio of moisture-flux over sqrt(TKE) / kg kg-1
real(kind=real_hmprec), parameter :: max_q_tot = 0.05
! Max allowed ratio of momentum-flux over sqrt(TKE) / m s-1
real(kind=real_hmprec), parameter :: max_wind = 20.0


contains


!----------------------------------------------------------------
! Subroutine to nullify the pointers in a turb_type structure
!----------------------------------------------------------------
subroutine turb_nullify( turb )
implicit none
type(turb_type), intent(in out) :: turb

turb % w_var => null()
turb % f_templ => null()
turb % f_q_tot => null()
turb % f_wind_u => null()
turb % f_wind_v => null()

turb % lengthscale => null()
turb % par_radius_amp => null()
turb % z_bl_top => null()

return
end subroutine turb_nullify


!----------------------------------------------------------------
! Subroutine to assign a list of pointers to the full 3-D fields
!----------------------------------------------------------------
! This is used to allow a super do-loop to be applied when
! we want to do exactly the same thing to all the fields.
! This routine also checks that all the required fields exist
! and have the right size and shape.
subroutine turb_list_make( turb )

use comorph_constants_mod, only: nx_full, ny_full, k_bot_conv, k_top_init,     &
                                 newline
use raise_error_mod, only: raise_fatal

implicit none

! turb object containing pointers to the turbulence fields
type(turb_type), intent(in out) :: turb

! Bounds of the input array pointers
integer :: lb(3), ub(3)

! Loop counter
integer :: i_turb

character(len=*), parameter :: routinename = "TURB_LIST_MAKE"


! Allocate list of pointers
allocate( turb%list(n_turb) )

! Assign pointers to fields
turb % list(i_w_var)%pt    => turb % w_var
turb % list(i_f_templ)%pt  => turb % f_templ
turb % list(i_f_q_tot)%pt  => turb % f_q_tot
turb % list(i_f_wind_u)%pt => turb % f_wind_u
turb % list(i_f_wind_v)%pt => turb % f_wind_v
! Note: can't include lengthscale or z_bl_top in the list because they
! aren't defined on half-levels like the others.

! Now loop over all the fields and check that the contained
! pointers all actually point at arrays with the required extent
do i_turb = 1, n_turb
  if ( .not. associated( turb % list(i_turb)%pt ) ) then
    call raise_fatal( routinename,                                             &
           "A required input turbulence field has not been "                 //&
           "assigned to any memory.  You must point:"               //newline//&
           "turb % " // trim(adjustl(turb_names(i_turb)))           //newline//&
           "at a target array before the call to comorph_ctl." )
  else
    lb = lbound( turb % list(i_turb)%pt )
    ub = ubound( turb % list(i_turb)%pt )
    ! Note: required upper bound for k is k_top_init+1, since the
    ! these fields are on half-levels, and the upper model-level interface
    ! of the top convective initiation level is needed.
    if ( .not. ( lb(1) <= 1 .and. ub(1) >= nx_full                             &
           .and. lb(2) <= 1 .and. ub(2) >= ny_full                             &
           .and. lb(3) <= k_bot_conv .and. ub(3) >= k_top_init+1 ) ) then
      call raise_fatal( routinename,                                           &
             "A required input turbulence field has "                        //&
             "insufficient extent / the wrong shape."               //newline//&
             "Field: turb % " // trim(adjustl(turb_names(i_turb))) )
    end if
  end if
end do

! Check turbulence length-scale
if ( .not. associated( turb % lengthscale) ) then
  call raise_fatal( routinename,                                               &
         "A required input turbulence field has not been "                 //  &
         "assigned to any memory.  You must point:"               //newline//  &
         "turb % lengthscale"                                     //newline//  &
         "at a target array before the call to comorph_ctl." )
else
  lb = lbound( turb % lengthscale )
  ub = ubound( turb % lengthscale )
  if ( .not. ( lb(1) <= 1 .and. ub(1) >= nx_full                               &
         .and. lb(2) <= 1 .and. ub(2) >= ny_full                               &
         .and. lb(3) <= k_bot_conv .and. ub(3) >= k_top_init ) ) then
    call raise_fatal( routinename,                                             &
           "A required input turbulence field has "                        //  &
           "insufficient extent / the wrong shape."               //newline//  &
           "Field: turb % lengthscale" )
  end if
end if

! Check radius amplification factor
if ( .not. associated( turb % par_radius_amp) ) then
  call raise_fatal( routinename,                                               &
         "A required input turbulence field has not been "                 //  &
         "assigned to any memory.  You must point:"               //newline//  &
         "turb % par_radius_amp"                                  //newline//  &
         "at a target array before the call to comorph_ctl." )
else
  lb(1:2) = lbound( turb % par_radius_amp )
  ub(1:2) = ubound( turb % par_radius_amp )
  if ( .not. ( lb(1) <= 1 .and. ub(1) >= nx_full                               &
         .and. lb(2) <= 1 .and. ub(2) >= ny_full ) ) then
    call raise_fatal( routinename,                                             &
           "A required input turbulence field has "                        //  &
           "insufficient extent / the wrong shape."               //newline//  &
           "Field: turb % par_radius_amp" )
  end if
end if

! Check for the 2-D BL-top height field
if ( .not. associated( turb % z_bl_top ) ) then
  call raise_fatal( routinename,                                               &
           "A required input turbulence field has not been "                 //&
           "assigned to any memory.  You must point:"               //newline//&
           "turb % z_bl_top"                                        //newline//&
           "at a target array before the call to comorph_ctl." )
else
  lb(1:2) = lbound( turb % z_bl_top )
  ub(1:2) = ubound( turb % z_bl_top )
  if ( .not. ( lb(1) <= 1 .and. ub(1) >= nx_full                               &
         .and. lb(2) <= 1 .and. ub(2) >= ny_full ) ) then
    call raise_fatal( routinename,                                             &
             "A required input turbulence field has "                        //&
             "insufficient extent / the wrong shape."               //newline//&
             "Field: turb % z_bl_top" )
  end if
end if


return
end subroutine turb_list_make


!----------------------------------------------------------------
! Subroutine to clear the above list of pointers
!----------------------------------------------------------------
subroutine turb_list_clear( turb )

implicit none

type(turb_type), intent(in out) :: turb

integer :: i_turb

! Nullify the pointers
do i_turb = 1, n_turb
  turb % list(i_turb)%pt => null()
end do

! Deallocate the list
deallocate( turb % list )

return
end subroutine turb_list_clear


!----------------------------------------------------------------
! Subroutine to check consistency between the turbulence fields
!----------------------------------------------------------------
subroutine turb_check_consistent( turb )

use comorph_constants_mod, only: nx_full, ny_full, k_bot_conv, k_top_init,     &
                                 newline, i_check_turb_consistent,             &
                                 i_check_bad_warn, i_check_bad_fatal
use raise_error_mod, only: raise_warning, raise_fatal

implicit none

! Input turb structure to be checked
type(turb_type), intent(in) :: turb

! Turbulent vertical velocity scale = sqrt(w_var)
real(kind=real_hmprec) :: w_pert

! Max length of string containing coordinates and values of
! any inconsistent values found
integer, parameter :: str_len = 500

! String containing coordinates and values of any bad values found
character(len=str_len) :: bad_values

! String containing coordinates and value for current bad value
character(len=100) :: current_bad

! Strings for the i,j,k coordinates
character(len=5) :: i_str
character(len=5) :: j_str
character(len=5) :: k_str
! String for the value
character(len=77) :: val_str

! Counter for current write position in bad_values string
integer :: n
! Number of characters occupied by current bad value string
integer :: dn

! Loop counters
integer :: i, j, k

character(len=*), parameter :: routinename = "TURB_CHECK_CONSISTENT"

! Initialise
n = 0
bad_values(1:str_len) = " "

! Turb fields defined on half-levels; used from lower interface
! of lowest full-level to upper interface of highest full-level
! used for convective initiation
outer: do k = k_bot_conv, k_top_init+1
  do j = 1, ny_full
    do i = 1, nx_full

      ! Find turbulent vertical velocity scale
      w_pert = sqrt(turb % w_var(i,j,k))

      ! Check fluxes don't exceed a max perturbation
      ! scale times the w scale.  If they do, this will imply
      ! implausibly large turbulent fluctuations in the scalars
      ! and winds, probably implying a problem with the input TKE.
      if ( abs(turb % f_templ(i,j,k)) > max_templ * w_pert .or.                &
           abs(turb % f_q_tot(i,j,k)) > max_q_tot * w_pert .or.                &
           abs(turb % f_wind_u(i,j,k)) > max_wind * w_pert .or.                &
           abs(turb % f_wind_v(i,j,k)) > max_wind * w_pert ) then

        ! If bad value found...

        ! Write coordinates and value to strings
        write(i_str,"(I5)") i
        write(j_str,"(I5)") j
        write(k_str,"(I5)") k
        write(val_str,"(7ES11.3)") turb%w_var(i,j,k),                          &
                                   turb%f_templ(i,j,k), turb%f_q_tot(i,j,k),   &
                                   turb%f_wind_u(i,j,k), turb%f_wind_v(i,j,k)

        ! Construct combined string with the coordinates and value
        current_bad = "(" // trim(adjustl(i_str)) // "," //                    &
                             trim(adjustl(j_str)) // "," //                    &
                             trim(adjustl(k_str)) // "): " //                  &
                             trim(adjustl(val_str))
        current_bad = adjustl(current_bad)

        ! Find number of characters needed to print it
        dn = len(trim(current_bad)) + 2  ! + 2 blanks to space out

        if ( n + dn <= str_len ) then
          ! If enough space to add this bad value to the error string,
          ! copy it in
          bad_values(n+1:n+dn) = current_bad(1:dn)
        else
          ! No more space in the string; exit
          exit outer
        end if

        ! Increment counter for space used up so far
        n = n + dn

      end if

    end do
  end do
end do outer

! If any inconsistency detected
if ( n > 0 ) then
  select case(i_check_turb_consistent)
  case (i_check_bad_fatal)
    ! Fatal error: kill the run.
    call raise_fatal( routinename,                                             &
                      "Input TKE is inconsistent with the input fluxes "    // &
                      "or diffusivities."                          //newline// &
                      "In at least the following: (i,j,k) "                 // &
                      "w_var k_mom k_sca f_templ f_q_tot f_wind_u f_wind_v : " &
                                                                   //newline// &
                      trim(adjustl(bad_values)) )
  case (i_check_bad_warn)
    ! Not fatal error: just print a warning.
    call raise_warning( routinename,                                           &
                      "Input TKE is inconsistent with the input fluxes "    // &
                      "or diffusivities."                          //newline// &
                      "In at least the following: (i,j,k) "                 // &
                      "w_var k_mom k_sca f_templ f_q_tot f_wind_u f_wind_v : " &
                                                                   //newline// &
                      trim(adjustl(bad_values)) )
  end select
end if  ! ( n > 0 )

return
end subroutine turb_check_consistent


end module turb_type_mod
