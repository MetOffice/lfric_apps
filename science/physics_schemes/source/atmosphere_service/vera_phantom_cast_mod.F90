! *********************************COPYRIGHT************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *********************************COPYRIGHT************************************
! Some of the content of this file has been produced with the assistance of
! Anthropic Claude Opus 5 (Claude Code).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: atmos_service_vera
!
! Part of the Vera scheme for diagnosing visual range.
!
! This module contains utility routines for casting phantom aerosol populations.
!

module vera_phantom_cast_mod

use vera_kind_mod,          only: wp => vera_real, wi => vera_integer

use vera_global_mod,        only: vera_constants          ,                    &
                                  vera_phantom            ,                    &
                                  vera_population_type    ,                    &
                                  vera_scaled             ,                    &
                                  vera_update_population

use vera_phantom_tools_mod, only: vera_lin_spacing        ,                    &
                                  vera_lin_intervals      ,                    &
                                  vera_log_spacing        ,                    &
                                  vera_log_normal         ,                    &
                                  vera_triangle           ,                    &
                                  vera_linear_flat_distn

implicit none

! Description:
!   This module contains utility subroutines for casting phantom
!   aerosol populations.
!
! Method:
!   The two suroutines in this module are:
!
!     vera_volume_scale
!       Scales the aerosol population to have the same total aerosol
!       volume as the old monodisperse MURK aerosol field.
!
!     vera_cast_rd_log_normal_b0_triangle
!       Casts a phantom aerosol population with a single log-normal size
!       distribution and a triangular hygroscopy distribution
!
!   For more detail, please refer to the Vera user guide.
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 standards.

! name of this module
character (len=*), parameter, private :: ModuleName='VERA_PHANTOM_CAST_MOD'

private

! make all the subroutines in this module Public,
! they are all called by vera_phantom_list_mod.F90
public :: vera_volume_scale                    ,                               &
          vera_cast_rd_b0                      ,                               &
          vera_cast_rd_log_normal_b0_triangle  ,                               &
          vera_cast_rd_log_normal_b0

contains

  !=============================================================================
  !
  ! vera_volume_scale
  !
  ! Scales the aerosol population to have the same total aerosol volume as the
  ! old monodisperse MURK aerosol field.
  !
  !=============================================================================

subroutine vera_volume_scale( cast_population )

  ! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!===========================================================================
! arguments for vera_volume_scale
!
!===========================================================================

! derived type aerosol population
type(vera_population_type), intent(in out) :: cast_population(:)

!===========================================================================
! local variables for vera_volume_scale
!===========================================================================

! ratio of the aerosol population volume to the volume of MURK aerosol
real (wp) :: volume_ratio

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_VOLUME_SCALE'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! executable code for vera_volume_scale
!===========================================================================

! compute the ration of the volume of the aerosol population
! to the volume of the MURK field
volume_ratio = sum( (cast_population%rd ** vera_constants%three) *             &
                     cast_population%nc ) / vera_scaled%nc_rd_cubed

! scale the aerosol number concentration to that the aerosol
! volume matches the MURK volume
cast_population%nc = cast_population%nc / volume_ratio

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_volume_scale

!=============================================================================
!
! vera_cast_rd_b0
!
! Casts a phantom aerosol population using the input values of Rd and B0.
!
! The phantom aerosol population generated is popped in the aerosol
! definition in the Vera global data module.
!
!=============================================================================

subroutine vera_cast_rd_b0( rd, b0, rd_pdf, b0_pdf )

  ! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!===========================================================================
! arguments for vera_cast_rd_b0
!
!===========================================================================

! values of Rd, the dry aerosol particle sizes, to use
real (wp),           intent(in) :: rd(:)

! values of B0, the dry aerosol particle hygroscopies, to use
real (wp),           intent(in) :: b0(:)

! pdf of the Rd values
real (wp), optional, intent(in) :: rd_pdf(:)

! pdf of the B0 values
real (wp), optional, intent(in) :: b0_pdf(:)

!===========================================================================
! local variables for vera_cast_rd_b0
!===========================================================================

! the cast aerosol population - there are three components:
!
!   %rd - list of particles sizes rd
!   %b0 - list of hygroscopy b0
!   %nc - list of number concentration, normalised so that sum(nc) = 1
!
! an aerosol species in the population is defined by the triple (rd, b0, nc)
!
type(vera_population_type), allocatable :: cast_population(:)

! number of input Rd
integer (wi) :: n_rd

! number of input B0
integer (wi) :: n_b0

! pdf of the Rd values to use in this subroutine
real    (wp) :: rd_pdf_use( size( rd ) )

! pdf of the B0 values to use in this subroutine
real    (wp) :: b0_pdf_use( size( b0 ) )

! counter to loop over the rd distribution
integer (wi) :: ii

! counter to loop over the b0 distribution
integer (wi) :: jj

! counter used to populate the serialised (rd*b0) population array
integer (wi) :: kk

! 2-d array of the joint rd and b0 pdf
real    (wp) :: rd_b0( 1:size( rd ), 1:size( b0 ) )

! 2-d array of the rd that map to the rd and b0 pdf
real    (wp) :: rd_rd( 1:size( rd ), 1:size( b0 ) )

! 2-d array of the b0 that map to the rd and b0 pdf
real    (wp) :: b0_b0( 1:size( rd ), 1:size( b0 ) )

! number of elements in the list of rd values, with pdf > 0
integer (wi) :: rd_elements

! number of elements in the list of b0 values, with pdf > 0
integer (wi) :: b0_elements

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_CAST_RD_B0'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! executable code for vera_cast_rd_b0
!===========================================================================

! number of input Rd
n_rd = size( rd )

! number of input B0
n_b0 = size( b0 )

!---------------------------------------------------------------------------
! gather the optional inputs to this subroutine

! set the rd_pdf to use
if (present(rd_pdf)) then

  ! use the rd_pdf passed to this subroutine
  rd_pdf_use = rd_pdf

else

  ! contruct a flat pdf, i.e. all particle sizes have the same number
  ! concentration
  rd_pdf_use = vera_constants%one / size( rd, kind = wp )

end if

! set the b0_pdf to use
if (present(b0_pdf)) then

  ! use the b0_pdf passed to this subroutine
  b0_pdf_use = b0_pdf

else

  ! contruct a flat pdf, i.e. all hygroscopies have the same number
  ! concentration
  b0_pdf_use = vera_constants%one / size( b0, kind = wp )

end if

! construct the joint rd and b0 pdf
! loop over the rd values
do ii = vera_constants%one_i, n_rd

  ! loop over the b0 values
  do jj = vera_constants%one_i, n_b0

    ! compute the joint pdf from the rd and b0 pdfs
    rd_b0(ii, jj) = rd_pdf_use(ii) * b0_pdf_use(jj)

    ! construct the 2-d array that maps the rd to the joint pdf
    rd_rd(ii, jj) = rd(ii)

    ! construct the 2-d array that maps the b0 to the joint pdf
    b0_b0(ii, jj) = b0(jj)

  end do

end do

! now to construct the derived type cast_population that defines the
! cast phantom population, weeding out the zero elements in the joint pdfs

! how many non zero elements in rd pdf?
do ii = vera_constants%one_i, size( rd_pdf_use )
  if ( rd_pdf_use(ii) <= vera_constants%zero ) then
    rd_pdf_use(ii) = - vera_constants%one
  end if
end do

rd_elements = n_rd - int( abs(sum(rd_pdf_use, mask =                           &
                                  (rd_pdf_use <= vera_constants%zero))) +      &
                          vera_constants%half )

!---------------------------------------------------------------------------
! how many non zero elements in b0 pdf?
do ii = vera_constants%one_i, size( b0_pdf_use )
  if ( b0_pdf_use(ii) <= vera_constants%zero ) then
    b0_pdf_use(ii) = - vera_constants%one
  end if
end do

b0_elements = n_b0 - int( abs(sum(b0_pdf_use, mask =                           &
                                  (b0_pdf_use <= vera_constants%zero))) +      &
                          vera_constants%half )

! total number of aersol species in cast phantom population
! is rd_elements*b0_elements
!
! form the derived type cast_population to hold the cast phantom population
if ( rd_elements*b0_elements > 0 ) then
  ! there is at least one aerosol species, so fill out cast_population
  allocate ( cast_population(rd_elements*b0_elements) )
  kk = vera_constants%zero_i
  ! Note loop ordering is intentionally "backwards" here to ensure
  ! desired ordering within the final type
  do jj = vera_constants%one_i, n_b0
    do ii = vera_constants%one_i, n_rd
      if ( rd_b0(ii, jj) > vera_constants%zero ) then
        kk = kk + vera_constants%one_i
        cast_population(kk)%rd = rd_rd(ii, jj)
        cast_population(kk)%b0 = b0_b0(ii, jj)
        cast_population(kk)%nc = rd_b0(ii, jj)
      end if
    end do
  end do

else

  ! there are no aerosol species in the population, so define
  ! cast_population as a zero element array
  allocate ( cast_population(rd_elements*b0_elements) )

end if

! update the aerosol definition in the Vera global data module
call vera_update_population( cast_population )

! tidy up - get rid of the temporary definition of the
!           cast aerosol population
if ( allocated( cast_population ) ) deallocate( cast_population )

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_cast_rd_b0

!=============================================================================
!
! vera_cast_rd_log_normal_b0_triangle
!
! Casts a phantom aerosol population with a single log-normal size
! distribution and a triangular hygroscopy distribution.
!
! The phantom aerosol population generated is popped in the aerosol
! definition in the Vera global data module.
!
!=============================================================================

subroutine vera_cast_rd_log_normal_b0_triangle( n_rd              ,            &
                                                rd_max , rd_min   ,            &
                                                rd_mode, rd_sigma ,            &
                                                n_b0              ,            &
                                                b0_max , b0_min   ,            &
                                                b0_peak             )

  ! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!===========================================================================
! arguments for vera_cast_rd_log_normal_b0_triangle
!
!===========================================================================

! number of rd values required
integer (wi), intent(in) :: n_rd

! maxmimum value of rd
real    (wp), intent(in) :: rd_max

! minmimum value of rd
real    (wp), intent(in) :: rd_min

! mode of the rd distribution
real    (wp), intent(in) :: rd_mode

! geometric standard deviation of the rd distribution,
! i.e. the width of the distribution
real    (wp), intent(in) :: rd_sigma

! number of b0 values required
integer (wi), intent(in) :: n_b0

! maxmimum value of b0 in the distribution
real    (wp), intent(in) :: b0_max

! minmimum value of b0 in the distribution
real    (wp), intent(in) :: b0_min

! b0 value corresponding to the peak of the triangular b0 distribution
! NOTE: x_peak must be in the range ]b0_min, b0_max[
real    (wp), intent(in) :: b0_peak

!===========================================================================
! local variables for vera_cast_rd_log_normal_b0_triangle
!===========================================================================

! computed logarithmically spaced rd values
real    (wp)             :: rd(1:n_rd)

! computed linearly spaced b0 values
real    (wp)             :: b0(1:n_b0)

! computed rd pdf
real    (wp)             :: rd_pdf(1:n_rd)

! computed b0 pdf
real    (wp)             :: b0_pdf(1:n_b0)

! pdf value to compute outside of the range
! b0 = [b0_min, b0_max], rd = [rd_min, rd_max]
! set this to zero, i.e. no chance of a particle outside of the
! specified ranges
real    (wp)             :: boundless = vera_constants%zero

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName =                                &
                                       'VERA_CAST_RD_LOG_NORMAL_B0_TRIANGLE'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! executable code for vera_cast_rd_log_normal_b0_triangle
!===========================================================================

! compute the rd values - these are logarithmically spaced
!                         across the range [rd_min, rd_max]
!                         with mode rd_mode and width rd_sigma
call vera_log_spacing( n_rd, rd_max, rd_min, rd )

! compute the log-normal pdf of the rd values,
! with mode rd_mode and width rd_sigma
call vera_log_normal( rd, rd_mode, rd_sigma     ,                              &
                      rd_max*vera_phantom%nudge ,                              &
                      rd_min/vera_phantom%nudge ,                              &
                      boundless, rd_pdf           )

! compute the b0 values - these are linearly spaced
!                         across the range [b0_min, b0_max]
call vera_lin_intervals( n_b0, b0_max, b0_min, b0 )

! compute the triangular pdf of the b0 values
call vera_triangle( b0, b0_peak, b0_max, b0_min, boundless, b0_pdf )

! cast an aerosol population with the computed rd, rd_pdf, b0, b0_pdf
! NOTE - this call will update the aerosol population definition
!        in the Vera global data module
call vera_cast_rd_b0( rd, b0, rd_pdf = rd_pdf, b0_pdf = b0_pdf )

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_cast_rd_log_normal_b0_triangle

!=============================================================================
!
! vera_cast_rd_log_normal_b0
!
! Casts a phantom aerosol population with a single log-normal size
! distribution and the input values of B0.
!
! The phantom aerosol population generated is popped in the aerosol
! definition in the Vera global data module.
!
!=============================================================================

subroutine vera_cast_rd_log_normal_b0( n_rd              ,                     &
                                       rd_max , rd_min   ,                     &
                                       rd_mode, rd_sigma ,                     &
                                       b0     , b0_pdf     )

  ! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!===========================================================================
! arguments for vera_cast_rd_log_normal_b0
!
!===========================================================================

! number of rd values required
integer (wi),           intent(in)  :: n_rd

! maxmimum value of rd
real    (wp),           intent(in)  :: rd_max

! minmimum value of rd
real    (wp),           intent(in)  :: rd_min

! mode of the rd distribution
real    (wp),           intent(in)  :: rd_mode

! geometric standard deviation of the rd distribution,
! i.e. the width of the distribution
real    (wp),           intent(in)  :: rd_sigma

! values of B0, the dry aerosol particle hygroscopies, to use
real    (wp),           intent(in)  :: b0(:)

! pdf of the B0 values
real    (wp), optional, intent (in) :: b0_pdf(:)

!===========================================================================
! local variables for vera_cast_rd_log_normal_b0
!===========================================================================

! computed logarithmically spaced rd values
real    (wp)                        :: rd(1:n_rd)

! computed rd pdf
real    (wp)                        :: rd_pdf(1:n_rd)

! pdf of the B0 values to use in this subroutine
real    (wp)                        :: b0_pdf_use( size( b0 ) )

! pdf value to compute outside of the range
! b0 = [b0_min, b0_max], rd = [rd_min, rd_max]
! set this to zero, i.e. no chance of a particle outside of the
! specified ranges
real    (wp)                        :: boundless = vera_constants%zero

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_CAST_RD_LOG_NORMAL_B0'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! executable code for vera_cast_rd_log_normal_b0
!===========================================================================

!---------------------------------------------------------------------------
! gather the optional inputs to this subroutine

! set the b0_pdf to use
if (present(b0_pdf)) then

  ! use the b0_pdf passed to this subroutine
  b0_pdf_use = b0_pdf

else

  ! contruct a flat pdf, i.e. all hygroscopies have the same number
  ! concentration
  b0_pdf_use = vera_constants%one / size( b0, kind = wp )

end if


! compute the rd values - these are logarithmically spaced
!                         across the range [rd_min, rd_max]
!                         with mode rd_mode and width rd_sigma
call vera_log_spacing( n_rd, rd_max, rd_min, rd )

! compute the log-normal pdf of the rd values,
! with mode rd_mode and width rd_sigma
call vera_log_normal( rd, rd_mode, rd_sigma     ,                              &
                      rd_max*vera_phantom%nudge ,                              &
                      rd_min/vera_phantom%nudge ,                              &
                      boundless, rd_pdf           )

! cast an aerosol population with the computed rd, rd_pdf, b0, b0_pdf
! NOTE - this call will update the aerosol population definition
!        in the Vera global data module
call vera_cast_rd_b0( rd, b0, rd_pdf = rd_pdf, b0_pdf = b0_pdf_use )

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_cast_rd_log_normal_b0

end module vera_phantom_cast_mod
