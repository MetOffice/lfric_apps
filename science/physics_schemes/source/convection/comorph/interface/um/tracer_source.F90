! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module tracer_source_mod

implicit none

contains

! Subroutine adds on in-parcel source or sink terms for tracer species
! (increments to in-parcel values over the current level-step).
!
! The rest of comorph treats tracers as generic passively-transported
! scalars; any calculations that are specific to the actual variables
! represented by the tracers in the UM should be done in here
! (e.g. scavenging of aerosol and chemistry fields by convective precip).
!
subroutine tracer_source( n_points, n_points_super,                            &
                          massflux_d, dq_prec, par_fields, par_tracers )

use comorph_constants_mod, only: real_cvprec, real_hmprec,                     &
                                 zero, n_tracers, gravity
use fields_type_mod, only: n_fields, i_q_cl

use ukca_option_mod,     only: l_ukca, l_ukca_plume_scav
use ukca_scavenging_mod, only: ukca_plume_scav

implicit none

! Number of points
integer, intent(in) :: n_points
! Array size for super-arrays (maybe larger than n_points where array
! is reused to save on having to deallocate / reallocate to smaller size)
integer, intent(in) :: n_points_super

! Parcel dry-mass flux / kg m-2
real(kind=real_cvprec), intent(in) :: massflux_d(n_points)

! In-parcel precip mixing-ratio production increment over the current
! level-step (summed over all precip species) / kg kg-1
real(kind=real_cvprec), intent(in) :: dq_prec(n_points)

! In-parcel primary fields (latest values)
real(kind=real_cvprec), intent(in) :: par_fields                               &
                                      ( n_points_super, n_fields )

! In-parcel tracer concentrations to be updated
real(kind=real_cvprec), intent(in out) :: par_tracers                          &
                                          ( n_points_super, n_tracers )

! Inputs to ukca_plume_scav, converted to host-model precision:
!  - In-parcel tracer concentrations to be updated
real(kind=real_hmprec) :: trapkp1( n_points, n_tracers )
!  - In-parcel cloud liquid water / kg kg-1
real(kind=real_hmprec) :: xpkp1( n_points )
!  - In-parcel precip production this level-step, converted to a grid-mean
!    precip flux contribution / kg m-2 s-1
real(kind=real_hmprec) :: prekp1( n_points )
! Mass-flux converted to Pa s-1
real(kind=real_hmprec) :: flxkp1( n_points )

! Loop counters
integer :: ic, i_field


if ( l_ukca .and. l_ukca_plume_scav ) then

  ! Convert inputs to host-model precision and do required unit conversions...

  do i_field = 1, n_tracers
    do ic = 1, n_points
      ! Convert the tracer fields
      trapkp1(ic,i_field) = real( par_tracers(ic,i_field), real_hmprec )
    end do
  end do

  do ic = 1, n_points
    ! Convert the liquid cloud water mixing-ratio
    xpkp1(ic) = real( par_fields(ic,i_q_cl), real_hmprec )

    ! Convert precip mixing-ratio increment to grid-mean precip flux:
    ! Flux = mixing-ratio increment * dry-mass flux / kg m-2 s-1
    ! Also ensuring prekp1 is not negative
    prekp1(ic) = real( max(dq_prec(ic),zero) * massflux_d(ic), real_hmprec )

    ! Convert mass-flux to Pa s-1
    flxkp1(ic) = real( massflux_d(ic) * gravity, real_hmprec )

    ! Note: these unit conversions are somewhat pointless;
    ! where ukca_plume_scav uses prekp1, it just divides it by mass-flux
    ! to get the in-parcel mixing-ratio increment back again,
    ! and where it uses the mass-flux it divides by g again to get it
    ! back in kg m-2 s-1 :P
  end do

  ! Call Plume scavenging routine
  call ukca_plume_scav( n_tracers, n_points, trapkp1, xpkp1, prekp1, flxkp1 )

  ! Convert the tracer fields back to comorph's precision
  ! Note: converting these to host-model precision and back again
  ! shouldn't change the values at points where no plume scavenging
  ! occurs, but only because we expect comorph to be running at
  ! single-precision whereas host-model runs at double-precision.
  ! If it was the other way round, we'd have truncation issues.
  do i_field = 1, n_tracers
    do ic = 1, n_points
       ! Convert the tracer fields
      par_tracers(ic,i_field) = real( trapkp1(ic,i_field), real_cvprec )
    end do
  end do

end if  ! ( l_ukca .and. l_ukca_plume_scav )


return
end subroutine tracer_source

end module tracer_source_mod
