! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module tracer_source_mod

implicit none

! Addresses of CASIM prognostic number concentrations in tracer super-array.
! These need to be set by the host-model before calling CoMorph.
integer :: i_tr_n_cl = 0
integer :: i_tr_n_rain = 0
integer :: i_tr_n_cf = 0
integer :: i_tr_n_snow = 0
integer :: i_tr_n_graup = 0

contains


! Subroutine adds on in-parcel source or sink terms for tracer species
! (increments to in-parcel values over the current level-step).
!
! The rest of comorph treats tracers as generic passively-transported
! scalars; any calculations that are specific to the actual variables
! represented by the tracers in the host-model should be done in here
! (e.g. scavenging of aerosol and chemistry fields by convective precip).
!
subroutine tracer_source( n_points, n_points_env, n_points_next, n_points_res, &
                          l_res_source, massflux_d, dq_prec,                   &
                          n_cond, flux_cond, dt_over_rhod_lz,                  &
                          env_fields, env_tracers,                             &
                          par_fields, par_tracers,                             &
                          res_source_tracers )

use comorph_constants_mod, only: real_cvprec, real_hmprec, min_float,          &
                                 zero, n_tracers, n_cond_species, gravity,     &
                                 i_cond_cl, i_cond_rain,                       &
                                 i_cond_cf, i_cond_snow, i_cond_graup
use fields_type_mod, only: n_fields, i_q_cl, i_qc_first

use ukca_option_mod,     only: l_ukca, l_ukca_plume_scav
use ukca_scavenging_mod, only: ukca_plume_scav, tracer_info
use comorph_um_namelist_mod, only: l_cv_numconcs

implicit none

! Number of points
integer, intent(in) :: n_points
! Array size for parcel super-arrays
integer, intent(in) :: n_points_next
! Array size for environment fields super-arrays
integer, intent(in) :: n_points_env
! Array size for resolved-scale source-term super-arrays
integer, intent(in) :: n_points_res

! Flag for whether to calculate resolved-scale source terms
logical, intent(in) :: l_res_source

! Parcel dry-mass flux / kg m-2
real(kind=real_cvprec), intent(in) :: massflux_d(n_points)

! In-parcel precip mixing-ratio production increment over the current
! level-step (summed over all precip species) / kg kg-1
real(kind=real_cvprec), intent(in) :: dq_prec(n_points)

! Condensate number concentrations
real(kind=real_cvprec), intent(in) :: n_cond( n_points, n_cond_species )

! Fall-flux of each hydrometeor species from parcel to environment
real(kind=real_cvprec), intent(in) :: flux_cond                                &
                                      ( n_points, n_cond_species )

! Factor dt / ( rho_dry lz ), used for converting between
! precip fall-fluxes and in-parcel mixing-ratio increments
real(kind=real_cvprec), intent(in) :: dt_over_rhod_lz(n_points)

! Environment primary fields
real(kind=real_cvprec), intent(in) :: env_fields                               &
                                      ( n_points_env, n_fields )
! Environment tracer concentrations
real(kind=real_cvprec), intent(in) :: env_tracers                              &
                                      ( n_points_env, n_tracers )

! In-parcel primary fields (latest values)
real(kind=real_cvprec), intent(in) :: par_fields                               &
                                      ( n_points_next, n_fields )
! In-parcel tracer concentrations to be updated
real(kind=real_cvprec), intent(in out) :: par_tracers                          &
                                          ( n_points_next, n_tracers )

! Resolved-scale source terms for tracers
real(kind=real_cvprec), optional, intent(in out) :: res_source_tracers         &
                                                    ( n_points_res, n_tracers )

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

! Address of each condensed water species number concentration in tracer array
integer :: i_tr_n_cond(n_cond_species)

! Loop counters
integer :: ic, i_field, i_cond


if ( l_ukca .and. l_ukca_plume_scav ) then

  ! Convert inputs to host-model precision and do required unit conversions...
  ! (only loop over UKCA tracers here)
  do i_field = tracer_info%i_ukca_first, tracer_info%i_ukca_last
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
  do i_field = tracer_info%i_ukca_first, tracer_info%i_ukca_last
    do ic = 1, n_points
       ! Convert the tracer fields
      par_tracers(ic,i_field) = real( trapkp1(ic,i_field), real_cvprec )
    end do
  end do

end if  ! ( l_ukca .AND. l_ukca_plume_scav )

if ( l_cv_numconcs ) then
  ! If running CoMorph with prognostic hydrometeor number concentrations,
  ! need to set the in-plume values, and the precip number
  ! passed to the environment by the fall-flux...

  ! Setup indices for addressing the number concentration for each condensate
  ! species in the tracer super-array
  if ( i_cond_cl > 0    )  i_tr_n_cond(i_cond_cl)    = i_tr_n_cl
  if ( i_cond_rain > 0  )  i_tr_n_cond(i_cond_rain)  = i_tr_n_rain
  if ( i_cond_cf > 0    )  i_tr_n_cond(i_cond_cf)    = i_tr_n_cf
  if ( i_cond_snow > 0  )  i_tr_n_cond(i_cond_snow)  = i_tr_n_snow
  if ( i_cond_graup > 0 )  i_tr_n_cond(i_cond_graup) = i_tr_n_graup

  ! Loop over active condensate species
  do i_cond = 1, n_cond_species
    if ( i_tr_n_cond(i_cond) > 0 ) then
      ! If this species has a number concentration in the tracer array...

      ! Set the in-plume number concentration equal to the value passed
      ! here from CoMorph's microphysics, so that we detrain that value
      ! to the grid-mean fields
      do ic = 1, n_points
        par_tracers(ic,i_tr_n_cond(i_cond)) = n_cond(ic,i_cond)
      end do

      if ( l_res_source .and. present(res_source_tracers) ) then
        ! If calculating resolved-scale source terms from precip fall...

        ! Estimate the precip fall-flux of number concentration.  This is
        ! just the fall-flux of mass but scaled by n_cond instead of q_cond
        ! (so same formula as in subroutine precip_res_source, but multiplied
        !  by n_cond and divided by q_cond).
        do ic = 1, n_points
          res_source_tracers(ic,i_tr_n_cond(i_cond))                           &
            = res_source_tracers(ic,i_tr_n_cond(i_cond))                       &
            + dt_over_rhod_lz(ic) * massflux_d(ic)                             &
            * par_tracers(ic,i_tr_n_cond(i_cond))                              &
            * ( flux_cond(ic,i_cond)                                           &
              / max( par_fields(ic,i_qc_first-1+i_cond), min_float ) )
        end do

      end if  ! ( l_res_source etc )

    end if  ! ( i_tr_n_cond(i_cond) > 0 )
  end do  ! i_cond = 1, n_cond_species

end if  ! ( l_cv_numconcs )

return
end subroutine tracer_source


! Subroutine to perform vertical homogenisation of resolved-scale source terms
! within the boundary-layer for tracers.  The homogenisation needs to be
! handled differently for different variables; this subroutine allows the
! method used to be taylored to the individual tracer fields.
subroutine tracer_homog_conv_bl( n_points_top, k_bl_top, n_fields_tot,         &
                                 cmpr, index_ic_top,                           &
                                 par_bl_top_fields, layer_mass_cmpr,           &
                                 fields_cmpr )

use comorph_constants_mod, only: real_cvprec, zero, one, sqrt_min_float,       &
                                 k_bot_conv, n_tracers
use fields_type_mod, only: n_fields
use cmpr_type_mod, only: cmpr_type

implicit none

! Number of points where convection crossed the BL-top at level k_bl_top
integer, intent(in) :: n_points_top

! Highest model-level below the BL-top
integer, intent(in) :: k_bl_top

! Number of transported fields (including tracers)
integer, intent(in) :: n_fields_tot

! Compression indices of points on each level where there are
! resolved-scale source terms in columns where convection hit
! the BL-top at k_bl_top
type(cmpr_type), intent(in) :: cmpr(k_bot_conv:k_bl_top)

! Indices of those points in the par_bl_top compression list
integer, intent(in) :: index_ic_top ( n_points_top, k_bot_conv:k_bl_top )

! Parcel mean properties at the first level above the boundary-layer
real(kind=real_cvprec), intent(in) :: par_bl_top_fields                        &
                                      ( n_points_top, n_fields_tot )

! Fraction of the BL-top mass-flux to be entrained from each level
real(kind=real_cvprec), intent(in) :: layer_mass_cmpr                          &
                                      ( n_points_top, k_bot_conv:k_bl_top )

! IN:  Compressed fields from each model-level
! OUT: The values entrained from each model-level
real(kind=real_cvprec), intent(in out) :: fields_cmpr                          &
                                          ( n_points_top,                      &
                                            n_fields_tot, k_bot_conv:k_bl_top )

! After this routine, the output entrained tracer values are scaled by the
! mass entrained from each level to obtain the vertically-homogenised
! resolved-scale source-terms.  This routine may modify the resulting
! source-terms by modifying the entrained values below...

! Vertical means of current tracer field over the boundary-layer
real(kind=real_cvprec) :: par_bl_mean_tracer(n_points_top)

! Perturbation factor
real(kind=real_cvprec) :: fac

! Loop counters
integer :: ic, ic2, k, i_field, i_tracer


! Loop over tracers
do i_tracer = 1, n_tracers

  ! Choose course of action based on the specific tracer variable...
  if ( i_tracer==i_tr_n_cl .or. i_tracer==i_tr_n_rain .or.                     &
       i_tracer==i_tr_n_cf .or. i_tracer==i_tr_n_snow .or.                     &
       i_tracer==i_tr_n_graup ) then

    ! For CASIM number concentrations, do nothing!
    ! (leave the entrained tracer values equal to grid-mean values as input).
    !
    ! The homogenisation code sets the entrained hydrometeor masses to the
    ! grid-mean values (and adjusts the entrained vapour and temperature to
    ! ensure conservation of water and moist enthalpy).  So we should
    ! entrain grid-mean values for the hydrometeor number concentrations
    ! too for consistency.  Since number is not conserved (any condensation
    ! or evaporation can create or destroy number), we don't need to
    ! calculate the entrained number such that it integrates to the value
    ! in the parcel at the BL-top.  Indeed, doing-so would be wrong if
    ! hydrometeors have been created in the updraft before it reached the
    ! BL-top...
    continue

  else  ! ( i_tracer )
    ! Assume other tracer species are conserved, so calculate entrained
    ! values such that the vertical integral of the source-terms within
    ! the BL equals the value in the parcel at the BL-top.

    ! Find fields super-array index
    i_field = n_fields + i_tracer

    ! Initialise vertical integral to zero
    do ic = 1, n_points_top
      par_bl_mean_tracer(ic) = zero
    end do

    do k = k_bot_conv, k_bl_top
      if ( cmpr(k) % n_points > 0 ) then

        ! Add up vertical integral of tracer
        do ic2 = 1, cmpr(k) % n_points
          ic = index_ic_top(ic2,k)
          par_bl_mean_tracer(ic) = par_bl_mean_tracer(ic)                      &
            + fields_cmpr(ic2,i_field,k) * layer_mass_cmpr(ic2,k)
        end do

      end if  ! ( cmpr(k) % n_points > 0 )
    end do  ! k = k_bot_conv, k_bl_top

    do k = k_bot_conv, k_bl_top
      if ( cmpr(k) % n_points > 0 ) then

        ! Add perturbations to entrained tracer so-as to
        ! scale the vertical integral of entrained
        ! tracer to the value in the parcel at BL-top
        do ic2 = 1, cmpr(k) % n_points
          ic = index_ic_top(ic2,k)

          ! Calc ratio of BL-top parcel value over mean value in
          ! the BL, with safety-check to avoid div-by-zero:
          if ( abs(par_bl_mean_tracer(ic)) > sqrt_min_float ) then
            fac = par_bl_top_fields(ic,i_field) / par_bl_mean_tracer(ic)
          else
            fac = one
          end if

          ! Apply correction factor
          fields_cmpr(ic2,i_field,k) = fields_cmpr(ic2,i_field,k) * fac

        end do

      end if  ! ( cmpr(k) % n_points > 0 )
    end do  ! k = k_bot_conv, k_bl_top

  end if  ! ( i_tracer )

end do  ! i_tracer = 1, n_tracers


return
end subroutine tracer_homog_conv_bl


! Subroutine to set parcel initial values of tracer fields
subroutine tracer_gen( n_points, n_points_super, n_fields_tot,                 &
                       fields_k, par_gen_mean, par_gen_core )

use comorph_constants_mod, only: real_cvprec, min_float,                       &
                                 n_tracers, n_cond_species, l_par_core,        &
                                 i_cond_cl, i_cond_rain,                       &
                                 i_cond_cf, i_cond_snow, i_cond_graup
use comorph_um_namelist_mod, only: l_cv_numconcs
use fields_type_mod, only: n_fields, i_qc_first

implicit none

! Number of points
integer, intent(in) :: n_points

! Size of input super-arrays (maybe larger than n_points,
! as arrays are dimensioned with the max size they could need
! on any level and reused for all levels)
integer, intent(in) :: n_points_super

! Total number of fields in the input super-arrays,
! including tracers if they are used
integer, intent(in) :: n_fields_tot

! Primary model-fields at level k
real(kind=real_cvprec), intent(in) :: fields_k                                 &
                                      ( n_points_super, n_fields_tot )

! In-parcel means of primary fields
real(kind=real_cvprec), intent(in out) :: par_gen_mean                         &
                                          ( n_points, n_fields_tot )
! Parcel core primary fields
real(kind=real_cvprec), intent(in out) :: par_gen_core                         &
                                          ( n_points, n_fields_tot )

! Address of each condensed water species number concentration in tracer array
integer :: i_tr_n_cond(n_cond_species)

! Loop counters
integer :: ic, i_field, i_cond


! For now, copy grid-mean tracer values into the parcel
do i_field = n_fields+1, n_fields_tot
  do ic = 1, n_points
    par_gen_mean(ic,i_field) = fields_k(ic,i_field)
  end do
end do

if ( l_cv_numconcs ) then
  ! If running CoMorph with prognostic hydrometeor number concentrations,
  ! need to set parcel initial number cosistent with mass where the
  ! masses aren't equal to the grid-mean...

  ! Setup indices for addressing the number concentration for each condensate
  ! species in the tracer super-array
  if ( i_cond_cl > 0    )  i_tr_n_cond(i_cond_cl)    = i_tr_n_cl
  if ( i_cond_rain > 0  )  i_tr_n_cond(i_cond_rain)  = i_tr_n_rain
  if ( i_cond_cf > 0    )  i_tr_n_cond(i_cond_cf)    = i_tr_n_cf
  if ( i_cond_snow > 0  )  i_tr_n_cond(i_cond_snow)  = i_tr_n_snow
  if ( i_cond_graup > 0 )  i_tr_n_cond(i_cond_graup) = i_tr_n_graup

  ! Loop over active condensate species
  do i_cond = 1, n_cond_species
    if ( i_tr_n_cond(i_cond) > 0 ) then
      ! If this species has a number concentration in the tracer array...
      ! Set fields array address of the number concentration (in tracer sect)
      i_field = n_fields + i_tr_n_cond(i_cond)
      ! Scale parcel number conc up by ratio of parcel initial condensate mass
      ! over grid-mean condensate mass.  This should account for when the
      ! parcel has triggered from within cloud and so should pick-up
      ! in-cloud hydrometeor number consistent with in-cloud mass.
      do ic = 1, n_points
        par_gen_mean(ic,i_field) = par_gen_mean(ic,i_field)                    &
           * ( par_gen_mean(ic,i_qc_first-1+i_cond)                            &
             / max( fields_k(ic,i_qc_first-1+i_cond), min_float ) )
      end do
    end if
  end do

end if  ! ( l_cv_numconcs )

if ( l_par_core ) then
  ! Set parcel core tracers equal to parcel mean tracers
  do i_field = n_fields+1, n_fields_tot
    do ic = 1, n_points
      par_gen_core(ic,i_field) = par_gen_mean(ic,i_field)
    end do
  end do
end if


return
end subroutine tracer_gen


end module tracer_source_mod
