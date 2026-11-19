! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

#if !defined(LFRIC)
module tracer_source_mod

implicit none

contains

! Subroutine adds on in-parcel source or sink terms for tracer species
! (increments to in-parcel values over the current level-step).
!
! The rest of comorph treats tracers as generic passively-transported
! scalars; any calculations that are specific to the actual variables
! represented by the tracers in the stand-alone test should be done in here
!
subroutine tracer_source( n_points, n_points_env, n_points_next, n_points_res, &
                          l_res_source, massflux_d, dq_prec,                   &
                          n_cond, flux_cond, dt_over_rhod_lz,                  &
                          env_fields, env_tracers,                             &
                          par_fields, par_tracers,                             &
                          res_source_tracers )

use comorph_constants_mod, only: real_cvprec, n_tracers, n_cond_species
use fields_type_mod, only: n_fields

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


! So far this is just a dummy routine that does nothing, since no
! source or sink terms are yet applied to tracers in the standalone test.


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

  ! Assume tracer species are conserved, so calculate entrained
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
        par_bl_mean_tracer(ic) = par_bl_mean_tracer(ic)                        &
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

end do  ! i_tracer = 1, n_tracers


return
end subroutine tracer_homog_conv_bl


! Subroutine to set parcel initial values of tracer fields
subroutine tracer_gen( n_points, n_points_super, n_fields_tot,                 &
                       fields_k, par_gen_mean, par_gen_core )

use comorph_constants_mod, only: real_cvprec, n_tracers, l_par_core
use fields_type_mod, only: n_fields

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

! Loop counters
integer :: ic, i_field

! For now, copy grid-mean tracer values into the parcel
do i_field = n_fields+1, n_fields_tot
  do ic = 1, n_points
    par_gen_mean(ic,i_field) = fields_k(ic,i_field)
  end do
end do

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
#endif
