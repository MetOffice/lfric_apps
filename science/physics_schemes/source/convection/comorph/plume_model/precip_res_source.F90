! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module precip_res_source_mod

implicit none

contains

! Subroutine to add on resolved-scale source terms from exchange
! of precipitation between the parcel and the environment.
! Precip is assumed to be transfered with the same temperature
! and wind velocity as its source air, and homogenise with its
! new surroundings on arrival.
subroutine precip_res_source( n_points, n_points_src_super, n_points_res_super,&
                              l_fall_in, dt_over_rhod_lz, massflux_d,          &
                              src_air_fields, flux_cond,                       &
                              res_source_fields )

use comorph_constants_mod, only: real_cvprec, cond_params, n_cond_species
use fields_type_mod, only: i_qc_first, i_qc_last, i_temperature,               &
                           i_wind_u, i_wind_w

implicit none

! Number of points
integer, intent(in) :: n_points

! Number of points in source air fields super-array
! (maybe larger than needed here, to save having to reallocate)
integer, intent(in) :: n_points_src_super
! Number of points in resolved-scale source term fields array
integer, intent(in) :: n_points_res_super

! Flag for precip falling into parcel from environment
logical, intent(in) :: l_fall_in

! Factor dt / ( rho_dry lz ), used for converting between
! precip fall-fluxes and in-parcel mixing-ratio increments
real(kind=real_cvprec), intent(in) :: dt_over_rhod_lz(n_points)

! Dry-mass flux in kg m-2 s-1 (per m2 of surface area)
real(kind=real_cvprec), intent(in) :: massflux_d(n_points)

! Super-array containing mean-fields in the source air
! (the air that the precip has fallen FROM)
! For precip falling out of the parcel and into the environment,
! this should be the mean parcel properties.
! But for precip falling from the environment into the parcel,
! this should be the mean environment properties.
real(kind=real_cvprec), intent(in) :: src_air_fields                           &
                              ( n_points_src_super, i_wind_u:i_temperature )

! Fall-flux of each hydrometeor species (in a super-array)
real(kind=real_cvprec), intent(in) :: flux_cond                                &
                                      ( n_points, n_cond_species )

! Resolved-scale source terms for primary fields
real(kind=real_cvprec), intent(in out) :: res_source_fields                    &
                              ( n_points_res_super, i_wind_u:i_qc_last )

! Work array:
! Mass-flux * parcel mixing-ratio increment over the level-step
! = environment condensate mass source term within the
!   model-level, in kg m-2 s-1
real(kind=real_cvprec) :: dmass_cond(n_points)

! Loop counters
integer :: ic, i_cond, i_field


! For each condensed water species...
do i_cond = 1, n_cond_species

  ! Change parcel precip fall-fluxes to resolved-scale flux
  ! of condensate onto the layer
  ! (find mass of condensate transfered to the environment per s)
  do ic = 1, n_points
    ! Convert the fall-flux:
    ! Scale by delta_t / ( rho_dry lz )
    !   converts to in-parcel increment of parcel
    !   mixing-ratio over the level-step.
    ! Scale by convective flux of dry-mass
    !   converts to flux of condensate (kg m-2 s-1)
    !   into the current model-level
    dmass_cond(ic) = flux_cond(ic,i_cond)                                      &
                     * dt_over_rhod_lz(ic) * massflux_d(ic)
  end do

  ! Change sign if precip is falling from the environment into the parcel
  if ( l_fall_in ) then
    do ic = 1, n_points
      dmass_cond(ic) = -dmass_cond(ic)
    end do
  end if

  ! Find field address of current species
  i_field = i_qc_first - 1 + i_cond
  ! Add on resolved-scale source term for the precip mixing ratio
  do ic = 1, n_points
    res_source_fields(ic,i_field) = res_source_fields(ic,i_field)              &
                                  + dmass_cond(ic)
  end do

  ! Add on source of enthalpy carried from the source air by the
  ! falling precip
  ! (using water species heat capacities stored in cond_params)
  do ic = 1, n_points
    res_source_fields(ic,i_temperature) = res_source_fields(ic,i_temperature)  &
      + dmass_cond(ic) * src_air_fields(ic,i_temperature)                      &
                         * cond_params(i_cond)%pt % cp
  end do

  ! Add on source of momentum carried from the source air by the
  ! falling precip
  do i_field = i_wind_u, i_wind_w
    do ic = 1, n_points
      res_source_fields(ic,i_field) = res_source_fields(ic,i_field)            &
        + dmass_cond(ic) * src_air_fields(ic,i_field)
    end do
  end do

end do  ! i_cond = 1, n_cond_species


return
end subroutine precip_res_source

end module precip_res_source_mod
