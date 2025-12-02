! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module homog_conv_bl_mod

implicit none

contains


! Subroutine to homogenize the convective
! resolved-scale source terms within the boundary-layer.
!
! The method used is to calculate a uniform profile of
! mass entrainment from each level below the BL-top such
! that the vertical integral of mass entrainment equals
! the parcel mass-flux at the BL-top.  Entrained air
! properties then also need to be worked out on each
! level below the BL-top such that their vertical integral
! equals the properties of the parcel at the BL-top.
! To do this while accurately adhering to conservation
! is fiddly; first the dry-mass and total-water entrained
! from each level need to be worked out; once these are
! known, the accurately conserving budgets for other fields
! can be calculated.
!
! Note that, when this routine is used for downdrafts which
! descend into the boundary-layer instead of updrafts rising
! out of it, the resulting resolved-scale sources correspond
! to a uniform profile of detrainment within the BL, instead
! of entrainment.
!
subroutine homog_conv_bl( n_conv_types, n_conv_layers,                         &
                          n_fields_tot, l_down,                                &
                          i_type, i_layr, k_bl_top,                            &
                          ij_first, ij_last, index_ic,                         &
                          grid, fields_np1, layer_mass,                        &
                          par_bl_top, res_source )

use comorph_constants_mod, only: real_cvprec, real_hmprec, zero, one,          &
                     nx_full, ny_full, k_bot_conv, k_top_conv,                 &
                     n_cond_species, n_cond_species_liq,                       &
                     L_con_0, L_sub_0, l_cv_cloudfrac, sqrt_min_float
use cmpr_type_mod, only: cmpr_type, cmpr_alloc
use grid_type_mod, only: grid_type
use fields_type_mod, only: fields_type, n_fields,                              &
                           i_qc_first, i_qc_last, i_q_vap,                     &
                           i_temperature, i_wind_u, i_wind_w,                  &
                           i_cf_liq, i_cf_bulk
use parcel_type_mod, only: parcel_type, i_massflux_d
use res_source_type_mod, only: res_source_type, i_ent, i_det
use compress_mod, only: compress
use calc_q_tot_mod, only: calc_q_tot
use set_cp_tot_mod, only: set_cp_tot
use dry_adiabat_mod, only: dry_adiabat
use calc_virt_temp_dry_mod, only: calc_virt_temp_dry


implicit none

! Number of convection types
integer, intent(in) :: n_conv_types

! Number of distinct convecting layers
integer, intent(in) :: n_conv_layers

! Number of transported fields (including tracers)
integer, intent(in) :: n_fields_tot

! Flag for downdraft versus updraft
logical, intent(in) :: l_down

! Indices of current convection type and layer
integer, intent(in) :: i_type
integer, intent(in) :: i_layr

! Highest model-level below the BL-top
integer, intent(in) :: k_bl_top

! ij indices of first and last point where convection hit BL-top
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last

! Array storing compression list indices of the resolved-scale
! source term arrays on a grid, to enable addressing res_source
! arrays from the par_bl_top compression list
integer, intent(in) :: index_ic                                                &
                       ( ij_first:ij_last, k_bot_conv:k_bl_top )

! Structure containing model grid-fields
type(grid_type), intent(in) :: grid

! Structure containing full 3-D primary fields
type(fields_type), intent(in) :: fields_np1

! Full 3-D array of dry-mass per unit surface area
! contained in each grid-cell
real(kind=real_hmprec), intent(in) :: layer_mass                               &
       ( nx_full, ny_full, k_bot_conv:k_top_conv )

! Structures containing parcel properties at the
! first level above the boundary-layer top, only for convection
! of the current type and layer which passed the BL-top at
! level k_bl_top
type(parcel_type), intent(in) :: par_bl_top

! Resolved-scale source terms to be modified
type(res_source_type), intent(in out) :: res_source                            &
          ( n_conv_types, n_conv_layers, k_bot_conv:k_top_conv )


! Total-water mixing ratio and liquid+ice-water enthalpy
! of the parcel at the BL-top
real(kind=real_cvprec) :: par_bl_top_q_tot                                     &
                          ( par_bl_top % cmpr % n_points )
real(kind=real_cvprec) :: par_bl_top_templ                                     &
                          ( par_bl_top % cmpr % n_points )

! Total-water mixing ratio and liquid+ice-water enthalpy
! of vertical-mean environment below the BL-top
real(kind=real_cvprec) :: par_bl_mean_q_tot                                    &
                          ( par_bl_top % cmpr % n_points )
real(kind=real_cvprec) :: par_bl_mean_templ                                    &
                          ( par_bl_top % cmpr % n_points )

! Compressed layer-masses
real(kind=real_cvprec) :: layer_mass_k                                         &
                          ( par_bl_top % cmpr % n_points )

! Compressed primary fields from one model-level
real(kind=real_cvprec) :: fields_k                                             &
             ( par_bl_top % cmpr % n_points, n_fields_tot )

! Vertical integral of layer-masses below the BL-top
real(kind=real_cvprec) :: layer_mass_bl                                        &
                          ( par_bl_top % cmpr % n_points )

! Vertical means of fields over the boundary-layer
real(kind=real_cvprec) :: par_bl_mean_fields                                   &
             ( par_bl_top % cmpr % n_points, n_fields_tot )

! Compressed environment total condensed water mixing-ratio
real(kind=real_cvprec) :: q_cond_tot                                           &
                          ( par_bl_top % cmpr % n_points )
! Compressed environment latent enthalpy term from condensate
real(kind=real_cvprec) :: lat_heat_term                                        &
                          ( par_bl_top % cmpr % n_points )

! Compressed pressure from current level and BL-top
real(kind=real_cvprec) :: pressure_k                                           &
                          ( par_bl_top % cmpr % n_points )
real(kind=real_cvprec) :: par_bl_top_pressure                                  &
                          ( par_bl_top % cmpr % n_points )

! Total heat capacity for conversion from temperature to enthalpy
real(kind=real_cvprec) :: cp_tot                                               &
                          ( par_bl_top % cmpr % n_points )
! Exner ratio for pressure change up to BL-top
real(kind=real_cvprec) :: exner_ratio                                          &
                          ( par_bl_top % cmpr % n_points )
! Total-water mixing ratio for conversion to momentum / dry-mass
real(kind=real_cvprec) :: q_tot                                                &
                          ( par_bl_top % cmpr % n_points )
! Dry virtual temperature T ( 1 + Rv/Rd q_vap )
real(kind=real_cvprec) :: virt_temp_dry                                        &
                          ( par_bl_top % cmpr % n_points )

! Latent heat term
real(kind=real_cvprec) :: Lc

! Perturbation factor
real(kind=real_cvprec) :: fac

! Compression indices of points on each level where there are
! resolved-scale source terms in columns where convection hit
! the BL-top at k_bl_top
type(cmpr_type), allocatable :: cmpr(:)
! Indices of those points in the par_bl_top and res-source
! compression lists
integer :: index_ic_top                                                        &
           ( par_bl_top % cmpr % n_points, k_bot_conv:k_bl_top )
integer :: index_ic_res                                                        &
           ( par_bl_top % cmpr % n_points, k_bot_conv:k_bl_top )

! Lower and upper bounds of full 3-D arrays
integer :: lb(3), ub(3)

! Loop counters
integer :: i, j, k, ic, ic2, ic3, i_cond, i_field



!----------------------------------------------------------------
! 1) Find points on each level where there are resolved-scale
!    source terms in columns where the parcel passed the BL-top
!    at model-level k_bl_top
!----------------------------------------------------------------

allocate( cmpr(k_bot_conv:k_bl_top) )

do k = k_bot_conv, k_bl_top

  ! Initialise number of points to zero
  cmpr(k) % n_points = 0

  do ic = 1, par_bl_top % cmpr % n_points
    ! Extract i,j coordinates of points which hit the BL-top
    i = par_bl_top % cmpr % index_i(ic)
    j = par_bl_top % cmpr % index_j(ic)
    ! Extract res-source compression index for those coords
    ic2 = index_ic( nx_full*(j-1)+i, k )
    if ( ic2 > 0 ) then
      ! If this is a res-source point, add it to the list.
      cmpr(k) % n_points = cmpr(k) % n_points + 1
      index_ic_res( cmpr(k) % n_points, k ) = ic2
      index_ic_top( cmpr(k) % n_points, k ) = ic
    end if
  end do

  ! If found any res_source points
  if ( cmpr(k) % n_points > 0 ) then
    ! Save their i,j, indices as well
    call cmpr_alloc( cmpr(k), cmpr(k) % n_points )
    do ic2 = 1, cmpr(k) % n_points
      ic = index_ic_top(ic2,k)
      cmpr(k) % index_i(ic2) = par_bl_top % cmpr % index_i(ic)
      cmpr(k) % index_j(ic2) = par_bl_top % cmpr % index_j(ic)
    end do
  end if

end do


!----------------------------------------------------------------
! 2) Calculate vertical mean environment dry-mass and water over
!    res_source points below the BL-top
!----------------------------------------------------------------

! Initialise vertical means over the BL to zero
do ic = 1, par_bl_top % cmpr % n_points
  layer_mass_bl(ic) = zero
end do
do i_field = i_q_vap, i_qc_last
  do ic = 1, par_bl_top % cmpr % n_points
    par_bl_mean_fields(ic,i_field) = zero
  end do
end do

do k = k_bot_conv, k_bl_top
  if ( cmpr(k) % n_points > 0 ) then

    ! Compress the layer-masses onto current points
    lb = lbound(layer_mass)
    ub = ubound(layer_mass)
    call compress( cmpr(k), lb(1:2), ub(1:2), layer_mass(:,:,k), layer_mass_k )

    ! Store layer-mass in the res_source mass array
    do ic2 = 1, cmpr(k) % n_points
      ic = index_ic_res(ic2,k)
      res_source(i_type,i_layr,k) % res_super(ic,i_ent)                        &
        = layer_mass_k(ic2)
    end do

    ! Add up vertical integral of layer-mass over the BL
    do ic2 = 1, cmpr(k) % n_points
      ic = index_ic_top(ic2,k)
      layer_mass_bl(ic) = layer_mass_bl(ic) + layer_mass_k(ic2)
    end do

    ! For each water-species...
    do i_field = i_q_vap, i_qc_last

      ! Compress the environment water onto current points
      lb = lbound(fields_np1%list(i_field)%pt)
      ub = ubound(fields_np1%list(i_field)%pt)
      call compress( cmpr(k), lb(1:2), ub(1:2),                                &
                     fields_np1%list(i_field)%pt(:,:,k), fields_k(:,i_field) )
      ! Remove spurious negative values from input data
      do ic2 = 1, cmpr(k) % n_points
        fields_k(ic2,i_field) = max( fields_k(ic2,i_field),                    &
                                     zero )
      end do

      ! We now overwrite the resolved-scale source terms
      ! with new entrained air properties and mass based
      ! on the mean environment properties
      do ic2 = 1, cmpr(k) % n_points
        ic = index_ic_res(ic2,k)
        res_source(i_type,i_layr,k) % fields_super(ic,i_field)                 &
          = fields_k(ic2,i_field)
      end do

      ! Add up vertical integral of water over the BL
      do ic2 = 1, cmpr(k) % n_points
        ic = index_ic_top(ic2,k)
        par_bl_mean_fields(ic,i_field)                                         &
          = par_bl_mean_fields(ic,i_field)                                     &
          + fields_k(ic2,i_field) * layer_mass_k(ic2)
      end do

    end do  ! i_field = i_q_vap, i_qc_last

  end if  ! ( cmpr(k) % n_points > 0 )
end do  ! k = k_bot_conv, k_bl_top

! Normalise the vertical means of water-masses
do i_field = i_q_vap, i_qc_last
  do ic = 1, par_bl_top % cmpr % n_points
    par_bl_mean_fields(ic,i_field)                                             &
      = par_bl_mean_fields(ic,i_field)                                         &
      / layer_mass_bl(ic)
  end do
end do


!----------------------------------------------------------------
! 3) Calculate total-water mixing-ratio of the parcel at the
!    BL-top and the mean environment below BL-top
!----------------------------------------------------------------

call calc_q_tot( par_bl_top % cmpr % n_points,                                 &
                 par_bl_top % cmpr % n_points,                                 &
                 par_bl_top % mean_super(:,i_q_vap),                           &
                 par_bl_top % mean_super(:,i_qc_first:i_qc_last),              &
                 par_bl_top_q_tot )
call calc_q_tot( par_bl_top % cmpr % n_points,                                 &
                 par_bl_top % cmpr % n_points,                                 &
                 par_bl_mean_fields(:,i_q_vap),                                &
                 par_bl_mean_fields(:,i_qc_first:i_qc_last),                   &
                 par_bl_mean_q_tot )


!----------------------------------------------------------------
! 4) Work out dry-mass and water-vapour to entrain from each
!    level below BL-top such that they integrate to the parcel
!    mass-flux and total-water at the BL-top
!----------------------------------------------------------------

do k = k_bot_conv, k_bl_top
  if ( cmpr(k) % n_points > 0 ) then

    ! Scale mass entrained from each layer such that
    ! its vertical integral equals the mass-flux of
    ! the parcel at the BL-top
    do ic2 = 1, cmpr(k) % n_points
      ic = index_ic_res(ic2,k)
      ic3 = index_ic_top(ic2,k)
      res_source(i_type,i_layr,k) % res_super(ic,i_ent)                        &
        = par_bl_top % par_super(ic3,i_massflux_d)                             &
          * ( res_source(i_type,i_layr,k) % res_super(ic,i_ent)                &
            / layer_mass_bl(ic3) )
    end do

    ! Calculate environment total condensed water using the
    ! condensate mixing-ratios we stashed in res_source earlier.
    do ic2 = 1, cmpr(k) % n_points
      q_cond_tot(ic2) = zero
    end do
    do i_field = i_qc_first, i_qc_last
      do ic2 = 1, cmpr(k) % n_points
        ic = index_ic_res(ic2,k)
        q_cond_tot(ic2) = q_cond_tot(ic2)                                      &
         + res_source(i_type,i_layr,k) % fields_super(ic,i_field)
      end do
    end do

    ! Add perturbations to entrained q_vap so-as to
    ! scale the vertical integral of entrained
    ! total-water to the value in the parcel at BL-top
    do ic2 = 1, cmpr(k) % n_points
      ic = index_ic_res(ic2,k)
      ic3 = index_ic_top(ic2,k)
      if ( par_bl_mean_q_tot(ic3) > zero ) then
        ! Avoid div-by-zero in the unlikely event that the
        ! BL water content has fallen to zero

        ! Fractional increase of total water required
        fac = par_bl_top_q_tot(ic3) / par_bl_mean_q_tot(ic3)

        ! Modification of vapour; we want:
        ! q_vap_new + q_cond = fac * ( q_vap_old + q_cond )
        ! Rearranging:
        ! q_vap_new = q_vap_old * fac + q_cond * (fac-1)
        res_source(i_type,i_layr,k) % fields_super(ic,i_q_vap)                 &
        = res_source(i_type,i_layr,k) % fields_super(ic,i_q_vap)               &
          * fac                                                                &
        + q_cond_tot(ic2) * ( fac - one )

      else

        res_source(i_type,i_layr,k) % fields_super(ic,i_q_vap)                 &
          = par_bl_top_q_tot(ic3)  ! Probably zero

      end if
    end do

  end if  ! ( cmpr(k) % n_points > 0 )
end do  ! k = k_bot_conv, k_bl_top


!----------------------------------------------------------------
! 5) Calculate vertical mean enthalpy per unit dry-mass over
!    res_source points below the BL-top
!----------------------------------------------------------------

! Initialise vertical integral of enthalpy
do ic = 1, par_bl_top % cmpr % n_points
  par_bl_mean_fields(ic,i_temperature) = zero
end do

do k = k_bot_conv, k_bl_top
  if ( cmpr(k) % n_points > 0 ) then

    ! Compress the environment temperature
    lb = lbound(fields_np1%list(i_temperature)%pt)
    ub = ubound(fields_np1%list(i_temperature)%pt)
    call compress( cmpr(k), lb(1:2), ub(1:2),                                  &
                   fields_np1%list(i_temperature)%pt(:,:,k),                   &
                   fields_k(:,i_temperature) )

    ! Extract the entrained water-species
    ! back from the res_source arrays
    do i_field = i_q_vap, i_qc_last
      do ic2 = 1, cmpr(k) % n_points
        ic = index_ic_res(ic2,k)
        fields_k(ic2,i_field)                                                  &
         = res_source(i_type,i_layr,k) % fields_super(ic,i_field)
      end do
    end do

    ! Convert environment temperature to enthalpy cp_tot * T
    ! (note we are using the already updated q_vap to set cp)
    call set_cp_tot( cmpr(k)%n_points,                                         &
                     par_bl_top % cmpr % n_points,                             &
                     fields_k(:,i_q_vap),                                      &
                     fields_k(:,i_qc_first:i_qc_last),                         &
                     cp_tot )
    do ic2 = 1, cmpr(k) % n_points
      fields_k(ic2,i_temperature) = fields_k(ic2,i_temperature)                &
                                  * cp_tot(ic2)
    end do

    ! Extract pressure from current level and BL-top
    lb = lbound(grid % pressure_full)
    ub = ubound(grid % pressure_full)
    call compress( cmpr(k), lb(1:2), ub(1:2),                                  &
                   grid % pressure_full(:,:,k), pressure_k )
    ! Note: k_bl_top is the index of the BL-top full-level; the BL-top
    ! parcel is defined at the next half-level up from this, so need
    ! to index pressure_half using k_bl_top+1
    lb = lbound(grid % pressure_half)
    ub = ubound(grid % pressure_half)
    call compress( cmpr(k), lb(1:2), ub(1:2),                                  &
                   grid % pressure_half(:,:,k_bl_top+1), par_bl_top_pressure )

    ! Calculate enthalpy of air if lifted to the BL-top
    call dry_adiabat( cmpr(k)%n_points,                                        &
                      par_bl_top % cmpr % n_points,                            &
                      pressure_k, par_bl_top_pressure,                         &
                      fields_k(:,i_q_vap),                                     &
                      fields_k(:,i_qc_first:i_qc_last),                        &
                      exner_ratio )
    do ic2 = 1, cmpr(k) % n_points
      fields_k(ic2,i_temperature) = fields_k(ic2,i_temperature)                &
                                  * exner_ratio(ic2)
    end do

    ! Copy this enthalpy into the resolved-scale source term
    ! arrays
    do ic2 = 1, cmpr(k) % n_points
      ic = index_ic_res(ic2,k)
      res_source(i_type,i_layr,k) %fields_super(ic,i_temperature)              &
          = fields_k(ic2,i_temperature)
    end do

    ! Integrate up vertical mean enthalpy
    do ic2 = 1, cmpr(k) % n_points
      ic = index_ic_res(ic2,k)
      ic3 = index_ic_top(ic2,k)
      par_bl_mean_fields(ic3,i_temperature)                                    &
        = par_bl_mean_fields(ic3,i_temperature)                                &
        + fields_k(ic2,i_temperature)                                          &
          * res_source(i_type,i_layr,k) % res_super(ic,i_ent)
    end do

  end if  ! ( cmpr(k) % n_points > 0 )
end do  ! k = k_bot_conv, k_bl_top

! Normalise the vertical mean enthalpy
do ic = 1, par_bl_top % cmpr % n_points
  par_bl_mean_fields(ic,i_temperature)                                         &
    = par_bl_mean_fields(ic,i_temperature)                                     &
    / par_bl_top % par_super(ic,i_massflux_d)
end do


!----------------------------------------------------------------
! 6) Calculate liquid+ice-water enthalpy of the parcel at the
!    BL-top and the mean environment below BL-top
!----------------------------------------------------------------

! Initialise to zero
do ic = 1, par_bl_top % cmpr % n_points
  par_bl_top_templ(ic) = zero
  par_bl_mean_templ(ic) = zero
end do

! Add latent heat contributions from condensed water species
do i_cond = 1, n_cond_species
  i_field = i_qc_first-1+i_cond
  ! Select the right latent heat, depending on water phase
  if ( i_cond <= n_cond_species_liq ) then
    Lc = L_con_0
  else
    Lc = L_sub_0
  end if
  do ic = 1, par_bl_top % cmpr % n_points
    ! Add latent enthalpy terms -Lc0 qc to cp_tot*T
    par_bl_top_templ(ic) = par_bl_top_templ(ic)                                &
                         - par_bl_top % mean_super(ic,i_field)*Lc
    par_bl_mean_templ(ic) = par_bl_mean_templ(ic)                              &
                          - par_bl_mean_fields(ic,i_field)*Lc
  end do
end do  ! i_cond = 1, n_cond_species

! Add on the rest of the enthalpy
do ic = 1, par_bl_top % cmpr % n_points
  par_bl_top_templ(ic) = par_bl_top_templ(ic)                                  &
                   + par_bl_top % mean_super(ic,i_temperature)
  par_bl_mean_templ(ic) = par_bl_mean_templ(ic)                                &
                   + par_bl_mean_fields(ic,i_temperature)
end do


!----------------------------------------------------------------
! 7) Work out enthalpy to entrain from each
!    level below BL-top such that it integrates to the parcel
!    enthalpy at the BL-top
!----------------------------------------------------------------

do k = k_bot_conv, k_bl_top
  if ( cmpr(k) % n_points > 0 ) then

    ! Extract the entrained enthalpy and water-species
    ! back from the res_source arrays
    do i_field = i_temperature, i_qc_last
      do ic2 = 1, cmpr(k) % n_points
        ic = index_ic_res(ic2,k)
        fields_k(ic2,i_field)                                                  &
         = res_source(i_type,i_layr,k) % fields_super(ic,i_field)
      end do
    end do

    ! Calculate environment latent enthalpy term due to
    ! condensed water
    do ic2 = 1, cmpr(k) % n_points
      lat_heat_term(ic2) = zero
    end do
    do i_cond = 1, n_cond_species
      i_field = i_qc_first-1+i_cond
      ! Select the right latent heat, depending on water phase
      if ( i_cond <= n_cond_species_liq ) then
        Lc = L_con_0
      else
        Lc = L_sub_0
      end if
      do ic2 = 1, cmpr(k) % n_points
        lat_heat_term(ic2) = lat_heat_term(ic2)                                &
         - fields_k(ic2,i_field) * Lc
      end do
    end do

    ! Add perturbations to entrained enthalpy so-as to
    ! scale the vertical integral of entrained
    ! liquid+ice-water enthalpy to the value in the
    ! parcel at BL-top
    do ic2 = 1, cmpr(k) % n_points
      ic = index_ic_top(ic2,k)

      ! Fractional increase of liquid+ice-water enthalpy required
      fac = par_bl_top_templ(ic) / par_bl_mean_templ(ic)

      ! Modification of enthalpy; we want:
      ! enth_new + lat_heat = fac * ( enth_old + lat_heat )
      ! Rearranging:
      ! enth_new = enth_old * fac + lat_heat * (fac-1)
      fields_k(ic2,i_temperature)                                              &
        = fields_k(ic2,i_temperature) * fac                                    &
        + lat_heat_term(ic2) * ( fac - one )
    end do

    ! Extract pressure from current level and BL-top
    lb = lbound(grid % pressure_full)
    ub = ubound(grid % pressure_full)
    call compress( cmpr(k), lb(1:2), ub(1:2),                                  &
                   grid % pressure_full(:,:,k), pressure_k )
    ! Note: k_bl_top is the index of the BL-top full-level; the BL-top
    ! parcel is defined at the next half-level up from this, so need
    ! to index pressure_half using k_bl_top+1
    lb = lbound(grid % pressure_half)
    ub = ubound(grid % pressure_half)
    call compress( cmpr(k), lb(1:2), ub(1:2),                                  &
                   grid % pressure_half(:,:,k_bl_top+1), par_bl_top_pressure )

    ! Calculate enthalpy of air after subsiding back from the
    ! BL-top down to the current level
    call dry_adiabat( cmpr(k)%n_points,                                        &
                      par_bl_top % cmpr % n_points,                            &
                      pressure_k, par_bl_top_pressure,                         &
                      fields_k(:,i_q_vap),                                     &
                      fields_k(:,i_qc_first:i_qc_last),                        &
                      exner_ratio )
    do ic2 = 1, cmpr(k) % n_points
      fields_k(ic2,i_temperature) = fields_k(ic2,i_temperature)                &
                                  / exner_ratio(ic2)
    end do

    ! Copy final enthalpy into the resolved-scale source term
    ! arrays
    do ic2 = 1, cmpr(k) % n_points
      ic = index_ic_res(ic2,k)
      res_source(i_type,i_layr,k) %fields_super(ic,i_temperature)              &
          = fields_k(ic2,i_temperature)
    end do

  end if  ! ( cmpr(k) % n_points > 0 )
end do  ! k = k_bot_conv, k_bl_top


!----------------------------------------------------------------
! 8) Calculate vertical mean momentum per unit dry-mass over
!    res_source points below the BL-top
!----------------------------------------------------------------

! Initialise vertical integrals
do i_field = i_wind_u, i_wind_w
  do ic = 1, par_bl_top % cmpr % n_points
    par_bl_mean_fields(ic,i_field) = zero
  end do
end do

do k = k_bot_conv, k_bl_top
  if ( cmpr(k) % n_points > 0 ) then

    ! Extract the entrained water-species
    ! back from the res_source arrays
    do i_field = i_q_vap, i_qc_last
      do ic2 = 1, cmpr(k) % n_points
        ic = index_ic_res(ic2,k)
        fields_k(ic2,i_field)                                                  &
         = res_source(i_type,i_layr,k) % fields_super(ic,i_field)
      end do
    end do

    ! Calculate total-water mixing-ratio
    call calc_q_tot( cmpr(k) % n_points,                                       &
                     par_bl_top % cmpr % n_points,                             &
                     fields_k(:,i_q_vap),                                      &
                     fields_k(:,i_qc_first:i_qc_last),                         &
                     q_tot )

    ! For each wind component:
    do i_field = i_wind_u, i_wind_w

      ! Compress the environment winds onto current points
      lb = lbound(fields_np1%list(i_field)%pt)
      ub = ubound(fields_np1%list(i_field)%pt)
      call compress( cmpr(k), lb(1:2), ub(1:2),                                &
                     fields_np1%list(i_field)%pt(:,:,k), fields_k(:,i_field) )

      ! Convert environment winds to momentum per unit dry-mass
      do ic2 = 1, cmpr(k) % n_points
        fields_k(ic2,i_field) = fields_k(ic2,i_field)                          &
          * ( one + q_tot(ic2) )
      end do

      ! Store environment momentum per unit dry-mass in the
      ! res_source arrays
      do ic2 = 1, cmpr(k) % n_points
        ic = index_ic_res(ic2,k)
        res_source(i_type,i_layr,k) % fields_super(ic,i_field)                 &
          = fields_k(ic2,i_field)
      end do

      ! Add up vertical integral of momentum
      do ic2 = 1, cmpr(k) % n_points
        ic = index_ic_res(ic2,k)
        ic3 = index_ic_top(ic2,k)
        par_bl_mean_fields(ic3,i_field)                                        &
          = par_bl_mean_fields(ic3,i_field)                                    &
          + fields_k(ic2,i_field)                                              &
            * res_source(i_type,i_layr,k) % res_super(ic,i_ent)
      end do

    end do  ! i_field = i_wind_u, i_wind_w

  end if  ! ( cmpr(k) % n_points > 0 )
end do  ! k = k_bot_conv, k_bl_top

! Normalise the vertical mean momentum
do i_field = i_wind_u, i_wind_w
  do ic = 1, par_bl_top % cmpr % n_points
    par_bl_mean_fields(ic,i_field)                                             &
      = par_bl_mean_fields(ic,i_field)                                         &
      / par_bl_top % par_super(ic,i_massflux_d)
  end do
end do


!----------------------------------------------------------------
! 9) Work out momentum to entrain from each
!    level below BL-top such that it integrates to the parcel
!    momentum at the BL-top
!----------------------------------------------------------------

do k = k_bot_conv, k_bl_top
  if ( cmpr(k) % n_points > 0 ) then

    ! Add perturbations to entrained momentum so-as to
    ! scale the vertical integral of entrained
    ! momentum to the value in the parcel at BL-top
    do i_field = i_wind_u, i_wind_w
      do ic2 = 1, cmpr(k) % n_points
        ic = index_ic_res(ic2,k)
        ic3 = index_ic_top(ic2,k)
        res_source(i_type,i_layr,k) % fields_super(ic,i_field)                 &
         = res_source(i_type,i_layr,k) % fields_super(ic,i_field)              &
         + par_bl_top % mean_super(ic3,i_field)                                &
         - par_bl_mean_fields(ic3,i_field)
      end do
    end do

  end if  ! ( cmpr(k) % n_points > 0 )
end do  ! k = k_bot_conv, k_bl_top


!----------------------------------------------------------------
! 10) Calculate entrained cloud-fractions
!----------------------------------------------------------------

if ( l_cv_cloudfrac ) then
  do k = k_bot_conv, k_bl_top
    if ( cmpr(k) % n_points > 0 ) then

      ! Extract the entrained water-species and enthalpy
      ! back from the res_source arrays
      do i_field = i_temperature, i_qc_last
        do ic2 = 1, cmpr(k) % n_points
          ic = index_ic_res(ic2,k)
          fields_k(ic2,i_field)                                                &
           = res_source(i_type,i_layr,k)%fields_super(ic,i_field)
        end do
      end do

      ! Convert entrained enthalpy back to temperature
      call set_cp_tot( cmpr(k)%n_points,                                       &
                       par_bl_top % cmpr % n_points,                           &
                       fields_k(:,i_q_vap),                                    &
                       fields_k(:,i_qc_first:i_qc_last),                       &
                       cp_tot )
      do ic2 = 1, cmpr(k) % n_points
        fields_k(ic2,i_temperature) = fields_k(ic2,i_temperature)              &
                                    / cp_tot(ic2)
      end do

      ! Calculate dry virtual temperature T ( 1 + Rv/Rd q_vap )
      call calc_virt_temp_dry( cmpr(k)%n_points,                               &
                               fields_k(:,i_temperature),                      &
                               fields_k(:,i_q_vap),                            &
                               virt_temp_dry )

      ! For each cloud-fraction field...
      do i_field = i_cf_liq, i_cf_bulk

        ! Compress environment cloud-fractions
        lb = lbound(fields_np1%list(i_field)%pt)
        ub = ubound(fields_np1%list(i_field)%pt)
        call compress( cmpr(k), lb(1:2), ub(1:2),                              &
                       fields_np1%list(i_field)%pt(:,:,k),                     &
                       fields_k(:,i_field) )
        ! Remove spurious negative values from input data
        do ic2 = 1, cmpr(k) % n_points
          fields_k(ic2,i_field) = max( fields_k(ic2,i_field),                  &
                                       zero )
        end do

        ! Convert to cloud volume per unit dry-mass, i.e.
        ! scale by dry virtual temperature T ( 1 + Rv/Rd q_vap )
        do ic2 = 1, cmpr(k) % n_points
          fields_k(ic2,i_field) = fields_k(ic2,i_field)                        &
                                * virt_temp_dry(ic2)
        end do

        ! Save entrained cloud volumes in the res_source array
        do ic2 = 1, cmpr(k) % n_points
          ic = index_ic_res(ic2,k)
          res_source(i_type,i_layr,k) % fields_super(ic,i_field)               &
            = fields_k(ic2,i_field)
        end do

      end do  ! i_field = i_cf_liq, i_cf_bulk

    end if  ! ( cmpr(k) % n_points > 0 )
  end do  ! k = k_bot_conv, k_bl_top
end if  ! ( l_cv_cloudfrac )


!----------------------------------------------------------------
! 11) Calculate entrained tracer concentrations
!----------------------------------------------------------------

if ( n_fields_tot > n_fields ) then

  ! Initialise vertical integrals of tracer in the BL
  do i_field = n_fields+1, n_fields_tot
    do ic = 1, par_bl_top % cmpr % n_points
      par_bl_mean_fields(ic,i_field) = zero
    end do
  end do

  do k = k_bot_conv, k_bl_top
    if ( cmpr(k) % n_points > 0 ) then

      ! For each tracer field
      do i_field = n_fields+1, n_fields_tot

        ! Compress the environment tracers onto current points
        lb = lbound(fields_np1%list(i_field)%pt)
        ub = ubound(fields_np1%list(i_field)%pt)
        call compress( cmpr(k), lb(1:2), ub(1:2),                              &
                       fields_np1%list(i_field)%pt(:,:,k),                     &
                       fields_k(:,i_field) )

        ! Store environment tracer in the res_source arrays
        do ic2 = 1, cmpr(k) % n_points
          ic = index_ic_res(ic2,k)
          res_source(i_type,i_layr,k) % fields_super(ic,i_field)               &
            = fields_k(ic2,i_field)
        end do

        ! Add up vertical integral of tracer
        do ic2 = 1, cmpr(k) % n_points
          ic = index_ic_res(ic2,k)
          ic3 = index_ic_top(ic2,k)
          par_bl_mean_fields(ic3,i_field)                                      &
            = par_bl_mean_fields(ic3,i_field)                                  &
            + fields_k(ic2,i_field)                                            &
              * res_source(i_type,i_layr,k) % res_super(ic,i_ent)
        end do

      end do  ! i_field = n_fields+1, n_fields_tot

    end if  ! ( cmpr(k) % n_points > 0 )
  end do  ! k = k_bot_conv, k_bl_top

  ! Normalise the vertical mean tracer
  do i_field = n_fields+1, n_fields_tot
    do ic = 1, par_bl_top % cmpr % n_points
      par_bl_mean_fields(ic,i_field)                                           &
        = par_bl_mean_fields(ic,i_field)                                       &
        / par_bl_top % par_super(ic,i_massflux_d)
    end do
  end do

  do k = k_bot_conv, k_bl_top
    if ( cmpr(k) % n_points > 0 ) then

      ! Add perturbations to entrained tracer so-as to
      ! scale the vertical integral of entrained
      ! tracer to the value in the parcel at BL-top
      do i_field = n_fields+1, n_fields_tot
        do ic2 = 1, cmpr(k) % n_points
          ic = index_ic_res(ic2,k)
          ic3 = index_ic_top(ic2,k)

          ! Calc ratio of BL-top parcel value over mean value in
          ! the BL, with safety-check to avoid div-by-zero:
          if ( abs(par_bl_mean_fields(ic3,i_field))                            &
               > sqrt_min_float ) then
            fac = par_bl_top % mean_super(ic3,i_field)                         &
                / par_bl_mean_fields(ic3,i_field)
          else
            fac = one
          end if

          ! Apply correction factor
          res_source(i_type,i_layr,k) % fields_super(ic,i_field)               &
           = res_source(i_type,i_layr,k)%fields_super(ic,i_field)              &
           * fac
        end do
      end do

    end if  ! ( cmpr(k) % n_points > 0 )
  end do  ! k = k_bot_conv, k_bl_top

end if  ! ( n_fields_tot > n_fields )


!----------------------------------------------------------------
! 12) Scale entrained properties by entrained mass to get
!     resolved-scale source terms
!----------------------------------------------------------------

! Note: convention is that entrainment implies a negative
! value for resolved-scale source terms, since it removes
! mass, water, enthalpy etc from the environment.

! Note that for downdrafts we actually have detrainment below
! the BL-top instead of entrainment

! For downdrafts...
if ( l_down ) then

  ! Calculated mass at each level is actually detrainment
  ! not entrainment; copy into detrainment field and
  ! zero the entrainment
  do k = k_bot_conv, k_bl_top
    if ( cmpr(k) % n_points > 0 ) then
      do ic2 = 1, cmpr(k) % n_points
        ic = index_ic_res(ic2,k)
        res_source(i_type,i_layr,k) % res_super(ic,i_det)                      &
          = res_source(i_type,i_layr,k) % res_super(ic,i_ent)
        res_source(i_type,i_layr,k) % res_super(ic,i_ent) = zero
      end do
    end if
  end do

  ! Multiply all detrained properties by detrained mass to
  ! get resolved-scale source term from detrainment
  do k = k_bot_conv, k_bl_top
    if ( cmpr(k) % n_points > 0 ) then
      do i_field = 1, n_fields_tot
        do ic2 = 1, cmpr(k) % n_points
          ic = index_ic_res(ic2,k)
          res_source(i_type,i_layr,k) % fields_super(ic,i_field)               &
          = res_source(i_type,i_layr,k)%fields_super(ic,i_field)               &
          * res_source(i_type,i_layr,k) % res_super(ic,i_det)
        end do
      end do
    end if
  end do

  ! For updrafts...
else  ! ( l_down )

  ! Zero the detrainment (only have uniform entrainment)
  do k = k_bot_conv, k_bl_top
    if ( cmpr(k) % n_points > 0 ) then
      do ic2 = 1, cmpr(k) % n_points
        ic = index_ic_res(ic2,k)
        res_source(i_type,i_layr,k) % res_super(ic,i_det) = zero
      end do
    end if
  end do

  ! Multiply all entrained properties by minus entrained mass to
  ! get resolved-scale sink term from entrainment
  do k = k_bot_conv, k_bl_top
    if ( cmpr(k) % n_points > 0 ) then
      do i_field = 1, n_fields_tot
        do ic2 = 1, cmpr(k) % n_points
          ic = index_ic_res(ic2,k)
          res_source(i_type,i_layr,k) % fields_super(ic,i_field)               &
          = res_source(i_type,i_layr,k)%fields_super(ic,i_field)               &
          * (-res_source(i_type,i_layr,k) % res_super(ic,i_ent))
        end do
      end do
    end if
  end do

end if  ! ( l_down )


! Deallocate compression indices
deallocate( cmpr )


return
end subroutine homog_conv_bl


end module homog_conv_bl_mod
