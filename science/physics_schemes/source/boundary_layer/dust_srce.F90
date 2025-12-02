! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Calculates flux of mineral dust entrained into atmosphere

! Method:
!   Calculates mineral dust flux as a function of friction velocity (U*),
!   soil moisture and soil porperties developed from description in:
!     "Modelling the atmospheric lifecycle..."
!     Woodward, JGR106, D16, pp18155-18166, 2001.
!   Treatment of soil moisture as in:
!     "Parametrization of the increase of..."
!     Fecan et al, Ann. Geophysicae 17, 149-157, 1999.!
!   NB: Currently only calculates flux from bare soil tiles

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dust

! Code description:
!  Language: Fortran 95
!  Programming standard : UMDP 3

! Subroutine Interface:
module dust_srce_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'DUST_SRCE_MOD'
contains

subroutine dust_srce(                                                          &
! in arguments
     land_pts, ntiles, tile_pts, tile_index,                                   &
     fland, tstar_tile, rhostar_land, soil_layer_moisture, snow_tile,          &
     u_s_std_tile, mrel_land, clay_land, sand_land, ho2r2_orog, emsc_land,     &
! out arguments
     dust_flux_tile, u_s_t_tile, u_s_t_dry_tile )

use conversions_mod, only: zerodegc
use dust_parameters_mod, only:                                                 &
     ! number of divisions that can be lifted from the surface and the
     ! number that can be blown horizontally along the surface and
     ! contribute to the lifting of the 1st NDIV divisions and the
     ! number of these that contain data in the dust soil ancillaries
     ndiv, ndivh, ndivl,                                                       &
     ! impact U*t derived from Bagnold (1941)
     ustd_bas,                                                                 &
     ! additional and multiplicative tunings to U* and U*t
     us_aa, us_am, ust_aa, ust_am,                                             &
     ! multiplicative tuning factor to level 1 soil moisture
     sm_corr,                                                                  &
     ! maximum clay fraction used in calculation of vertical flux
     clay_max, l_dust_clay_as_max,                                             &
     ! various limits used in diagnosis of horizontal flux
     u_s_min, snowmin, h_orog_limit, fland_lim,                                &
     ! constants used in horiz. and vert. fluxes (see module for details)
     horiz_c, horiz_d, vert_a, vert_b, vert_c,                                 &
     ! switch to diagnose vertical flux using a fixed size distribution
     ! if off the vertical flux is proportional to the horizontal flux
     l_fix_size_dist,                                                          &
     ! proportion of flux from each bin if using the fixed distribution
     size_dist,                                                                &
     ! Two-bin dust flag
     l_twobin_dust,                                                            &
     ! dust emission scaling
     l_dust_emp_sc
use jules_soil_mod, only: sm_levels
use jules_surface_types_mod, only: ntype, soil
use planet_constants_mod, only: g
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Subroutine arguments
integer, intent(in) :: land_pts       ! No. land points in whole grid
integer, intent(in) :: ntiles         ! No. land-surface tile types
integer, intent(in) :: tile_pts(ntype)! Total number of tiles
integer, intent(in) :: tile_index(land_pts,ntype)
                                    ! Index of tiles on landpts
real(kind=real_umphys), intent(in) :: fland(land_pts)
                                      ! Land fraction on land points
real(kind=real_umphys), intent(in) :: tstar_tile(land_pts,ntiles)
                                    ! Surface temp on tiles
real(kind=real_umphys), intent(in) :: rhostar_land(land_pts)
                                    ! Surface air density on land pts
real(kind=real_umphys), intent(in) :: soil_layer_moisture(land_pts,sm_levels)
                                    ! Soil moisture (kg m-2)
real(kind=real_umphys), intent(in) :: snow_tile(land_pts,ntiles)
                                    ! Lying snow on tiles (kg m-2)
real(kind=real_umphys), intent(in) :: u_s_std_tile(land_pts,ntiles)
                                    ! Friction velocity on tiles
real(kind=real_umphys), intent(in) :: mrel_land(land_pts,ndivl)
                                    ! Relative soil mass per size div
real(kind=real_umphys), intent(in) :: clay_land(land_pts)
                                    ! Soil clay fraction
real(kind=real_umphys), intent(in) :: sand_land(land_pts)
                                    ! Soil sand fraction
real(kind=real_umphys), intent(in) :: ho2r2_orog(land_pts)
                                    ! Half peak to trough ht of orog
real(kind=real_umphys), intent(in) :: emsc_land(land_pts)
                                    ! emission scaling factor
! Outputs
real(kind=real_umphys), intent(out):: dust_flux_tile(land_pts,ntiles,ndiv)
                                    ! Dust flux (kg m-2 s-1)
real(kind=real_umphys), intent(out):: u_s_t_tile(land_pts,ntiles,ndivh)
                                    ! Thresh. friction vel. per tile
                                    ! (all 9 divisions)
real(kind=real_umphys), intent(out):: u_s_t_dry_tile(land_pts,ntiles,ndivh)
                                    ! Thresh. frict. vel. per tile
                                    ! excluding soil moisture effects
! Local variables
integer :: l                          !index of land pt
integer :: m                          !loop counter, tile types
integer :: n                          !loop counter, tile points
integer :: idiv                       !loop counter, dust divisions

real(kind=real_umphys) :: ratio       ! u_s_t/u_s_std
real(kind=real_umphys) :: mrel7(land_pts)
                                      ! mrel for 31.6 - 100 um radius
real(kind=real_umphys) :: mrel8(land_pts)               ! mrel for 100 - 316 um
real(kind=real_umphys) :: mrel9(land_pts)               ! mrel for 316 - 1000 um
real(kind=real_umphys) :: mrel_land_all(land_pts,ndivh) !all the mrels

real(kind=real_umphys) :: horiz_flux_tile(land_pts,ntiles,ndivh)
                                    !horizontal flux per tile
real(kind=real_umphys) :: tot_horiz_flux_tile(land_pts,ntiles)
                                    !total over all divs
real(kind=real_umphys) :: horiz_flux_789(land_pts,ntiles)
                                    !total over div 7,8,9
real(kind=real_umphys) :: horiz_flux_1to6(land_pts,ntiles)
                                    !total over divs 1 to 6
real(kind=real_umphys) :: us          ! "tuned" U* on soil tiles
real(kind=real_umphys) :: smt(land_pts)
                                      ! Soil moisture term in U*t calc
real(kind=real_umphys) :: smt_work    ! Temp variable in smt calc

integer :: soil_tile                     ! Tile on which to calculate
                                         ! bare soil dust flux

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='DUST_SRCE'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! parameters in block below - ndivh, ndivl
!$OMP PARALLEL DEFAULT(none) private(idiv,m,n,l)                               &
!$OMP SHARED(ndiv,ntiles,tile_pts,tile_index,dust_flux_tile,                   &
!$OMP   tot_horiz_flux_tile,horiz_flux_789,horiz_flux_1to6,                    &
!$OMP   u_s_t_tile,u_s_t_dry_tile,horiz_flux_tile,soil,soil_tile,              &
!$OMP   mrel7,mrel8,mrel9,mrel_land_all,sand_land,land_pts,mrel_land)
! Initialisation
!$OMP do SCHEDULE(STATIC) COLLAPSE(2)
do idiv = 1, ndiv
  do m = 1, ntiles
    do n = 1, tile_pts(m)
      l = tile_index(n,m)
      dust_flux_tile(l,m,idiv) = 0.0
    end do !tile_pts
  end do !ntiles
end do !ndiv
!$OMP end do NOWAIT

!$OMP do SCHEDULE(STATIC)
do m = 1, ntiles
  do n = 1, tile_pts(m)
    l = tile_index(n,m)
    tot_horiz_flux_tile(l,m) = 0.0
    horiz_flux_789(l,m) = 0.0
    horiz_flux_1to6(l,m) = 0.0
  end do
end do
!$OMP end do NOWAIT

!$OMP do SCHEDULE(STATIC) COLLAPSE(2)
do idiv = 1, ndivh
  do m = 1, ntiles
    do n = 1, tile_pts(m)
      l = tile_index(n,m)
      u_s_t_tile(l,m,idiv) = 0.0
      u_s_t_dry_tile(l,m,idiv) = 0.0
      horiz_flux_tile(l,m,idiv) = 0.0
    end do
  end do
end do
!$OMP end do

!$OMP SINGLE
! In 1 bin scheme, bare soil dust flux is calculated on single aggregated tile
if (ntiles > 1) then
  soil_tile = soil
else
  soil_tile = 1
end if
!$OMP end SINGLE

! Calculate relative mass for all divs and put into single array
!$OMP do SCHEDULE(STATIC)
do l = 1, land_pts
  mrel7(l)=sand_land(l)*0.312
  mrel8(l)=sand_land(l)*0.312
  mrel9(l)=sand_land(l)*0.312
  mrel_land_all(l,7)=mrel7(l)
  mrel_land_all(l,8)=mrel8(l)
  mrel_land_all(l,9)=mrel9(l)
end do !land_pts
!$OMP end do NOWAIT

! This is explicitly specified as ndivl=6, even for the two-bin scheme
!$OMP do SCHEDULE(STATIC) COLLAPSE(2)
do idiv = 1, ndivl
  do l = 1, land_pts
    mrel_land_all(l,idiv)=mrel_land(l,idiv)
  end do
end do
!$OMP end do
!$OMP end PARALLEL


! Increase U*t for "tuned" soil moisture (w) > a threshold w' (fecan'99)
! smt_work = w-w'
! w is % volmetric sm, but approximate as same as s_l_m in top JULES level
! ( w = s_l_m * 100(for %age) * (100cm/10cm) / rho=998 ~ s_l_m)
do l = 1, land_pts
  smt_work = soil_layer_moisture(l,1) * sm_corr -                              &
       14.0*clay_land(l)*clay_land(l) - 17.0*clay_land(l)
  if (smt_work > 0.0) then
    smt(l) =  ((1.0 + 1.21* smt_work **.68 )**.5)
  else
    smt(l) = 1.0
  end if
end do

! Loop over points in each div calculating dust flux, if any
do idiv = 1, ndivh
  do m = 1, ntiles
    do n = 1, tile_pts(m)
      l = tile_index(n,m)

      ! Calculate "tuned" U* and/or U*t for this tile

      u_s_t_dry_tile(l,m,idiv)=ustd_bas(idiv)*ust_am+ust_aa
      u_s_t_tile(l,m,idiv)=ustd_bas(idiv)*smt(l)*ust_am+ust_aa

      ! use the U* on bare soil, regardless of whether this is a bare
      ! soil tile or not - we're looking at the U* of any bare soil in
      ! the tile, the fraction of which is determined by
      ! dust_calc_emiss_frac
      ! For the 1 tile case, u_s_std_tile is calculated as bare soil only
      us=u_s_std_tile(l,soil_tile)*us_am+us_aa

      ! Horizontal flux

      if ( (u_s_t_tile(l,m,idiv) < us) .and.                                   &
           (fland(l) >  fland_lim) .and.                                       &
           (tstar_tile(l,m) > zerodegc) .and.                                  &
           (snow_tile(l,m)  <   snowmin) .and.                                 &
           (mrel_land_all(l,idiv)  >   0.0) .and.                              &
           (ho2r2_orog(l)  <=  h_orog_limit) .and.                             &
           ((soil_layer_moisture(l,1)*sm_corr <                                &
           (clay_land(l)+.12)/.03)) .and.                                      &
           (us > u_s_min) ) then

        ratio = u_s_t_tile(l,m,idiv) / us
        horiz_flux_tile(l,m,idiv) = horiz_c * rhostar_land(l) *                &
             horiz_d * us**3.0 * (1.0 + ratio) *                               &
             (1.0 - ratio*ratio) * mrel_land_all(l,idiv) / g
        tot_horiz_flux_tile(l,m)=                                              &
             tot_horiz_flux_tile(l,m)+horiz_flux_tile(l,m,idiv)

      end if
    end do !tile_pts
  end do ! tiles
end do !ndiv

do m = 1, ntiles
  do n = 1, tile_pts(m)
    l = tile_index(n,m)
    horiz_flux_789(l,m)=horiz_flux_tile(l,m,7)+                                &
         horiz_flux_tile(l,m,8)+horiz_flux_tile(l,m,9)
    horiz_flux_1to6(l,m)=horiz_flux_tile(l,m,1)+horiz_flux_tile(l,m,2)         &
        +horiz_flux_tile(l,m,3)+horiz_flux_tile(l,m,4)+                        &
         horiz_flux_tile(l,m,5)+horiz_flux_tile(l,m,6)
  end do
end do

! parameters in block below - vert_a, clay_max, vert_b, vert_c
!$OMP PARALLEL DEFAULT(none) private(idiv,m,n,l) SHARED(l_fix_size_dist,       &
!$OMP  l_twobin_dust,l_dust_clay_as_max,l_dust_emp_sc,                         &
!$OMP  ndiv,ntiles,tile_pts,tile_index,snow_tile,snowmin,                      &
!$OMP  dust_flux_tile,size_dist,horiz_flux_1to6,                               &
!$OMP  horiz_flux_789,clay_land,mrel_land,emsc_land,horiz_flux_tile,fland)
! Vertical flux
! Again, this calculation is for bare soil only, but could be
! extended to other tile types
if (l_fix_size_dist .or. l_twobin_dust) then
  !   Lift dust using fixed size distribution
!$OMP do SCHEDULE(STATIC) COLLAPSE(2)
  do idiv = 1, ndiv
    do m = 1, ntiles
      do n = 1, tile_pts(m)
        l = tile_index(n,m)
        if (snow_tile(l,m) < snowmin) then
          ! Calculate proportion of vertical flux to come from each bin
          dust_flux_tile(l,m,idiv) = vert_c * size_dist(idiv) *                &
               (horiz_flux_1to6(l,m)+horiz_flux_789(l,m))
          !           Calculate vertical flux
          if (l_dust_clay_as_max) then
            dust_flux_tile(l,m,idiv) =  dust_flux_tile(l,m,idiv) *             &
               10**(vert_a * clay_max + vert_b)
          else
            dust_flux_tile(l,m,idiv) =  dust_flux_tile(l,m,idiv) *             &
               10**(vert_a * min(clay_land(l),clay_max) + vert_b)
          end if
        end if
      end do
    end do
  end do
!$OMP end do
else
  !   Lift dust using diagnosed size distribution
!$OMP do SCHEDULE(STATIC) COLLAPSE(2)
  do idiv = 1, ndiv
    do m = 1, ntiles
      do n = 1, tile_pts(m)
        l = tile_index(n,m)
        !         Calculate proportion of vertical flux to come from each bin
        if ((snow_tile(l,m) < snowmin) .and.                                   &
            (mrel_land(l,idiv) > 0.0) .and.                                    &
            (horiz_flux_1to6(l,m) > 0.0)) then
          dust_flux_tile(l,m,idiv) = vert_c *                                  &
               (1.0 + horiz_flux_789(l,m)/horiz_flux_1to6(l,m)) *              &
               horiz_flux_tile(l,m,idiv)
        end if
        !         Calculate vertical flux
        if (l_dust_clay_as_max) then
          dust_flux_tile(l,m,idiv) =  dust_flux_tile(l,m,idiv) *               &
               10**(vert_a * clay_max + vert_b)
        else
          dust_flux_tile(l,m,idiv) =  dust_flux_tile(l,m,idiv) *               &
              10**(vert_a * min(clay_land(l),clay_max) + vert_b)
        end if
      end do
    end do
  end do
!$OMP end do
end if

! Multiply by scalings: Correct fluxes for land fraction - for coastal points
! and, if applied, the emission scaling field
if (l_dust_emp_sc) then
!$OMP do SCHEDULE(STATIC) COLLAPSE(2)
  do idiv = 1, ndiv
    do m = 1, ntiles
      do n = 1, tile_pts(m)
        l = tile_index(n,m)
        ! Horizontal flux not yet output as a diagnositic
        ! If added, do the following calc (but for 9 divisions)
        !     horiz_flux_tile(l,m,idiv)=horiz_flux_tile(l,m,idiv)*fland(l)
        dust_flux_tile(l,m,idiv) = dust_flux_tile(l,m,idiv) * fland(l) *       &
                                      emsc_land(l)
      end do !tile_pts
    end do !ntiles
  end do !ndiv
!$OMP end do
else
!$OMP do SCHEDULE(STATIC) COLLAPSE(2)
  do idiv = 1, ndiv
    do m = 1, ntiles
      do n = 1, tile_pts(m)
        l = tile_index(n,m)
        ! Horizontal flux not yet output as a diagnositic
        ! If added, do the following calc (but for 9 divisions)
        !     horiz_flux_tile(l,m,idiv)=horiz_flux_tile(l,m,idiv)*fland(l)
        dust_flux_tile(l,m,idiv) = dust_flux_tile(l,m,idiv) * fland(l)
      end do !tile_pts
    end do !ntiles
  end do !ndiv
!$OMP end do
end if ! l_dust_emp_sc

!$OMP end PARALLEL

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine dust_srce


subroutine dust_emp_scaling(land_pts, ntiles, tile_pts, tile_index, frac,      &
                            mrel_land, clay_land, sand_land, psti_land,        &
                            emsc_land)
! Calculates the dust emission scaling factor using the preferential source
! term input as a dust emission potential - this is compared to the current
! model physics and ancillary input data's dust emission under idealised
! conditions...

use jules_surface_types_mod, only: ntype
use dust_parameters_mod, only: ndiv, ndivl, ndivh

implicit none

! Subroutine arguments
integer, intent(in) :: land_pts       ! No. land points in whole grid
integer, intent(in) :: ntiles         ! No. land-surface tile types
integer, intent(in) :: tile_pts(ntype)! Total number of tiles
integer, intent(in) :: tile_index(land_pts,ntype)
                                    ! Index of tiles on landpts
real(kind=real_umphys), intent(in) :: frac(land_pts,ntype)
                                    ! in Fractions of surface types.
real(kind=real_umphys), intent(in) :: mrel_land(land_pts,ndivl)
                                    ! Relative soil mass per size div
real(kind=real_umphys), intent(in) :: clay_land(land_pts)
                                    ! Soil clay fraction
real(kind=real_umphys), intent(in) :: sand_land(land_pts)
                                    ! Soil sand fraction
real(kind=real_umphys), intent(in) :: psti_land(land_pts)
                                    ! preferential source term on land points
! Outputs: dust emission scaling (on land points)
real(kind=real_umphys), intent(out):: emsc_land(land_pts)

! Local variables
integer :: l                          !index of land pt
integer :: m                          !loop counter, tile types
integer :: n                          !loop counter, tile points
integer :: idiv                       !loop counter, dust divisions
real(kind=real_umphys):: dust_flux_tile(land_pts,ntiles,ndiv)
real(kind=real_umphys):: u_s_t_tile(land_pts,ntiles,ndivh)
real(kind=real_umphys):: u_s_t_dry_tile(land_pts,ntiles,ndivh)
real(kind=real_umphys):: dust_flux_tot(land_pts)

real(kind=real_umphys):: dust_flux_min ! min value of dust flux

! points where the total idealised dust emission is below this will
! be treated as zero. Outputting/calculating this offline from output
! emissions suggests a value of:
dust_flux_min = 1e-9 ! 1 microgram per m2 per s

! Calculate the idealised dust emission:
call dust_srce_ideal(land_pts, ntiles, tile_pts, tile_index, mrel_land,        &
                     clay_land, sand_land,                                     &
                     dust_flux_tile, u_s_t_tile, u_s_t_dry_tile)

! Initialise dust_flux_tot
do l = 1, land_pts
  dust_flux_tot(l) = 0.0
end do ! land_pts initialisation
! Sum up dust_flux_tot over the different size bins
do idiv = 1, ndiv
  do m = 1, ntiles
    do n = 1, tile_pts(m)
      l = tile_index(n,m)
      dust_flux_tot(l) = dust_flux_tot(l) + dust_flux_tile(l,m,idiv) * frac(l,m)
    end do !TILE_PTS
  end do !NTILES
end do !NDIV
do l = 1, land_pts
  if (dust_flux_tot(l) < dust_flux_min) then
    ! check for tiny numbers first:
    emsc_land(l) = 0.0
  else
    emsc_land(l) = psti_land(l) / dust_flux_tot(l)
  end if
end do ! land_pts

return
end subroutine dust_emp_scaling


subroutine dust_srce_ideal(                                                    &
! in arguments
     land_pts, ntiles, tile_pts, tile_index, mrel_land, clay_land, sand_land,  &
! out arguments
     dust_flux_tile, u_s_t_tile, u_s_t_dry_tile )
! A wrapper around dust_srce() that calls it with idealised meteorological
! inputs - can be used to compare against dust emission potential maps to
! produce emission scaling factor maps.
!
! Directly set in the dust emission potential definition:
! u_s_std_tile(land_pts,ntiles) - set all to 0.69
! soil_layer_moisture(land_pts,sm_levels) - set all to 0.0
! snow_tile(land_pts,ntiles) - set all to 0.0
! Not directly set, but assumed less imoprtant:
! tstar_tile(land_pts,ntiles) - set all to 300.0 K
! fland(land_pts) - set to 1.0
! ho2r2_orog() - set to 0.0
! emsc_land(land_pts) - all 1.0

use dust_parameters_mod, only:  us_aa, us_am, ndiv, ndivl, ndivh
use jules_surface_types_mod, only: ntype
use jules_soil_mod, only: sm_levels

implicit none

! Subroutine arguments
integer, intent(in) :: land_pts       ! No. land points in whole grid
integer, intent(in) :: ntiles         ! No. land-surface tile types
integer, intent(in) :: tile_pts(ntype)! Total number of tiles
integer, intent(in) :: tile_index(land_pts,ntype)
                                    ! Index of tiles on landpts
real(kind=real_umphys), intent(in) :: mrel_land(land_pts,ndivl)
                                    ! Relative soil mass per size div
real(kind=real_umphys), intent(in) :: clay_land(land_pts)
                                    ! Soil clay fraction
real(kind=real_umphys), intent(in) :: sand_land(land_pts)
                                    ! Soil sand fraction
! Outputs
real(kind=real_umphys), intent(out):: dust_flux_tile(land_pts,ntiles,ndiv)
                                    ! Dust flux (kg m-2 s-1)
real(kind=real_umphys), intent(out):: u_s_t_tile(land_pts,ntiles,ndivh)
                                    ! Thresh. friction vel. per tile
                                    ! (all 9 divisions)
real(kind=real_umphys), intent(out):: u_s_t_dry_tile(land_pts,ntiles,ndivh)
                                    ! Thresh. frict. vel. per tile
                                    ! excluding soil moisture effects


! Local arguments
real(kind=real_umphys) :: fland_idl(land_pts)
                                    ! Land fraction on land points
real(kind=real_umphys) :: rhostar_land_idl(land_pts)
                                    ! Air density of surface level
real(kind=real_umphys) :: tstar_tile_idl(land_pts,ntiles)
                                    ! Surface temp on tiles
real(kind=real_umphys) :: soil_layer_moisture_idl(land_pts,sm_levels)
                                    ! Soil moisture (kg m-2)
real(kind=real_umphys) :: snow_tile_idl(land_pts,ntiles)
                                    ! Lying snow on tiles (kg m-2)
real(kind=real_umphys) :: u_s_std_tile_idl(land_pts,ntiles)
                                    ! Friction velocity on tiles
real(kind=real_umphys) :: emsc_land_idl(land_pts)
                                    ! Dust emission sclaing factor
real(kind=real_umphys) :: ho2r2_orog_idl(land_pts)
                                    ! Half peak to trough ht of orog

! Now create the idealised values:
! for ustar, apply the tuning parameters in reverse:
u_s_std_tile_idl = (0.69 - us_aa) / us_am
! the other parameters are defined to match dust emission potential
! psuedo obs:
soil_layer_moisture_idl = 0.0 ! dry soil
snow_tile_idl = 0.0 ! no snow
rhostar_land_idl = 1.19586906 ! ACAO std atmosphere at 250m (typical height
                              ! of the in-situ wind tunnel obs)
tstar_tile_idl = 300.0 ! a representative temperature
fland_idl = 1.0 ! all land
ho2r2_orog_idl = 0.0 ! flat

emsc_land_idl = 1.0 ! unscaled

! Now the dust_srce call, but with idealised inputs:
call dust_srce(land_pts, ntiles, tile_pts, tile_index, fland_idl,              &
               tstar_tile_idl, rhostar_land_idl, soil_layer_moisture_idl,      &
               snow_tile_idl, u_s_std_tile_idl, mrel_land,                     &
               clay_land, sand_land, ho2r2_orog_idl, emsc_land_idl,            &
               dust_flux_tile, u_s_t_tile, u_s_t_dry_tile )

return
end subroutine dust_srce_ideal
end module dust_srce_mod
