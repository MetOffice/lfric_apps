! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Description:
!   Calculates the dust emission fraction for all land surface tiles, or
!   provides an aggregate value in the case of a single tile scheme.

! Method:
!   For a given input tile it calculates the bare soil emission fraction
!   for that tile, using the required inputs. The different methods for
!   doing this calculation are dependent on the dust_veg_emiss switch:
!   0: Do not allow dust emission from tiles other than the bare soil tile
!   1: Allow emission from all vegetated tiles, according to their bare soil
!      radiative fraction (used to calculate the surface albedo), scaled by
!      the dust_veg_sc constant for each land surface type.
!   2 & 3: Allow emission from all vegetated tiles, with the fraction
!      increasing linearly from 0.0 at some threshold LAI, to 1.0 at LAI=0.0.
!      LAIt=0.1 for option 2, and LAIt=0.3 for option 3.
!      Option 2: Yoshika et al., 2007 (doi:10.1175/JCLI4056.1) and
!        Mahowald et al., 2006 (doi:10.1029/2005JD006653)
!      Option 3: Mahowald et al., 2011 (doi:10.5194/bg-8-387-2011)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dust

! Code Description:
!   Language: Fortran 95
!   Programming standard : UMDP 3
module dust_calc_emiss_frac_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'DUST_CALC_EMISS_FRAC_MOD'
contains

subroutine dust_calc_emiss_frac(land_pts, ntiles, tile_pts, tile_index,        &
    frac, lai_ft, dust_emiss_frac)

use parkind1, only: jpim, jprb       !DrHook
use yomhook,  only: lhook, dr_hook   !DrHook
use jules_surface_types_mod, only: npft, ntype, soil
use pftparm, only: dust_veg_scj
use dust_parameters_mod, only: dust_veg_emiss

implicit none

! Subroutine arguments
! in arguments
integer, intent(in) :: land_pts       ! No. land points in whole grid
integer, intent(in) :: ntiles         ! No. land-surface tile types
integer, intent(in) :: tile_pts(ntype)! Total number of tiles
integer, intent(in) :: tile_index(land_pts,ntype)
                                    ! Index of tiles on landpts
real(kind=real_umphys),    intent(in) :: frac(land_pts,ntype)
                                    ! in Fractions of surface types.
real(kind=real_umphys),    intent(in) :: lai_ft(land_pts,npft)
                                 ! in Leaf area index
! out arguments
real(kind=real_umphys), intent(out):: dust_emiss_frac(land_pts,ntiles)

! Local variables
integer :: l                          !index of land pt
integer :: m                          !loop counter, tiles or tile types
integer :: n                          !loop counter, tile points
real(kind=real_umphys)    :: rfracbs(land_pts,npft)
                                      !raditive fraction of bare soil
real(kind=real_umphys)    :: dust_veg_sc(npft)
                                      !scaling factor for each PFT
real(kind=real_umphys)    :: dust_lai_thresh            !LAI threshold value

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='DUST_CALC_EMISS_FRAC'

! End of header

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! explicitly initialise the dust_emiss_frac array as all zeroes as it will be
! calculated as a sum over tiles
!$OMP PARALLEL DEFAULT(none) private(l,m)                                      &
!$OMP SHARED(ntiles,land_pts,dust_emiss_frac,npft,rfracbs,dust_veg_emiss)
!$OMP do SCHEDULE(STATIC) COLLAPSE(2)
do m=1,ntiles
  do l=1,land_pts
    dust_emiss_frac(l,m) = 0.0
  end do
end do
!$OMP end do
! and test to see if this does anything...
if (dust_veg_emiss /= 0) then
!$OMP do SCHEDULE(STATIC) COLLAPSE(2)
  do m=1,npft
    do l=1,land_pts
      rfracbs(l,m) = 0.0
    end do
  end do
!$OMP end do
end if
!$OMP end PARALLEL

! start with the null option, where the emission is only from the bare soil
! tile, and none from vegetation
if (dust_veg_emiss == 0) then

  if (ntiles == 1) then
    m = soil
    do n = 1, tile_pts(m)
      l = tile_index(n,m)
      dust_emiss_frac(l,1) = frac(l,soil)
    end do

  else
    m = soil
    do n = 1, tile_pts(m)
      l = tile_index(n,m)
      dust_emiss_frac(l,m) = frac(l,soil)
    end do

  end if

else if (dust_veg_emiss == 1) then
  ! Option 1 is to include the bare soil radiative fraciton of bare soil
  ! of all vegetated tiles, using the radiative fraction from the leaf
  ! area index:

  ! Get the radiative fraction of bare soil within each of the vegetated tiles:
  ! assuming that the vegetated ones are always first (which they should be)
  ! get the scaling parameters for each PFT, for JULES or non-JULES:
  do m = 1, npft
    dust_veg_sc(m)=dust_veg_scj(m)
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none) private(l,n)                  &
!$OMP SHARED(tile_pts,tile_index,rfracbs,lai_ft,m)
    do n = 1, tile_pts(m)
      l = tile_index(n,m)
      rfracbs(l,m) = exp(-lai_ft(l,m)/2.0)
    end do
!$OMP end PARALLEL do
  end do


  ! For 1 tile model, the output is an aggregated emission fraction of all tiles
  if (ntiles == 1) then
    ! start with the bare soil fraction, where there are tile_pts for that
    ! fraction:
!$OMP PARALLEL DEFAULT(none) private(l,m,n)                                    &
!$OMP SHARED(soil,tile_pts,tile_index,frac,dust_emiss_frac,dust_veg_sc,        &
!$OMP npft,rfracbs)
    m = soil
!$OMP do SCHEDULE(STATIC)
    do n = 1, tile_pts(m)
      l = tile_index(n,m)
      dust_emiss_frac(l,1) = frac(l,m)
    end do
!$OMP end do

    ! Add on the radiative fractions of bare soil over the vegetated tiles
    do m = 1, npft
!$OMP do SCHEDULE(STATIC)
      do n = 1, tile_pts(m)
        l = tile_index(n,m)
          ! output on 1 tile, so (l,1)
        dust_emiss_frac(l,1) = (dust_veg_sc(m)*frac(l,m)*rfracbs(l,m))         &
                             + dust_emiss_frac(l,1)
      end do
!$OMP end do
    end do
!$OMP end PARALLEL
    ! ends the 1 tile version
  else
    ! For the multiple tile code, basically the same but each tile has an array:

    ! start with the vegetated tiles (assumes they are first in the tile arrays)
!$OMP PARALLEL DEFAULT(none) private(l,m,n)                                    &
!$OMP SHARED(npft,dust_emiss_frac,dust_veg_sc,tile_pts,tile_index,rfracbs,frac,&
!$OMP soil)
    do m = 1, npft
!$OMP do SCHEDULE(STATIC)
      do n = 1, tile_pts(m)
        l = tile_index(n,m)
        dust_emiss_frac(l,m) = dust_veg_sc(m)*frac(l,m)*rfracbs(l,m)
      end do
!$OMP end do
    end do

    ! then do the bare soil tile
    m = soil
!$OMP do SCHEDULE(STATIC)
    do n = 1, tile_pts(m)
      l = tile_index(n,m)
      dust_emiss_frac(l,m) = frac(l,m)
    end do
!$OMP end do
!$OMP end PARALLEL

    ! other tiles (urban, lake, glacier) are left as zero, so do not emit dust

  end if ! ends 1 tile or multiple tile test

else if (dust_veg_emiss == 2 .or. dust_veg_emiss == 3) then
  ! Get the fraction of bare soil within each of the vegetated tiles, for dust
  ! emission, which increases linearly, as LAI decreases, below a threshold.
  !
  ! which is the same scheme with an LAI threshold value of 0.3
  ! (dust_veg_emiss == 3).
  if (dust_veg_emiss == 2) then
    dust_lai_thresh=0.1
  else if (dust_veg_emiss == 3) then
    dust_lai_thresh=0.3
  end if

  do m=1,npft
    do n=1,tile_pts(m)
      l = tile_index(n,m)
      rfracbs(l,m) = max(0.0, (dust_lai_thresh-lai_ft(l,m))/dust_lai_thresh)
    end do
  end do

  ! get the sclaing parameters for each PFT, for JULES or non-JULES:
  dust_veg_sc(:)=dust_veg_scj(:)

  if (ntiles == 1) then
    ! start with the bare soil fraction, where there are tile_pts for that
    ! fraction:
    m = soil
    do n = 1, tile_pts(m)
      l = tile_index(n,m)
      dust_emiss_frac(l,1) = frac(l,m)
    end do
    ! Add on the radiative fractions of bare soil over the vegetated tiles
    do m = 1, npft
      do n = 1, tile_pts(m)
        l = tile_index(n,m)
        ! output on 1 tile, so (l,1)
        dust_emiss_frac(l,1) = (dust_veg_sc(m)*frac(l,m)*rfracbs(l,m)) +       &
                               dust_emiss_frac(l,1)
      end do
    end do
    ! ends the 1 tile version
  else
    ! For the multiple tile code, basically the same but each tile has an array:
    ! start with the vegetated tiles (assumes they are first in the tile arrays)
    do m = 1, npft
      do n = 1, tile_pts(m)
        l = tile_index(n,m)
        dust_emiss_frac(l,m) = dust_veg_sc(m)*frac(l,m)*rfracbs(l,m)
      end do
    end do
    ! then do the bare soil tile
    m = soil
    do n = 1, tile_pts(m)
      l = tile_index(n,m)
      dust_emiss_frac(l,m) = frac(l,m)
    end do
    ! other tiles (urban, lake, glacier) are left as zero, so do not emit dust
  end if ! ends 1 tile or multiple tile test

end if ! Ends test on dust_veg_emiss switch (currently 0,1,2,3).

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine dust_calc_emiss_frac
end module dust_calc_emiss_frac_mod
