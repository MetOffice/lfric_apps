! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Function to calculate cloud droplet number concentration.

! Purpose:
!   Cloud droplet number concentration is calculated from aerosol
!   concentration or else fixed values are assigned.

! Method:
!   Sulphate aerosol mass concentration is converted to number
!   concentration by assuming a log-normal size distribution.
!   Sea-salt and/or biomass-burning aerosols may then
!   be added if required. The total is then converted to cloud
!   droplet concentration following the parametrization of
!   Jones et al. (1994) and lower limits are imposed.
!   Alternatively, fixed droplet values are assigned if the
!   parametrization is not required.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: atmos_service_number_droplet

! Description of Code:
!   FORTRAN 90

!- ---------------------------------------------------------------------
module number_droplet_mod

use um_types, only: real_umphys

implicit none

private
public :: number_droplet, min_cdnc_sea_ice, min_cdnc_land, snow_depth_thresh,  &
          land_frac_thresh

! Variables for this module only
logical :: nd_init = .false.         ! Are particle volumes initialised?
real(kind=real_umphys)    :: particle_volume_nh42so4
                                     ! Mean volume of (NH4)2SO4
real(kind=real_umphys)    :: particle_volume_biomass
                                     ! Mean volume of biomass smoke
real(kind=real_umphys)    :: particle_volume_biogenic
                                     ! Mean volume of biogenic aerosol
real(kind=real_umphys)    :: particle_volume_ocff
                                     ! Mean volume of OCFF aerosol
real(kind=real_umphys)    :: particle_volume_nitrate
                                     ! Mean volume of nitrate aerosol


!             Median radius of log-normal distribution for (NH4)2SO4
real(kind=real_umphys), parameter :: radius_0_nh42so4  = 9.5e-8
!             Geometric standard deviation of same
real(kind=real_umphys), parameter :: sigma_0_nh42so4   = 1.4
!             Density of ammonium sulphate aerosol
real(kind=real_umphys), parameter :: density_nh42so4   = 1.769e+03
!             Median radius of log-normal distribution for biomass smoke
real(kind=real_umphys), parameter :: radius_0_biomass  = 1.2e-07
!             Geometric standard deviation of same
real(kind=real_umphys), parameter :: sigma_0_biomass   = 1.30
!             Density of biomass smoke aerosol
real(kind=real_umphys), parameter :: density_biomass   = 1.35e+03
!             Median radius of log-normal dist. for biogenic aerosol
real(kind=real_umphys), parameter :: radius_0_biogenic = 9.5e-08
!             Geometric standard deviation of same
real(kind=real_umphys), parameter :: sigma_0_biogenic  = 1.50
!             Density of biogenic aerosol
real(kind=real_umphys), parameter :: density_biogenic  = 1.3e+03
!             Median radius of log-normal dist. for OCFF aerosol
real(kind=real_umphys), parameter :: radius_0_ocff     = 0.12e-06
!             Geometric standard deviation of same
real(kind=real_umphys), parameter :: sigma_0_ocff      = 1.30
!             Density of OCFF aerosol
real(kind=real_umphys), parameter :: density_ocff      = 1350.0
!             Median radius of log-normal dist. for nitrate aerosol
real(kind=real_umphys), parameter :: radius_0_nitrate  = 9.5e-8
!             Geometric standard deviation of same
real(kind=real_umphys), parameter :: sigma_0_nitrate   = 1.4
!             Density of nitrate aerosol
real(kind=real_umphys), parameter :: density_nitrate   = 1.725e+03

!             Minimum droplet numbers [m-3]
real(kind=real_umphys), parameter :: min_cdnc_sea_ice  = 5.0e+06
real(kind=real_umphys), parameter :: min_cdnc_land     = 35.0e+06
!             Depth of snow threshold [m]
real(kind=real_umphys), parameter :: snow_depth_thresh = 5000.0
!             Land fraction threshold []
real(kind=real_umphys), parameter :: land_frac_thresh  = 0.2


character(len=*), parameter, private :: ModuleName='NUMBER_DROPLET_MOD'

contains
!- ---------------------------------------------------------------------
subroutine number_droplet_init()
use conversions_mod, only: pi
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
implicit none

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='NUMBER_DROPLET_INIT'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  ! Initialize particle volume constants
particle_volume_nh42so4 =                                                      &
    (4.0e+00*pi/3.0e+00)*radius_0_nh42so4**3                                   &
    *exp(4.5e+00*(log(sigma_0_nh42so4))**2)

particle_volume_biomass =                                                      &
    (4.0e+00*pi/3.0e+00)*radius_0_biomass**3                                   &
    *exp(4.5e+00*(log(sigma_0_biomass))**2)

particle_volume_biogenic =                                                     &
    (4.0e+00*pi/3.0e+00)*radius_0_biogenic**3                                  &
    *exp(4.5e+00*(log(sigma_0_biogenic))**2)

particle_volume_ocff =                                                         &
    (4.0e+00*pi/3.0e+00)*radius_0_ocff**3                                      &
    *exp(4.5e+00*(log(sigma_0_ocff))**2)

particle_volume_nitrate =                                                      &
    (4.0e+00*pi/3.0e+00)*radius_0_nitrate**3                                   &
    *exp(4.5e+00*(log(sigma_0_nitrate))**2)

! Now volumes are initialised
nd_init = .true.

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine number_droplet_init

!- ---------------------------------------------------------------------

subroutine number_droplet(                                                     &
     i_start, i_end, j_start, j_end,                                           &
     k_start, k_end, k_loop_start, k_loop_end,                                 &
     l_aerosol_droplet, l_nh42so4,                                             &
     accum_sulphate, diss_sulphate,                                            &
     l_seasalt_ccn, sea_salt_film, sea_salt_jet,                               &
     l_biogenic_ccn, biogenic,                                                 &
     l_biomass_ccn, biomass_aged, biomass_cloud,                               &
     l_ocff_ccn, ocff_aged, ocff_cloud,                                        &
     l_nitrate_ccn,                                                            &
     density_air,                                                              &
     snow_depth,                                                               &
     land_fract,                                                               &
     n_drop,                                                                   &
     nitr_acc, nitr_diss, aerosol                                              &
   )

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
use mphys_constants_mod, only: ntot_land, ntot_sea, m0_murk, n0_murk
use lsp_autoc_consts_mod, only: power_murk
use murk_inputs_mod, only: l_murk_vis

implicit none

integer, intent(in) :: i_start      ! i indexing start
integer, intent(in) :: i_end        ! i indexing end
integer, intent(in) :: j_start      ! j indexing start
integer, intent(in) :: j_end        ! j indexing end
integer, intent(in) :: k_start      ! k indexing start
integer, intent(in) :: k_end        ! k indexing end
integer, intent(in) :: k_loop_start ! k looping start
integer, intent(in) :: k_loop_end   ! k looping end

!             Flag to use aerosols to find droplet number
logical, intent(in) :: l_aerosol_droplet
!             Is the input "sulphate" aerosol in the form of
!             ammonium sulphate (T) or just sulphur (F)? Used in the
!             same way for the nitrate aerosol, in form of ammonium
!             nitrate (T) or just nitrogen (F).
logical, intent(in) :: l_nh42so4
!             Is sea-salt aerosol to be used?
logical, intent(in) :: l_seasalt_ccn
!             Is biomass smoke aerosol to be used?
logical, intent(in) :: l_biomass_ccn
!             Is biogenic aerosol to be used?
logical, intent(in) :: l_biogenic_ccn
!             Is fossil-fuel organic carbon aerosol to be used?
logical, intent(in) :: l_ocff_ccn
!             Is ammonium nitrate aerosol to be used?
logical, intent(in) :: l_nitrate_ccn

!             Mixing ratio of accumulation-mode sulphate aerosol
real(kind=real_umphys), intent(in) ::                                          &
                   accum_sulphate ( i_start : i_end,                           &
                                    j_start : j_end,                           &
                                    k_start : k_end )

!             Mixing ratio of dissolved sulphate aerosol
real(kind=real_umphys), intent(in) ::                                          &
                    diss_sulphate ( i_start : i_end,                           &
                                    j_start : j_end,                           &
                                    k_start : k_end )

!             Mixing ratio of aged biomass smoke
real(kind=real_umphys), intent(in) ::                                          &
                     biomass_aged ( i_start : i_end,                           &
                                    j_start : j_end,                           &
                                    k_start : k_end )

!             Mixing ratio of in-cloud biomass smoke
real(kind=real_umphys), intent(in) ::                                          &
                     biomass_cloud( i_start : i_end,                           &
                                    j_start : j_end,                           &
                                    k_start : k_end )

!             Number concentration of film-mode sea salt aerosol (m-3)
real(kind=real_umphys), intent(in) ::                                          &
                     sea_salt_film( i_start : i_end,                           &
                                    j_start : j_end,                           &
                                    k_start : k_end )

!             Number concentration of jet-mode sea salt aerosol (m-3)
real(kind=real_umphys), intent(in) ::                                          &
                     sea_salt_jet ( i_start : i_end,                           &
                                    j_start : j_end,                           &
                                    k_start : k_end )

!             Mixing ratio of biogenic aerosol
real(kind=real_umphys), intent(in) ::                                          &
                     biogenic     ( i_start : i_end,                           &
                                    j_start : j_end,                           &
                                    k_start : k_end )

!             Mixing ratio of aged fossil-fuel organic carbon
real(kind=real_umphys), intent(in) ::                                          &
                     ocff_aged    ( i_start : i_end,                           &
                                    j_start : j_end,                           &
                                    k_start : k_end )

!             Mixing ratio of in-cloud fossil-fuel organic carbon
real(kind=real_umphys), intent(in) ::                                          &
                     ocff_cloud   ( i_start : i_end,                           &
                                    j_start : j_end,                           &
                                    k_start : k_end )

!             Density of air (kg m-3)
real(kind=real_umphys), intent(in) ::                                          &
                     density_air  ( i_start : i_end,                           &
                                    j_start : j_end,                           &
                                    k_start : k_end )

!             Snow depth (m; >5000 is flag for ice-sheets)
real(kind=real_umphys), intent(in) ::                                          &
                     snow_depth   ( i_start : i_end,                           &
                                    j_start : j_end)

!             Land fraction
real(kind=real_umphys), intent(in) ::                                          &
                     land_fract   ( i_start : i_end,                           &
                                    j_start : j_end)


!             Returned number concentration of cloud droplets (m-3)
real(kind=real_umphys), intent(out) ::                                         &
                             n_drop( i_start : i_end,                          &
                                     j_start : j_end,                          &
                                     k_start : k_end )

#if defined(IFORT_VERSION) && (IFORT_VERSION < 17000007)
! Avoid an Intel compiler bug.
! Affects the specific case where:
!  * pointer actual argument is passed through to an optional dummy argument;
!  * the dummy argument is not assumed shape and is non-pointer;
!  * the flag `-fpscomp logicals` is set  (implied by -standard-semantics).
! Workaround is to make the dummy argument assumed shape.

!             Mixing ratio of accumulation-mode nitrate aerosol
real(kind=real_umphys), intent(in), optional ::                                &
                     nitr_acc     ( i_start : , j_start : , k_start : )

!             Mixing ratio of dissolved-mode nitrate aerosol
real(kind=real_umphys), intent(in), optional ::                                &
                     nitr_diss    ( i_start : , j_start : , k_start : )

#else

!             Mixing ratio of accumulation-mode nitrate aerosol
real(kind=real_umphys), intent(in), optional ::                                &
                     nitr_acc     ( i_start : i_end,                           &
                                    j_start : j_end,                           &
                                    k_start : k_end )

!             Mixing ratio of dissolved-mode nitrate aerosol
real(kind=real_umphys), intent(in), optional ::                                &
                     nitr_diss    ( i_start : i_end,                           &
                                    j_start : j_end,                           &
                                    k_start : k_end )
#endif

! "Murk" aerosol mass mixing ratio for use in visibility calculations
real(kind=real_umphys), intent(in out), optional :: aerosol(i_start : i_end,   &
                                                            j_start : j_end,   &
                                                            k_start : k_end )
!     Local variables:

integer :: i       ! Looper
integer :: j       ! Looper
integer :: k       ! Looper

!  Number density of CCN
real(kind=real_umphys) :: n_ccn ( i_start : i_end,                             &
                                  j_start : j_end,                             &
                                  k_start : k_end )


integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='NUMBER_DROPLET'

!---------------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Note the following call needs to be run by a single thread at a time
! to avoid race/deadlock conditions as this routine is sometimes called in
! a parallel context.
!$OMP CRITICAL
if (.not. nd_init) call number_droplet_init()
!$OMP end CRITICAL

if (l_aerosol_droplet) then
!$OMP PARALLEL DEFAULT(SHARED) private(i, j, k)
!$OMP do SCHEDULE(STATIC)
  do k = k_loop_start, k_loop_end
    do j = j_start, j_end
      do i = i_start, i_end

        !        If active, aerosol concentrations are used to calculate the
        !        number of CCN, which is then used to determine the number
        !        concentration of cloud droplets (m-3).

        if (l_nh42so4) then
          !           Input data have already been converted to ammonium sulphate.
          n_ccn(i,j,k)=(accum_sulphate(i,j,k)+diss_sulphate(i,j,k))            &
             *density_air(i,j,k)/(density_nh42so4*particle_volume_nh42so4)
        else
          !           Convert m.m.r. of sulphur to ammonium sulphate by
          !           multiplying by ratio of molecular weights:
          n_ccn(i,j,k)=(accum_sulphate(i,j,k)+diss_sulphate(i,j,k))*4.125      &
             *density_air(i,j,k)/(density_nh42so4*particle_volume_nh42so4)
        end if

        if (l_seasalt_ccn) then
          n_ccn(i,j,k)=n_ccn(i,j,k)+sea_salt_film(i,j,k)+sea_salt_jet(i,j,k)
        end if

        if (l_biomass_ccn) then
          n_ccn(i,j,k)=n_ccn(i,j,k)+((biomass_aged(i,j,k)+biomass_cloud(i,j,k))&
                       *density_air(i,j,k)                                     &
                       /(density_biomass*particle_volume_biomass))
        end if

        if (l_biogenic_ccn) then
          n_ccn(i,j,k)=n_ccn(i,j,k)+(biogenic(i,j,k)*density_air(i,j,k)        &
                       /(density_biogenic*particle_volume_biogenic))
        end if

        if (l_ocff_ccn) then
          n_ccn(i,j,k)=n_ccn(i,j,k)+((ocff_aged(i,j,k)+ocff_cloud(i,j,k))      &
                      *density_air(i,j,k)/(density_ocff*particle_volume_ocff))
        end if

        if (l_nitrate_ccn  .and. present(nitr_acc)                             &
                           .and. present(nitr_diss)) then
          if (l_nh42so4) then
            !             Input data have already been converted to ammonium nitrate
            n_ccn(i,j,k)=n_ccn(i,j,k) + ((nitr_acc(i,j,k) + nitr_diss(i,j,k))  &
              * density_air(i,j,k)/(density_nitrate*particle_volume_nitrate))
          else
            !             Convert m.m.r. of nitrogen to ammonium nitrate by
            !             multiplying by ratio of molecular weights
            n_ccn(i,j,k)=n_ccn(i,j,k) + ((nitr_acc(i,j,k) + nitr_diss(i,j,k))  &
                         * 5.714 * density_air(i,j,k)                          &
                         / (density_nitrate*particle_volume_nitrate))
          end if
        end if

        !        Apply relation of Jones et al. (1994) to get droplet number
        !        and apply minimum value (equivalent to 5 cm-3):

        n_drop(i,j,k)=3.75e+08*(1.0e+00-exp(-2.5e-9*n_ccn(i,j,k)))

        if ( n_drop(i,j,k)  <   min_cdnc_sea_ice ) then
          n_drop(i,j,k) = min_cdnc_sea_ice
        end if

        ! If gridbox is more than 20% land and this land is not covered
        ! by an ice-sheet, use larger minimum droplet number (=35 cm-3):

        if (land_fract(i,j)  >   land_frac_thresh  .and.                       &
           snow_depth(i,j)  <   snow_depth_thresh .and.                        &
           n_drop(i,j,k)  <  min_cdnc_land ) then
          n_drop(i,j,k) = min_cdnc_land
        end if

      end do ! i
    end do ! j
  end do ! k
!$OMP end do NOWAIT

  if (l_murk_vis .and. present(aerosol)) then
!$OMP do SCHEDULE(STATIC)
    do k = k_loop_start, k_loop_end
      do j = j_start, j_end
        do i = i_start, i_end
          ! Invert relationship from visibility scheme between aerosol mass
          ! and aerosol number (n_ccn)
          ! The factor 2 is a scaling factor to get realistic aerosol masses
          ! in a high resolution model from a low resolution climatology
          aerosol(i,j,k) = ((2.0*n_ccn(i,j,k) / n0_murk)**(1.0/power_murk))    &
                         * m0_murk * 1.0e9
        end do ! i
      end do ! j
    end do ! k
!$OMP end do
  end if
!$OMP end PARALLEL

else  ! l_aerosol_droplet

  !        Without aerosols, the number of droplets is fixed; a simple
  !        50% criterion is used for land or sea in this case.
  do k = k_loop_start, k_loop_end
    do j = j_start, j_end
      do i = i_start, i_end

        if (land_fract(i,j)  >=  0.5) then
          n_drop(i,j,k)=ntot_land
        else
          n_drop(i,j,k)=ntot_sea
        end if
      end do
    end do
  end do

end if ! l_aerosol_droplet

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine number_droplet

end module number_droplet_mod
