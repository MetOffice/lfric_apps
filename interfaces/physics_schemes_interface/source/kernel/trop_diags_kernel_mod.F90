!-----------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Tropopause aviation diagnostics

module trop_diags_kernel_mod

use argument_mod,      only : arg_type,                        &
                              GH_FIELD, GH_REAL, GH_SCALAR,    &
                              GH_READ, GH_WRITE, CELL_COLUMN,  &
                              ANY_DISCONTINUOUS_SPACE_1
use fs_continuity_mod, only : Wtheta
use constants_mod,     only : r_def, i_def, l_def
use kernel_mod,        only : kernel_type

implicit none

private

! layered input fields before flat output fields for nlayers > 1
type, public, extends(kernel_type) :: trop_diags_kernel_type
  private
  type(arg_type) :: meta_args(10) = (/                                  &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    & ! theta_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    & ! exner_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    & ! height_wth
       arg_type(GH_SCALAR, GH_REAL, GH_READ),                            & ! g_over_r
       arg_type(GH_SCALAR, GH_REAL, GH_READ),                            & ! p_zero
       arg_type(GH_SCALAR, GH_REAL, GH_READ),                            & ! kappa
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! trop_pres
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! trop_temp
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! trop_height
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1)  & ! trop_icao_height
    /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: trop_diags_code
end type

public :: trop_diags_code

contains

!> @brief Calculate tropopause aviation diagnostics (pressure, temperature and height)
!> @details
!> Find interval (k to k+1) in which tropopause lies by checking the lapse
!> rate of surrounding intervals, and then calculate the corresponding
!> height, pressure and temperature at the tropopause by using the lapse rates
!> above and below the interval identified.
!>
!> If an input field is not empty, other input fields upon which it depends
!> are guaranteed non-empty by the caller.
!>
!> If the tropopause is not found, outputs are set to missing data.
!>
!> Method:
!> 1) Loop over model_levels
!>
!>    i)  Calculate lapse rate field for interval k~k+1
!>        and the lapse rate for the interval below k-1~k
!>
!>    ii) Where the lapse rate of k~k+1 < 2C and the lapse rate of the interval
!>         below is +ve  -- we may have found the tropopause interval -- search
!>         upwards for the next model level that at least 2km above k (k2km)
!>         When (if found) check that the lapse rate k~k2km < 2C.
!>         If it does then k is the bottom model level of the WMO tropopause.
!>         Otherwise we cycle around again looking at the next k level.
!>
!> 2) Having found the model level k that is close to the bottom of the
!>    tropopause, we shall then use the lapse rates of the intervals above and
!>    below the k level to calculate the tropopause height at the intersection
!>    of these lapse rate lines.
!>
!>       Equation of first line:
!>         in form y=mx+c where m=-lapselwr and c = templwr + lapselwr*heightlwr
!>       i)  TropTemp = - lapselwr*Tropheight + (templwr + lapselwr*heightlwr)
!>       Equation of second line:
!>         in form y=mx+c where m=-LapseUpr and c = TempUpr + LapseUpr*HeightUpr
!>       ii) TropTemp = - lapseupr*Tropheight + (tempupr + lapseupr*heightupr)
!>
!>       Point of intersection is given by:
!>         - lapseupr*tropheight + (tempupr + lapseupr*heightupr) =
!>                         - lapselwr*tropheight + (templwr + lapselwr*heightlwr)
!>       => (lapselwr - lapseypr)*tropheight =
!>               (templwr + lapselwr*heightlwr) - (tempupr + lapseupr*heightupr)
!>       => TropHeight = (templwr + lapselwr*heightlwr) -
!>               (tempupr + lapseupr*heightupr)/(lapselwr - lapseupr)
!>
!>     iii) Calculate tropopause temperature by substuting trop height into
!>          equation of first line:
!>          Troptemp = templwr - lapselwr * (Tropheight - heightlwr)
!>
!>     iv) Calculate tropopause pressure:
!>          TropPress = Presslwr * (TropTemp/templwr)**(g_over_r/lapselwr)
!>
!> @param[in]  nlayers                Number of cells in a column
!> @param[in]  theta_wth              Potential temperature
!> @param[in]  exner_wth              Exner pressure in wth space
!> @param[in]  height_wth             Height of wth levels above surface
!> @param[in]  g_over_r               Gravity / specific dry air gas constant
!> @param[in]  p_zero                 Reference surface pressure
!> @param[in]  kappa                  Ratio R/cp (dry air)
!> @param[out] trop_pres              Pressure at the tropopause
!> @param[out] trop_temp              Temperature at the tropopause
!> @param[out] trop_height            Height of the tropopause above surface
!> @param[out] trop_icao_height       ICAO height of the tropopause
!> @param[in]  ndf_wth                No. DOFs per cell for wth space
!> @param[in]  undf_wth               No. unique DOFs for wth space
!> @param[in]  map_wth                Dofmap for wth space column base cell
!> @param[in]  ndf_2d                 No. DOFs per cell for 2D space
!> @param[in]  undf_2d                No. unique DOFs for 2D space
!> @param[in]  map_2d                 Dofmap for 2D space column base cell
subroutine trop_diags_code(nlayers,                    &
                           theta_wth,                  &
                           exner_wth,                  &
                           height_wth,                 &
                           g_over_r,                   &
                           p_zero,                     &
                           kappa,                      &
                           trop_pres,                  &
                           trop_temp,                  &
                           trop_height,                &
                           trop_icao_height,           &
                           ndf_wth, undf_wth, map_wth, &
                           ndf_2d, undf_2d, map_2d)

  use missing_data_mod,         only : rmdi, imdi
  use icao_heights_kernel_mod,  only : icao_heights_kernel_code
  use empty_data_mod,           only : empty_real_data

  implicit none

  ! Arguments (kernel)
  integer(i_def), intent(in) :: nlayers
  integer(i_def), intent(in) :: ndf_2d, ndf_wth
  integer(i_def), intent(in) :: undf_2d, undf_wth

  integer(i_def), dimension(ndf_2d),  intent(in) :: map_2d
  integer(i_def), dimension(ndf_wth), intent(in) :: map_wth

  ! Arguments (algorithm)
  real(r_def), dimension(undf_wth), intent(in) :: theta_wth
  real(r_def), dimension(undf_wth), intent(in) :: exner_wth
  real(r_def), dimension(undf_wth), intent(in) :: height_wth
  real(r_def), intent(in)                      :: g_over_r
  real(r_def), intent(in)                      :: p_zero
  real(r_def), intent(in)                      :: kappa

  real(r_def), pointer, dimension(:), intent(inout) :: trop_pres
  real(r_def), pointer, dimension(:), intent(inout) :: trop_temp
  real(r_def), pointer, dimension(:), intent(inout) :: trop_height
  real(r_def), pointer, dimension(:), intent(inout) :: trop_icao_height


  ! Local variables
  integer(i_def) :: k, k2km
  integer(i_def) :: trop_level
  real(r_def) :: t_wth(nlayers)
  real(r_def) :: lapse, lapse_below, lapse_above, lapse_2km
  real(r_def) :: delta_lapse, pres_wth

  ! Parameters for WMO tropopause definition
  real(r_def), parameter :: lapse_trop = 0.002_r_def   ! WMO 2K/km threshold, in K/m
  real(r_def), parameter :: dz_trop = 2000.0_r_def     ! m

  ! Small number used to protect divisions from near-zero denominators
  real(r_def), parameter :: vsmall = 1.0e-9_r_def

  ! cut off limits to be used in tropopause calculations.
  ! todo: put this in a constants module? in the um it was in pws_diags_mod
  ! arbritary limits for high and low trop levels for search
  real(r_def), parameter :: heightcut_top = 22000.0_r_def
  real(r_def), parameter :: heightcut_bot = 4500.0_r_def
  ! max temp allowed for tropopause
  real(r_def), parameter :: tempcut = 243.0_r_def

  logical(l_def) :: do_pres, do_temp, do_height, do_icao_height


  ! A field has non-empty data if it was requested or is needed as a dependency.
  do_pres        = .not. associated(trop_pres,        empty_real_data)
  do_temp        = .not. associated(trop_temp,        empty_real_data)
  do_height      = .not. associated(trop_height,      empty_real_data)
  do_icao_height = .not. associated(trop_icao_height, empty_real_data)

  ! Initialise all outputs to missing data. They are only meaningful when a
  ! lapse-rate (WMO) tropopause has been located.
  ! todo: do this to fields which are initialised (not just the requested)
  if (do_pres)        trop_pres(map_2d(1))        = rmdi
  if (do_temp)        trop_temp(map_2d(1))        = rmdi
  if (do_height)      trop_height(map_2d(1))      = rmdi
  if (do_icao_height) trop_icao_height(map_2d(1)) = rmdi

  ! Temperature on wth levels
  ! todo: put into following loop
  do k=1, nlayers
    t_wth(k) = theta_wth(map_wth(1)+k) * exner_wth(map_wth(1)+k)
  end do


  ! Locate the lapse-rate (WMO) tropopause
  trop_level = imdi
  ! Note: The UM looped from 1 to nlayers, explaining that it wouldn't go
  !       out of bounds due to the height range search criteria.
  !       Here, we limit the loop for robustness, e.g for edge-case testing.
  do k=3, nlayers-2

    if (t_wth(k) < tempcut .and.             &
        height_wth(map_wth(1)+k) > heightcut_bot .and.  &
        height_wth(map_wth(1)+k) < heightcut_top) then

      lapse = (t_wth(k) - t_wth(k+1)) / (height_wth(map_wth(1)+k+1) - height_wth(map_wth(1)+k))
      lapse_below = (t_wth(k-1) - t_wth(k)) / (height_wth(map_wth(1)+k) - height_wth(map_wth(1)+k-1))

      if (lapse < lapse_trop .and. lapse_below > 0.0_r_def) then
        ! Lapse rate has dropped below the threshold. If this is maintained
        ! for 2km above then the WMO criteria for the tropopause has been met.
        do k2km=k, nlayers
          if (height_wth(map_wth(1)+k2km) > heightcut_top) exit

          if ((height_wth(map_wth(1)+k2km)-height_wth(map_wth(1)+k)) >= dz_trop ) then
            lapse_2km = (t_wth(k) - t_wth(k2km)) / (height_wth(map_wth(1)+k2km) - height_wth(map_wth(1)+k))
            ! if 2km interval also < 2 then we have the tropopause level
            if (lapse_2km < lapse_trop) then
              trop_level = k
              exit
            end if
          end if

        end do  ! looking upwards for 2km
      end if  ! lapse values in range

    end if

    ! did we find the tropopause?
    if (trop_level /= imdi) exit

  end do


  ! if level found
  if (trop_level /= imdi) then
    k = trop_level

    ! todo: lapse_below was already calculated in the loop above
    lapse_above = (t_wth(k+1) - t_wth(k+2)) / (height_wth(map_wth(1)+k+2) - height_wth(map_wth(1)+k+1))
    lapse_below = (t_wth(k-1) - t_wth(k)) / (height_wth(map_wth(1)+k) - height_wth(map_wth(1)+k-1))
    delta_lapse = lapse_below - lapse_above
    if ( abs(delta_lapse) < vsmall ) then
      if ( delta_lapse >= 0.0_r_def ) delta_lapse =  vsmall
      if ( delta_lapse <  0.0_r_def ) delta_lapse = -vsmall
    end if

    ! height of tropopause between k and k+1
    if (do_height) then
      trop_height(map_2d(1)) = (                                     &
          (t_wth(k)   + (lapse_below * height_wth(map_wth(1)+k)))  &
        - (t_wth(k+1) + (lapse_above * height_wth(map_wth(1)+k+1)))    &
        ) / delta_lapse

      ! ensure trop level doesn't undershoot
      if (trop_height(map_2d(1)) < height_wth(map_wth(1)+k)) then
        trop_height(map_2d(1)) = height_wth(map_wth(1)+k)
      end if
      ! or overshoot
      if (trop_height(map_2d(1)) > height_wth(map_wth(1)+k+1)) then
        trop_height(map_2d(1)) = height_wth(map_wth(1)+k+1)
      end if

    end if

    ! temperature at tropopause
    if (do_temp) then
      trop_temp(map_2d(1)) = t_wth(k) -                                     &
        lapse_below * (trop_height(map_2d(1)) - height_wth(map_wth(1)+k-1))
    end if

    ! pressure at tropopause is derived from the hydrostatic equation
    if ( abs(lapse_below) < vsmall ) then
      if ( lapse_below >= 0.0_r_def ) lapse_below =  vsmall
      if ( lapse_below <  0.0_r_def ) lapse_below = -vsmall
    end if
    if (do_pres) then
      pres_wth = p_zero * exner_wth(map_wth(1)+k) ** (1.0_r_def / kappa)
      trop_pres(map_2d(1)) = pres_wth *                           &
        (trop_temp(map_2d(1)) / t_wth(k)) ** (g_over_r / lapse_below)
    end if

    ! ICAO height of the tropopause
    if (do_icao_height) then
      call icao_heights_kernel_code( &
        nlayers, trop_icao_height, trop_pres, g_over_r, &
        ndf_2d, undf_2d, map_2d)
    end if

  end if

end subroutine trop_diags_code

end module trop_diags_kernel_mod

