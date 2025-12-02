!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Process ocean and sea ice data coming through the coupler

module process_o2a_kernel_mod

use kernel_mod,                only : kernel_type
use argument_mod,              only : arg_type,              &
                                      GH_FIELD, GH_REAL,     &
                                      GH_INTEGER,            &
                                      GH_READ, GH_READWRITE, &
                                      GH_SCALAR,             &
                                      GH_WRITE, ANY_DISCONTINUOUS_SPACE_1, &
                                      ANY_DISCONTINUOUS_SPACE_2, &
                                      CELL_COLUMN
use constants_mod,             only : r_def, i_def
use sci_constants_mod,         only : zero_C_in_K
use section_choice_config_mod, only : iau_sst
use derived_config_mod,        only : l_couple_ocean, l_couple_sea_ice

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: process_o2a_kernel_type
  private
  type(arg_type) :: meta_args(14) = (/                    &
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), &
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), &
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), &
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), &
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), &
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), &
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), &
       arg_type(GH_FIELD, GH_INTEGER, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                             &
       arg_type(GH_SCALAR, GH_REAL, GH_READ),                                &
       arg_type(GH_SCALAR, GH_REAL, GH_READ),                                &
       arg_type(GH_SCALAR, GH_REAL, GH_READ)                                 &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: process_o2a_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: process_o2a_code
contains

!> @brief Process ocean and sea ice data coming through the coupler
!! @param[in] nlayers Number of layers
!! @param[in,out] sea_surf_temp        the sea surface temperature field
!! @param[in]     sea_surf_temp_pert   the sea surface temperature perturbations
!! @param[in,out] sea_ice_fraction     the sea ice fraction field
!! @param[in,out] sea_ice_thickness    the sea ice thickness field
!! @param[in,out] sea_ice_layer_t      the sea ice layer temperature field
!! @param[in,out] sea_ice_conductivity the sea ice thermal conductivity
!! @param[in,out] sea_ice_snow_depth   the amount of snow on sea ice
!! @param[in,out] melt_pond_fraction   the fraction of melt ponds on the sea ice
!! @param[in,out] melt_pond_depth      the depth of the melt ponds on the sea ice
!! @param[in,out] ocn_cpl_point        the mask to determine whether coupled
!! @param[in] n_sea_ice_tile Number of sea ice tiles
!! @param[in] T_freeze_h2o_sea Temperature at which sea water freezes, [K]
!! @param[in] therm_cond_sice Thermal conductivity of sea-ice (W/m/K)
!! @param[in] therm_cond_sice_snow Thermal conductivity of snow on zero layer sea-ice (W/m/K)
!! @param[in] ndf_ocn Number of degrees of freedom per cell for the ocean surface
!! @param[in] undf_ocn Number of unique degrees of freedom  for the ocean surface
!! @param[in] map_ocn Dofmap for the cell at the base of the column for the ocean surface
!! @param[in] ndf_ice Number of degrees of freedom per cell for the sea ice
!! @param[in] undf_ice Number of unique degrees of freedom  for the sea ice
!! @param[in] map_ice Dofmap for the cell at the base of the column for the sea ice


subroutine process_o2a_code(nlayers, sea_surf_temp, sea_surf_temp_pert, &
                              sea_ice_fraction, sea_ice_thickness,      &
                              sea_ice_layer_t, sea_ice_conductivity,    &
                              sea_ice_snow_depth,                       &
                              melt_pond_fraction, melt_pond_depth,      &
                              ocn_cpl_point,                            &
                              n_sea_ice_tile,                           &
                              T_freeze_h2o_sea,                         &
                              therm_cond_sice,                          &
                              therm_cond_sice_snow,                     &
                              ndf_ocn, undf_ocn, map_ocn,               &
                              ndf_ice, undf_ice, map_ice)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_ocn
  integer(kind=i_def), intent(in) :: undf_ocn
  integer(kind=i_def), dimension(ndf_ocn),  intent(in) :: map_ocn
  integer(kind=i_def), intent(in) :: ndf_ice
  integer(kind=i_def), intent(in) :: undf_ice
  integer(kind=i_def), dimension(ndf_ice),  intent(in) :: map_ice
  real(kind=r_def), dimension(undf_ocn), intent(inout) :: sea_surf_temp
  real(kind=r_def), dimension(undf_ocn), intent(in)    :: sea_surf_temp_pert
  real(kind=r_def), dimension(undf_ice), intent(inout) :: sea_ice_fraction
  real(kind=r_def), dimension(undf_ice), intent(inout) :: sea_ice_thickness
  real(kind=r_def), dimension(undf_ice), intent(inout) :: sea_ice_layer_t
  real(kind=r_def), dimension(undf_ice), intent(inout) :: sea_ice_conductivity
  real(kind=r_def), dimension(undf_ice), intent(inout) :: sea_ice_snow_depth
  real(kind=r_def), dimension(undf_ice), intent(inout) :: melt_pond_fraction
  real(kind=r_def), dimension(undf_ice), intent(inout) :: melt_pond_depth
  integer(kind=i_def), dimension(undf_ocn), intent(inout) :: ocn_cpl_point
  integer(kind=i_def), intent(in) :: n_sea_ice_tile
  real(kind=r_def), intent(in) :: T_freeze_h2o_sea
  real(kind=r_def), intent(in) :: therm_cond_sice
  real(kind=r_def), intent(in) :: therm_cond_sice_snow

  ! Parameters
  ! Minimum sea ice thickness allowed through the coupler
  real(kind=r_def), parameter :: hi_min = 0.05_r_def
  ! Minimum sea ice fraction allowed through the coupler
  real(kind=r_def), parameter :: aicenmin = 1.0e-02_r_def
  ! Maximum sea ice layer temperature allowed through the coupler
  real(kind=r_def), parameter :: ti_max = zero_C_in_K
  ! Minimum sea ice layer temperature allowed through the coupler
  real(kind=r_def), parameter :: ti_min = 200.0_r_def

  ! Local variables
  integer(kind=i_def) :: i
  real(kind=r_def)    :: total_sea_ice_fraction
  real(kind=r_def)    :: recip_sea_ice_fraction

  if (l_couple_ocean) then
    if (sea_surf_temp(map_ocn(1))>1.0_r_def) then
      ocn_cpl_point(map_ocn(1))=1_i_def
    else
      ocn_cpl_point(map_ocn(1))=0_i_def
    end if
  end if

  total_sea_ice_fraction = 0.0_r_def

  if (l_couple_sea_ice) then
    do i = 0, n_sea_ice_tile - 1

        ! Make sure sea ice is thicker than hi_min
        if ( sea_ice_fraction(map_ice(1) + i) > 0.0_r_def ) then
          if ( sea_ice_thickness(map_ice(1) + i)/sea_ice_fraction(map_ice(1)+i)  &
                                                                  < hi_min ) then
            sea_ice_fraction(map_ice(1) + i) = 0.0_r_def
          end if
        end if

        if ( sea_ice_fraction(map_ice(1) + i)  < aicenmin ) then
          sea_ice_fraction(map_ice(1) + i) = 0.0_r_def
          sea_ice_thickness(map_ice(1) + i) = 0.0_r_def
          sea_ice_snow_depth(map_ice(1) + i) = 0.0_r_def
          sea_ice_layer_t(map_ice(1) + i) = 0.0_r_def
          sea_ice_conductivity(map_ice(1) + i) = 0.0_r_def
          melt_pond_fraction(map_ice(1) + i) = 0.0_r_def
          melt_pond_depth(map_ice(1) + i) = 0.0_r_def
        end if

        total_sea_ice_fraction = total_sea_ice_fraction                          &
                                              + sea_ice_fraction(map_ice(1) + i)

    end do
  end if

  ! Add perturbations to sea surface temperature if they have been set
  ! in the IAU
  if ( iau_sst ) then
    sea_surf_temp(map_ocn(1)) = sea_surf_temp(map_ocn(1)) &
                                       + sea_surf_temp_pert(map_ocn(1))
  end if

  ! Set the sea surface temperature to be at the freezing point of water
  ! when there is sea ice.
  if (l_couple_ocean) then
    if ( total_sea_ice_fraction > 0.0_r_def ) then
      sea_surf_temp(map_ocn(1)) = T_freeze_h2o_sea
    end if
  end if

  if (l_couple_sea_ice) then
    do i = 0, n_sea_ice_tile - 1

        ! Reduce category ice fractions if aggregate > 1
        if (total_sea_ice_fraction > 1.0_r_def) then
          sea_ice_fraction(map_ice(1) + i) = sea_ice_fraction(map_ice(1) + i)    &
                                                        / total_sea_ice_fraction
        end if

        ! Increase sea ice thickness due to overlying snow
        sea_ice_thickness(map_ice(1) + i) = sea_ice_thickness(map_ice(1) + i) +  &
                                        sea_ice_snow_depth(map_ice(1) + i)*      &
                                        (therm_cond_sice/therm_cond_sice_snow)

        ! Undo sea ice fraction scaling done on ocean side of coupler
        if ( sea_ice_fraction(map_ice(1) + i) > 0.0_r_def ) then
          recip_sea_ice_fraction = 1.0_r_def / sea_ice_fraction(map_ice(1) + i)
          sea_ice_thickness(map_ice(1) + i) = sea_ice_thickness(map_ice(1) + i)  &
                                              * recip_sea_ice_fraction
          sea_ice_snow_depth(map_ice(1)+i) = sea_ice_snow_depth(map_ice(1) + i)  &
                                              * recip_sea_ice_fraction
          sea_ice_layer_t(map_ice(1) + i) = sea_ice_layer_t(map_ice(1) + i)      &
                                              * recip_sea_ice_fraction
          sea_ice_conductivity(map_ice(1)+i) = sea_ice_conductivity(map_ice(1)+i)&
                                              * recip_sea_ice_fraction
        end if

        ! Apply bounds to sea ice layer temperatures
        if ( sea_ice_layer_t(map_ice(1) + i) > ti_max ) then
          sea_ice_layer_t(map_ice(1) + i) = ti_max
        end if
        if ( sea_ice_layer_t(map_ice(1) + i) < ti_min .and.                      &
                             sea_ice_layer_t(map_ice(1) + i) /= 0.0_r_def ) then
          sea_ice_layer_t(map_ice(1) + i) = ti_min
        end if

        ! Apply bounds to melt ponds
        if ( melt_pond_fraction(map_ice(1) + i) <= 0.0_r_def ) then
          melt_pond_fraction(map_ice(1) + i) = 0.0_r_def
          melt_pond_depth(map_ice(1) + i) = 0.0_r_def
        end if
        if ( melt_pond_fraction(map_ice(1) + i) > 1.0_r_def ) then
          melt_pond_fraction(map_ice(1) + i) = 1.0_r_def
        end if
        if ( melt_pond_depth(map_ice(1) + i) < 0.0_r_def ) then
          melt_pond_depth(map_ice(1) + i) = 0.0_r_def
        end if
    end do
  end if

end subroutine process_o2a_code

end module process_o2a_kernel_mod
