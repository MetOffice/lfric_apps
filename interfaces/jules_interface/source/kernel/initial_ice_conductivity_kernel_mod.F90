!-------------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Initialise Jules surface fields on tiles
!> @details Non-standard Surface fields (pseudo-levels) aren't as yet not
!>  implemented in LFRic. As an interim measure Higher-order W3 fields have
!>  been used to mimic psuedo-level field behaviour. This code is written
!>  based on this interim measure and will need to be updated when
!>  suitable infrastructure is available (Ticket #2081)
module initial_ice_conductivity_kernel_mod
  use argument_mod,  only: arg_type,                       &
                           GH_FIELD, GH_REAL, GH_INTEGER,  &
                           GH_WRITE, CELL_COLUMN,          &
                           ANY_DISCONTINUOUS_SPACE_1,      &
                           ANY_DISCONTINUOUS_SPACE_2
  use constants_mod, only: r_def, i_def
  use kernel_mod,    only: kernel_type
  use jules_control_init_mod, only: n_sea_ice_tile, first_sea_ice_tile
  use jules_sea_seaice_config_mod,     only: therm_cond_sice => kappai
  implicit none
  private
  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: initial_ice_conductivity_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                    &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1) &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: initial_ice_conductivity_code
  end type initial_ice_conductivity_kernel_type
  public :: initial_ice_conductivity_code
contains
  !> @param[in]     nlayers              The number of layers
  !> @param[in,out] sea_ice_fraction     Sea Ice Fractions on categories
  !> @param[in,out] sea_ice_thickness    Sea Ice Thickness on categories
  !> @param[in,out] sea_ice_conductivity Sea Ice Conductivity on categories
  !> @param[in]     ndf_ice              Number of DOFs per cell for tiles
  !> @param[in]     undf_ice             Number of total DOFs for tiles
  !> @param[in]     map_ice              Dofmap for cell for surface tiles
  subroutine initial_ice_conductivity_code(nlayers,                           &
                                           sea_ice_fraction,                  &
                                           sea_ice_thickness,                 &
                                           sea_ice_conductivity,              &
                                           ndf_ice, undf_ice, map_ice)
    implicit none
    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_ice, undf_ice
    integer(kind=i_def), intent(in) :: map_ice(ndf_ice)
    real(kind=r_def), intent(inout) :: sea_ice_fraction(undf_ice)
    real(kind=r_def), intent(inout) :: sea_ice_thickness(undf_ice)
    real(kind=r_def), intent(inout) :: sea_ice_conductivity(undf_ice)
    ! Internal variables
    integer(kind=i_def) :: i
    real(kind=r_def) :: min_ice_thick, max_ice_cond
    !Taken from UM recon value
    max_ice_cond = 25.0_r_def
    min_ice_thick = 8.0_r_def * therm_cond_sice/max_ice_cond
    do i=0,n_sea_ice_tile-1
        if (sea_ice_thickness(map_ice(1)+i) >= min_ice_thick) then
            sea_ice_conductivity(map_ice(1)+i) = (8.0_r_def * therm_cond_sice) &
                                              / sea_ice_thickness(map_ice(1)+i)
        else
            sea_ice_conductivity(map_ice(1)+i) = max_ice_cond
        endif
    end do
  end subroutine initial_ice_conductivity_code
end module initial_ice_conductivity_kernel_mod
