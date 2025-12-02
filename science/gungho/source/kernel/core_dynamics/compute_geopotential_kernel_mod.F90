!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Computes the geopotential field for shallow (constant gradient with height)
!> and deep (varying gradient with height) geometries. The geopotential is
!> normalised so that it is zero at the domain surface (radius = planet_radius).
!> For shallow geometries, this gives geopotential = gravity*(radius - planet_radius)
!> For deep geometries, this gives geopotential = gravity*planet_radius/radius*(radius - planet_radius)
!> and for Cartesian geometries the simpler form geopotential = gravity*height
!> is used.
!>
module compute_geopotential_kernel_mod

  use argument_mod,              only: arg_type, func_type,    &
                                       GH_FIELD, GH_REAL,      &
                                       GH_INTEGER, GH_LOGICAL, &
                                       GH_READ, GH_WRITE,      &
                                       GH_SCALAR,              &
                                       ANY_SPACE_9, GH_BASIS,  &
                                       CELL_COLUMN, GH_EVALUATOR
  use constants_mod,             only: r_def, i_def, l_def
  use coord_transform_mod,       only: xyz2llr
  use fs_continuity_mod,         only: W3
  use kernel_mod,                only: kernel_type

  ! Configuration modules
  use finite_element_config_mod, only: COORD_SYSTEM_XYZ
  use base_mesh_config_mod,      only: GEOMETRY_SPHERICAL

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel.
  !> Contains the metadata needed by the Psy layer
  type, public, extends(kernel_type) :: compute_geopotential_kernel_type
    private
    type(arg_type) :: meta_args(7) = (/                        &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, W3),          &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ,  ANY_SPACE_9), &

         ! Configuration
         arg_type(GH_SCALAR,  GH_REAL, GH_READ),    & ! gravity
         arg_type(GH_SCALAR,  GH_REAL, GH_READ),    & ! planet_radius
         arg_type(GH_SCALAR,  GH_INTEGER, GH_READ), & ! coord_system
         arg_type(GH_SCALAR,  GH_LOGICAL, GH_READ), & ! shallow
         arg_type(GH_SCALAR,  GH_INTEGER, GH_READ)  & ! geometry

         /)
    type(func_type) :: meta_funcs(1) = (/ &
         func_type(ANY_SPACE_9, GH_BASIS) &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: compute_geopotential_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: compute_geopotential_code

contains

!! @param[in] nlayers Number of layers
!! @param[in,out] phi Geopotential array
!! @param[in] chi_1 1st physical coordinate field
!! @param[in] chi_2 2nd physical coordinate field
!! @param[in] chi_3 3rd physical coordinate field
!! @param[in] gravity       Planet gravity
!! @param[in] planet_radius Planet radius
!! @param[in] coord_system  Enumeration for coordinate system
!! @param[in] shallow       Logical to run in shallow atmosphere mode
!! @param[in] geometry      Enumeration for the mesh geometry
!! @param[in] ndf_w3 Number of degrees of freedom per cell for w3
!! @param[in] undf_w3 Number of unique degrees of freedom for w3
!! @param[in] map_w3 Dofmap for the cell at the base of the column for w3
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi
!! @param[in] undf_chi Number of unique degrees of freedom for chi
!! @param[in] map_chi Dofmap for the cell at the base of the column for chi
!! @param[in] chi_basis Chi basis functions evaluated at w3 nodes
subroutine compute_geopotential_code(nlayers, phi,               &
                                     chi_1, chi_2, chi_3,        &
                                     gravity,                    &
                                     planet_radius,              &
                                     coord_system,               &
                                     shallow,                    &
                                     geometry,                   &
                                     ndf_w3, undf_w3, map_w3,    &
                                     ndf_chi, undf_chi, map_chi, &
                                     chi_basis)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                           :: nlayers
  integer(kind=i_def), intent(in)                           :: ndf_w3
  integer(kind=i_def), intent(in)                           :: undf_w3
  integer(kind=i_def), intent(in)                           :: ndf_chi
  integer(kind=i_def), intent(in)                           :: undf_chi
  integer(kind=i_def), dimension(ndf_w3), intent(in)        :: map_w3
  integer(kind=i_def), dimension(ndf_chi), intent(in)       :: map_chi
  real(kind=r_def), dimension(undf_w3), intent(inout)       :: phi
  real(kind=r_def), dimension(undf_chi), intent(in)         :: chi_1
  real(kind=r_def), dimension(undf_chi), intent(in)         :: chi_2
  real(kind=r_def), dimension(undf_chi), intent(in)         :: chi_3
  real(kind=r_def), dimension(1,ndf_chi,ndf_w3), intent(in) :: chi_basis

  ! Configuration variables
  real(kind=r_def),    intent(in) :: gravity
  real(kind=r_def),    intent(in) :: planet_radius
  integer(kind=i_def), intent(in) :: coord_system
  logical(kind=l_def), intent(in) :: shallow
  integer(kind=i_def), intent(in) :: geometry


  ! Internal variables
  integer(kind=i_def) :: df, dfc, k
  real(kind=r_def)    :: coord(3)
  real(kind=r_def)    :: lat, lon, radius, shallow_switch, height

  real(kind=r_def), dimension(ndf_chi) :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def) :: phi_shallow, phi_deep

  ! If geometry is spherical then calculate geopotential using radius
  if ( geometry == GEOMETRY_SPHERICAL ) then
    ! We introduce a shallow_switch, which controls whether we assume
    ! a constant geopotential with radius or whether to use inverse square law
    if ( shallow ) then
      shallow_switch = 1.0_r_def
    else
      shallow_switch = 0.0_r_def
    end if
    ! For Cartesian coordinate system, obtain radius from (X,Y,Z)
    if (coord_system == COORD_SYSTEM_XYZ) then
      do k = 0, nlayers-1
        do dfc = 1, ndf_chi
          chi_1_e(dfc) = chi_1( map_chi(dfc) + k)
          chi_2_e(dfc) = chi_2( map_chi(dfc) + k)
          chi_3_e(dfc) = chi_3( map_chi(dfc) + k)
        end do

        do df = 1, ndf_w3
          coord(:) = 0.0_r_def
          do dfc = 1, ndf_chi
            coord(1) = coord(1) + chi_1_e(dfc)*chi_basis(1,dfc,df)
            coord(2) = coord(2) + chi_2_e(dfc)*chi_basis(1,dfc,df)
            coord(3) = coord(3) + chi_3_e(dfc)*chi_basis(1,dfc,df)
          end do
          call xyz2llr(coord(1), coord(2), coord(3), lon, lat, radius)

            phi_shallow = gravity*(radius - planet_radius)
            phi_deep    = -gravity*planet_radius*(planet_radius/radius - 1.0_r_def)
            phi(map_w3(df) + k) = shallow_switch*phi_shallow &
                                + (1.0_r_def-shallow_switch)*phi_deep

        end do
      end do
    ! The spherical coordinate system already has radius as chi_3
    else
      do k = 0, nlayers-1
        do df = 1, ndf_w3
          radius = planet_radius
          do dfc = 1, ndf_chi
            radius = radius + chi_3( map_chi(dfc) + k )*chi_basis(1,dfc,df)
          end do

          phi_shallow = gravity*(radius - planet_radius)
          phi_deep    = -gravity*planet_radius*(planet_radius/radius - 1.0_r_def)
          phi(map_w3(df) + k) = shallow_switch*phi_shallow &
                              + (1.0_r_def-shallow_switch)*phi_deep
        end do
      end do
    end if

  ! Otherwise domain is planar Cartesian and chi_3 is the height
  else

    do k = 0, nlayers-1
      do dfc = 1, ndf_chi
        chi_3_e(dfc) = chi_3( map_chi(dfc) + k )
      end do

      do df = 1, ndf_w3
        height = 0.0_r_def
        do dfc = 1, ndf_chi
          height = height + chi_3_e(dfc)*chi_basis(1,dfc,df)
        end do

        phi(map_w3(df) + k) =  gravity*height

      end do
    end do
  end if

end subroutine compute_geopotential_code

end module compute_geopotential_kernel_mod
