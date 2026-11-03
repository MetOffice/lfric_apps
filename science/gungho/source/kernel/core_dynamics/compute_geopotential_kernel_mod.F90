!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Computes the geopotential field for shallow (constant gradient with
!!        height) and deep (varying gradient with height) gravity assumptions.
!> @details The geopotential is normalised so that it is zero at the domain
!!          surface (radius = planet_radius).
!!          For shallow geometries, this gives:
!!            geopotential = g*(radius - planet_radius)
!!          For deep geometries, this gives:
!!            geopotential = g*planet_radius/radius*(radius - planet_radius)
module compute_geopotential_kernel_mod

  use argument_mod,              only: arg_type, func_type,                    &
                                       GH_FIELD, GH_SCALAR,                    &
                                       GH_REAL,GH_INTEGER, GH_LOGICAL,         &
                                       GH_READ, GH_WRITE,                      &
                                       CELL_COLUMN
  use constants_mod,             only: r_def, i_def, l_def
  use fs_continuity_mod,         only: W3
  use kernel_mod,                only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel.
  !> Contains the metadata needed by the Psy layer
  type, public, extends(kernel_type) :: compute_geopotential_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                        &
        arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W3),                         & ! geopotential
        arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),                         & ! height
        ! Configuration
        arg_type(GH_SCALAR, GH_REAL,    GH_READ),                              & ! gravity
        arg_type(GH_SCALAR, GH_REAL,    GH_READ),                              & ! planet_radius
        arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                               & ! shallow

    /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: compute_geopotential_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: compute_geopotential_code

contains

!> @param[in]     nlayers       Number of layers in mesh
!> @param[in,out] phi           Geopotential field to compute
!> @param[in]     height        Field containing height of W3 DoFs
!> @param[in]     gravity       Planet gravity at surface
!> @param[in]     planet_radius Planet radius
!> @param[in]     shallow       Logical to run in shallow gravity mode
!> @param[in]     ndf_w3        Num DoFs per cell for W3
!> @param[in]     undf_w3       Num DoFs in this partition for W3
!> @param[in]     map_w3        DoFmap for W3
subroutine compute_geopotential_code(nlayers, phi,                             &
                                     height,                                   &
                                     gravity,                                  &
                                     planet_radius,                            &
                                     shallow,                                  &
                                     ndf_w3, undf_w3, map_w3)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_w3
  integer(kind=i_def), intent(in)    :: undf_w3
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
  real(kind=r_def),    intent(inout) :: phi(undf_w3)
  real(kind=r_def),    intent(in)    :: height(undf_w3)
  ! Configuration variables
  real(kind=r_def),    intent(in)    :: gravity
  real(kind=r_def),    intent(in)    :: planet_radius
  logical(kind=l_def), intent(in)    :: shallow


  ! Internal variables
  integer(kind=i_def) :: df, w3_idx_b, w3_idx_t
  real(kind=r_def)    :: shallow_switch
  real(kind=r_def)    :: phi_shallow(nlayers)
  real(kind=r_def)    :: phi_deep(nlayers)

  ! We introduce a shallow_switch, which controls whether we assume
  ! a constant geopotential with radius or whether to use inverse square law
  if ( shallow ) then
    shallow_switch = 1.0_r_def
  else
    shallow_switch = 0.0_r_def
  end if

  do df = 1, ndf_w3
    w3_idx_b = map_w3(df)
    w3_idx_t = map_w3(df) + nlayers - 1
    phi_shallow(:) = gravity*height(w3_idx_b:w3_idx_t)
    phi_deep(:) =                                                              &
        gravity*planet_radius*height(w3_idx_b:w3_idx_t)                        &
        / (planet_radius + height(w3_idx_b:w3_idx_t))

    phi(w3_idx_b:w3_idx_t) =                                                   &
        shallow_switch*phi_shallow + (1.0_r_def-shallow_switch)*phi_deep
  end do

end subroutine compute_geopotential_code

end module compute_geopotential_kernel_mod
