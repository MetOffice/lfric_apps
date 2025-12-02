!-------------------------------------------------------------------------------
!(c) Crown copyright 2023 Met Office. All rights reserved.
!The file LICENCE, distributed with this code, contains details of the terms
!under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Compute the SKEB2 divergent increments
module skeb_dudiv_kernel_mod

  use argument_mod,      only: arg_type, GH_FIELD,          &
                               GH_REAL, GH_WRITE, GH_READ,  &
                               CELL_COLUMN, GH_INTEGER,     &
                               GH_SCALAR, STENCIL, CROSS,   &
                               ANY_DISCONTINUOUS_SPACE_1
  use fs_continuity_mod, only: W3, W2
  use constants_mod,     only: r_def, i_def
  use kernel_mod,        only: kernel_type

  implicit none

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: skeb_dudiv_kernel_type
    private
    type(arg_type) :: meta_args(6) = (/                              &
    arg_type(GH_FIELD, GH_REAL, GH_WRITE, W2),                       & ! du_div_w2
    arg_type(GH_FIELD, GH_REAL, GH_READ, W3, STENCIL(CROSS)),        & ! psif_hat
    arg_type(GH_FIELD, GH_REAL, GH_READ, W2),                        & ! dx_in_w2
    arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1), & ! latitude
    arg_type(GH_SCALAR, GH_INTEGER, GH_READ ),                       & ! skeb_level_bottom
    arg_type(GH_SCALAR, GH_INTEGER, GH_READ )                        & ! skeb_level_top
    /)
    integer :: operates_on = CELL_COLUMN

  contains
    procedure, nopass :: skeb_dudiv_code
  end type skeb_dudiv_kernel_type

  public skeb_dudiv_code
contains

  !> @brief Calculate SKEB2 divergent increments
  !> @param[in]    nlayers     The number of layers
  !> @param[out]   du_div_w2   Divergent wind increments
  !> @param[in]    psif_hat    Velocity potential
  !> @param[in]    map_w3_sten_size Size of the stencil for w3 fields
  !> @param[in]    map_w3_sten      Stencil map for w3 fields
  !> @param[in]    dx_in_w2    delta-x and w2 dofs
  !> @param[in]    latitude    latitude of w2 dofs
  !> @param[in]    skeb_level_bottom Bottom SKEB level
  !> @param[in]    skeb_level_top    Top SKEB level
  !> @param[in]    ndf_w2      Number of DOFs per cell for w2 space
  !> @param[in]    undf_w2     Number of unique DOFs  for w2 space
  !> @param[in]    map_w2      dofmap for the cell at the base of the column for w2 space
  !> @param[in]    ndf_w3      Number of DOFs per cell for density space
  !> @param[in]    undf_w3     Number of unique DOFs  for density space
  !> @param[in]    map_w3      dofmap for the cell at the base of the column for density space
  !> @param[in]    ndf_2d      Number of degrees of freedom per cell for 2D space
  !> @param[in]    undf_2d     Number of unique degrees of freedom for 2D space
  !> @param[in]    map_2d      Dofmap for the cell at the base of the column for 2D space
  subroutine skeb_dudiv_code(nlayers,           &
                             du_div_w2,         &
                             psif_hat,          &
                             map_w3_sten_size,  &
                             map_w3_sten,       &
                             dx_at_w2,          &
                             latitude,          &
                             skeb_level_bottom, &
                             skeb_level_top,    &
                             ndf_w2,            &
                             undf_w2,           &
                             map_w2,            &
                             ndf_w3,            &
                             undf_w3,           &
                             map_w3,            &
                             ndf_2d,            &
                             undf_2d,           &
                             map_2d)

    implicit none

    !Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_w3, ndf_w2, ndf_2d
    integer(kind=i_def), intent(in) :: undf_w3, undf_w2, undf_2d
    integer(kind=i_def), intent(in) :: map_w3_sten_size
    integer(kind=i_def), intent(in),  dimension(ndf_w3)   :: map_w3
    integer(kind=i_def), intent(in),  dimension(ndf_w2)   :: map_w2
    integer(kind=i_def), intent(in),  dimension(ndf_2d)   :: map_2d
    integer(kind=i_def), intent(in),  dimension(ndf_w3,map_w3_sten_size) :: map_w3_sten

    ! Fields
    real(kind=r_def),    intent(in),    dimension(undf_w3)  :: psif_hat
    real(kind=r_def),    intent(in),    dimension(undf_w2)  :: dx_at_w2
    real(kind=r_def),    intent(in),    dimension(undf_2d)  :: latitude
    real(kind=r_def),    intent(inout), dimension(undf_w2)  :: du_div_w2

    ! Scalars
    integer(kind=i_def), intent(in) :: skeb_level_bottom, skeb_level_top

    integer(kind=i_def) :: k, df, fac

!
! Here we're calculating:
!
! (du, dv) = (-dpsi/dx, -dpsi/dy)*sgn(phi)
!
! Calculated as a finite difference using the following:
! x-comp = i - i+1
! dof 3 = 1 - 4 = 1 - df+1
! dof 1 = 2 - 1 = df+1 - 1
! y-comp = j - j+1 but because beta = -y we need j+1 - j
! dof 4 = 5 - 1 = df+1 - 1
! dof 2 = 1 - 3 = 1 - df+1
!
    do df = 1,4

      if (df == 1 .or. df == 4) then
        fac = -1.0_r_def
      else
        fac = 1.0_r_def
      end if

      do k = skeb_level_bottom, skeb_level_top

        du_div_w2(map_w2(df)+k-1) = fac*sign(1.0_r_def, latitude(map_2d(df))) &
                                    * (psif_hat(map_w3_sten(1,1)+k-1) - &
                                       psif_hat(map_w3_sten(1,df+1)+k-1)) / &
                                       dx_at_w2(map_w2(df)+k-1)

      end do

    end do

  end subroutine skeb_dudiv_code

end module skeb_dudiv_kernel_mod
