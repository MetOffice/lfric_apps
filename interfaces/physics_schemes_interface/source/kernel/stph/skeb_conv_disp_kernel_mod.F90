    !-------------------------------------------------------------------------------
    !(c) Crown copyright 2022 Met Office. All rights reserved.
    !The file LICENCE, distributed with this code, contains details of the terms
    !under which the code may be used.
    !-------------------------------------------------------------------------------
    !> @brief Compute the SKEB2 CAPE CONVECTION MASK
    module skeb_conv_disp_kernel_mod

    use argument_mod,      only: arg_type, GH_FIELD,          &
                                 GH_REAL, GH_WRITE, GH_READ,  &
                                 ANY_DISCONTINUOUS_SPACE_1,   &
                                 CELL_COLUMN, GH_INTEGER,     &
                                 GH_SCALAR
    use fs_continuity_mod, only: W3, Wtheta
    use constants_mod,     only: r_def, i_def
    use kernel_mod,        only: kernel_type

    implicit none

    !> Kernel metadata for Psyclone
    type, public, extends(kernel_type) :: skeb_conv_disp_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                                 &
    arg_type(GH_FIELD, GH_REAL, GH_WRITE, W3),                          & ! cdisp
    arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                          & ! rho
    arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                      & ! diff_flux
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1),   & ! cape
    arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                      & ! height
    arg_type(GH_SCALAR, GH_INTEGER, GH_READ ),                          & ! skeb_level_bottom
    arg_type(GH_SCALAR, GH_INTEGER, GH_READ ),                          & ! skeb_level_top
    arg_type(GH_SCALAR, GH_REAL, GH_READ )                              & ! gravity
    /)
    integer :: operates_on = CELL_COLUMN

    contains
    procedure, nopass :: skeb_conv_disp_code
    end type skeb_conv_disp_kernel_type

    public skeb_conv_disp_code
    contains

    !> @brief Buid the convective CAPE mask
    !> @param[in]    nlayers     The number of layers
    !> @param[inout] cdisp       Convective dissipation
    !> @param[in]    diff_flux   Massflux difference (up-down)
    !> @param[in]    rho         Dry density in its native (w3) space
    !> @param[in]    cape        CAPE  (2D field)
    !> @param[in]    skeb_level_bottom Bottom SKEB level
    !> @param[in]    skeb_level_top    Top SKEB level
    !> @param[in]    gravity     Gravity
    !> @param[in]    rdz_w3      SpacingHeight of potential temperature space levels above surface
    !> @param[in]    ndf_w3      Number of DOFs per cell for density space
    !> @param[in]    undf_w3     Number of unique DOFs  for density space
    !> @param[in]    map_w3      dofmap for the cell at the base of the column for density space
    !> @param[in]    ndf_wth     Number of DOFs per cell for potential temperature space
    !> @param[in]    undf_wth    Number of unique DOFs for potential temperature space
    !> @param[in]    map_wth     dofmap for the cell at the base of the column for potential temperature space
    !> @param[in]    ndf_2d      Number of DOFs per cell for 2D fields
    !> @param[in]    undf_2d     Number of unique DOFs  for 2D fields
    !> @param[in]    map_2d      dofmap for the cell at the base of the column for 2D fields

    subroutine skeb_conv_disp_code(nlayers,    &
                                   cdisp,      &
                                   rho,        &
                                   diff_flux,  &
                                   cape,       &
                                   height_wth, &
                                   skeb_level_bottom, &
                                   skeb_level_top,    &
                                   gravity,    &
                                   ndf_w3,     &
                                   undf_w3,    &
                                   map_w3,     &
                                   ndf_wth,    &
                                   undf_wth,   &
                                   map_wth,    &
                                   ndf_2d,     &
                                   undf_2d,    &
                                   map_2d)

    implicit none

    !Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth, ndf_w3, ndf_2d
    integer(kind=i_def), intent(in) :: undf_wth, undf_w3, undf_2d
    integer(kind=i_def), intent(in),  dimension(ndf_wth)  :: map_wth
    integer(kind=i_def), intent(in),  dimension(ndf_w3)   :: map_w3
    integer(kind=i_def), intent(in),  dimension(ndf_2d)   :: map_2d

    ! Fields
    real(kind=r_def),    intent(inout), dimension(undf_w3)  :: cdisp
    real(kind=r_def),    intent(in),    dimension(undf_w3)  :: rho
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: diff_flux
    real(kind=r_def),    intent(in),    dimension(undf_2d)  :: cape
    real(kind=r_def),    intent(in),    dimension(undf_wth)  :: height_wth

    ! Scalars
    integer(kind=i_def), intent(in) :: skeb_level_bottom, skeb_level_top
    real(kind=r_def), intent(in) :: gravity

    ! CAPE out with values below 0.01 removed
    real(kind=r_def) :: cape_kk

    integer(kind=i_def) :: k
    ! ---------------------------------------------------------------

    ! Remove negative/noise values in CAPE array
    cape_kk = cape(map_2d(1))
    if (cape_kk < 0.01_r_def) cape_kk=0.0_r_def

    ! ---------------------------------------------------------------
    ! Calculate cdisp =  vertical profile from centred derivative
    !                     of dIF f_flux
    !                     Horizontal (xy) field from cape
    ! --------------------------------------------------------------------
    !   Need to divide by (g*rho*dz) to convert Pa.s^-1 to s^-1
    !   Then: multiply by CAPE (J.kg^-1) => J.kg^-1.s^-1  or  m^2.s^-3

    do k=skeb_level_bottom, skeb_level_top

      cdisp(map_w3(1) + k-1)=ABS(diff_flux(map_wth(1) + k+1) - diff_flux(map_wth(1) + k-1) ) *cape_kk &
                          /( gravity * rho(map_w3(1) + k-1) * &
                        (height_wth(map_wth(1)+k+1)-height_wth(map_wth(1)+k-1)))
    end do

    end subroutine skeb_conv_disp_code

    end module skeb_conv_disp_kernel_mod
