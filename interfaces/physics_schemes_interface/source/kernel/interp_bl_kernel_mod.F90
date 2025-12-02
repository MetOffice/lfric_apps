!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Interpolates BL momentum variables from W3 to W2 dofs
!> @details Takes all the variables required for BL momentum mixing and
!>          interpolates them from their lowest order W3 dof to W2 dofs
!>          so that wind increments can be calculated in their native space

module interp_bl_kernel_mod

  use kernel_mod,               only: kernel_type
  use argument_mod,             only: arg_type,                 &
                                      GH_FIELD, GH_REAL,        &
                                      GH_INC, GH_READ,          &
                                      ANY_SPACE_1, CELL_COLUMN, &
                                      ANY_DISCONTINUOUS_SPACE_1

  use constants_mod,            only: r_def, i_def
  use fs_continuity_mod,        only: W2, W3, WTHETA
  use kernel_mod,               only: kernel_type

  use jules_control_init_mod,   only: n_surf_interp

  implicit none

  private

  !----------------------------------------------------------------------------
  ! Public types
  !----------------------------------------------------------------------------
  !> Kernel metadata type.
  type, public, extends(kernel_type) :: interp_bl_kernel_type
    private
    type(arg_type) :: meta_args(9) = (/                                  &
         arg_type(GH_FIELD, GH_REAL, GH_READ, WTHETA),                    &! rhokm_bl
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1), &! surf_interp
         arg_type(GH_FIELD, GH_REAL, GH_READ, WTHETA),                    &! ngstress_bl
         arg_type(GH_FIELD, GH_REAL, GH_READ, W3),                        &! wetrho_in_w3
         arg_type(GH_FIELD, GH_REAL, GH_INC,  W2),                        &! rhokm_w2
         arg_type(GH_FIELD, GH_REAL, GH_INC,  ANY_SPACE_1),               &! surf_interp_w2
         arg_type(GH_FIELD, GH_REAL, GH_INC,  W2),                        &! ngstress_w2
         arg_type(GH_FIELD, GH_REAL, GH_INC,  W2),                        &! wetrho_in_w2
         arg_type(GH_FIELD, GH_REAL, GH_READ, W2)                         &! w2_rmultiplicity
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: interp_bl_code
  end type interp_bl_kernel_type

  !----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !----------------------------------------------------------------------------
  public :: interp_bl_code

contains

  !> @brief Subroutine to do the re-mapping
  !> @param[in]     nlayers       Number of layers
  !> @param[in]     rhokm_bl      Momentum eddy diffusivity on BL levels
  !> @param[in]     surf_interp   Surface variables for interpolation
  !> @param[in]     ngstress_bl   Non-gradient stress function on BL levels
  !> @param[in]     wetrho_in_w3  Wet density in W3 space
  !> @param[in,out] rhokm_w2      Momentum eddy diffusivity mapped to cell faces
  !> @param[in,out] surf_interp_w2 Surface variables on cell faces
  !> @param[in,out] ngstress_w2   NG stress function mapped to cell faces
  !> @param[in,out] wetrho_in_w2  Wet density in W2 space
  !> @param[in]     w2_rmultiplicity Reciprocal of multiplicity for w2
  !> @param[in]     ndf_wth       Number of DOFs per cell for potential temperature space
  !> @param[in]     undf_wth      Number of unique DOFs for potential temperature space
  !> @param[in]     map_wth       dofmap for the cell at the base of the column for potential temperature space
  !> @param[in]     ndf_surf      Number of DOFs per cell for surface variables
  !> @param[in]     undf_surf     Number of unique DOFs for surface variables
  !> @param[in]     map_surf      dofmap for the cell at the base of the column for surface variables
  !> @param[in]     ndf_w3        Number of DOFs per cell for density space
  !> @param[in]     undf_w3       Number of unique DOFs for density space
  !> @param[in]     map_w3        dofmap for the cell at the base of the column for density space
  !> @param[in]     ndf_w2        Number of DOFs per cell for W2 space
  !> @param[in]     undf_w2       Number of unique DOFs for W2 space
  !> @param[in]     map_w2        dofmap for the cell at the base of the column for W2 space
  !> @param[in]     ndf_w2_2d     Number of DOFs per cell for W2 surface space
  !> @param[in]     undf_w2_2d    Number of unique DOFs for W2 surface space
  !> @param[in]     map_w2_2d     dofmap for the cell at the base of the column for W2 surface space
  subroutine interp_bl_code(nlayers,          &
                            rhokm_bl,         &
                            surf_interp,      &
                            ngstress_bl,      &
                            wetrho_in_w3,     &
                            rhokm_w2,         &
                            surf_interp_w2,   &
                            ngstress_w2,      &
                            wetrho_in_w2,     &
                            w2_rmultiplicity, &
                            ndf_wth,          &
                            undf_wth,         &
                            map_wth,          &
                            ndf_surf,         &
                            undf_surf,        &
                            map_surf,         &
                            ndf_w3,           &
                            undf_w3,          &
                            map_w3,           &
                            ndf_w2,           &
                            undf_w2,          &
                            map_w2,           &
                            ndf_w2_2d,        &
                            undf_w2_2d,       &
                            map_w2_2d)

    !---------------------------------------
    ! UM modules containing switches or global constants
    !---------------------------------------
    use nlsizes_namelist_mod, only: bl_levels

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers

    integer(kind=i_def), intent(in) :: ndf_wth, ndf_w3, ndf_w2, ndf_w2_2d
    integer(kind=i_def), intent(in) :: undf_wth, undf_w3, undf_w2, undf_w2_2d
    integer(kind=i_def), intent(in) :: map_wth(ndf_wth)
    integer(kind=i_def), intent(in) :: map_w3(ndf_w3)
    integer(kind=i_def), intent(in) :: map_w2(ndf_w2)
    integer(kind=i_def), intent(in) :: map_w2_2d(ndf_w2_2d)

    integer(kind=i_def), intent(in) :: ndf_surf, undf_surf
    integer(kind=i_def), intent(in) :: map_surf(ndf_surf)

    real(kind=r_def), dimension(undf_wth), intent(in) :: rhokm_bl,           &
                                                         ngstress_bl

    real(kind=r_def), dimension(undf_w3),  intent(in) :: wetrho_in_w3

    real(kind=r_def), dimension(undf_surf), intent(in)  :: surf_interp

    real(kind=r_def), dimension(undf_w2),  intent(in) :: w2_rmultiplicity

    real(kind=r_def), dimension(undf_w2), intent(inout) :: rhokm_w2,           &
                                                           ngstress_w2,        &
                                                           wetrho_in_w2

    real(kind=r_def), dimension(undf_w2_2d), intent(inout) :: surf_interp_w2

    ! Internal variables
    integer :: k, df

    ! Map the BL variables from w3 to w2 space and
    ! wtheta to a vertically shifted w2 space

    do df = 1,4
      do k = 0, n_surf_interp-1
        ! Single level variables which need mapping
        surf_interp_w2(map_w2_2d(df)+k) = surf_interp_w2(map_w2_2d(df)+k) +    &
                                          w2_rmultiplicity(map_w2(df)) *       &
                                          surf_interp(map_surf(1)+k)
      end do

      do k = 0, bl_levels-1
        ! Strictly speaking this is not a w2 field as it is the normal to the
        ! cell face at the top/bottom of the cell
        rhokm_w2(map_w2(df) + k) = rhokm_w2(map_w2(df) + k) +                  &
                                   w2_rmultiplicity(map_w2(df) + k) *          &
                                   rhokm_bl(map_wth(1) + k)
        wetrho_in_w2(map_w2(df) + k) = wetrho_in_w2(map_w2(df) + k) +          &
                                       w2_rmultiplicity(map_w2(df) + k) *      &
                                       wetrho_in_w3(map_w3(1) + k)
      end do

      do k = 1, bl_levels-1
        ! Strictly speaking this is not a w2 field as it is the normal to the
        ! cell face at the top/bottom of the cell
        ngstress_w2(map_w2(df) + k) = ngstress_w2(map_w2(df) + k) +            &
                                      w2_rmultiplicity(map_w2(df) + k) *       &
                                      ngstress_bl(map_wth(1) + k)
      end do
    end do

  end subroutine interp_bl_code

end module interp_bl_kernel_mod
