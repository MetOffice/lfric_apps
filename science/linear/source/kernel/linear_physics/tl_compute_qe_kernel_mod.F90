!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes coefficients Q and E used in the TLM boundary layer scheme for momentum.
module tl_compute_qe_kernel_mod

  use argument_mod,      only : arg_type, func_type,                          &
                                GH_FIELD, GH_OPERATOR,                        &
                                GH_SCALAR, GH_INTEGER,                        &
                                GH_READ, GH_INC, GH_READWRITE,                &
                                GH_REAL, CELL_COLUMN, GH_BASIS, GH_EVALUATOR, &
                                ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,     only : r_def, i_def, r_um
  use fs_continuity_mod, only : W1, W2, W3, Wtheta
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  type, public, extends(kernel_type) :: tl_compute_qe_kernel_type
    private
    type(arg_type) :: meta_args(14) = (/                                   &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, W3),                    & ! Q
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, W3),                    & ! E
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                        & ! ls_rho
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                        & ! height_w3
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    & ! height_wth
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1 ), & ! ls_land_fraction
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ ),                        & ! log_layer
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ ),                        & ! Blevs_m
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ ),                        & ! e_folding_levs_m
         arg_type(GH_SCALAR, GH_REAL,    GH_READ ),                        & ! u_land_m
         arg_type(GH_SCALAR, GH_REAL,    GH_READ ),                        & ! u_sea_m
         arg_type(GH_SCALAR, GH_REAL,    GH_READ ),                        & ! z_land_m
         arg_type(GH_SCALAR, GH_REAL,    GH_READ ),                        & ! z_sea_m
         arg_type(GH_SCALAR, GH_REAL,    GH_READ )                         & ! L_0_m
     /)
   integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: tl_compute_qe_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: tl_compute_qe_code

contains

!> @brief Computes coefficients Q and E used in the TLM boundary layer scheme for momentum.
!> @details The TLM BL scheme is described in Var Scientific Documentation Paper 55,
!>          https://wwwspice/~frva/VAR/view/var-2025.03.0/doc/VSDP55_1A.pdf
!! @param[in]     nlayers          Number of layers
!! @param[in,out] Q                Coefficient for TLM boundary layer
!! @param[in,out] E                Coefficient for TLM boundary layer
!! @param[in]     ls_rho           Linearisation state for density
!! @param[in]     height_w3        Height of w3 space levels above the surface
!! @param[in]     height_wth       Height of wth space levels above surface
!! @param[in]     ls_land_fraction Land area fraction
!! @param[in]     log_layer        Number of levels in log layer
!! @param[in]     Blevs_m          Number of levels in momentum boundary layer
!! @param[in]     e_folding_levs_m e-folding level in exchange coefficient for momemtum
!! @param[in]     z_land_m         Effective roughness length for momemtum over land
!! @param[in]     z_sea_m          Effective roughness length for momemtum over sea
!! @param[in]     u_land_m         Friction velocity for momemtum over land
!! @param[in]     u_sea_m          Friction velocity for momemtum over sea
!! @param[in]     L_0_m            Mixing length for momemtum
!! @param[in]     ndf_w3           Number of degrees of freedom per cell for w3 space
!! @param[in]     undf_w3          Number of unique degrees of freedom for w3 space
!! @param[in]     map_w3           Dofmap for the cell at column base for w3
!! @param[in]     ndf_wtheta       Number of degrees of freedom per cell for wtheta space
!! @param[in]     undf_wtheta      Number of unique degrees of freedom for wtheta space
!! @param[in]     map_wtheta       Dofmap for the cell at the base of the column for wtheta
!! @param[in]     ndf_2d           Number of degrees of freedom per cell for 2D field
!! @param[in]     undf_2d          Number of unique degrees of freedom for 2D field
!! @param[in]     map_2d           Dofmap for 2D field
subroutine tl_compute_qe_code( nlayers,                             &
                               Q,                                   &
                               E,                                   &
                               ls_rho,                              &
                               height_w3,                           &
                               height_wth,                          &
                               ls_land_fraction,                    &
                               log_layer,                           &
                               Blevs_m,                             &
                               e_folding_levs_m,                    &
                               u_land_m,                            &
                               u_sea_m,                             &
                               z_land_m,                            &
                               z_sea_m,                             &
                               L_0_m,                               &
                               ndf_w3, undf_w3, map_w3,             &
                               ndf_wtheta, undf_wtheta, map_wtheta, &
                               ndf_2d, undf_2d, map_2d )

  implicit none

  ! Arguments
  integer(kind=i_def),                         intent(in)    :: nlayers
  integer(kind=i_def),                         intent(in)    :: undf_w3, ndf_w3
  integer(kind=i_def), dimension(ndf_w3),      intent(in)    :: map_w3
  integer(kind=i_def),                         intent(in)    :: undf_wtheta, ndf_wtheta
  integer(kind=i_def), dimension(ndf_wtheta),  intent(in)    :: map_wtheta
  integer(kind=i_def),                         intent(in)    :: undf_2d, ndf_2d
  integer(kind=i_def), dimension(ndf_2d),      intent(in)    :: map_2d
  real(kind=r_def),    dimension(undf_w3),     intent(inout) :: Q ! (0:BLevs_m)
  real(kind=r_def),    dimension(undf_w3),     intent(inout) :: E ! (BLevs_m)
  real(kind=r_def),    dimension(undf_w3),     intent(in)    :: ls_rho
  real(kind=r_def),    dimension(undf_w3),     intent(in)    :: height_w3
  real(kind=r_def),    dimension(undf_wtheta), intent(in)    :: height_wth
  real(kind=r_def),    dimension(undf_2d),     intent(in)    :: ls_land_fraction
  integer(kind=i_def),                         intent(in)    :: log_layer,Blevs_m,e_folding_levs_m
  real(kind=r_def),                            intent(in)    :: u_land_m, u_sea_m, z_land_m, z_sea_m, L_0_m

  ! Internal variables
  integer(kind=i_def)         :: df, k
  real(kind=r_def)            :: roughness_length_m
  real(kind=r_def), parameter :: Von_Karman = 0.4_r_def
  real(kind=r_def)            :: L_diff_m(1:BLevs_m)
  real(kind=r_def)            :: u1
  real(kind=r_def)            :: num
  real(kind=r_def)            :: denom

  df = 1 ! map_w3 only has range 1:1, map_wtheta has range 1:2 with 1 being lower face of cell
  roughness_length_m = z_land_m * ls_land_fraction(map_2d(df)) + z_sea_m * (1.0_r_def - ls_land_fraction(map_2d(df)))

  ! Define the mixing length, L_diff
  ! L_diff is on centre of horizontal faces - for convenience indexed by centre of cell above.
  ! Vertical numbering (k index) matches that in old Var schemePF_bdy_lyr.f90 (see link above),
  ! in which centre of lowest cell is rho level 1 (allowing for the fact that the PF model is for New Dynamics,
  ! which has an extra rho level which means that W3 k=0 in LFRic is equivalent to rho levels k=1 in VAR).
  ! L_diff_m(k) defined for 1 <= k <= BLevs_m.
  do k = 1, BLevs_m
    if (k <= Log_layer) then
      ! num and denom to take log of
      num = height_w3(map_w3(df) + k) - height_wth(map_wtheta(df)) + roughness_length_m
      denom = height_w3(map_w3(df) + k - 1) - height_wth(map_wtheta(df)) + roughness_length_m
      ! num and denom of final expression
      denom = log(num / denom) + Von_Karman * (height_w3(map_w3(df) + k) - height_w3(map_w3(df) + k - 1)) / L_0_m
      num = (height_w3(map_w3(df) + k) - height_w3(map_w3(df) + k - 1))
    else
      num = height_wth(map_wtheta(df) + k) - height_wth(map_wtheta(df)) + roughness_length_m
      denom = 1.0_r_def + (height_wth(map_wtheta(df) + k) - height_wth(map_wtheta(df)) + roughness_length_m) / L_0_m
    end if
    L_diff_m(k) = Von_Karman * num / denom
  end do

  ! Define Q and E
  ! E is on cell centres
  ! Q is on centre of horizontal faces - for convenience indexed by centre of cell above
  ! Q(k) defined for 0 <= k <= BLevs_m
  ! E(k) defined for 1 <= k <= BLevs_m
  u1 = u_land_m * ls_land_fraction(map_2d(df)) + u_sea_m * (1.0_r_def - ls_land_fraction(map_2d(df)))
  do k = 0, BLevs_m
    if (k == 0) then
      ! num and denom to take log of
      num = height_w3(map_w3(df) + 0) - height_wth(map_wtheta(df) + 0) + roughness_length_m
      denom = roughness_length_m
      ! num and denom of final expression
      denom = log(num / denom)
      num = Von_Karman * u1
      Q(map_w3(df) + k) = num / denom
    else ! i.e., k >= 1
      num = height_wth(map_wtheta(df)) - height_wth(map_wtheta(df) + k)
      denom = height_wth(map_wtheta(df) + e_folding_levs_m) - height_wth(map_wtheta(df))
      Q(map_w3(df) + k) = L_diff_m(k) * u1 * exp(num / denom)
      E(map_w3(df) + k) = ls_rho(map_w3(df) + k - 1) * (height_wth(map_wtheta(df) + k) - height_wth(map_wtheta(df) + k - 1))
    end if
  end do

end subroutine tl_compute_qe_code

end module tl_compute_qe_kernel_mod
