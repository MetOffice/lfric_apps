!-------------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Transform CLASSIC 2-bin dust variables to UKCA 2-mode dust variables.

module dust_transform_kernel_mod

use argument_mod,      only: arg_type,           &
                             GH_FIELD, GH_REAL,  &
                             GH_SCALAR, GH_READ, &
                             GH_WRITE, CELL_COLUMN

use fs_continuity_mod, only: WTHETA
use kernel_mod,        only: kernel_type
use constants_mod,     only: r_def, i_def, pi

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel.
!> Contains the metadata needed by the PSy layer.

type, public, extends(kernel_type) :: dust_transform_kernel_type
  private
  type(arg_type) :: meta_args(8) = (/ &
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE,  WTHETA), & ! n_acc_ins
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE,  WTHETA), & ! acc_ins_du
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE,  WTHETA), & ! n_cor_ins
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE,  WTHETA), & ! cor_ins_du
       arg_type(GH_FIELD,  GH_REAL, GH_READ,   WTHETA), & ! dust1_mmr
       arg_type(GH_FIELD,  GH_REAL, GH_READ,   WTHETA), & ! dust2_mmr
       arg_type(GH_SCALAR, GH_REAL, GH_READ),           & ! rhop
       arg_type(GH_SCALAR, GH_REAL, GH_READ)            & ! rgas
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: dust_transform_code
end type dust_transform_kernel_type

public :: dust_transform_code

contains

  !> @details Transform the UM CLASSIC 2-bin dust fields into the LFRIC UKCA
  !!          2-mode dust variables. This transformation is a simple scheme
  !!          combining two transforms already in the UM and UKCA code,
  !!          including:
  !!             1. The UM's reconfiguration transforms between CLASSIC 6 bin
  !!                and 2 bin schemes (unified model documentation paper 020,
  !!                section 6.7),
  !!             2. The UKCA's dust emissions scheme transforms from CLASSIC 6
  !!                bins to UKCA 2 modes (unified model documentation paper 084,
  !!                section 10.2.6).
  !!          Thus the first step is to transform from CLASSIC 2 bin to 6 bin
  !!          mass concentrations (lines 172 - 176). The second step transforms
  !!          from 6 bin CLASSIC mass concentrations to 2 mode UKCA number and
  !!          mass concentrations (lines 180 - 186). Note that AOD at 550nm is
  !!          conserved in this transformation.
  !> @param[in]      nlayers             Number of layers
  !> @param[in,out]  n_acc_ins           UKCA accumulation mode number/volume
  !> @param[in,out]  acc_ins_du          UKCA accumulation mode mass
  !!                                     concentration
  !> @param[in,out]  n_cor_ins           UKCA coarse insoluble mode
  !!                                     number/volume
  !> @param[in,out]  cor_ins_du          UKCA coarse insoluble mode mass
  !!                                     concentration
  !> @param[in]      dust1_mmr           CLASSIC dust bin 1 mass concentration
  !> @param[in]      dust2_mmr           CLASSIC dust bin 2 mass concentration
  !> @param[in]      rhop                Density of a dust particle (quartz)
  !!                                     (= 2.65e+3)
  !> @param[in]      rgas                Specific gas constant (= 287.05 J/kg/K)
  !> @param[in]      ndf_wtheta          Number of degrees of freedom per cell
  !!                                     for potential temperature space
  !> @param[in]      undf_wtheta         Number of unique degrees of freedom in
  !!                                     segment for potential temperature
  !!                                     space
  !> @param[in]      map_wtheta          Dofmap for a segment of any potential
  !!                                     temperature space field

  subroutine dust_transform_code(nlayers,     &
                                 n_acc_ins,   &
                                 acc_ins_du,  &
                                 n_cor_ins,   &
                                 cor_ins_du,  &
                                 dust1_mmr,   &
                                 dust2_mmr,   &
                                 rhop,        &
                                 rgas,        &
                                 ndf_wtheta,  &
                                 undf_wtheta, &
                                 map_wtheta)

    !---------------------------------------
    ! UM modules containing global constants
    !---------------------------------------
    use chemistry_constants_mod, only: boltzmann

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers, ndf_wtheta, undf_wtheta
    integer(kind=i_def), intent(in), dimension(ndf_wtheta) :: map_wtheta
    real(kind=r_def), intent(inout), dimension(undf_wtheta) :: n_acc_ins
    real(kind=r_def), intent(inout), dimension(undf_wtheta) :: acc_ins_du
    real(kind=r_def), intent(inout), dimension(undf_wtheta) :: n_cor_ins
    real(kind=r_def), intent(inout), dimension(undf_wtheta) :: cor_ins_du
    real(kind=r_def), intent(in), dimension(undf_wtheta) :: dust1_mmr
    real(kind=r_def), intent(in), dimension(undf_wtheta) :: dust2_mmr
    real(kind=r_def), intent(in) :: rhop
    real(kind=r_def), intent(in) :: rgas

    ! Internal variables
    integer(kind=i_def) :: k, df
    real(kind=r_def) :: volconst

    ! Number of dust bins in the 6-bin and 2-bin UM schemes
    integer(kind=i_def), parameter :: n_6bins = 6_i_def
    integer(kind=i_def), parameter :: n_2bins = 2_i_def

    ! Specific extinction coefficients (m2.kg-1)
    real(kind=r_def), parameter :: kext6(n_6bins) = (/ 652.797_r_def, &
                                                       3626.7_r_def,  &
                                                       976.290_r_def, &
                                                       260.675_r_def, &
                                                       77.6338_r_def, &
                                                       23.9075_r_def /)
    real(kind=r_def), parameter :: kext2(n_2bins) = (/ 700.367_r_def, &
                                                       141.453_r_def /)

    ! Cube of the mean diameters of the 6-bin dust scheme (in m3)
    real(kind=r_def), parameter :: drep6c(n_6bins) = (/ 1.42261347e-21_r_def, &
                                                        4.49873508e-20_r_def, &
                                                        1.42261347e-18_r_def, &
                                                        4.49873508e-17_r_def, &
                                                        1.42261347e-15_r_def, &
                                                        4.49873508e-14_r_def /)

    ! Mass concentration equivalents for the 6-bin UM scheme
    real(kind=r_def) :: bin6_2, bin6_3, bin6_4, bin6_5

    ! Factors to convert mass concentrations from the 2-bin to 6-bin schemes
    real(kind=r_def) :: p1(n_6bins)
    real(kind=r_def) :: p2(n_6bins)

    ! Initialise p1 and p2
    p1 = (/ 0.0_r_def, 0.035_r_def, 0.22_r_def, 0.745_r_def, &
            0.0_r_def, 0.0_r_def /)
    p2 = (/ 0.0_r_def, 0.0_r_def, 0.0_r_def, 0.219_r_def, &
            0.781_r_def, 0.0_r_def /)

    ! Because we want to preserve AOD at 550nm when moving from 2 bins
    ! to 6 bins (rather than total mass) the proportions must be
    ! adjusted to account for extinction efficiency
    p1 = p1 * kext2(1) / sum( kext6 * p1 )
    p2 = p2 * kext2(2) / sum( kext6 * p2 )

    ! Factor to translate from mass to number concentration
    volconst = 6_i_def * boltzmann / (rhop * rgas * pi)

    ! Perform the conversions
    do k = 0, nlayers - 1
      do df = 1, ndf_wtheta

        ! Transform into equivalent mass concentrations in the 6-bin scheme
        ! following UM code: dustbin_conversion_mod.F90
        bin6_2 = dust1_mmr(map_wtheta(df) + k) * p1(2)
        bin6_3 = dust1_mmr(map_wtheta(df) + k) * p1(3)
        bin6_4 = dust1_mmr(map_wtheta(df) + k) * p1(4) &
                 + dust2_mmr(map_wtheta(df) + k) * p2(4)
        bin6_5 = dust2_mmr(map_wtheta(df) + k) * p2(5)

        ! Transform to 2-mode UKCA variables following GLOMAP/UKCA code:
        ! ukca_prim_du_mod.F90
        n_acc_ins(map_wtheta(df) + k) = volconst * ( bin6_2 / drep6c(2) &
                                        + (bin6_3 * 0.5_r_def) / drep6c(3) )
        n_cor_ins(map_wtheta(df) + k) = volconst * ( ( bin6_3 * 0.5_r_def ) / &
                                        drep6c(3) + bin6_4 / drep6c(4) &
                                        + bin6_5 / drep6c(5) )
        acc_ins_du(map_wtheta(df) + k) = bin6_2 + bin6_3 * 0.5_r_def
        cor_ins_du(map_wtheta(df) + k) = bin6_3 * 0.5_r_def + bin6_4 + bin6_5

      end do
    end do

  end subroutine dust_transform_code

end module dust_transform_kernel_mod
