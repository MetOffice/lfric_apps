!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Implements Watkins algorithm to adjust vertical wind to avoid
!!        breaching Lipschitz conditions

module watkins_kernel_mod

use argument_mod,                only : arg_type, GH_SCALAR,       &
                                        GH_FIELD, GH_REAL,         &
                                        GH_WRITE, GH_READ,         &
                                        ANY_DISCONTINUOUS_SPACE_2, &
                                        GH_INTEGER, CELL_COLUMN
use fs_continuity_mod,           only : W3, W2v, W2
use constants_mod,               only : r_tran, i_def
use kernel_mod,                  only : kernel_type
use reference_element_mod,       only : E, S, N, W, B, T

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: watkins_kernel_type
  private
  type(arg_type) :: meta_args(6) = (/                                          &
       arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W2v),                         & ! first_v_wind
       arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_2),   & ! lipschitz_max_field
       arg_type(GH_FIELD,  GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2),   & ! watkins_failures
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2),                          & ! wind
       arg_type(GH_SCALAR, GH_REAL,    GH_READ),                               & ! dt
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3)                           & ! detj
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: watkins_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: watkins_code

contains

!> @brief Implements Watkins algorithm to adjust vertical wind to avoid
!!        breaching Lipschitz conditions
!> @param[in]     nlayers               The number of layers in the mesh
!> @param[in,out] first_v_wind          Wind for first vertical step
!> @param[in,out] lipschitz_max_field   Field to store maximum Lipschitz number
!> @param[in,out] watkins_failures      Integer field set to 1 if Watkins fails
!> @param[in]     wind                  3D wind
!> @param[in]     dt                    Transport time step
!> @param[in]     detj                  Det(J) at W3: the volume of cells
!> @param[in]     ndf_w2v               Num of DoFs per cell for W2V
!> @param[in]     undf_w2v              Num of W2V DoFs for this partition
!> @param[in]     map_w2v               Map of lowest-cell W2V DoFs
!> @param[in]     ndf_w3_2d             Num of DoFs per cell for W3_2d
!> @param[in]     undf_w3_2d            Num of W3_2d DoFs for this partition
!> @param[in]     map_w3_2d             Map of W3_2d DoFs
!> @param[in]     ndf_w2                Num of DoFs per cell for W2
!> @param[in]     undf_w2               Num of W2 DoFs for this partition
!> @param[in]     map_w2                Map of lowest-cell W2 DoFs
!> @param[in]     ndf_w3                Number of DoFs per cell for W3
!> @param[in]     undf_w3               Num of W3 DoFs for this partition
!> @param[in]     map_w3                Map of lowest-cell W3 DoFs
subroutine watkins_code( nlayers,             &
                         first_v_wind,        &
                         lipschitz_max_field, &
                         watkins_failures,    &
                         wind,                &
                         dt,                  &
                         detj,                &
                         ndf_w2v,             &
                         undf_w2v,            &
                         map_w2v,             &
                         ndf_w3_2d,           &
                         undf_w3_2d,          &
                         map_w3_2d,           &
                         ndf_w2,              &
                         undf_w2,             &
                         map_w2,              &
                         ndf_w3,              &
                         undf_w3,             &
                         map_w3 )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_w2v, ndf_w3, ndf_w2, ndf_w3_2d
  integer(kind=i_def), intent(in)    :: undf_w2v, undf_w3, undf_w2, undf_w3_2d
  integer(kind=i_def), intent(in)    :: map_w2(ndf_w2)
  integer(kind=i_def), intent(in)    :: map_w2v(ndf_w2v)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in)    :: map_w3_2d(ndf_w3_2d)
  real(kind=r_tran),   intent(in)    :: dt
  real(kind=r_tran),   intent(in)    :: detj(undf_w3)
  real(kind=r_tran),   intent(in)    :: wind(undf_w2)
  real(kind=r_tran),   intent(inout) :: first_v_wind(undf_w2v)
  real(kind=r_tran),   intent(inout) :: lipschitz_max_field(undf_w3_2d)
  integer(kind=i_def), intent(inout) :: watkins_failures(undf_w3_2d)

  integer(kind=i_def) :: k
  real(kind=r_tran)   :: lip_z_half, lip_hori
  real(kind=r_tran)   :: lip_hori_km1, lip_z_half_km1
  real(kind=r_tran)   :: lip_x_kp1, lip_y_kp1, lip_hori_kp1, lip_z_half_kp1
  real(kind=r_tran)   :: lip_lim_km1, lip_lim, lip_lim_kp1, slack_km1, remainder
  real(kind=r_tran)   :: max_lip_init, max_lip_new

  ! The threshold sets the maximum allowed Lipschitz number under this scheme,
  ! while the tolerance is used for checking whether the algorithm has failed,
  ! and this is slacker than machine precision to avoid machine precision errors
  ! being detected as failures.
  real(kind=r_tran), parameter :: threshold = 0.9_r_tran
  real(kind=r_tran), parameter :: tolerance = 1.0E-3_r_tran

  ! Set failures to zero
  watkins_failures(map_w3_2d(1)) = 0_i_def

  ! -------------------------------------------------------------------------- !
  ! Fill first V wind with Strang split values
  ! -------------------------------------------------------------------------- !

  ! Set the bottom value
  first_v_wind(map_w2v(1)) = 0.0_r_tran

  do k = 1, nlayers - 1
    ! Simply set vertical wind from original wind
    ! Factor of 0.5 for Strang splitting
    first_v_wind(map_w2v(1)+k) = 0.5_r_tran*wind(map_w2(B)+k)
  end do

  ! Set the top values
  first_v_wind(map_w2v(1)+nlayers) = 0.0_r_tran

  ! -------------------------------------------------------------------------- !
  ! Adjustment to avoid breaching Lipschitz conditions
  ! -------------------------------------------------------------------------- !

  slack_km1 = 0.0_r_tran

  ! Compute Lipschitz numbers for this cell
  lip_x_kp1 = dt*(wind(map_w2(E)) - wind(map_w2(W)))/detj(map_w3(1))
  lip_y_kp1 = dt*(wind(map_w2(S)) - wind(map_w2(N)))/detj(map_w3(1))
  lip_z_half_kp1 = dt*(first_v_wind(map_w2v(2)) - first_v_wind(map_w2v(1)))/detj(map_w3(1))
  lip_hori_kp1 = MAX(lip_x_kp1, lip_y_kp1, lip_x_kp1 + lip_y_kp1)
  lip_lim_kp1 = MAX(lip_z_half_kp1, lip_z_half_kp1 + lip_hori_kp1)
  max_lip_init = lip_lim_kp1
  max_lip_new = 0.0_r_tran

  do k = 0, nlayers - 1

    ! ------------------------------------------------------------------------ !
    ! Compute Lipschitz numbers
    ! ------------------------------------------------------------------------ !

    ! Recompute Lipschitz numbers for cell below -------------------------------
    if (k > 0) then

      lip_hori_km1 = lip_hori
      lip_z_half_km1 = dt*(first_v_wind(map_w2v(2)+k-1) - first_v_wind(map_w2v(1)+k-1))/detj(map_w3(1)+k-1)
      lip_lim_km1 = MAX(lip_z_half_km1, lip_z_half_km1 + lip_hori_km1)

      ! Calculate amount of slack to Lipschitz conditions for cell below
      slack_km1 = (threshold - lip_lim_km1) * detj(map_w3(1)+k-1) / dt
    end if

    ! Lipschitz numbers for this cell ------------------------------------------
    lip_hori = lip_hori_kp1
    lip_z_half = dt*(first_v_wind(map_w2v(2)+k) - first_v_wind(map_w2v(1)+k))/detj(map_w3(1)+k)
    lip_lim = MAX(lip_z_half, lip_z_half + lip_hori)

    ! Lipschitz numbers for cell above -----------------------------------------
    if (k < nlayers - 1) then
      lip_x_kp1 = dt*(wind(map_w2(E)+k+1) - wind(map_w2(W)+k+1))/detj(map_w3(1)+k+1)
      lip_y_kp1 = dt*(wind(map_w2(S)+k+1) - wind(map_w2(N)+k+1))/detj(map_w3(1)+k+1)
      lip_z_half_kp1 = dt*(first_v_wind(map_w2v(2)+k+1) - first_v_wind(map_w2v(1)+k+1))/detj(map_w3(1)+k+1)
      lip_hori_kp1 = MAX(lip_x_kp1, lip_y_kp1, lip_x_kp1 + lip_y_kp1)
      lip_lim_kp1 = MAX(lip_z_half_kp1, lip_z_half_kp1 + lip_hori_kp1)
      max_lip_init = MAX(max_lip_init, lip_lim_kp1)
    end if

    ! ------------------------------------------------------------------------ !
    ! Adjust wind if a Lipschitz condition is breached
    ! ------------------------------------------------------------------------ !

    if (lip_lim > threshold) then  ! Wind needs adjusting

      ! First try to adjust wind below this cell -------------------------------
      ! If there is enough slack from cell below, then Lipschitz conditions
      ! in this cell can be satisfied by adjusting the lower wind
      if (slack_km1 > (lip_lim - threshold) * detj(map_w3(1)+k) / dt) then

        ! Check whether cell above might need some of this slack too
        if (lip_lim_kp1 > threshold .AND. k < nlayers - 1) then
          ! TODO: don't need to use all slack here! Arbitrarily using half
          first_v_wind(map_w2v(1)+k) = first_v_wind(map_w2v(1)+k)              &
            + (lip_lim - threshold) * detj(map_w3(1)+k) / dt                   &
            + 0.5_r_tran * (slack_km1 - (lip_lim - threshold) * detj(map_w3(1)+k) / dt)

        else  ! Don't take cell above into account
          first_v_wind(map_w2v(1)+k) = first_v_wind(map_w2v(1)+k)              &
            + (lip_lim - threshold) * detj(map_w3(1)+k) / dt
        end if

      ! Need to also adjust wind above this cell -------------------------------
      else
        ! First, use up all of slack to adjust wind below
        first_v_wind(map_w2v(1)+k) = first_v_wind(map_w2v(1)+k) + slack_km1
        remainder = (lip_lim - threshold) * detj(map_w3(1)+k) / dt - slack_km1

        ! Adjust wind above (only when not top cell)
        if (k < nlayers - 1) then
          ! TODO: this is a naive method for adjusting wind above. Can we do better?
          first_v_wind(map_w2v(2)+k) = first_v_wind(map_w2v(2)+k) - remainder
        end if

      end if  ! Just adjusting wind below cell, or both winds

      ! Calculate Lipschitz number from cell below (which now won't change)
      if (k > 0) then
        lip_z_half_km1 = dt*(first_v_wind(map_w2v(2)+k-1) - first_v_wind(map_w2v(1)+k-1))/detj(map_w3(1)+k-1)
        lip_lim_km1 = MAX(lip_z_half_km1, lip_z_half_km1 + lip_hori_km1)
        max_lip_new = MAX(max_lip_new, lip_lim_km1)
      end if

    else
      ! Update running maximum Lipschitz number
      ! As nothing was changed, existing lip_lim_km1 is correct
      if (k > 0) then
        max_lip_new = MAX(max_lip_new, lip_lim_km1)
      end if

    end if  ! Lipschitz condition was breached

  end do  ! Loop through layers

  ! Update running maximum Lipschitz number with top layer
  k = nlayers - 1
  lip_z_half = dt*(first_v_wind(map_w2v(2)+k) - first_v_wind(map_w2v(1)+k))/detj(map_w3(1)+k)
  lip_lim = MAX(lip_z_half, lip_z_half + lip_hori)
  max_lip_new = MAX(max_lip_new, lip_lim)

  ! -------------------------------------------------------------------------- !
  ! Check if algorithm has failed, and if so take original wind
  ! -------------------------------------------------------------------------- !

  if (max_lip_new > max_lip_init + tolerance .AND. &
      max_lip_new > threshold + tolerance) then
    watkins_failures(map_w3_2d(1)) = 1_i_def

    ! Set the bottom value
    first_v_wind(map_w2v(1)) = 0.0_r_tran

    do k = 1, nlayers - 1
      ! Factor of 0.5 for Strang splitting
      first_v_wind(map_w2v(1)+k) = 0.5_r_tran*wind(map_w2(B)+k)
    end do

    ! Set the top values
    first_v_wind(map_w2v(1)+nlayers) = 0.0_r_tran

  else if (max_lip_new > threshold + tolerance) then
    watkins_failures(map_w3_2d(1)) = 1_i_def
  end if

  ! As we have calculated max Lipschitz numbers, store it for logging
  lipschitz_max_field(map_w3_2d(1)) = MIN(max_lip_new, max_lip_init)

end subroutine watkins_code

end module watkins_kernel_mod
