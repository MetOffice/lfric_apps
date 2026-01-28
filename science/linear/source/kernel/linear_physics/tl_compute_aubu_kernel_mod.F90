!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module tl_compute_aubu_kernel_mod

  use argument_mod,      only : arg_type, func_type,   &
                                GH_FIELD, GH_OPERATOR, &
                                GH_SCALAR, GH_INTEGER, &
                                GH_READ, GH_INC, GH_READWRITE,      &
                                GH_REAL, CELL_COLUMN, GH_BASIS, GH_EVALUATOR, &
                                ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,     only : r_def, i_def, r_um
  use fs_continuity_mod, only : W1, W2, W3
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  type, public, extends(kernel_type) :: tl_compute_aubu_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                   &
         arg_type(GH_FIELD, GH_REAL, GH_INC, W2),      &  ! Auv
         arg_type(GH_FIELD, GH_REAL, GH_INC, W2),      &  ! Buv_inv
         arg_type(GH_FIELD, GH_REAL, GH_READ, W3),     &  ! Q
         arg_type(GH_FIELD, GH_REAL, GH_READ, W3),     &  ! E
         arg_type(GH_FIELD, GH_REAL, GH_READ, W2),     &  ! height_w2
         arg_type(GH_FIELD, GH_REAL, GH_READ, W2),     &  ! w2_rmultiplicity
         arg_type(GH_SCALAR, GH_REAL,    GH_READ ),  & ! dt
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ )   & ! Blevs_m
     /)

   integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: tl_compute_aubu_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: tl_compute_aubu_code

contains

!> @brief Computes coefficients Auv and Buv_inv
!! @param[in] nlayers Number of layers
!! @param[in] ndf_w2 Number of degrees of freedom per cell for the output field
!! @param[in] undf_w2 Unique number of degrees of freedom  for the output field
!! @param[in] map_w2 Dofmap for the cell at the base of the column for the output field
subroutine       tl_compute_aubu_code(nlayers,             &
                                     Auv,                  &
                                     Buv_inv,              &
                                     Q,                    &
                                     E,                    &
                                     height_w2,            &
                                     w2_rmultiplicity,     &
                                     dt,                   &
                                     Blevs_m,                 &
                                     ndf_w2, undf_w2, map_w2, &
                                     ndf_w3, undf_w3, map_w3)


  implicit none

  ! Arguments
  integer(kind=i_def),                     intent(in) ::  nlayers

  integer(kind=i_def),                     intent(in) :: undf_w3, ndf_w3
  integer(kind=i_def), dimension(ndf_w3),  intent(in) :: map_w3

  integer(kind=i_def),                     intent(in) :: undf_w2, ndf_w2
  integer(kind=i_def), dimension(ndf_w2),  intent(in) :: map_w2

    real(kind=r_def), dimension(undf_w2),  intent(in) :: w2_rmultiplicity

   REAL(kind=r_def), dimension(undf_w2), intent(inout) :: Auv  !! (0:BLevs_m)
   REAL(kind=r_def), dimension(undf_w2), intent(inout) :: Buv_inv  !!  Use inverse of Buv as this is what is averaged

   REAL(kind=r_def), dimension(undf_w3), intent(in) :: Q  !! (0:BLevs_m)
   REAL(kind=r_def), dimension(undf_w3), intent(in) :: E  !! (BLevs_m)

  real(kind=r_def), dimension(undf_w2),   intent(in) :: height_w2

  real(kind=r_def),                      intent(in)  :: dt
  integer(kind=i_def),                   intent(in)  :: Blevs_m

  ! Internal variables
   integer(kind=i_def) :: df, df3, k

df3=1

    do df = 1,4
      do k = 0, BLevs_m

      IF (k == 0) THEN
        Auv(map_w2(df)+k) = Auv(map_w2(df)+k)+         &
                  w2_rmultiplicity(map_w2(df) ) *      &
                  Q(map_w3(df3) )

      ELSE ! 1 <= k <= BLevs_m

        Auv(map_w2(df)+k) = Auv(map_w2(df)+k)+         &
                  w2_rmultiplicity(map_w2(df)+k ) *      &
                  Q(map_w3(df3) + k ) / ( height_w2(map_w2(df)+k) -  height_w2(map_w2(df)+k-1) )

        Buv_inv(map_w2(df)+k) = Buv_inv(map_w2(df)+k)+       &
                  ( w2_rmultiplicity(map_w2(df)+k ) * E( map_w3(df3) + k) ) / dt

      END IF

    end do ! k = 0, BLevs_m
  end do ! df = 1,4

end subroutine tl_compute_aubu_code

end module  tl_compute_aubu_kernel_mod
