!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
module atl_bl_inc_kernel_mod

  use argument_mod,      only : arg_type,              &
                                GH_FIELD, GH_OPERATOR, &
                                GH_SCALAR, GH_INTEGER, &
                                GH_READ, GH_INC,       &
                                GH_REAL, CELL_COLUMN,  &
                                ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,     only : r_def, i_def, r_um
  use fs_continuity_mod, only : W1, W2, W3
  use kernel_mod,        only : kernel_type
  use reference_element_mod,    only: N

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  type, public, extends(kernel_type) :: atl_bl_inc_kernel_type
  PRIVATE
  TYPE(arg_type) :: meta_args(7) = (/ &
  arg_type(GH_FIELD, GH_REAL, GH_INC, W2), &
  arg_type(GH_FIELD, GH_REAL, GH_READ, W2), &
  arg_type(GH_FIELD, GH_REAL, GH_READ, W2), &
  arg_type(GH_FIELD, GH_REAL, GH_READ, W2), &
  arg_type(GH_FIELD, GH_INTEGER, GH_READ, ANY_DISCONTINUOUS_SPACE_1), &
  arg_type(GH_FIELD, GH_INTEGER, GH_READ, ANY_DISCONTINUOUS_SPACE_1), &
  arg_type(GH_SCALAR, GH_INTEGER, GH_READ) &
  /)
  INTEGER :: operates_on = CELL_COLUMN
  CONTAINS
  PROCEDURE, NOPASS :: atl_bl_inc_code
  END TYPE

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: atl_bl_inc_code

  contains

  !> @brief (Adjoint of) computes boundary layer u inc
  !! @param[in] nlayers Number of layers
  !! @param[in,out] u_inc Output
  !! @param[in] ndf_w2 Number of degrees of freedom per cell for the output field
  !! @param[in] undf_w2 Unique number of degrees of freedom  for the output field
  !! @param[in] map_w2 Dofmap for the cell at the base of the column for the output field
  subroutine atl_bl_inc_code(        nlayers,              &
                                     u_inc,                &
                                     u,                    &
                                     Auv,                  &
                                     Buv_inv,              &
                                     face_selector_ew,     &
                                     face_selector_ns,     &
                                     Blevs_m,                &
                                     ndf_w2, undf_w2, map_w2,         &
                                     ndf_w3_2d, undf_w3_2d, map_w3_2d )

  implicit none

  ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_w2
    integer(kind=i_def), intent(in) :: ndf_w2
    integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
    real(kind=r_def), dimension(undf_w2), intent(inout) :: u
    real(kind=r_def), dimension(undf_w2), intent(in) :: auv
    real(kind=r_def), dimension(undf_w2), intent(in) :: buv_inv
    real(kind=r_def), dimension(undf_w2), intent(inout) :: u_inc
    integer(kind=i_def), intent(in) :: ndf_w3_2d
    integer(kind=i_def), intent(in) :: undf_w3_2d
    integer(kind=i_def), dimension(ndf_w3_2d), intent(in) :: map_w3_2d
    integer(kind=i_def), dimension(undf_w3_2d), intent(in) :: face_selector_ew
    integer(kind=i_def), dimension(undf_w3_2d), intent(in) :: face_selector_ns
    integer(kind=i_def), intent(in) :: blevs_m

 ! Internal variables
     integer(kind=i_def) :: df
    integer(kind=i_def) :: k
    integer(kind=i_def) :: j
    real(kind=r_def), dimension(blevs_m) :: a0
    real(kind=r_def), dimension(blevs_m) :: a1
    real(kind=r_def), dimension(blevs_m) :: a2
    real(kind=r_def), dimension(blevs_m) :: u_rhs
    real(kind=r_def), dimension(blevs_m) :: u_out
    real(kind=r_def), dimension(blevs_m) :: factor_u
    integer :: idx
    integer :: idx_1

    do j = face_selector_ew(map_w3_2d(1)) + face_selector_ns(map_w3_2d(1)), 1, -1
      df = j
      if (j == 3 .AND. face_selector_ns(map_w3_2d(1)) == 2 .AND. face_selector_ew(map_w3_2d(1)) == 1) then
        df = n
      end if

      u_out = 0.0_r_def
      u_rhs = 0.0_r_def
      a0(:) = 0.0_r_def
      a1(:) = 0.0_r_def
      a2(:) = 0.0_r_def
      factor_u(:) = 0.0_r_def

! ===================== set up coeffs a0,a1,a2,factor_u  =====================

        a0(1)    = 1.0_r_def+(Auv(map_w2(df)+1)+Auv(map_w2(df)+0))/Buv_inv(map_w2(df)+1)
        a1(1)    = -Auv(map_w2(df)+1)/Buv_inv(map_w2(df)+1)

         DO k=2,BLevs_m-1
        a0(k)    = 1.0_r_def+(Auv(map_w2(df)+k)+Auv(map_w2(df)+k-1))/Buv_inv(map_w2(df)+k)
        a2(k)    = -Auv(map_w2(df)+k-1)/Buv_inv(map_w2(df)+k)
        a1(k)    = -Auv(map_w2(df)+k)/Buv_inv(map_w2(df)+k)
          end do

        a0(BLevs_m)    = 1.0_r_def+Auv(map_w2(df)+BLevs_m-1)/Buv_inv(map_w2(df)+BLevs_m)
        a2(BLevs_m)    = -Auv(map_w2(df)+BLevs_m-1)/Buv_inv(map_w2(df)+BLevs_m)

        a0(1) = 1.0_r_def/a0(1)

      do k=2,BLevs_m
        factor_u(k) = a2(k)*a0(k-1)
        a0(k)       = 1.0_r_def/(a0(k)-factor_u(k)*a1(k-1))
       END DO

! ================  Adj of Solve for u_inc and transform to upper triangular form =============

      do k = 1, blevs_m - 1

        u_out(k) = u_out(k) + u_inc(map_w2(df) + k - 1)
        u_inc(map_w2(df) + k - 1) = 0.0_r_def

        u_out(k + 1) = u_out(k + 1) + (-a0(k) * a1(k) * u_out(k))
        u_rhs(k) = u_rhs(k) + a0(k) * u_out(k)
        u_out(k) = 0.0_r_def
      enddo

      u_out(blevs_m) = u_out(blevs_m) + u_inc(map_w2(df) + blevs_m - 1)
      u_inc(map_w2(df) + blevs_m - 1) = 0.0_r_def

      u_rhs(blevs_m) = u_rhs(blevs_m) + a0(blevs_m) * u_out(blevs_m)
      u_out(blevs_m) = 0.0_r_def

      do k = blevs_m, 2, -1
        u_rhs(k - 1) = u_rhs(k - 1) + (-factor_u(k) * u_rhs(k))
      enddo

      u(blevs_m + map_w2(df) - 2) = u(blevs_m + map_w2(df) - 2) + &
      auv(blevs_m + map_w2(df) - 1) * u_rhs(blevs_m) / buv_inv(blevs_m + map_w2(df))
      u(blevs_m + map_w2(df) - 1) = u(blevs_m + map_w2(df) - 1) - &
      auv(blevs_m + map_w2(df) - 1) * u_rhs(blevs_m) / buv_inv(blevs_m + map_w2(df))
      u_rhs(blevs_m) = 0.0_r_def

      do k = blevs_m - 1, 2, -1
        u(k + map_w2(df)) = u(k + map_w2(df)) + auv(k + map_w2(df)) * u_rhs(k) / buv_inv(k + map_w2(df))
        u(k + map_w2(df) - 1) = u(k + map_w2(df) - 1) - auv(k + map_w2(df)) * u_rhs(k) / buv_inv(k + map_w2(df))
        u(k + map_w2(df) - 2) = u(k + map_w2(df) - 2) + auv(k + map_w2(df) - 1) * u_rhs(k) / buv_inv(k + map_w2(df))
        u(k + map_w2(df) - 1) = u(k + map_w2(df) - 1) - auv(k + map_w2(df) - 1) * u_rhs(k) / buv_inv(k + map_w2(df))
        u_rhs(k) = 0.0_r_def
      enddo

      u(map_w2(df) + 1) = u(map_w2(df) + 1) + auv(map_w2(df) + 1) * u_rhs(1) / buv_inv(map_w2(df) + 1)
      u(map_w2(df)) = u(map_w2(df)) - auv(map_w2(df) + 1) * u_rhs(1) / buv_inv(map_w2(df) + 1)
      u(map_w2(df)) = u(map_w2(df)) - auv(map_w2(df)) * u_rhs(1) / buv_inv(map_w2(df) + 1)
      u_rhs(1) = 0.0_r_def

      do idx_1 = blevs_m, 1, -1
        u_out(idx_1) = 0.0_r_def
      enddo
      do idx = blevs_m, 1, -1
        u_rhs(idx) = 0.0_r_def
      enddo

    enddo

  end subroutine atl_bl_inc_code

end module atl_bl_inc_kernel_mod
