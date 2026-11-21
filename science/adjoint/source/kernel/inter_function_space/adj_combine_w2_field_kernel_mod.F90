!-----------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Combine W2h and W2v fields into a single W2 field
module adj_combine_w2_field_kernel_mod

  use, intrinsic :: iso_fortran_env, only: real32, real64

  use argument_mod,                  only : arg_type, GH_INTEGER,      &
                                            GH_FIELD, GH_REAL,         &
                                            GH_READWRITE, GH_INC,      &
                                            GH_READ,                   &
                                            ANY_DISCONTINUOUS_SPACE_2, &
                                            ANY_DISCONTINUOUS_SPACE_3, &
                                            CELL_COLUMN
  use constants_mod,                 only : i_def, l_def
  use fs_continuity_mod,             only : W2, W2h, W2v
  use kernel_mod,                    only : kernel_type
  use sci_face_selector_support_mod, only : face_from_face_selector

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: adj_combine_w2_field_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                           &
         ! NB: This is to be used to write to a continuous W2 field, but using
         ! a discontinuous data pattern, so use discontinuous metadata
         arg_type(GH_FIELD, GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL,    GH_INC,       W2h),                       &
         arg_type(GH_FIELD, GH_REAL,    GH_READWRITE, W2v),                       &
         arg_type(GH_FIELD, GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_3), &
         arg_type(GH_FIELD, GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_3)  &
         /)
    integer :: operates_on = CELL_COLUMN
  end type
  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: adj_combine_w2_field_code

  ! Generic interface for real32 and real64 types
  interface adj_combine_w2_field_code
    module procedure  &
      adj_combine_w2_field_code_real32, &
      adj_combine_w2_field_code_real64
  end interface

contains

!> @brief Combine W2h and W2v fields into a single W2 field
!> @param[in]     nlayers           Number of layers in the mesh
!> @param[in,out] uvw               Full W2 field
!> @param[in]     uv                Input W2h field
!> @param[in]     w                 Input W2v field
!> @param[in]     face_selector_ew  2D field indicating which W/E faces
!!                                  to loop over in this column
!> @param[in]     face_selector_ns  2D field indicating which N/S faces
!!                                  to loop over in this column
!> @param[in]     ndf_w2            Num of DoFs for W2 per cell
!> @param[in]     undf_w2           Num of DoFs for this partition for W2
!> @param[in]     map_w2            Map for W2
!> @param[in]     ndf_w2h           Num of DoFs for W2h per cell
!> @param[in]     undf_w2h          Num of DoFs for this partition for W2h
!> @param[in]     map_w2h           Map for W2h
!> @param[in]     ndf_w2v           Num of DoFs for W2v per cell
!> @param[in]     undf_w2v          Num of DoFs for this partition for W2v
!> @param[in]     map_w2v           Map for W2v
!> @param[in]     ndf_w3_2d         Num of DoFs for 2D W3 per cell
!> @param[in]     undf_w3_2d        Num of DoFs for this partition for 2D W3
!> @param[in]     map_w3_2d         Map for 2D W3

! REAL64 PRECISION
! ==================
subroutine adj_combine_w2_field_code_real64( nlayers,                         &
                                             uvw, uv, w,                      &
                                             face_selector_ew,                &
                                             face_selector_ns,                &
                                             ndf_w2,  undf_w2,  map_w2,       &
                                             ndf_w2h, undf_w2h, map_w2h,      &
                                             ndf_w2v, undf_w2v, map_w2v,      &
                                             ndf_w3_2d, undf_w3_2d, map_w3_2d )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w2, ndf_w2h, ndf_w2v
  integer(kind=i_def), intent(in) :: undf_w2, undf_w2h, undf_w2v
  integer(kind=i_def), intent(in) :: ndf_w3_2d, undf_w3_2d
  integer(kind=i_def), dimension(ndf_w2),     intent(in)    :: map_w2
  integer(kind=i_def), dimension(ndf_w2h),    intent(in)    :: map_w2h
  integer(kind=i_def), dimension(ndf_w2v),    intent(in)    :: map_w2v
  integer(kind=i_def), dimension(ndf_w3_2d),  intent(in)    :: map_w3_2d
  real(kind=real64),   dimension(undf_w2),    intent(inout) :: uvw
  real(kind=real64),   dimension(undf_w2h),   intent(inout) :: uv
  real(kind=real64),   dimension(undf_w2v),   intent(inout) :: w
  integer(kind=i_def), dimension(undf_w3_2d), intent(in)    :: face_selector_ew
  integer(kind=i_def), dimension(undf_w3_2d), intent(in)    :: face_selector_ns

  ! Internal variables
  integer(kind=i_def) :: df, k, j
  integer(kind=i_def) :: hori_dofs_to_do
  logical(kind=l_def) :: lowest_order

  if (ndf_w2 == 6) then
    lowest_order = .true.
  else
    lowest_order = .false.
  end if

  if (lowest_order) then
    hori_dofs_to_do = ABS(face_selector_ew(map_w3_2d(1))) + ABS(face_selector_ns(map_w3_2d(1)))
  else
    hori_dofs_to_do = ndf_w2h
  end if

  ! Loop over vertical W2 DoFs
  if (ndf_w2v == 2) then
    ! Optimal looping for lowest-order space
    do k = nlayers, 0, -1
      w(map_w2v(1) + k) = w(map_w2v(1) + k) &
                        + uvw(map_w2(ndf_w2h+1) + k)
      uvw(map_w2(ndf_w2h+1) + k) = 0.0_real64
    end do

  else
    ! General looping for not lowest-order
    do df = ndf_w2v, 1, -1
      do k = nlayers-1, 0, -1
        w(map_w2v(df) + k) = w(map_w2v(df) + k) &
                          + uvw(map_w2(ndf_w2h+df) + k)
        uvw(map_w2(ndf_w2h+df) + k) = 0.0_real64
      end do
    end do
  end if

  ! Loop over horizontal W2 DoFs
  do j = hori_dofs_to_do, 1, -1
    df = face_from_face_selector(j, face_selector_ew(map_w3_2d(1)), face_selector_ns(map_w3_2d(1)))
    if (.not. lowest_order) df = j

    ! Loop over layers
    do k = nlayers-1, 0, -1
      uv(map_w2h(df) + k) = uv(map_w2h(df) + k) &
                          + uvw(map_w2(df) + k)
      uvw(map_w2(df) + k) = 0.0_real64
    end do
  end do

end subroutine adj_combine_w2_field_code_real64

! REAL32 PRECISION
! ==================
subroutine adj_combine_w2_field_code_real32( nlayers,                         &
                                             uvw, uv, w,                      &
                                             face_selector_ew,                &
                                             face_selector_ns,                &
                                             ndf_w2,  undf_w2,  map_w2,       &
                                             ndf_w2h, undf_w2h, map_w2h,      &
                                             ndf_w2v, undf_w2v, map_w2v,      &
                                             ndf_w3_2d, undf_w3_2d, map_w3_2d )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w2, ndf_w2h, ndf_w2v
  integer(kind=i_def), intent(in) :: undf_w2, undf_w2h, undf_w2v
  integer(kind=i_def), intent(in) :: ndf_w3_2d, undf_w3_2d
  integer(kind=i_def), dimension(ndf_w2),     intent(in)    :: map_w2
  integer(kind=i_def), dimension(ndf_w2h),    intent(in)    :: map_w2h
  integer(kind=i_def), dimension(ndf_w2v),    intent(in)    :: map_w2v
  integer(kind=i_def), dimension(ndf_w3_2d),  intent(in)    :: map_w3_2d
  real(kind=real32),   dimension(undf_w2),    intent(inout) :: uvw
  real(kind=real32),   dimension(undf_w2h),   intent(inout) :: uv
  real(kind=real32),   dimension(undf_w2v),   intent(inout) :: w
  integer(kind=i_def), dimension(undf_w3_2d), intent(in)    :: face_selector_ew
  integer(kind=i_def), dimension(undf_w3_2d), intent(in)    :: face_selector_ns

  ! Internal variables
  integer(kind=i_def) :: df, k, j
  integer(kind=i_def) :: hori_dofs_to_do
  logical(kind=l_def) :: lowest_order

  if (ndf_w2 == 6) then
    lowest_order = .true.
  else
    lowest_order = .false.
  end if

  if (lowest_order) then
    hori_dofs_to_do = ABS(face_selector_ew(map_w3_2d(1))) + ABS(face_selector_ns(map_w3_2d(1)))
  else
    hori_dofs_to_do = ndf_w2h
  end if

  ! Loop over vertical W2 DoFs
  if (ndf_w2v == 2) then
    ! Optimal looping for lowest-order space
    do k = nlayers, 0, -1
      w(map_w2v(1) + k) = w(map_w2v(1) + k) &
                        + uvw(map_w2(ndf_w2h+1) + k)
      uvw(map_w2(ndf_w2h+1) + k) = 0.0_real32
    end do

  else
    ! General looping for not lowest-order
    do df = ndf_w2v, 1, -1
      do k = nlayers-1, 0, -1
        w(map_w2v(df) + k) = w(map_w2v(df) + k) &
                           + uvw(map_w2(ndf_w2h+df) + k)
        uvw(map_w2(ndf_w2h+df) + k) = 0.0_real32
      end do
    end do
  end if

  ! Loop over horizontal W2 DoFs
  do j = hori_dofs_to_do, 1, -1
    df = face_from_face_selector(j, face_selector_ew(map_w3_2d(1)), face_selector_ns(map_w3_2d(1)))
    if (.not. lowest_order) df = j

    ! Loop over layers
    do k = nlayers-1, 0, -1
      uv(map_w2h(df) + k) = uv(map_w2h(df) + k) &
                          + uvw(map_w2(df) + k)
      uvw(map_w2(df) + k) = 0.0_real32
    end do
  end do

end subroutine adj_combine_w2_field_code_real32


end module adj_combine_w2_field_kernel_mod
