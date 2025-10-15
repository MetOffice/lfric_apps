!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes masks which indicate halo cells on each side of a partition
module halo_mask_xy_kernel_mod

  use argument_mod,          only : arg_type, func_type,         &
                                    GH_FIELD, GH_REAL, GH_WRITE, &
                                    GH_READ, GH_INTEGER,         &
                                    ANY_DISCONTINUOUS_SPACE_3,   &
                                    STENCIL, CROSS2D,            &
                                    HALO_CELL_COLUMN

  use constants_mod,         only : r_def, i_def, l_def, IMDI
  use kernel_mod,            only : kernel_type
  use reference_element_mod, only : W, S, N, E

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: halo_mask_xy_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                        &
        arg_type(GH_FIELD, GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_3),   &
        arg_type(GH_FIELD, GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_3),   &
        arg_type(GH_FIELD, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3,    &
                                                            STENCIL(CROSS2D))  &
    /)
    integer :: operates_on = HALO_CELL_COLUMN
  contains
    procedure, nopass :: halo_mask_xy_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: halo_mask_xy_code

contains

!> @brief Computes the distance of columns from the edges of cubed sphere panels
!> @param[in]     nlayers           Number of layers in mesh
!> @param[in,out] halo_mask_x       2D field to indicate whether this column is
!!                                  in the halos on the x side of the partition
!> @param[in,out] halo_mask_y       2D field to indicate whether this column is
!!                                  in the halos on the y side of the partition
!> @param[in]     halo_mask         2D field indicating whether this column is
!!                                  a halo column
!> @param[in]     stencil_sizes     Size of each branch of the cross stencil
!> @param[in]     stencil_map       Map for W3 cross stencil
!> @param[in]     max_stencil_size  Maximum size of a cross stencil branch
!> @param[in]     ndf_w3            Num DoFs per cell for W3
!> @param[in]     undf_w3           Num DoFs for this partition for W3
!> @param[in]     map_w3            DoFmap for W3
subroutine halo_mask_xy_code( nlayers,                                         &
                              halo_mask_x,                                     &
                              halo_mask_y,                                     &
                              halo_mask,                                       &
                              stencil_sizes,                                   &
                              max_stencil_size,                                &
                              stencil_map,                                     &
                              ndf_w3, undf_w3, map_w3                          &
                            )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_w3
  integer(kind=i_def), intent(in)    :: undf_w3
  integer(kind=i_def), intent(in)    :: max_stencil_size
  integer(kind=i_def), intent(in)    :: stencil_sizes(4)

  integer(kind=i_def), intent(inout) :: halo_mask_x(undf_w3)
  integer(kind=i_def), intent(inout) :: halo_mask_y(undf_w3)
  integer(kind=i_def), intent(in)    :: halo_mask(undf_w3)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in)    :: stencil_map(ndf_w3,max_stencil_size,4)

  ! Local variables
  integer(kind=i_def) :: i, sum
  logical(kind=l_def) :: halo_x, halo_y

  ! Not in the halos (kernel shouldn't be called here anyway)
  if (halo_mask(map_w3(1)) == 0) then
    halo_mask_x(map_w3(1)) = 0
    halo_mask_y(map_w3(1)) = 0

  else
    ! Check X side of halos ----------------------------------------------------
    halo_x = .false.

    ! Check W side of halos
    ! We can't just check the end W cell, since we could be skipping across
    ! the owned cells in the centre of the partition. To determine this,
    ! sum value of all cells in the W direction
    sum = 0
    do i = 1, stencil_sizes(W)
      sum = sum + halo_mask(stencil_map(1, i, W))
    end do
    if (sum /= stencil_sizes(W)) then
      halo_x = .true.
    end if

    ! Check E side of halos
    ! We can't just check the end E cell, since we could be skipping across
    ! the owned cells in the centre of the partition. To determine this,
    ! sum value of all cells in the E direction
    sum = 0
    do i = 1, stencil_sizes(E)
      sum = sum + halo_mask(stencil_map(1, i, E))
    end do
    if (sum /= stencil_sizes(E)) then
      halo_x = .true.
    end if

    ! Check Y side of halos ----------------------------------------------------
    halo_y = .false.

    ! Check S side of halos
    ! We can't just check the end S cell, since we could be skipping across
    ! the owned cells in the centre of the partition. To determine this,
    ! sum value of all cells in the S direction
    sum = 0
    do i = 1, stencil_sizes(S)
      sum = sum + halo_mask(stencil_map(1, i, S))
    end do
    if (sum /= stencil_sizes(S)) then
      halo_y = .true.
    end if

    ! Check N side of halos
    ! We can't just check the end N cell, since we could be skipping across
    ! the owned cells in the centre of the partition. To determine this,
    ! sum value of all cells in the N direction
    sum = 0
    do i = 1, stencil_sizes(N)
      sum = sum + halo_mask(stencil_map(1, i, N))
    end do
    if (sum /= stencil_sizes(N)) then
      halo_y = .true.
    end if

    ! Set halo masks ----------------------------------------------------------
    if (.not. halo_x .and. .not. halo_y) then
      ! Corners of domain
      halo_mask_x(map_w3(1)) = 1
      halo_mask_y(map_w3(1)) = 1

    else if (halo_x .and. .not. halo_y) then
      ! X side of domain
      halo_mask_x(map_w3(1)) = 1
      halo_mask_y(map_w3(1)) = 0

    else if (.not. halo_x .and. halo_y) then
      ! Y side of domain
      halo_mask_x(map_w3(1)) = 0
      halo_mask_y(map_w3(1)) = 1

    else if (halo_x .and. halo_y) then
      ! Something has gone wrong, set value to be missing data
      halo_mask_x(map_w3(1)) = IMDI
      halo_mask_y(map_w3(1)) = IMDI
    end if
  end if

end subroutine halo_mask_xy_code

end module halo_mask_xy_kernel_mod
