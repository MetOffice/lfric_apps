! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module defining and managing structures to store EasyAerosol
! distributions.
!
! Method:
!
!  Provide structure and memory allocation/deallocation routines.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: radiation_control
!
! Code description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.
!
! Contained subroutines in this module:
!   allocate_easyaerosol_rad           (Public)
!   allocate_easyaerosol_cdnc          (Public)
!   deallocate_easyaerosol_rad         (Public)
!   deallocate_easyaerosol_cdnc        (Public)
!
! --------------------------------------------------------------------------
module def_easyaerosol

use um_types, only: real_umphys

implicit none

!
! Optical properties: 3D distributions with also a dependence on
! spectral waveband.
!
type :: t_easyaerosol_rad
  integer :: dim1
  integer :: dim2
  integer :: dim3
  integer :: dim4
  real(kind=real_umphys), allocatable :: extinction(:,:,:,:)
  real(kind=real_umphys), allocatable :: absorption(:,:,:,:)
  real(kind=real_umphys), allocatable :: asymmetry(:,:,:,:)
end type t_easyaerosol_rad

!
! Cloud droplet number concentrations: 3D distributions
!
type :: t_easyaerosol_cdnc
  integer :: dim1
  integer :: dim2
  integer :: dim3
  real(kind=real_umphys), allocatable :: cdnc(:,:,:)
end type t_easyaerosol_cdnc

integer :: i, j, k, l

character(len=*), parameter, private :: ModuleName='DEF_EASYAEROSOL'

contains

!
! Memory allocation routines
!
subroutine allocate_easyaerosol_rad(rad, row_length, rows, model_levels,       &
                                    n_wavebands)

use yomhook,         only: lhook, dr_hook
use parkind1,        only: jprb, jpim

implicit none

type (t_easyaerosol_rad), intent(in out) :: rad
integer, intent(in) :: row_length
integer, intent(in) :: rows
integer, intent(in) :: model_levels
integer, intent(in) :: n_wavebands

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character (len=*), parameter :: RoutineName = 'ALLOCATE_EASYAEROSOL_RAD'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,                          &
                        zhook_in, zhook_handle)

rad%dim1 = row_length
rad%dim2 = rows
rad%dim3 = model_levels
rad%dim4 = n_wavebands

if (.not. allocated(rad%extinction)) then
  allocate(rad%extinction(rad%dim1, rad%dim2, rad%dim3, rad%dim4))
end if
!$OMP PARALLEL DEFAULT(none)                                                   &
!$OMP private(i,j,k,l)                                                         &
!$OMP SHARED(rad)
do l=1, rad%dim4
!$OMP do SCHEDULE(STATIC)
  do k=1, rad%dim3
    do j=1, rad%dim2
      do i=1, rad%dim1
        rad%extinction(i,j,k,l) = 0.0
      end do
    end do
  end do
!$OMP end do NOWAIT
end do
!$OMP end PARALLEL

if (.not. allocated(rad%absorption)) then
  allocate(rad%absorption(rad%dim1, rad%dim2, rad%dim3, rad%dim4))
end if
!$OMP PARALLEL DEFAULT(none)                                                   &
!$OMP private(i,j,k,l)                                                         &
!$OMP SHARED(rad)
do l=1, rad%dim4
!$OMP do SCHEDULE(STATIC)
  do k=1, rad%dim3
    do j=1, rad%dim2
      do i=1, rad%dim1
        rad%absorption(i,j,k,l) = 0.0
      end do
    end do
  end do
!$OMP end do NOWAIT
end do
!$OMP end PARALLEL

if (.not. allocated(rad%asymmetry)) then
  allocate(rad%asymmetry(rad%dim1, rad%dim2, rad%dim3, rad%dim4))
end if
!$OMP PARALLEL DEFAULT(none)                                                   &
!$OMP private(i,j,k,l)                                                         &
!$OMP SHARED(rad)
do l=1, rad%dim4
!$OMP do SCHEDULE(STATIC)
  do k=1, rad%dim3
    do j=1, rad%dim2
      do i=1, rad%dim1
        rad%asymmetry(i,j,k,l) = 0.0
      end do
    end do
  end do
!$OMP end do NOWAIT
end do
!$OMP end PARALLEL


if (lhook) call dr_hook(ModuleName//':'//RoutineName,                          &
                        zhook_out, zhook_handle)
return

end subroutine allocate_easyaerosol_rad

subroutine allocate_easyaerosol_cdnc(easy_cdnc, row_length, rows, model_levels)

use yomhook,         only: lhook, dr_hook
use parkind1,        only: jprb, jpim

implicit none

type (t_easyaerosol_cdnc), intent(in out) :: easy_cdnc
integer, intent(in) :: row_length
integer, intent(in) :: rows
integer, intent(in) :: model_levels

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character (len=*), parameter :: RoutineName = 'ALLOCATE_EASYAEROSOL_CDNC'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,                          &
                        zhook_in, zhook_handle)

easy_cdnc%dim1 = row_length
easy_cdnc%dim2 = rows
easy_cdnc%dim3 = model_levels

if (.not. allocated(easy_cdnc%cdnc)) then
  allocate(easy_cdnc%cdnc(easy_cdnc%dim1, easy_cdnc%dim2, easy_cdnc%dim3))
end if
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(easy_cdnc)
do k=1, easy_cdnc%dim3
  do j=1, easy_cdnc%dim2
    do i=1, easy_cdnc%dim1
      easy_cdnc%cdnc(i,j,k) = 0.0
    end do
  end do
end do
!$OMP end PARALLEL do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,                          &
                        zhook_out, zhook_handle)
return

end subroutine allocate_easyaerosol_cdnc

!
! Memory deallocation routines
!
subroutine deallocate_easyaerosol_rad(rad)

implicit none

type (t_easyaerosol_rad), intent(in out) :: rad

if (allocated(rad%extinction)) deallocate(rad%extinction)
if (allocated(rad%absorption)) deallocate(rad%absorption)
if (allocated(rad%asymmetry))  deallocate(rad%asymmetry)

return

end subroutine deallocate_easyaerosol_rad

subroutine deallocate_easyaerosol_cdnc(easy_cdnc)

implicit none

type (t_easyaerosol_cdnc), intent(in out) :: easy_cdnc

if (allocated(easy_cdnc%cdnc)) deallocate(easy_cdnc%cdnc)

return

end subroutine deallocate_easyaerosol_cdnc

end module def_easyaerosol
