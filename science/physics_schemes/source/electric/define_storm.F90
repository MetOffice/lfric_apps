! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Purpose: Within the electric scheme, this routine defines whether
! a storm exists within each model column based on a given definition
! of 'a storm'. It then sets a 2D logical (storm_field) with the
! .true. condition being used to determine a storm point.

! The current definition for a storm is based entirely
! on graupel water path (gwp) but the code is flexible to add
! other definitions in future.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Electric

module define_storm_mod

use atm_fields_bounds_mod,  only: tdims
use electric_inputs_mod,    only: storm_definition, graupel_only,              &
                                  graupel_and_ice, gwp_thresh, tiwp_thresh

! Dr Hook modules
use yomhook,                only: lhook, dr_hook
use parkind1,               only: jprb, jpim

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName='DEFINE_STORM_MOD'

contains

subroutine define_storm( gwp, tiwp, storm_field, nspts )

implicit none

real(kind=real_umphys), intent(in)  :: gwp( tdims%i_start : tdims%i_end,       &
                                            tdims%j_start : tdims%j_end )
! Graupel water path [kg m-3]

real(kind=real_umphys), intent(in)  :: tiwp( tdims%i_start : tdims%i_end,      &
                                             tdims%j_start : tdims%j_end )
! Total ice water path [kg m-3]

logical, intent(in out) :: storm_field ( tdims%i_start : tdims%i_end,          &
                                         tdims%j_start : tdims%j_end )

integer, intent(out)   :: nspts ! number of storm points

!==============================================================
! Local variables
!==============================================================

integer :: i, j

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='DEFINE_STORM'

!==================================================================
! Start the subroutine
!==================================================================

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! This is a very simple definition so far ...

nspts = 0

select case (storm_definition)

case (graupel_only)

  ! This is the original version of the scheme, which defines a storm solely
  ! based on whether there is graupel in the model above a given threshold.
  ! However, as some schemes (e.g. McCaul) also use ice water path to
  ! generate lightning, this could mean that some lightning is missed.

!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none) private(i,j)                  &
!$OMP REDUCTION(+:nspts) SHARED(tdims,gwp,gwp_thresh,storm_field)
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      if ( gwp(i,j) > gwp_thresh ) then
        storm_field(i,j) = .true.
        nspts = nspts + 1
      end if
    end do ! i
  end do  ! j
!$OMP end PARALLEL do

case (graupel_and_ice)

  ! This newer version defines a storm when either there is graupel water path
  ! above a set threshold or ice water path above a threshold. It is thought
  ! to be more robust with the McCaul scheme and avoids issues with only storm
  ! cores producing lightning.

!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none) private(i,j)                  &
!$OMP REDUCTION(+:nspts) SHARED(tdims,gwp,gwp_thresh,tiwp,tiwp_thresh,         &
!$OMP storm_field)
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      if ( gwp(i,j) > gwp_thresh .or. tiwp(i,j) > tiwp_thresh ) then
        storm_field(i,j) = .true.
        nspts = nspts + 1
      end if
    end do ! i
  end do  ! j
!$OMP end PARALLEL do

end select

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return

end subroutine define_storm

end module define_storm_mod
