! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate the effect of convection upon the large-scale cloud fractions.
!
module pc2_environ_mod

use um_types, only: real_umphys

implicit none

! Description:
! Calculate the effect of convection upon the large-scale cloud fractions.
!
! Method:
!   See UM Documentation paper No 27
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_cloud
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.

character(len=*), parameter, private :: ModuleName='PC2_ENVIRON_MOD'

contains

! Subroutine Interface:
subroutine pc2_environ (   npnts,                                              &
                    qclek, qclekp1, qcfek, qcfekp1,                            &
                    cflek, cflekp1,  cffek,  cffekp1,                          &
                    qclpk,          qcfpk,                                     &
                    dqclek, dqcfek,                                            &
                    dqclekp1, dqcfekp1,                                        &
                    l_q_interact,                                              &
                    bterm,                                                     &
                    ! Out
                    dcflek, dcffek, dbcfek,                                    &
                    dcflekp1, dcffekp1, dbcfekp1,                              &
                    !Indirect indexing
                    idx,ni)

use cv_run_mod, only: eff_dcfl, eff_dcff

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------

integer,intent(in) :: npnts         ! Number of points

real(kind=real_umphys),intent(in) :: qclek(npnts) ! Env. qcl in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: qclekp1(npnts)
                                    ! Env. qcl in layer k+1 (kg/kg)
real(kind=real_umphys),intent(in) :: qcfek(npnts) ! Env. qcf in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: qcfekp1(npnts)
                                    ! Env. qcf in layer k+1 (kg/kg)
real(kind=real_umphys),intent(in) :: cflek(npnts)
                                    ! Env. liquid cloud volume fraction
                                    ! in layer k
real(kind=real_umphys),intent(in) :: cflekp1(npnts)
                                    ! Env. liquid cloud volume fraction
                                    ! in layer k+1
real(kind=real_umphys),intent(in) :: cffek(npnts)
                                    ! Env. frozen cloud volume fraction
                                    ! in layer k
real(kind=real_umphys),intent(in) :: cffekp1(npnts)
                                    ! Env. frozen cloud volume fraction
                                    ! in layer k+1
real(kind=real_umphys),intent(in) :: qclpk(npnts) ! Par. qcl in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: qcfpk(npnts) ! Par. qcf in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: dqclek(npnts)
                                    ! Increment to qcl in layer k (kg/kg/s)
real(kind=real_umphys),intent(in) :: dqcfek(npnts)
                                    ! Increment to qcf in layer k (kg/kg/s)
real(kind=real_umphys),intent(in) :: dqclekp1(npnts)
                                    ! Increment to qcl in layer k+1 (kg/kg/s)
real(kind=real_umphys),intent(in) :: dqcfekp1(npnts)
                                    ! Increment to qcf in layer k+1 (kg/kg/s)

logical,intent(in) :: l_q_interact  ! True if PC2 is switched on
logical,intent(in) :: bterm(npnts)  ! Mask for parcels which terminate
                                    ! in layer k+1

!----------------------------------------------------------------------
! Variables which are output
!----------------------------------------------------------------------
! Convection increments to model fields at level k
real(kind=real_umphys),intent(in out) :: dcflek(npnts)
                                    ! Increment to liquid cloud volume fraction
                                    ! in layer k
real(kind=real_umphys),intent(in out) :: dcffek(npnts)
                                    ! Increment to frozen cloud volume fraction
                                    ! in layer k
real(kind=real_umphys),intent(in out) :: dbcfek(npnts)
                                    ! Increment to total cloud volume fraction
                                    ! in layer k
! Convection increments to model fields at level k+1
real(kind=real_umphys),intent(in out) :: dcflekp1(npnts)
                                    ! Increment to liquid cloud volume fraction
                                    ! in layer k+1
real(kind=real_umphys),intent(in out) :: dcffekp1(npnts)
                                    ! Increment to frozen cloud volume fraction
                                    ! in layer k+1
real(kind=real_umphys),intent(in out) :: dbcfekp1(npnts)
                                    ! Increment to total cloud volume fraction
                                    ! in layer k+1


! Indices of work points
integer, intent(in)  :: idx(npnts)
integer, intent(in)  :: ni

!-----------------------------------------------------------------------
! Variables that are defined locally
!-----------------------------------------------------------------------

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='PC2_ENVIRON'

integer :: i,m              ! loop counter

real(kind=real_umphys) :: denom ! Denominator in cloud increment calculation

! Parameters
real(kind=real_umphys), parameter :: ls0 = 5.0e-5  ! Minimum value for ls - l

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

if (l_q_interact) then
!DIR$ IVDEP
  do m=1, ni
    i = idx(m)
    !-----------------------------------------------------------------------
    !   Increment to Liquid Cloud Volume Fraction
    !-----------------------------------------------------------------------
        ! Increment at level k
    denom = max( (qclpk(i)-qclek(i) ), ls0)

    dcflek(i) = eff_dcfl * dqclek(i)*(1.0 - cflek(i)) / denom

    ! Increment at level k+1
    ! For efficiency only calc at termination. This test could be removed as
    ! the k+1 increments will always be overwritten by the k'th increments
    ! on the next iteration of the level loop except at termination.
    if (bterm(i)) then
      ! The use of qclpk at level k+1 is intended because qclpkp1 will
      ! always be 0.0 since convection terminated at level k.
      denom = max( (qclpk(i)-qclekp1(i) ), ls0)

      dcflekp1(i) = eff_dcfl * dqclekp1(i)*(1.0 - cflekp1(i)) / denom
    else
      dcflekp1(i) = 0.0
    end if

    !-----------------------------------------------------------------------
    !   Increment to Frozen Cloud Volume Fraction
    !-----------------------------------------------------------------------
        ! Increment at level k
    denom = max( (qcfpk(i)-qcfek(i) ), ls0)

    dcffek(i) = eff_dcff * dqcfek(i)*(1.0 - cffek(i)) / denom

    ! Increment at level k+1
    if (bterm(i)) then
      ! The use of qcfpk at level k+1 is intended because qcfpkp1 will
      ! always be 0.0 since convection terminated at level k.
      denom = max( (qcfpk(i)-qcfekp1(i) ), ls0)

      dcffekp1(i) = eff_dcff * dqcfekp1(i)*(1.0 - cffekp1(i)) / denom
    else
      dcffekp1(i) = 0.0
    end if

    !
    !-----------------------------------------------------------------------
    !   Increment to Total Cloud Volume Fraction
    !-----------------------------------------------------------------------
        ! Increment at level k
    dbcfek(i) = dcflek(i) + dcffek(i)

    ! Increment at level k+1
    if (bterm(i)) then
      dbcfekp1(i) = dcflekp1(i) + dcffekp1(i)
    else
      dbcfekp1(i) = 0.0
    end if

  end do  !loop over npnts

else
  !-----------------------------------------------------------------------
  !   Not PC2 so set all cloud fraction increments to zero
  !-----------------------------------------------------------------------
!DIR$ IVDEP
  do m=1, ni
    i = idx(m)
    dcflek(i)     = 0.0
    dcflekp1(i)   = 0.0
    dcffek(i)     = 0.0
    dcffekp1(i)   = 0.0
    dbcfek(i)     = 0.0
    dbcfekp1(i)   = 0.0
  end do

end if !l_q_interact

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine pc2_environ
end module pc2_environ_mod
