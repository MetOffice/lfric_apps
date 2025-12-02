! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine to calculate the new effective precipitation fraction
! when combining 2 precipitation masses with different fractions.
! Assumes the two masses are maximally overlapped.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

module lsp_combine_precfrac_mod

implicit none

character(len=*), parameter, private :: ModuleName='LSP_COMBINE_PRECFRAC_MOD'

contains

subroutine lsp_combine_precfrac( points,                                       &
                                 mass1, mass2, frac1, frac2 )

use lsprec_mod,       only: zero, small_number
use um_types,         only: real_lsprec
use mphys_inputs_mod, only: i_update_precfrac, i_homog_areas, i_sg_correl

! Dr Hook Modules
use yomhook,          only: lhook, dr_hook
use parkind1,         only: jprb, jpim

implicit none

! Number of points
integer, intent(in) :: points

! Masses of the 2 precip masses being combined
real(kind=real_lsprec), intent(in) :: mass1(points)
real(kind=real_lsprec), intent(in) :: mass2(points)

! Area fractions of the 2 precip masses being combined
! (frac1 gets overwritten with the updated combined fraction)
real(kind=real_lsprec), intent(in out) :: frac1(points)
real(kind=real_lsprec), intent(in) :: frac2(points)

! Temporary storage for re-used term in formula
real(kind=real_lsprec) :: tmp, tmp1, tmp2

! Loop counter
integer :: i


! Stuff for DrHook
integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_COMBINE_PRECFRAC'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


! A representative mean mass within the combined precipitation fraction fc
! can be written as a mass-weighted mean mass:
!
! <m> / fc = <m*m> / <m>
!
! where < > denotes the area-weighted mean.
! Rearranging, we can derive an expression for an effective combined
! fraction, so-as to yield the desired representative in-cloud value:
!
! fc = <m>^2 / <m^2>    (1)
!
! We have 2 fractions f1, f2 each with an associated mass m1, m2,
! and we assume the two are maximally correlated / overlapped.

if ( i_update_precfrac == i_homog_areas ) then
  ! Assume the two masses are homogeneous within their respective areas

  do i = 1, points
    if ( mass2(i) > zero ) then
      ! Only update frac1 if mass2 is nonzero
      if ( mass1(i) > zero ) then
        ! Both of the masses are non-zero...

        ! If fraction f2 > f1, then we have a region of area f2-f1 which
        ! only contains the mass associated with f2, and a region f1
        ! which contains the sum of the two masses.  Therefore we can write
        ! the new area-weighted mean mass as:
        !
        ! <m> = (f2 - f1) (m2/f2) + f1 (m1/f1 + m2/f2)
        !     = f2 (m2/f2)        + f1 (m1/f2)
        !     = m1 + m2    (2)
        !
        ! We can write the area-weighted mean of m^2 as:
        !
        ! <m^2> = (f2 - f1) (m2/f2)^2 + f1 (m1/f1 + m2/f2)^2
        !       = f2 (m2/f2)^2        + f1 (m1/f1)^2   + 2 f1 (m1/f1) (m2/f2)
        !       = 1/f2 m2^2  +  1/f1 m1^2  +  1/f2 2 m1 m2
        !       = 1/f2 (m1 + m2)^2 + (1/f1 - 1/f2) m1^2    (3)
        !
        ! Substituting (2) and (3) into (1), we have:
        !
        ! fc = (m1 + m2)^2 / ( 1/f2 (m1 + m2)^2 + (1/f1 - 1/f2) m1^2 )
        !    = f1 f2 / ( f1 + (f2 - f1) ( m1/(m1 + m2) )^2 )    (4)
        !
        ! This equation (4) is used to set the new effective combined
        ! fraction when f2 > f1.  In the other case (f2 < f1), the formula
        ! is equivalent but with the labels 1 and 2 swapped.

        ! Note: the brackets in the equations below are important, as
        ! the order in which the terms are calculated matters;
        ! floating point underflow can occasionally result otherwise,
        ! wrongly yielding frac1=0.0
        if ( frac2(i) > frac1(i) ) then
          ! Case where frac1 is fully contained within frac2
          tmp = mass1(i) / ( mass1(i) + mass2(i) )
          frac1(i) = frac2(i) * ( frac1(i)                                     &
            / max( frac1(i) + ( frac2(i) - frac1(i) ) * tmp*tmp,               &
                   small_number ) )
        else
          ! Case where frac2 is fully contained within frac1
          tmp = mass2(i) / ( mass1(i) + mass2(i) )
          frac1(i) = frac1(i) * ( frac2(i)                                     &
            / max( frac2(i) + ( frac1(i) - frac2(i) ) * tmp*tmp,               &
                   small_number ) )
        end if

      else  ! mass1(i) == zero
        frac1(i) = frac2(i)
      end if  ! mass1(i) == zero
    else  ! mass2(i) == zero
      ! Reset frac to zero if both masses are zero
      if ( .not. mass1(i) > zero )  frac1(i) = zero
    end if  ! mass2(i) == zero
  end do  ! i = 1, points

else if ( i_update_precfrac == i_sg_correl ) then
  ! Account for inhomgeneity of masses, with spatial correlation

  ! Recall our formula for the effective fraction, as a function of
  ! some mass which varies over a spatial dimension x:
  ! fc = <m(x)>^2 / <m(x)^2>    (1)
  !
  ! The combined grid-mean of m is trivially
  ! <m> = < m1(x) + m2(x) >
  ! <m> = <m1> + <m2>        (2)
  !
  ! The grid-mean of m^2 is then:
  ! <m^2> = < ( m1(x) + m2(x) )^2 >
  !       = < m1(x)^2 + 2 m1(x) m2(x) + m2(x)^2 >
  !       = <m1^2> + 2 <m1 m2> + <m2^2>             (3)
  !
  ! From (1) we have
  ! <m1^2> = <m1>^2 / f1
  ! <m2^2> = <m2>^2 / f2
  ! But the term <m1 m2> depends on the shapes of the sub-grid spatial
  ! distributions of m1(x) and m2(x), and how they are correlated.
  !
  ! If they are uncorrelated, then we will have:
  ! <m1 m2>  =  <m1> <m2>
  !
  ! If they are maximally correlated, then based on dimensional and symmetry
  ! arguments, we parameterise it as:
  !
  ! <m1 m2>  =  <m1> <m2> / sqrt(f1 f2)  =  sqrt( <m1^2> <m2^2> )
  !
  ! Therefore (3) becomes:
  !
  ! <m^2> = <m1>^2 / f1  +  2 <m1> <m2> / sqrt(f1 f2)  +  <m2>^2 / f2
  !
  ! Note that our choice for <m1 m2> conveniently allows this to be factorised:
  !
  ! <m^2> = ( <m1> / sqrt(f1) + <m2> / sqrt(f2) )^2    (4)
  !
  ! Substituting (2) and (4) into (1), we have:
  ! fc = ( <m1> + <m2> )^2 / ( <m1> / sqrt(f1) + <m2> / sqrt(f2) )^2
  !    = f1 f2 / ( m1/(m1+m2) sqrt(f2) + m2/(m1+m2) sqrt(f1) )^2      (5)
  !
  ! Note this corresponds to just taking the mass-weighted mean of 1/sqrt(f).

  do i = 1, points
    if ( mass2(i) > zero ) then
      ! Only update frac1 if mass2 is nonzero
      if ( mass1(i) > zero ) then
        ! Both of the masses are non-zero...

        ! From (5) above:
        ! fc = f1 f2 / ( m1/(m1+m2) sqrt(f2) + m2/(m1+m2) sqrt(f1) )^2
        tmp1 = mass1(i) / ( mass1(i) + mass2(i) )
        tmp2 = mass2(i) / ( mass1(i) + mass2(i) )
        tmp = tmp1 * sqrt(frac2(i)) + tmp2 * sqrt(frac1(i))
        if ( frac2(i) > frac1(i) ) then
          frac1(i) = frac2(i) * ( frac1(i) / max( tmp*tmp, small_number ) )
        else
          frac1(i) = frac1(i) * ( frac2(i) / max( tmp*tmp, small_number ) )
        end if

      else  ! mass1(i) == zero
        frac1(i) = frac2(i)
      end if  ! mass1(i) == zero
    else  ! mass2(i) == zero
      ! Reset frac to zero if both masses are zero
      if ( .not. mass1(i) > zero )  frac1(i) = zero
    end if  ! mass2(i) == zero
  end do  ! i = 1, points

end if  ! ( i_update_precfrac )


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_combine_precfrac


end module lsp_combine_precfrac_mod
