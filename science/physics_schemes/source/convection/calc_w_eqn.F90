! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Calculates parcel vertical velocity diagnostic (w)

module calc_w_eqn_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName='CALC_W_EQN_MOD'

contains

subroutine calc_w_eqn ( npnts, bterm, blowst                                   &
                      , ekp14, ekp34, zk, zkp12, zkp1                          &
                      , thek,  thekp1,  thpk,  thpkp1                          &
                      , qek,   qekp1,   qpk,   qpkp1                           &
                      , qclek, qclekp1, qclpk, qclpkp1                         &
                      , qcfek, qcfekp1, qcfpk, qcfpkp1                         &
                      , w2p_k, w2p_kp1, idx, ni)

use planet_constants_mod, only: c_virtual, g

use cv_run_mod, only: cnv_wat_load_opt
use cv_param_mod,                                                              &
    only: gamma_in_w_eqn, watload_opt, w2pi                                    &
        , cumulus_r, drag_coeff, k2_const, fix_alpha, alpha_opt                &
        , wCalcMethod, NegBuoyOpt, gamma_b, NegBuoyMinW

use yomhook,             only: lhook, dr_hook
use parkind1,            only: jprb, jpim

implicit none

!----------------------------------------------------------------------
! Description:
!   Calculates parcel vertical velocity (w) based on
!   [Simpson & Wiggert 1969]. Equation is discretised so that
!   a forward difference is used for the first layer in the parcel
!   profile. The remaining layers used centred difference.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.0 programming standards.
!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------

! Vector lengths and loop counters

integer, intent(in) :: npnts        ! Vector length (nconv)

real(kind=real_umphys), intent(in)    :: ekp14(npnts)
                                    ! Entrainment rate on k+1/4 (-)
real(kind=real_umphys), intent(in)    :: ekp34(npnts)
                                    ! Entrainment rate on k+3/4 (-)


! Masks for points where ...
logical, intent(in) :: bterm(npnts) ! ...parcels terminate in layer k+1
logical, intent(in) :: blowst(npnts)! ...stability is low enough for
                                    !    convection to occur


! Layer variables
real(kind=real_umphys), intent(in) :: zk    (npnts) ! Height of theta level k
real(kind=real_umphys), intent(in) :: zkp12 (npnts)
                                  ! Height of theta level k+1/2
real(kind=real_umphys), intent(in) :: zkp1  (npnts) ! Height of theta level k+1

! Potential Temperature (Theta) on ...
real(kind=real_umphys), intent(in) :: thek  (npnts) ! ...layer centre k   (env)
real(kind=real_umphys), intent(in) :: thekp1(npnts) ! ...layer centre k+1 (env)
real(kind=real_umphys), intent(in) :: thpk  (npnts) ! ...layer centre k   (par)
real(kind=real_umphys), intent(in) :: thpkp1(npnts) ! ...layer centre k+1 (par)

! Specific humidity on ...
real(kind=real_umphys), intent(in) :: qek   (npnts) ! ...layer centre k   (env)
real(kind=real_umphys), intent(in) :: qekp1 (npnts) ! ...layer centre k+1 (env)
real(kind=real_umphys), intent(in) :: qpk   (npnts) ! ...layer centre k   (par)
real(kind=real_umphys), intent(in) :: qpkp1 (npnts) ! ...layer centre k+1 (par)

! Liquid water content on ...
real(kind=real_umphys), intent(in) :: qclek   (npnts)
                                    ! ...layer centre k   (env)
real(kind=real_umphys), intent(in) :: qclekp1 (npnts)
                                    ! ...layer centre k+1 (env)
real(kind=real_umphys), intent(in) :: qclpk   (npnts)
                                    ! ...layer centre k   (par)
real(kind=real_umphys), intent(in) :: qclpkp1 (npnts)
                                    ! ...layer centre k+1 (par)

! Ice water content on ...
real(kind=real_umphys), intent(in) :: qcfek   (npnts)
                                    ! ...layer centre k   (env)
real(kind=real_umphys), intent(in) :: qcfekp1 (npnts)
                                    ! ...layer centre k+1 (env)
real(kind=real_umphys), intent(in) :: qcfpk   (npnts)
                                    ! ...layer centre k   (par)
real(kind=real_umphys), intent(in) :: qcfpkp1 (npnts)
                                    ! ...layer centre k+1 (par)



! Output variables
! (Parcel vertical velocity)^2 on ...
real(kind=real_umphys), intent(in out) :: w2p_k  (npnts)
                                       ! ...layer centre k   [(m/s)^2]
real(kind=real_umphys), intent(in out) :: w2p_kp1(npnts)
                                       ! ...layer centre k+1 [(m/s)^2]


! Indices of work points
integer, intent(in)  :: idx(npnts)
integer, intent(in)  :: ni


!=====================================================================

  ! Local variables
  !===================
real(kind=real_umphys) :: alpha  (npnts)
real(kind=real_umphys) :: coeff1
real(kind=real_umphys) :: dz
real(kind=real_umphys) :: thvpk
real(kind=real_umphys) :: thvek
real(kind=real_umphys) :: thvpkp1
real(kind=real_umphys) :: thvekp1
real(kind=real_umphys) :: term2
real(kind=real_umphys) :: BuoyExcess

integer :: i, m


! Local parameters so as not to use magic numbers
! (development)
!===================================================================
!
! For method use to calculate w or w^2
integer, parameter :: method_simpsonwiggert = 1
integer, parameter :: method_buoyscale      = 2

! For setting calcuation of alpha
integer, parameter :: alpha_simpsonwiggert  = 1
integer, parameter :: alpha_um              = 2
integer, parameter :: alpha_const           = 3

! For treatment of negative buoyancy occurance
integer, parameter :: NegBuoyLimit = 1
integer, parameter :: NegBuoyZero  = 2

! For treatment of water loading
integer, parameter :: watload_simpsonwiggert = 1
integer, parameter :: watload_um             = 2

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CALC_W_EQN'

!---------------------------------------------------------------------
! End of declarations
!---------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  ! Current development notes (Subject to further testing)
  ! =======================================================
  ! The following notes apply to some sensitivity tests
  ! on using a fixing value of alpha in the SCM. These are
  ! not intended as general advice or recommendations
  ! =======================================================
  ! alpha_opt = um
  ! fix_alpha = 0.01   ! This alpha dies
  ! fix_alpha = 0.005  ! This alpha seems noise creeping in, effectively dead
  ! fix_alpha = 0.002  ! Works okay
  ! fix_alpha = 0.001  ! Nice and smooth some buoyancy violations
  ! fix_alpha = 0.0001 ! Smooth, similar to none, peaks higher up, possibly
                       ! too low

w2p_kp1(:) = 0.0

select case (wCalcMethod)

case (method_simpsonwiggert)

  select case (alpha_opt)
  case (alpha_SimpsonWiggert)
    alpha(:) = (0.75*k2_const+drag_coeff)*(3.0/(8.0*cumulus_r))

  case (alpha_um)
    ! Use convection scheme entrainment coefficents
    alpha(:) =  0.5*(  ekp14(:)/(zkp12(:) - zk(:)   )                          &
                     + ekp34(:)/(zkp1(:)  - zkp12(:)) )
  case (alpha_const)
    alpha(:) = fix_alpha

  case DEFAULT
    alpha(:) = 0.0
  end select

  coeff1 = 2.0*g/(1+gamma_in_w_eqn)

  ! Note: Could use module to get the layer thicknesses,
  !       though the height in the modules have halos etc, and these
  !       pointare compressed... not a problem if the indexes are available
  !       though how to relate the halos.?
  !
  ! For now pass down height levels like the rest of the convection scheme.

  ! Give parcel at surface an initial kick to start it off?

  ! Calculate limit of Radius, to small and the velocity squared goes -ve
  ! which is not possible. i.e. the drag term (for the given seetings)
  ! is too large.

  ! Calculate min. radius required to prevent -ve w^2 for a given level

!DIR$ IVDEP
  do m=1, ni
    i = idx(m)

    ! Calculate w2p if this point is marked as convecting
    if ( .not. bterm(i) ) then

      ! Calculate Parcel/Environment Virtual Theta for level k
      !-------------------------------------------------------
      if (cnv_wat_load_opt == 1) then

        ! How to apply water loading?
        if (watload_opt == watload_simpsonwiggert) then
          thvpk = (thpk(i)*(1.0 + c_virtual*qpk(i)))*(1-qclpk(i)-qcfpk(i))
          thvek = (thek(i)*(1.0 + c_virtual*qek(i)))*(1-qclek(i)-qcfek(i))
        else if (watload_opt == watload_um) then
          thvpk =  thpk(i)*(1.0 + c_virtual*qpk(i)-qclpk(i)-qcfpk(i))
          thvek =  thek(i)*(1.0 + c_virtual*qek(i)-qclek(i)-qcfek(i))
        end if

      else
        ! Do not apply water loading
        thvpk = thpk(i) * (1.0 + c_virtual*qpk(i))
        thvek = thek(i) * (1.0 + c_virtual*qek(i))
      end if


      ! Take into account of water loading in parcel
      ! Assume no water loading for now
      term2 = (thvpk/thvek) - 1.0


      ! This is the first level of a convective layer?
      ! If so, assume initial w=w2pi for this convecting layer
      if (blowst(i)) w2p_k(i) = w2pi


      ! Calculate w2p on level k+1 using forward difference
      dz = zkp1(i) - zk(i)
      w2p_kp1(i) = dz*(coeff1*term2 - alpha(i)*w2p_k(i)) + w2p_k(i)

    end if ! bterm test
  end do

case (method_buoyscale)

  ! Calculate w via scaled buoyancy on level k+1
  ! This can only be done here because the call to
  ! calc_w_eqn is at the end of convec2 when the input
  ! parcel properties from convec2 at level k+1 are now
  ! known

  !=======================================================
  ! NOTE: May need to add a check for the initiation of
  !       convection so as to apply scaling to initial
  !       buoyancy excess
  !=======================================================
  !
!DIR$ IVDEP
  do m=1, ni
    i = idx(m)

    ! Calculate w2p if this point is marked as convecting
    if ( .not. bterm(i) ) then

      ! Calculate Parcel/Environment Virtual Theta for level k+1
      !---------------------------------------------------------
      if (cnv_wat_load_opt == 1) then
        thvpkp1 = thpkp1(i) * (1.0+c_virtual*qpkp1(i)-qclpkp1(i)-qcfpkp1(i))
        thvekp1 = thekp1(i) * (1.0+c_virtual*qekp1(i)-qclekp1(i)-qcfekp1(i))
      else
        thvpkp1 = thpkp1(i) * (1.0+c_virtual*qpkp1(i))
        thvekp1 = thekp1(i) * (1.0+c_virtual*qekp1(i))
      end if

      BuoyExcess = g * (thvpkp1 - thvekp1) / thvekp1

      w2p_kp1(i) = (gamma_b*BuoyExcess) * (gamma_b*BuoyExcess)

      if (BuoyExcess < 0.0) then
        select case(NegBuoyOpt)
        case (NegBuoyLimit)
          w2p_kp1(i) = (gamma_b*NegBuoyMinW) * (gamma_b*NegBuoyMinW)
        case (NegBuoyZero)
          w2p_kp1(i) = 0.0
        end select
      end if

    end if
  end do

end select

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine calc_w_eqn

end module calc_w_eqn_mod
