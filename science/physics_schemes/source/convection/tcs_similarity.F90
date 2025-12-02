! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Module to calculate similarity functions for warm rain tcs.
!
module tcs_similarity


use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
use ereport_mod, only: ereport
use um_types, only: real_umphys

implicit none
!
! Description:
! Module to calculate similarity functions for warm rain tcs.
!
! Method:
!   Calculate the similarity profiles using standard functions.
!   Parameters for these functions are determined from a
!   set predefined in tcs_parameters_warm based on the type of
!   convection that has been diagnosed (conv_type).
!   For specific detail see:
!   <reference to documentation to go here, once available>
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.1 programming standards.
!


! Working note: This needs to be tidied up a little bit and to make
! sure that the <= and < refer exactly to in cloud for rho and theta
! levels.

! Working note: Need to check whether these functions are on rho or
! theta levels and that they are used consistently in the rest of the
! code.

! Working note: Might be nice to get rid of the mask if possible.
! (Just need to make sure bounds and initialisations are set
! correctly in the rest of the code)



character(len=*), parameter, private :: ModuleName='TCS_SIMILARITY'

contains

subroutine calc_similarity(eta_theta, eta_rho, conv_type, sim)

use tcs_parameters_warm,       only:                                           &
   similarity_coefficient, sim_coeff_nonp_sh,                                  &
   sim_coeff_warm_cg

use tcs_class_similarity,      only:                                           &
   similarity

use tcs_common_warm,           only:                                           &
   scales

use errormessagelength_mod, only: errormessagelength

implicit none
!-------------------------------------------------------------------
! Subroutine Arguments
!-------------------------------------------------------------------
real(kind=real_umphys), intent(in) :: eta_theta(:,:)
                    ! non-dimensional height of cloud levels
real(kind=real_umphys), intent(in) :: eta_rho(:,:)
                    ! non-dimensional height of cloud levels on rho levels
integer, intent(in) :: conv_type(:)
                    ! Indicator of type of convection
                    !    1=non-precipitating shallow
                    !    2=drizzling shallow
                    !    3=warm congestus

! Sim is now declared as intent(inout) so that memory allocation
! (and deallocation) can be done in the calling routine.
type(similarity), intent(in out):: sim

!-----------------------------------------------------------------
! Variables defined locally
! Note that sim%n_xx and sim%nlev are dimensions of similarity
! functions and are used as shorthand instead of size(sim%...,n)
!-----------------------------------------------------------------
real(kind=real_umphys) :: wexpG(sim%n_xx,sim%nlev)
real(kind=real_umphys) :: wexpK(sim%n_xx,sim%nlev)
real(kind=real_umphys) :: wexpF(sim%n_xx,sim%nlev)
real(kind=real_umphys) :: wexpB(sim%n_xx,sim%nlev)

real(kind=real_umphys) :: arr1(sim%n_xx)   ! An array expression
real(kind=real_umphys) :: arr2(sim%n_xx)   ! An array expression

real(kind=real_umphys) :: zi(sim%n_xx,sim%nlev)     ! Shifted height

integer :: i,k ! loop counter

type(similarity_coefficient), pointer :: sc
                    ! coefficients for similarity functions

!
! Error reporting variables
!
integer :: ErrorStatus

character (len=errormessagelength) :: Cmessage

character (len=*), parameter :: RoutineName='CALC_SIMILARITY'

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
do i=1,sim%n_xx

  ! Select the appropriate coefficients
  select case (conv_type(i))
  case (1:2)
    ! shallow convection (NB eventually may want to be split
    !          into precipitating(1) and non-precipitating(2))
    sc=>sim_coeff_nonp_sh
  case (3:4)
    ! congestus convection (NB eventually may want to be split
    !                       into warm(3) and ice(4))
    sc=>sim_coeff_warm_cg
  case DEFAULT
    ! Should be set to one of the preceding options
    ErrorStatus = -1
    write(cmessage,'(A6,I2)') 'conv_type not recognized: ', conv_type(i)

    call Ereport(RoutineName, ErrorStatus, Cmessage )
  end select

  ! Working note: a1 and a2 need to be rationalized and the
  ! constants put into tcs_parameters_warm
  arr1(i) = sc%f0(1)*(0.2*scales%wsc_o_mb(i) - 2.0)                            &
     + 0.1*scales%wsc_o_mb(i)
  arr2(i) = (sc%f0(1)/sc%f0(2))*(0.2*scales%wsc_o_mb(i) - 2.0)

  !-----------------------------------------------------------------------
  !
  !  Functions on rho levels
  !
  !-----------------------------------------------------------------------

  ! Calculate functions
  wexpK(i,:) = 1.0 - exp(-sc%k(2)*eta_rho(i,:))
  wexpF(i,:) = 1.0 - exp(-sc%Fng(3)*eta_rho(i,:))
  wexpB(i,:) = 1.0 - exp(-sc%b(3)*eta_rho(i,:))
  wexpG(i,:) = exp(-sc%g(5)*eta_rho(i,:))

  do k=1,sim%nlev
    if (eta_rho(i,k) < 1.0 - spacing(eta_rho(i,k))) then
      sim%k_func_rho(i,k)   = sc%k(1)*(1.0 + eta_rho(i,k))*wexpK(i,k)
      sim%fng_func_rho(i,k) = 1.0 - sc%Fng(1)*(1.0+sc%Fng(2)*eta_rho(i,k))     &
         *wexpF(i,k)
      sim%b_func_rho(i,k)   = sc%b(1)*(1.0 + eta_rho(i,k)                      &
         - sc%b(2)*eta_rho(i,k)*eta_rho(i,k))                                  &
         *wexpB(i,k)
      sim%g_func_rho(i,k)   = sc%g(1)*(eta_rho(i,k)                            &
         - sc%g(2)*eta_rho(i,k)*eta_rho(i,k)                                   &
         + (sc%g(3) - sc%g(4)*eta_rho(i,k))                                    &
         * wexpG(i,k) - sc%g(3))

      sim%cmask(i,k)=1.0
      sim%cmask_rho(i,k)=1.0

      if (eta_rho(i,k) < sc%pc2_d(1)) then
        sim%pc2_detr(i,k) = sc%pc2_d(2)
      else
        sim%pc2_detr(i,k) = sc%pc2_d(2) +                                      &
           (eta_rho(i,k)-sc%pc2_d(1))/(1.0-sc%pc2_d(1))*sc%pc2_d(3)
      end if

    end if
  end do

  zi(i,:)= (eta_rho(i,:) - sc%fw(2))/sc%fw(3)
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(k)                                                               &
!$OMP SHARED(i,sim,eta_rho,sc,zi)
  do k=1,sim%nlev
    if (eta_rho(i,k) <= 1.0 + spacing(eta_rho(i,k))) then
      sim%fql_func_rho(i,k) = (1.0+eta_rho(i,k)                                &
         - sc%fql(1)*eta_rho(i,k)*eta_rho(i,k))                                &
         *( 1.0- sc%fql(2)*exp(-sc%fql(3)*eta_rho(i,k)))
      sim%gql_func_rho(i,k) = 1.0 - sc%gql(1)*exp(-sc%gql(2)*eta_rho(i,k))
      sim%fw_func_rho(i,k)  = sc%fw(1)*exp(-zi(i,k)*zi(i,k)/2.0)               &
         + sc%fw(4) + sc%fw(5)*eta_rho(i,k)                                    &
         + sc%fw(6)*eta_rho(i,k)*eta_rho(i,k)
    end if
  end do
!$OMP end PARALLEL do

  !-----------------------------------------------------------------------
  !
  !  Functions on theta levels
  !
  !-----------------------------------------------------------------------



  wexpK(i,:) = 1.0 - exp(-sc%k(2)*eta_theta(i,:))
  wexpG(i,:) = exp(-sc%g(5)*eta_theta(i,:))

  zi(i,:)=(eta_theta(i,:) - sc%fw(2))/sc%fw(3)

!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(k)                                                               &
!$OMP SHARED(i,sim,eta_theta,sc,wexpG,wexpK,arr1,arr2,zi)
  do k=1,sim%nlev
    if (eta_theta(i,k) <= 1.0 + spacing(eta_theta(i,k))) then
      sim%g_func(i,k)  = sc%g(1)*(eta_theta(i,k)                               &
         - sc%g(2)*eta_theta(i,k)*eta_theta(i,k)                               &
         + (sc%g(3) - sc%g(4)*eta_theta(i,k))*wexpG(i,k)                       &
         - sc%g(3))

      sim%k_func(i,k)  = sc%k(1)*(1.0 + eta_theta(i,k))*wexpK(i,k)

      sim%f0_func(i,k) = 1.0 - arr1(i)*eta_theta(i,k)                          &
         + arr2(i)*(1.0-exp(-sc%f0(2)*eta_theta(i,k)))

      sim%f1_func(i,k) = sc%f1(1)*(eta_theta(i,k)**sc%f1(2))

      sim%ftheta_func(i,k) = sc%fth(1)                                         &
         - sc%fth(2)*exp(-sc%fth(3)*eta_theta(i,k))
      sim%fw_func(i,k) = sc%fw(1)*exp(-zi(i,k)*zi(i,k)/2.0)                    &
         + sc%fw(4) + sc%fw(5)*eta_theta(i,k)                                  &
         + sc%fw(6)*eta_theta(i,k)*eta_theta(i,k)
    end if
  end do
!$OMP end PARALLEL do

end do
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine calc_similarity

end module tcs_similarity
