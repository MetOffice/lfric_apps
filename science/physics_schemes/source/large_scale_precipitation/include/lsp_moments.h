! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Conversion between moments of PSD

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Purpose:
!   Obtains a specified moment of the in-cloud  ice particle size
!   distribution given grid box mean input ice water content.

! Method:
!   Follows the generic ice particle size distribution
!   calculation methods, based on the method of Field et al 2007. Snow Size
!   Distribution Parameterization for Midlatitude and Tropical Ice Clouds.
!   J. Atmos. Sci., 64, pp 4346-4365.

!   In the midlatitudes, method b) should relax to method a), however
!   for comparisons and completeness, both methods have been included here.
!   The extra flexibility will allow selection of the appropriate routine,
!   but method b) is currently always used operationally.

! Description of Code:
!   Fortran 90.
!   This code is written to UMDP3 programming standards.

!   Documentation: UMDP 26.

!   This subroutine calculates the quantity p_moment_out = M(n_out)
!   given the quantity (rho qcf / cfice)=ai*M(bi) where M(bi) is the
!   bi moment of the particle size distribution (PSD), corresponding
!   to the mass diameter relationship m(D) = ai D**bi.

!   It first calculates M(2) given the input quantities rho, qcf and ai
!   and then uses M(2) to calculate the output quantity p_out*M(out),
!   which can be used directly in the transfer calculations.
!   The conversion is from Field et al and gives
!   M(n)=a(n,Tc)*M(2)**b(n,Tc) where a and b are specified parameters.

!   This subroutine can be used to convert directly from M(x) to M(y)
!   by inputting ai=1, bi=x, n_out=y, rho=[array of ones], qcf=M(x),
!   cficei=[array of ones] and outputting moment_out=M(y).


! Subroutine Arguments

integer, intent(in) ::                                                         &
points
! Number of points to calculate

real (kind=prec), intent(in) ::                                                &
n_out,                                                                         &
! Moment of the PSD to be output
rho(points),                                                                   &
! Air density / kg m-3
t(points),                                                                     &
! Temperature / K
qcf(points),                                                                   &
! Ice water content / kg kg-1
cficei(points)
! 1 / ice cloud fraction

real (kind=prec), intent(out) ::                                               &
moment_out(points)
! n_out moment of the in-cloud
! particle size distn

! Local Variables

integer ::                                                                     &
i, j
! Loop counter for points

real (kind=prec) ::                                                            &
one_over_ai,                                                                   &
! 1 / ai
a_tot_1, b_tot_1, c_tot_1, a_tot_2, b_tot_2, c_tot_2,                          &
! Totals of conversion factors for global
! version
temp_in(points), tc_1_in(points) , tc_2_in(points),                            &
temp_out(points), tc_1_out(points) , tc_2_out(points),                         &
! temporaries for vector computations
ai,bi
! Local versions of ai and bi at the correct precision

! The following values are from Field et al (2007)
! and represent the conversion parameters between moments.
real (kind=prec), parameter::                                                  &
  gl_a(3) = [13.6078_prec, -7.76062_prec, 0.478694_prec ]
real (kind=prec), parameter::                                                  &
  gl_b(3) = [-0.0360722_prec, 0.0150830_prec, 0.00149453_prec ]
real (kind=prec), parameter::                                                  &
  gl_c(3) = [0.806856_prec, 0.00581022_prec, 0.0456723_prec ]

real (kind=prec), parameter :: smallnum = 2.2e-14_prec
! Small value used in if tests


integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

! Start the subroutine

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!Get correct precision versions of ai and bi
ai = real(ai_mod, kind=prec)
bi = real(bi_mod, kind=prec)

one_over_ai = 1.0_prec / ai

! calculate loop invariant quantities

a_tot_1 = gl_a(1)+gl_a(2)*bi+gl_a(3)*bi**2
b_tot_1 = gl_b(1)+gl_b(2)*bi+gl_b(3)*bi**2
c_tot_1 = gl_c(1)+gl_c(2)*bi+gl_c(3)*bi**2

! Now do the exponential of a_tot (save variables)
a_tot_1 = exp(a_tot_1)
c_tot_1 = 1.0_prec/c_tot_1

a_tot_2 = gl_a(1)+gl_a(2)*n_out+gl_a(3)*n_out**2
b_tot_2 = gl_b(1)+gl_b(2)*n_out+gl_b(3)*n_out**2
c_tot_2 = gl_c(1)+gl_c(2)*n_out+gl_c(3)*n_out**2

a_tot_2 = exp(a_tot_2)

!--------------------------------------------------
! Use the global version of the particle size distn
! Field et al (2007) J. Atmos. Sci.
!--------------------------------------------------

!-----------------------------------------------
! Form the bi moment of the ice particle size
! distribution from the ice water content.
!-----------------------------------------------
j = 0
do i = 1, points
  ! pack into vector
  if (qcf(i) > smallnum) then
    j = j + 1
    tc_1_in(j) = (t(i) - zerodegc)*b_tot_1
    tc_2_in(j) = (t(i) - zerodegc)*b_tot_2
    temp_in(j) = rho(i) * qcf(i) * cficei(i) * one_over_ai
    temp_in(j) = max(temp_in(j),0.0_prec)
  end if
end do

call exp_v(j, tc_1_in, tc_1_out)
call exp_v(j, tc_2_in, tc_2_out)

tc_1_in(1:j) = a_tot_1*tc_1_out(1:j)
call oneover_v(j, tc_1_in, tc_1_out)
temp_in(1:j) = temp_in(1:j)*tc_1_out(1:j)

!-----------------------------------------------
! Calculate the second moment of the PSD
!-----------------------------------------------
call powr_v(j,temp_in, c_tot_1, temp_out)
temp_in(1:j) = temp_out(1:j)
!-----------------------------------------------
! Calculate the n_out moment of the PSD
!-----------------------------------------------
call powr_v(j,temp_in, c_tot_2, temp_out)

temp_out(1:j) = a_tot_2*tc_2_out(1:j)*temp_out(1:j)
j = 0

do i = 1, points
  ! unpack
  if (qcf(i) > smallnum) then
    j = j + 1
    moment_out(i)=temp_out(j)
  else  ! qcf > 0
    ! Set the output moment to zero
    moment_out(i) = 0.0_prec
  end if  ! qcf > 0
end do  ! i

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
