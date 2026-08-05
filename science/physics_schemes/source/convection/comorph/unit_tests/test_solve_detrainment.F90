! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

#if !defined(LFRIC)
! Standalone unit test for the iterative solution of the detrainment rate
! used in comorph.
program test_solve_detrainment

use comorph_constants_mod, only: real_cvprec, one
use solve_detrainment_mod, only: solve_detrainment

implicit none

! Number of points in list
integer, parameter :: n_points = 7

! Inputs and outputs for the tested routine, solve_detrainment
! (see solve_detrainment.F90 for more detail).
integer :: nc
integer :: index_ic(n_points)

real(kind=real_cvprec) :: mean_buoy(n_points)
real(kind=real_cvprec) :: core_buoy(n_points)
real(kind=real_cvprec) :: power(n_points)
real(kind=real_cvprec) :: delta_tv_k(n_points)

real(kind=real_cvprec) :: x_edge(n_points)
real(kind=real_cvprec) :: frac(n_points)

! Loop counter
integer :: ic


! Set list of indices of work points to just list all the points
nc = n_points
do ic = 1, n_points
  index_ic(ic) = ic
end do

! Below are the input data from various "problem points", where the iteration
! failed to converge when using previous versions of the iterative method.
! This unit test retains and tests these, to ensure future changes won't
! break these cases.

mean_buoy(1)  = 2.576858e-03
core_buoy(1)  = 2.637893e-03
power(1)      = 4.000000e+00
delta_tv_k(1) = 9.849160e-01

mean_buoy(2)  = 4.278708e-03
core_buoy(2)  = 4.400779e-03
power(2)      = 4.000000e+00
delta_tv_k(2) = 1.108564e-01

mean_buoy(3)  = 8.896179e-01
core_buoy(3)  = 8.938599e-01
power(3)      = 1.050316e+00
delta_tv_k(3) = 1.004896e+00

mean_buoy(4)  = 5.099487e-02
core_buoy(4)  = 5.108643e-02
power(4)      = 4.000000e+00
delta_tv_k(4) = 2.843389e-01

mean_buoy(5)  = 1.980591e-02
core_buoy(5)  = 1.983643e-02
power(5)      = 4.000000e+00
delta_tv_k(5) = 4.354659e-02

mean_buoy(6)  = 4.059449e-02
core_buoy(6)  = 4.068604e-02
power(6)      = 2.917735e+00
delta_tv_k(6) = 1.199562e-01

mean_buoy(7)  = 7.882690e-02
core_buoy(7)  = 7.885742e-02
power(7)      = 0.000000e+00
delta_tv_k(7) = 3.343286e+00

! Initialise non-detrained fraction and associated edge position to 1.
do ic = 1, n_points
  x_edge(ic) = one
  frac(ic)   = one
end do

! Call comorph's detrainment solve routine.
call solve_detrainment( n_points, nc, index_ic,                                &
                        mean_buoy, core_buoy,                                  &
                        power, delta_tv_k,                                     &
                        x_edge, frac )


end program test_solve_detrainment
#endif
