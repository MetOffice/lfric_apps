! *********************************COPYRIGHT************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *********************************COPYRIGHT************************************
! Some of the content of this file has been produced with the assistance of
! Anthropic Claude Opus 5 (Claude Code).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: atmos_service_vera
!
! Part of the Vera scheme for diagnosing visual range.
!
! This module uses log-linear interpolation to compute the extinction
! effeciency, Qext, from a lookup table as a function of particle radius.
!

module vera_mie_mod

use vera_kind_mod, only: wp => vera_real, wi => vera_integer

implicit none

! Description:
!   This module is used by the Vera visibility scheme, and computes the
!   extinction efficiency Qext using a lookup table of values computed
!   using simple Mie scattering.
!
! Method:
!   This module comprises a lookup table of Qext values and particle
!   radii. The table covers particle radii over the range
!   radius=[ 0.001, 10]microns, i.e. four orders of magnitude. The
!   values of the radii tabulated are logarithmically spaced, with
!   100 values per decade, i.e. there are a total of 401 radii
!   entries in the lookup table, together with their corresponding
!   Qext values.
!
!
!   There is a single subroutine that reads this table:
!
!     vera_blumel
!       Computes the Blumel approximation to Qext for very small
!       particles. The size parameter, xx, is defined as
!
!         xx = 2.pi.r. / lambda
!
!       where r is the particle radius and lambda is the wavelength of the
!       light. The Blumel approximation is valid for xx < 1, which for
!       visibile light (lambda=550nm), corresponds to r < 0.088 micrometres.
!
!     vera_mie_lookup
!       Reads from the lookup table, using log-linear interpolation,
!       i.e. logarithmic in particle radius and linear in Qext.
!
!       The table tabulates values of Qext for particles within the
!       radius range radius=[ 0.001, 10]microns. Outside of this range
!       the returned Qext are set to:
!
!         radius < 0.001 microns, Qext computed using Blumel approximation
!
!         radius > 10    microns, large particle limit Qext = 2
!
!   For more detail, please refer to the Vera user guide.
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 standards.

! name of this module
character(len=*), parameter, private :: ModuleName='VERA_MIE_MOD'

private

! make the subroutines in this module Public,
! they are all called by vera_scheme_mod.F90
public :: vera_mie_lookup, vera_blumel

! number of values in the Qext lookup table
!
! this is hardwired and is dependant on how many entries
! are in the lookup table
integer (wi), parameter :: n_mie = 401_wi

! lookup table values of particle raddi
integer (wi)            :: mie_radii_index
real    (wp)            :: mie_radii(1:n_mie)

data (mie_radii(mie_radii_index), mie_radii_index =    1,   40 )   /           &
 1.00000e-09_wp,  1.02329e-09_wp,  1.04713e-09_wp,  1.07152e-09_wp,            &
 1.09648e-09_wp,  1.12202e-09_wp,  1.14815e-09_wp,  1.17490e-09_wp,            &
 1.20226e-09_wp,  1.23027e-09_wp,  1.25893e-09_wp,  1.28825e-09_wp,            &
 1.31826e-09_wp,  1.34896e-09_wp,  1.38038e-09_wp,  1.41254e-09_wp,            &
 1.44544e-09_wp,  1.47911e-09_wp,  1.51356e-09_wp,  1.54882e-09_wp,            &
 1.58489e-09_wp,  1.62181e-09_wp,  1.65959e-09_wp,  1.69824e-09_wp,            &
 1.73780e-09_wp,  1.77828e-09_wp,  1.81970e-09_wp,  1.86209e-09_wp,            &
 1.90546e-09_wp,  1.94984e-09_wp,  1.99526e-09_wp,  2.04174e-09_wp,            &
 2.08930e-09_wp,  2.13796e-09_wp,  2.18776e-09_wp,  2.23872e-09_wp,            &
 2.29087e-09_wp,  2.34423e-09_wp,  2.39883e-09_wp,  2.45471e-09_wp /
data (mie_radii(mie_radii_index), mie_radii_index =   41,   80 )   /           &
 2.51189e-09_wp,  2.57040e-09_wp,  2.63027e-09_wp,  2.69153e-09_wp,            &
 2.75423e-09_wp,  2.81838e-09_wp,  2.88403e-09_wp,  2.95121e-09_wp,            &
 3.01995e-09_wp,  3.09030e-09_wp,  3.16228e-09_wp,  3.23594e-09_wp,            &
 3.31131e-09_wp,  3.38844e-09_wp,  3.46737e-09_wp,  3.54813e-09_wp,            &
 3.63078e-09_wp,  3.71535e-09_wp,  3.80189e-09_wp,  3.89045e-09_wp,            &
 3.98107e-09_wp,  4.07380e-09_wp,  4.16869e-09_wp,  4.26580e-09_wp,            &
 4.36516e-09_wp,  4.46684e-09_wp,  4.57088e-09_wp,  4.67735e-09_wp,            &
 4.78630e-09_wp,  4.89779e-09_wp,  5.01187e-09_wp,  5.12861e-09_wp,            &
 5.24807e-09_wp,  5.37032e-09_wp,  5.49541e-09_wp,  5.62341e-09_wp,            &
 5.75440e-09_wp,  5.88844e-09_wp,  6.02560e-09_wp,  6.16595e-09_wp /
data (mie_radii(mie_radii_index), mie_radii_index =   81,  120 )   /           &
 6.30957e-09_wp,  6.45654e-09_wp,  6.60693e-09_wp,  6.76083e-09_wp,            &
 6.91831e-09_wp,  7.07946e-09_wp,  7.24436e-09_wp,  7.41310e-09_wp,            &
 7.58578e-09_wp,  7.76247e-09_wp,  7.94328e-09_wp,  8.12831e-09_wp,            &
 8.31764e-09_wp,  8.51138e-09_wp,  8.70964e-09_wp,  8.91251e-09_wp,            &
 9.12011e-09_wp,  9.33254e-09_wp,  9.54993e-09_wp,  9.77237e-09_wp,            &
 1.00000e-08_wp,  1.02329e-08_wp,  1.04713e-08_wp,  1.07152e-08_wp,            &
 1.09648e-08_wp,  1.12202e-08_wp,  1.14815e-08_wp,  1.17490e-08_wp,            &
 1.20226e-08_wp,  1.23027e-08_wp,  1.25893e-08_wp,  1.28825e-08_wp,            &
 1.31826e-08_wp,  1.34896e-08_wp,  1.38038e-08_wp,  1.41254e-08_wp,            &
 1.44544e-08_wp,  1.47911e-08_wp,  1.51356e-08_wp,  1.54882e-08_wp /
data (mie_radii(mie_radii_index), mie_radii_index =  121,  160 )   /           &
 1.58489e-08_wp,  1.62181e-08_wp,  1.65959e-08_wp,  1.69824e-08_wp,            &
 1.73780e-08_wp,  1.77828e-08_wp,  1.81970e-08_wp,  1.86209e-08_wp,            &
 1.90546e-08_wp,  1.94984e-08_wp,  1.99526e-08_wp,  2.04174e-08_wp,            &
 2.08930e-08_wp,  2.13796e-08_wp,  2.18776e-08_wp,  2.23872e-08_wp,            &
 2.29087e-08_wp,  2.34423e-08_wp,  2.39883e-08_wp,  2.45471e-08_wp,            &
 2.51189e-08_wp,  2.57040e-08_wp,  2.63027e-08_wp,  2.69153e-08_wp,            &
 2.75423e-08_wp,  2.81838e-08_wp,  2.88403e-08_wp,  2.95121e-08_wp,            &
 3.01995e-08_wp,  3.09030e-08_wp,  3.16228e-08_wp,  3.23594e-08_wp,            &
 3.31131e-08_wp,  3.38844e-08_wp,  3.46737e-08_wp,  3.54813e-08_wp,            &
 3.63078e-08_wp,  3.71535e-08_wp,  3.80189e-08_wp,  3.89045e-08_wp /
data (mie_radii(mie_radii_index), mie_radii_index =  161,  200 )   /           &
 3.98107e-08_wp,  4.07380e-08_wp,  4.16869e-08_wp,  4.26580e-08_wp,            &
 4.36516e-08_wp,  4.46684e-08_wp,  4.57088e-08_wp,  4.67735e-08_wp,            &
 4.78630e-08_wp,  4.89779e-08_wp,  5.01187e-08_wp,  5.12861e-08_wp,            &
 5.24807e-08_wp,  5.37032e-08_wp,  5.49541e-08_wp,  5.62341e-08_wp,            &
 5.75440e-08_wp,  5.88844e-08_wp,  6.02560e-08_wp,  6.16595e-08_wp,            &
 6.30957e-08_wp,  6.45654e-08_wp,  6.60693e-08_wp,  6.76083e-08_wp,            &
 6.91831e-08_wp,  7.07946e-08_wp,  7.24436e-08_wp,  7.41310e-08_wp,            &
 7.58578e-08_wp,  7.76247e-08_wp,  7.94328e-08_wp,  8.12831e-08_wp,            &
 8.31764e-08_wp,  8.51138e-08_wp,  8.70964e-08_wp,  8.91251e-08_wp,            &
 9.12011e-08_wp,  9.33254e-08_wp,  9.54993e-08_wp,  9.77237e-08_wp /
data (mie_radii(mie_radii_index), mie_radii_index =  201,  240 )   /           &
 1.00000e-07_wp,  1.02329e-07_wp,  1.04713e-07_wp,  1.07152e-07_wp,            &
 1.09648e-07_wp,  1.12202e-07_wp,  1.14815e-07_wp,  1.17490e-07_wp,            &
 1.20226e-07_wp,  1.23027e-07_wp,  1.25893e-07_wp,  1.28825e-07_wp,            &
 1.31826e-07_wp,  1.34896e-07_wp,  1.38038e-07_wp,  1.41254e-07_wp,            &
 1.44544e-07_wp,  1.47911e-07_wp,  1.51356e-07_wp,  1.54882e-07_wp,            &
 1.58489e-07_wp,  1.62181e-07_wp,  1.65959e-07_wp,  1.69824e-07_wp,            &
 1.73780e-07_wp,  1.77828e-07_wp,  1.81970e-07_wp,  1.86209e-07_wp,            &
 1.90546e-07_wp,  1.94984e-07_wp,  1.99526e-07_wp,  2.04174e-07_wp,            &
 2.08930e-07_wp,  2.13796e-07_wp,  2.18776e-07_wp,  2.23872e-07_wp,            &
 2.29087e-07_wp,  2.34423e-07_wp,  2.39883e-07_wp,  2.45471e-07_wp /
data (mie_radii(mie_radii_index), mie_radii_index =  241,  280 )   /           &
 2.51189e-07_wp,  2.57040e-07_wp,  2.63027e-07_wp,  2.69153e-07_wp,            &
 2.75423e-07_wp,  2.81838e-07_wp,  2.88403e-07_wp,  2.95121e-07_wp,            &
 3.01995e-07_wp,  3.09030e-07_wp,  3.16228e-07_wp,  3.23594e-07_wp,            &
 3.31131e-07_wp,  3.38844e-07_wp,  3.46737e-07_wp,  3.54813e-07_wp,            &
 3.63078e-07_wp,  3.71535e-07_wp,  3.80189e-07_wp,  3.89045e-07_wp,            &
 3.98107e-07_wp,  4.07380e-07_wp,  4.16869e-07_wp,  4.26580e-07_wp,            &
 4.36516e-07_wp,  4.46684e-07_wp,  4.57088e-07_wp,  4.67735e-07_wp,            &
 4.78630e-07_wp,  4.89779e-07_wp,  5.01187e-07_wp,  5.12861e-07_wp,            &
 5.24807e-07_wp,  5.37032e-07_wp,  5.49541e-07_wp,  5.62341e-07_wp,            &
 5.75440e-07_wp,  5.88844e-07_wp,  6.02560e-07_wp,  6.16595e-07_wp /
data (mie_radii(mie_radii_index), mie_radii_index =  281,  320 )   /           &
 6.30957e-07_wp,  6.45654e-07_wp,  6.60693e-07_wp,  6.76083e-07_wp,            &
 6.91831e-07_wp,  7.07946e-07_wp,  7.24436e-07_wp,  7.41310e-07_wp,            &
 7.58578e-07_wp,  7.76247e-07_wp,  7.94328e-07_wp,  8.12831e-07_wp,            &
 8.31764e-07_wp,  8.51138e-07_wp,  8.70964e-07_wp,  8.91251e-07_wp,            &
 9.12011e-07_wp,  9.33254e-07_wp,  9.54993e-07_wp,  9.77237e-07_wp,            &
 1.00000e-06_wp,  1.02329e-06_wp,  1.04713e-06_wp,  1.07152e-06_wp,            &
 1.09648e-06_wp,  1.12202e-06_wp,  1.14815e-06_wp,  1.17490e-06_wp,            &
 1.20226e-06_wp,  1.23027e-06_wp,  1.25893e-06_wp,  1.28825e-06_wp,            &
 1.31826e-06_wp,  1.34896e-06_wp,  1.38038e-06_wp,  1.41254e-06_wp,            &
 1.44544e-06_wp,  1.47911e-06_wp,  1.51356e-06_wp,  1.54882e-06_wp /
data (mie_radii(mie_radii_index), mie_radii_index =  321,  360 )   /           &
 1.58489e-06_wp,  1.62181e-06_wp,  1.65959e-06_wp,  1.69824e-06_wp,            &
 1.73780e-06_wp,  1.77828e-06_wp,  1.81970e-06_wp,  1.86209e-06_wp,            &
 1.90546e-06_wp,  1.94984e-06_wp,  1.99526e-06_wp,  2.04174e-06_wp,            &
 2.08930e-06_wp,  2.13796e-06_wp,  2.18776e-06_wp,  2.23872e-06_wp,            &
 2.29087e-06_wp,  2.34423e-06_wp,  2.39883e-06_wp,  2.45471e-06_wp,            &
 2.51189e-06_wp,  2.57040e-06_wp,  2.63027e-06_wp,  2.69153e-06_wp,            &
 2.75423e-06_wp,  2.81838e-06_wp,  2.88403e-06_wp,  2.95121e-06_wp,            &
 3.01995e-06_wp,  3.09030e-06_wp,  3.16228e-06_wp,  3.23594e-06_wp,            &
 3.31131e-06_wp,  3.38844e-06_wp,  3.46737e-06_wp,  3.54813e-06_wp,            &
 3.63078e-06_wp,  3.71535e-06_wp,  3.80189e-06_wp,  3.89045e-06_wp /
data (mie_radii(mie_radii_index), mie_radii_index =  361,  400 )   /           &
 3.98107e-06_wp,  4.07380e-06_wp,  4.16869e-06_wp,  4.26580e-06_wp,            &
 4.36516e-06_wp,  4.46684e-06_wp,  4.57088e-06_wp,  4.67735e-06_wp,            &
 4.78630e-06_wp,  4.89779e-06_wp,  5.01187e-06_wp,  5.12861e-06_wp,            &
 5.24807e-06_wp,  5.37032e-06_wp,  5.49541e-06_wp,  5.62341e-06_wp,            &
 5.75440e-06_wp,  5.88844e-06_wp,  6.02560e-06_wp,  6.16595e-06_wp,            &
 6.30957e-06_wp,  6.45654e-06_wp,  6.60693e-06_wp,  6.76083e-06_wp,            &
 6.91831e-06_wp,  7.07946e-06_wp,  7.24436e-06_wp,  7.41310e-06_wp,            &
 7.58578e-06_wp,  7.76247e-06_wp,  7.94328e-06_wp,  8.12831e-06_wp,            &
 8.31764e-06_wp,  8.51138e-06_wp,  8.70964e-06_wp,  8.91251e-06_wp,            &
 9.12011e-06_wp,  9.33254e-06_wp,  9.54993e-06_wp,  9.77237e-06_wp /
data (mie_radii(mie_radii_index), mie_radii_index =  401,  401 )   /           &
 1.00000e-05_wp /

! lookup table values of extinction efficiency, Qext
integer (wi) :: mie_qext_index
real    (wp) :: mie_qext(1:n_mie)

data (mie_qext(mie_qext_index), mie_qext_index =    1,   40 )      /           &
 1.29779e-04_wp,  1.32803e-04_wp,  1.35897e-04_wp,  1.39064e-04_wp,            &
 1.42304e-04_wp,  1.45620e-04_wp,  1.49013e-04_wp,  1.52485e-04_wp,            &
 1.56038e-04_wp,  1.59674e-04_wp,  1.63395e-04_wp,  1.67203e-04_wp,            &
 1.71099e-04_wp,  1.75087e-04_wp,  1.79167e-04_wp,  1.83343e-04_wp,            &
 1.87616e-04_wp,  1.91989e-04_wp,  1.96465e-04_wp,  2.01044e-04_wp,            &
 2.05731e-04_wp,  2.10527e-04_wp,  2.15435e-04_wp,  2.20457e-04_wp,            &
 2.25597e-04_wp,  2.30858e-04_wp,  2.36241e-04_wp,  2.41750e-04_wp,            &
 2.47388e-04_wp,  2.53158e-04_wp,  2.59063e-04_wp,  2.65106e-04_wp,            &
 2.71291e-04_wp,  2.77621e-04_wp,  2.84099e-04_wp,  2.90728e-04_wp,            &
 2.97514e-04_wp,  3.04458e-04_wp,  3.11566e-04_wp,  3.18840e-04_wp /
data (mie_qext(mie_qext_index), mie_qext_index =   41,   80 )      /           &
 3.26286e-04_wp,  3.33906e-04_wp,  3.41706e-04_wp,  3.49689e-04_wp,            &
 3.57861e-04_wp,  3.66225e-04_wp,  3.74786e-04_wp,  3.83550e-04_wp,            &
 3.92520e-04_wp,  4.01703e-04_wp,  4.11104e-04_wp,  4.20727e-04_wp,            &
 4.30578e-04_wp,  4.40663e-04_wp,  4.50989e-04_wp,  4.61560e-04_wp,            &
 4.72383e-04_wp,  4.83465e-04_wp,  4.94812e-04_wp,  5.06431e-04_wp,            &
 5.18330e-04_wp,  5.30510e-04_wp,  5.42992e-04_wp,  5.55768e-04_wp,            &
 5.68872e-04_wp,  5.82277e-04_wp,  5.96005e-04_wp,  6.10061e-04_wp,            &
 6.24465e-04_wp,  6.39244e-04_wp,  6.54385e-04_wp,  6.69871e-04_wp,            &
 6.85806e-04_wp,  7.02086e-04_wp,  7.18775e-04_wp,  7.35858e-04_wp,            &
 7.53400e-04_wp,  7.71434e-04_wp,  7.89859e-04_wp,  8.08798e-04_wp /
data (mie_qext(mie_qext_index), mie_qext_index =   81,  120 )      /           &
 8.28172e-04_wp,  8.48133e-04_wp,  8.68577e-04_wp,  8.89551e-04_wp,            &
 9.11162e-04_wp,  9.33218e-04_wp,  9.55979e-04_wp,  9.79249e-04_wp,            &
 1.00324e-03_wp,  1.02791e-03_wp,  1.05318e-03_wp,  1.07925e-03_wp,            &
 1.10612e-03_wp,  1.13368e-03_wp,  1.16212e-03_wp,  1.19131e-03_wp,            &
 1.22142e-03_wp,  1.25246e-03_wp,  1.28446e-03_wp,  1.31747e-03_wp,            &
 1.35159e-03_wp,  1.38675e-03_wp,  1.42318e-03_wp,  1.46075e-03_wp,            &
 1.49960e-03_wp,  1.53981e-03_wp,  1.58151e-03_wp,  1.62469e-03_wp,            &
 1.66949e-03_wp,  1.71593e-03_wp,  1.76417e-03_wp,  1.81438e-03_wp,            &
 1.86658e-03_wp,  1.92091e-03_wp,  1.97753e-03_wp,  2.03662e-03_wp,            &
 2.09830e-03_wp,  2.16278e-03_wp,  2.23016e-03_wp,  2.30081e-03_wp /
data (mie_qext(mie_qext_index), mie_qext_index =  121,  160 )      /           &
 2.37480e-03_wp,  2.45244e-03_wp,  2.53404e-03_wp,  2.61982e-03_wp,            &
 2.71013e-03_wp,  2.80533e-03_wp,  2.90578e-03_wp,  3.01190e-03_wp,            &
 3.12412e-03_wp,  3.24294e-03_wp,  3.36892e-03_wp,  3.50259e-03_wp,            &
 3.64464e-03_wp,  3.79572e-03_wp,  3.95663e-03_wp,  4.12814e-03_wp,            &
 4.31120e-03_wp,  4.50679e-03_wp,  4.71600e-03_wp,  4.93998e-03_wp,            &
 5.18005e-03_wp,  5.43761e-03_wp,  5.71419e-03_wp,  6.01150e-03_wp,            &
 6.33138e-03_wp,  6.67582e-03_wp,  7.04704e-03_wp,  7.44746e-03_wp,            &
 7.87970e-03_wp,  8.34665e-03_wp,  8.85145e-03_wp,  9.39756e-03_wp,            &
 9.98872e-03_wp,  1.06291e-02_wp,  1.13230e-02_wp,  1.20756e-02_wp,            &
 1.28921e-02_wp,  1.37784e-02_wp,  1.47408e-02_wp,  1.57863e-02_wp /
data (mie_qext(mie_qext_index), mie_qext_index =  161,  200 )      /           &
 1.69225e-02_wp,  1.81577e-02_wp,  1.95008e-02_wp,  2.09616e-02_wp,            &
 2.25509e-02_wp,  2.42800e-02_wp,  2.61617e-02_wp,  2.82093e-02_wp,            &
 3.04377e-02_wp,  3.28627e-02_wp,  3.55015e-02_wp,  3.83723e-02_wp,            &
 4.14952e-02_wp,  4.48911e-02_wp,  4.85829e-02_wp,  5.25947e-02_wp,            &
 5.69522e-02_wp,  6.16825e-02_wp,  6.68142e-02_wp,  7.23775e-02_wp,            &
 7.84036e-02_wp,  8.49254e-02_wp,  9.19764e-02_wp,  9.95915e-02_wp,            &
 1.07806e-01_wp,  1.16656e-01_wp,  1.26178e-01_wp,  1.36408e-01_wp,            &
 1.47384e-01_wp,  1.59142e-01_wp,  1.71719e-01_wp,  1.85154e-01_wp,            &
 1.99486e-01_wp,  2.14756e-01_wp,  2.31011e-01_wp,  2.48300e-01_wp,            &
 2.66684e-01_wp,  2.86232e-01_wp,  3.07029e-01_wp,  3.29174e-01_wp /
data (mie_qext(mie_qext_index), mie_qext_index =  201,  240 )      /           &
 3.52785e-01_wp,  3.78002e-01_wp,  4.04976e-01_wp,  4.33872e-01_wp,            &
 4.64851e-01_wp,  4.98055e-01_wp,  5.33595e-01_wp,  5.71528e-01_wp,            &
 6.11850e-01_wp,  6.54499e-01_wp,  6.99358e-01_wp,  7.46280e-01_wp,            &
 7.95109e-01_wp,  8.45700e-01_wp,  8.97949e-01_wp,  9.51811e-01_wp,            &
 1.00731e+00_wp,  1.06458e+00_wp,  1.12382e+00_wp,  1.18534e+00_wp,            &
 1.24945e+00_wp,  1.31642e+00_wp,  1.38635e+00_wp,  1.45908e+00_wp,            &
 1.53422e+00_wp,  1.61116e+00_wp,  1.68921e+00_wp,  1.76772e+00_wp,            &
 1.84612e+00_wp,  1.92410e+00_wp,  2.00161e+00_wp,  2.07899e+00_wp,            &
 2.15678e+00_wp,  2.23545e+00_wp,  2.31507e+00_wp,  2.39533e+00_wp,            &
 2.47567e+00_wp,  2.55557e+00_wp,  2.63463e+00_wp,  2.71273e+00_wp /
data (mie_qext(mie_qext_index), mie_qext_index =  241,  280 )      /           &
 2.79013e+00_wp,  2.86739e+00_wp,  2.94464e+00_wp,  3.02099e+00_wp,            &
 3.09481e+00_wp,  3.16437e+00_wp,  3.22827e+00_wp,  3.28600e+00_wp,            &
 3.33825e+00_wp,  3.38645e+00_wp,  3.43118e+00_wp,  3.47232e+00_wp,            &
 3.50934e+00_wp,  3.54154e+00_wp,  3.56860e+00_wp,  3.59076e+00_wp,            &
 3.60709e+00_wp,  3.61610e+00_wp,  3.61673e+00_wp,  3.60874e+00_wp,            &
 3.59382e+00_wp,  3.57403e+00_wp,  3.55035e+00_wp,  3.52322e+00_wp,            &
 3.49147e+00_wp,  3.45459e+00_wp,  3.41174e+00_wp,  3.36328e+00_wp,            &
 3.31025e+00_wp,  3.25394e+00_wp,  3.19684e+00_wp,  3.14051e+00_wp,            &
 3.08499e+00_wp,  3.02790e+00_wp,  2.96902e+00_wp,  2.91032e+00_wp,            &
 2.85339e+00_wp,  2.79934e+00_wp,  2.74929e+00_wp,  2.70282e+00_wp /
data (mie_qext(mie_qext_index), mie_qext_index =  281,  320 )      /           &
 2.65797e+00_wp,  2.61498e+00_wp,  2.57555e+00_wp,  2.54097e+00_wp,            &
 2.51076e+00_wp,  2.48357e+00_wp,  2.45924e+00_wp,  2.43730e+00_wp,            &
 2.41936e+00_wp,  2.40551e+00_wp,  2.39423e+00_wp,  2.38521e+00_wp,            &
 2.37774e+00_wp,  2.37399e+00_wp,  2.37218e+00_wp,  2.37228e+00_wp,            &
 2.37248e+00_wp,  2.37423e+00_wp,  2.37708e+00_wp,  2.37982e+00_wp,            &
 2.38106e+00_wp,  2.38195e+00_wp,  2.38256e+00_wp,  2.38095e+00_wp,            &
 2.37764e+00_wp,  2.37391e+00_wp,  2.36900e+00_wp,  2.36286e+00_wp,            &
 2.35720e+00_wp,  2.35146e+00_wp,  2.34595e+00_wp,  2.34171e+00_wp,            &
 2.33773e+00_wp,  2.33447e+00_wp,  2.33182e+00_wp,  2.32910e+00_wp,            &
 2.32657e+00_wp,  2.32329e+00_wp,  2.31966e+00_wp,  2.31535e+00_wp /
data (mie_qext(mie_qext_index), mie_qext_index =  321,  360 )      /           &
 2.30987e+00_wp,  2.30392e+00_wp,  2.29673e+00_wp,  2.28892e+00_wp,            &
 2.28050e+00_wp,  2.27199e+00_wp,  2.26392e+00_wp,  2.25675e+00_wp,            &
 2.25081e+00_wp,  2.24625e+00_wp,  2.24272e+00_wp,  2.24002e+00_wp,            &
 2.23744e+00_wp,  2.23502e+00_wp,  2.23236e+00_wp,  2.22977e+00_wp,            &
 2.22713e+00_wp,  2.22443e+00_wp,  2.22163e+00_wp,  2.21855e+00_wp,            &
 2.21530e+00_wp,  2.21182e+00_wp,  2.20834e+00_wp,  2.20483e+00_wp,            &
 2.20106e+00_wp,  2.19721e+00_wp,  2.19348e+00_wp,  2.18999e+00_wp,            &
 2.18679e+00_wp,  2.18393e+00_wp,  2.18142e+00_wp,  2.17908e+00_wp,            &
 2.17659e+00_wp,  2.17365e+00_wp,  2.17034e+00_wp,  2.16721e+00_wp,            &
 2.16475e+00_wp,  2.16289e+00_wp,  2.16115e+00_wp,  2.15912e+00_wp /
data (mie_qext(mie_qext_index), mie_qext_index =  361,  400 )      /           &
 2.15662e+00_wp,  2.15365e+00_wp,  2.15044e+00_wp,  2.14759e+00_wp,            &
 2.14536e+00_wp,  2.14344e+00_wp,  2.14152e+00_wp,  2.13960e+00_wp,            &
 2.13761e+00_wp,  2.13544e+00_wp,  2.13321e+00_wp,  2.13099e+00_wp,            &
 2.12888e+00_wp,  2.12698e+00_wp,  2.12516e+00_wp,  2.12322e+00_wp,            &
 2.12118e+00_wp,  2.11916e+00_wp,  2.11744e+00_wp,  2.11600e+00_wp,            &
 2.11431e+00_wp,  2.11225e+00_wp,  2.11034e+00_wp,  2.10875e+00_wp,            &
 2.10713e+00_wp,  2.10551e+00_wp,  2.10397e+00_wp,  2.10249e+00_wp,            &
 2.10076e+00_wp,  2.09909e+00_wp,  2.09768e+00_wp,  2.09627e+00_wp,            &
 2.09473e+00_wp,  2.09328e+00_wp,  2.09191e+00_wp,  2.09049e+00_wp,            &
 2.08908e+00_wp,  2.08768e+00_wp,  2.08638e+00_wp,  2.08517e+00_wp /
data (mie_qext(mie_qext_index), mie_qext_index =  401,  401 )      /           &
 2.08380e+00_wp /

contains

  !=============================================================================
  !
  ! vera_blumel
  !
  ! This subroutine computes the Blumel approximation to the Mie scattering
  ! Qext for very small particles.
  !
  ! The size parameter, xx, is defined as
  !
  !         xx = 2.pi.r. / lambda
  !
  ! where r is the particle radius and lambda is the wavelength of the
  ! light. The Blumel approximation is valid for xx < 1, which for
  ! visibile light (lambda=550nm), corresponds to r < 0.08753522 micrometres.
  !
  !=============================================================================

subroutine vera_blumel( radii, qext )

use vera_global_mod, only: vera_mie       ,                                    &
                           vera_constants

! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!===========================================================================
! arguments for vera_blumel
!
!===========================================================================

! aerosol particle radii
real    (wp), intent(in)  :: radii(:)

! computed extinction efficiency, Qext
real    (wp), intent(out) :: qext(1:size(radii))

!===========================================================================
! local variables for vera_blumel
!
!===========================================================================

! size parameter  = 2.pi.r / lambda
real    (wp) :: xx(1:size(radii))

! square of the refractive index
complex (wp) :: m_squared

! parameter n
complex (wp) :: n(1:size(radii))

! parameter d
complex (wp) :: d(1:size(radii))

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_BLUMEL'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! start of executable code for vera_blumel
!
!===========================================================================

! compute the size parameter xx = 2.pi.r / lambda
xx = vera_constants%two_pi * radii / vera_mie%wavelength

! compute the square of the refractive index
m_squared = vera_mie%refractive_index ** vera_constants%two

! compute the parameter n
n = xx**(vera_constants%three) * ( vera_constants%one - m_squared )

! compute the parameter d
d = n + ( vera_constants%complex_i * vera_constants%three_halves *             &
          (vera_constants%two + m_squared) )

! compute the Blumel approximation to Qext
qext = real(( vera_constants%six / (xx ** vera_constants%two) ) *              &
       real( n * conjg(d) )                                /                   &
       ( d * conjg(d) ))

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_blumel

!=============================================================================
!
! vera_mie_lookup
!
! This subroutine uses log-linear interpoloation to compute the extinction
! efficiency Qext as a function of the particle radius, using a lookup table.
!
!=============================================================================

subroutine vera_mie_lookup( radii, qext )

use vera_global_mod, only: vera_mie       ,                                    &
                           vera_constants ,                                    &
                           vera_visbty

! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!===========================================================================
! arguments for vera_mie_lookup
!
!===========================================================================

! aerosol particle radii
real    (wp), intent(in)  :: radii(:)

! computed extinction efficiency, Qext
real    (wp), intent(out) :: qext(1:size(radii))

!===========================================================================
! local variables for vera_mie_lookup
!
!===========================================================================

! number of input particle radii
integer (wi)              :: n_radii

! counter to loop over the particle radii
integer (wi)              :: loop_radii

! index for upper bound of particle size in the lookup table
integer (wi)              :: upper_radius

! range of Qext between adjacent entries in the lookup table
real    (wp)              :: delta_qext

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_MIE_LOOKUP'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! start of executable code for vera_mie_lookup
!
!===========================================================================

! check to see if geometric scattering is preferred rather than Mie
! scattering, and check if the old VISBTY scattering scheme is required
if (   (vera_mie%geometric_scattering == vera_mie%geometric_scattering_on)     &
   .or. (vera_mie%visbty_scattering == vera_mie%visbty_scattering_on) ) then

  ! use geometric scattering and set Qext to 2, the large particle limit
  qext = vera_mie%qext_upper_limit

  ! check to see if the VISBTY scattering scheme is required
  if ( vera_mie%visbty_scattering == vera_mie%visbty_scattering_on ) then

    ! compute the effective extinction efficiency
    ! Qeff = eta . Qext(large particle limit)
    qext = vera_visbty%eta * vera_mie%qext_upper_limit

  end if

else
  ! else use Mie scattering to compute Qext

  ! how many input particle radii are there?
  n_radii = size(radii)

  ! loop over the particle radii
  do loop_radii = vera_constants%one_i, n_radii

    ! check to see if the particle radius is outside the scope of the
    ! lookup table, i.e. particle is too small or too large
    if ( ( radii(loop_radii) > maxval(mie_radii) ) .or.                        &
         ( radii(loop_radii) < minval(mie_radii) )      ) then

      ! large particle, so use large particle limit for Qext
      if ( radii(loop_radii) > maxval(mie_radii) ) then
        qext(loop_radii) = vera_mie%qext_upper_limit
      end if

      ! small particle, so compute Qext using Blumel approximation
      if ( radii(loop_radii) < minval(mie_radii) ) then

        if ( vera_mie%switch_blumel == vera_mie%switch_blumel_on ) then
          ! use Blumel's approximation to Qext
          call vera_blumel( [radii(loop_radii)], qext(loop_radii) )
        else
          ! use the lower limit for Qext
          qext(loop_radii) = vera_mie%qext_lower_limit
        end if
      end if

      ! else particle is in the range of the lookup table
    else

      ! find the index of the upper particle size in the lookup table
      upper_radius = minloc( mie_radii, dim = 1, mask =                        &
             ( ( mie_radii - radii(loop_radii) ) >= vera_constants%zero ) )

      ! compute the range of Qext spanning the input particle radius
      delta_qext   = mie_qext(upper_radius                       ) -           &
                     mie_qext(upper_radius - vera_constants%one_i)

      ! compute the required Qext as a log-linear interpolation
      ! i.e. log in the radius, linear in Qext
      qext(loop_radii) = ( (( log(radii(loop_radii)                     /      &
                      mie_radii(upper_radius - vera_constants%one_i)) ) /      &
                      ( log(mie_radii(upper_radius)                     /      &
                      mie_radii(upper_radius - vera_constants%one_i)) ) )      &
                      * delta_qext )                                    +      &
                      mie_qext(upper_radius - vera_constants%one_i)

      ! end of checking the particle size
    end if

    ! end of looping over the input particle radii
  end do

  ! end of checking for geometric scattering
end if

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_mie_lookup

end module vera_mie_mod
