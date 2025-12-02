! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! A module containing constants/parameters for Sulphur Cycle Chemistry
!
module c_sulchm_mod

use um_types, only: real_umphys

implicit none

! Description:
!   This module contains constants for the
!   Sulphur Cycle Chemistry
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Aerosols
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 standards.
!----------------------------------------------------------------------

! Timescale for dissolved SO4 to evaporate (5 mins)
real(kind=real_umphys), parameter :: evaptau = 300.0 ! (s)

! Timescale for accumulation mode particles to nucleate once they enter a cloud
real(kind=real_umphys), parameter :: nuctau = 30.0 ! (s)

! Sulphur/2*Nitrogen
real(kind=real_umphys), parameter :: relm_s_2n = 3.206 / 2.80

! Reaction rate for DMS+OH  (cc/mcl/s)
real(kind=real_umphys), parameter :: k_dms_oh = 9.1e-12

! Rate coeff for CH3SO2+O3 -> CH3SO3+O2 (cc/mcl/s)
real(kind=real_umphys), parameter :: k4_ch3so2_o3 = 1.0e-14

! Rate coeff for CH3SO3+HO2 -> MSA+O2
real(kind=real_umphys), parameter :: k5_ch3so3_ho2 = 4.0e-11

! High pressure reaction rate limit (cc/mcl/s) from STOCHEM model
real(kind=real_umphys), parameter :: k_so2oh_hi = 2.0e-12

! Branching ratio
real(kind=real_umphys), parameter :: brat_so2 = 0.9
                                             ! SO2 in DMS oxidn
real(kind=real_umphys), parameter :: brat_msa = 1.0 - brat_so2
                                             ! MSA in DMS oxidn

! Relative atomic mass
real(kind=real_umphys), parameter :: relm_s_h2o2 = 3.206 / 3.40
                                              ! Sulphur/RMM_H2O2


! Power of temp dependence of K_SO2OH_LO
real(kind=real_umphys), parameter :: parh = 3.3

! Parameters for calcn of low pressure reaction rate of SO2 with OH
real(kind=real_umphys), parameter :: k1 = 4.0e-31   ! (cc/mcl)2/s
real(kind=real_umphys), parameter :: t1 = 300.0     ! K

! Parameters for calcn of HO2 + HO2 reaction rate
real(kind=real_umphys), parameter :: k2 = 2.2e-13
real(kind=real_umphys), parameter :: t2 = 600.0     ! K
real(kind=real_umphys), parameter :: k3 = 1.9e-33
real(kind=real_umphys), parameter :: t3 = 890.0     ! K
real(kind=real_umphys), parameter :: k4 = 1.4e-21
real(kind=real_umphys), parameter :: t4 = 2200.0    ! K


! Parameters for interpolation between LO and HI reaction rate limits
real(kind=real_umphys), parameter :: fc = 0.45
real(kind=real_umphys), parameter :: fac1 = 1.1904 ! ( 0.75-1.27*log10(FC)

! Air parcel lifetime in cloud (3 hours)
real(kind=real_umphys), parameter :: cloudtau = 1.08e4 ! s

! Chem lifetime in cloud before oxidn (15 mins)
real(kind=real_umphys), parameter :: chemtau = 9.0e2 ! s

! Min mmr of O3 required for oxidn (kg/kg, equiv. 10ppbv)
real(kind=real_umphys), parameter :: o3_min = 1.6e-8

! Threshold for cloud liquid water (kg/kg)
real(kind=real_umphys), parameter :: thold = 1.0e-8

! Threshold concn of accu mode particles below which PSI=1 (m-3)
real(kind=real_umphys), parameter :: num_star = 1.0e6

! Geometric standard devn of the Aitken mode distn
real(kind=real_umphys), parameter :: sigma_ait = 1.30

! Mean radius of particles (m)
real(kind=real_umphys), parameter :: rad_ait = 6.5e-9  ! Aitken mode
real(kind=real_umphys), parameter :: rad_acc = 95.0e-9 ! Acccumulation mode

! Mole fraction of S in particle
real(kind=real_umphys), parameter :: chi = 32.0 / 132.0

! Standard devn of particle size distn for accum mode
real(kind=real_umphys), parameter :: sigma = 1.4

end module c_sulchm_mod

