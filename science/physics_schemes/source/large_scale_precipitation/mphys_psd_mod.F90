! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module mphys_psd_mod

! Description:
! Holds particle-size distribution constants required by the large-scale
! precipitation scheme
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

use um_types, only: real_umphys

implicit none


!---------------------------------------------------------------------
! CALCULATION OF FALL SPEEDS
!---------------------------------------------------------------------

      ! Do we wish to calculate the ice fall velocities?
      ! true if calculate speeds, false if specify speeds

logical, parameter :: l_calcfall = .true.

!---------------------------------------------------------------------
! RAIN PARAMETERS
!---------------------------------------------------------------------


      ! Drop size distribution for rain: N(D) =  N0 D^m exp(-lambda D)
      ! where N0 = X1R lambda^X2R  and m = X4R

!     real, parameter :: x1r is set in the run_precip namelist
!     real, parameter :: x2r is set in the run_precip namelist
!     real, parameter :: x4r is set in the run_precip namelist

      ! Fall speed diameter relationship for rain: vt(D) = CR D^DR

real(kind=real_umphys), parameter :: cr = 386.8
real(kind=real_umphys), parameter :: dr = 0.67

 ! Abel & Shipway (2007) fall speed diameter relationships for rain:
 !     vt(D) = C1R D^D1R exp(-h1r D) + C2R D^D2R exp(-h2r D)

real(kind=real_umphys), parameter :: c1r = 4845.1
real(kind=real_umphys), parameter :: d1r = 1.0
real(kind=real_umphys), parameter :: h1r = 195.0
real(kind=real_umphys), parameter :: c2r = -446.009
real(kind=real_umphys), parameter :: d2r = 0.782127
real(kind=real_umphys), parameter :: h2r = 4085.35

!---------------------------------------------------------------------
! ICE AGGREGATE PARAMETERS
!---------------------------------------------------------------------

      ! Particle size distribution for ice aggregates:
      ! N(D) = N0 D^m exp(-lambda D)
      ! where N0 = X1I TCG lambda^X2I,
      ! m = X4I and TCG = exp(- X3I T[deg C])

!     real, parameter :: x1i is set in the run_precip namelist

real(kind=real_umphys), parameter :: x2i = 0.0
real(kind=real_umphys), parameter :: x3i = 0.1222
real(kind=real_umphys), parameter :: x4i = 0.0

 ! Mass diameter relationship for ice:  m(D) = AI D^BI
 ! These are set in the run_precip namelist.
 ! Recommended values for the generic particle size distribution are
 ! from Brown and Francis and are
 ! AI = AIC = 1.85E-2, BI = BIC = 1.9.
 ! If l_calcfall is changed from .true. then the generic psd values
 ! should be set below to ci = cic = 8.203 and di = dic = 0.2888

 ! real, parameter :: ai   )
 ! real, parameter :: bi   ) Values are set in the run_precip namelist
 ! real, parameter :: aic  )
 ! real, parameter :: bic  ) Values are copied in lspcon module

 ! The area diameter relationships are only used if
 ! L_CALCFALL = .true.
 ! Area diameter relationship for ice:  Area(D) = RI D^SI

real(kind=real_umphys), parameter :: ri = 0.131
real(kind=real_umphys), parameter :: si = 1.88


 ! The Best/Reynolds relationships are only used if
 ! L_CALCFALL = .true.
 ! Relationship between Best number and Reynolds number:
 ! Re(D)  = LSP_EI(C) Be^LSP_FI(C)

 ! real, parameter :: lsp_ei  )
 ! real, parameter :: lsp_fi  ) Set in mphys_constants_mod module

 ! The fall speeds of ice particles are only used if
 ! L_CALCFALL = .false.
 ! Fall speed diameter relationships for ice:
 ! vt(D) = CI D^DI

real(kind=real_umphys), parameter :: ci0 = 14.3
real(kind=real_umphys), parameter :: di0 = 0.416

!---------------------------------------------------------------------
! ICE CRYSTAL PARAMETERS
!---------------------------------------------------------------------

      ! Particle size distribution for ice crystals:

      ! N(D) = N0 D^m exp(-lambda D)
      ! where N0 = X1I TCG lambda^X2I,
      ! m = X4I and TCG = exp(- X3I T[deg C])

!    real, parameter :: x1ic is set in the run_precip namelist

real(kind=real_umphys), parameter :: x2ic = 0.0

real(kind=real_umphys), parameter :: x3ic = 0.1222
real(kind=real_umphys), parameter :: x4ic = 0.0

 ! The area diameter relationships are only used if
 ! L_CALCFALL = .true.
 ! Area diameter relationship for ice:  Area(D) = RI D^SI


real(kind=real_umphys), parameter :: ric = 0.131
real(kind=real_umphys), parameter :: sic = 1.88

 ! The Best/Reynolds relationships are only used if
 ! L_CALCFALL = .true.
 ! Relationship between Best number and Reynolds number:
 ! Re(D)  = LSP_EI(C) Be^LSP_FI(C)

 ! real, parameter :: lsp_eic  )
 ! real, parameter :: lsp_fic  ) Set in mphys_constants_mod module

 ! The fall speeds of ice particles are only used if
 ! L_CALCFALL = .false.
 ! Fall speed diameter relationships for ice:
 ! vt(D) = CI D^DI

real(kind=real_umphys), parameter :: cic0 = 74.5
real(kind=real_umphys), parameter :: dic0 = 0.640

!---------------------------------------------------------------------
! GRAUPEL PARAMETERS
!---------------------------------------------------------------------

! Drop size distribution for graupel: N(D) =  N0 D^m exp(-lambda D)
! where N0 = X1G lambda^X2G and m = X4G.
! The following values are saved in this module and set in
! lspcon and hence are not parameters
real(kind=real_umphys) :: x1g
real(kind=real_umphys) :: x2g
real(kind=real_umphys) :: x4g

! Original values from LEM
real(kind=real_umphys), parameter :: x1gl = 5.0e25
real(kind=real_umphys), parameter :: x2gl = -4.0
real(kind=real_umphys), parameter :: x4gl = 2.5

! Field et al new PSD values
real(kind=real_umphys), parameter :: x1gf = 7.9e9
real(kind=real_umphys), parameter :: x2gf = -2.58
real(kind=real_umphys), parameter :: x4gf = 0.0

 ! Mass diameter relationship for graupel:  m(D) = AG D^BG

real(kind=real_umphys), parameter :: ag  = 261.8  ! 500kg/m3
real(kind=real_umphys), parameter :: bg  = 3.0

 ! Fall speed diameter relationship for graupel: vt(D) = CG D^DG

real(kind=real_umphys), parameter :: cg  = 253.0
real(kind=real_umphys), parameter :: dg  = 0.734


end module mphys_psd_mod
