! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates the total PM10 and PM2.5 concentrations and for each aerosol type

module calc_pm_diags_mod

use um_types, only: real_umphys

implicit none

! Global parameters available to other routines.
! Note that all variables whose name starts with 'denom_' are parameters
! related to the volume distribution and were calculated offline by the
! standalone program calc_pm_params.f90.
!
! Parameters for sulphate
! What is advected is sulphur. Therefore, it is necessary to convert m.m.r
! of sulphur to m.m.r of ammonium sulphate by multiplying by the ratio of
! their molecular weights: [(NH4)2]SO4 / S = 132 / 32 =  4.125
real(kind=real_umphys), parameter :: s_conv_fac      = 4.125
real(kind=real_umphys), parameter :: denom_su_ait    = 0.832154093803605184e-20
real(kind=real_umphys), parameter :: denom_su_acc    = 0.317222975403968212e-16
!
! Parameters for black carbon (soot)
real(kind=real_umphys), parameter :: denom_bc_fr     = 0.132771498349419138e-16
real(kind=real_umphys), parameter :: denom_bc_ag     = 0.132771498349419138e-16
!
! Parameters for OCFF
real(kind=real_umphys), parameter :: denom_ocff_fr   = 0.231243545609857872e-16
real(kind=real_umphys), parameter :: denom_ocff_ag   = 0.399588846813834411e-16
!
! Parameters for biogenic aerosol (SOA)
real(kind=real_umphys), parameter :: denom_soa       = 0.293507300762206618e-16

! Paremeter for unit conversion, from kg to micrograms
real(kind=real_umphys), parameter :: kg_to_micg        = 1.0e9

character(len=*), parameter, private :: ModuleName='CALC_PM_DIAGS_MOD'


end module calc_pm_diags_mod

