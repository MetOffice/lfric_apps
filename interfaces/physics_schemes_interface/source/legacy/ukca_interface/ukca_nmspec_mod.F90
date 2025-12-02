! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: To initialize the nmspec array
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided by The University of Cambridge,
!  by, University of Leeds, University of Oxford, and
!  The Met. Office.  See www.ukca.ac.uk
!
!  Called from addres
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_UM
!
!
! Code description:
!   Language: Fortran
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------
!
module ukca_nmspec_mod

use ukca_tracer_stash,    only: a_max_ukcavars, a_ukca_first, a_ukca_last

implicit none

private
public ::  nm_spec, ukca_set_nmspec, nm_spec_active, nmspec_len,               &
           ukca_name2index

! Names for tracers (aka species)
! Once the  conflicting tracers from the RAQ chemistry are moved
! this should become a parameter array
integer, parameter :: nmspec_len=10
character(len=nmspec_len), save :: nm_spec(a_max_ukcavars)
! Names for active UKCA tracers - set at run time
character(len=nmspec_len), allocatable, save :: nm_spec_active(:)

character(len=*), parameter, private :: ModuleName='UKCA_NMSPEC_MOD'

contains

subroutine ukca_set_nmspec()

use ukca_option_mod,      only: l_ukca_raq, l_ukca_raqaero, tr_ukca_a
use UM_ParCore,           only: mype
use umPrintMgr,           only: umPrint, umMessage, PrintStatus, PrStatus_Oper
use yomhook,              only: lhook, dr_hook
use parkind1,             only: jprb, jpim

implicit none

integer :: i, j

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='UKCA_SET_NMSPEC'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!  This list must agree with the STASHmaster file.

! This is the generic list of tracers; values for RAQ and RAQ+aerosol chemistry
! are overwritten later.
! Tracers 98,99 & 100 are for lumped Nitrogen, Br and Cl for stratospheric
! chemistry, but can only be renamed in STASHmaster file, not in advt or
! nm_spec.
nm_spec(1:a_max_ukcavars) = [                                                  &
'O3        ','NO        ','NO3       ','NO2       ','N2O5      ',              &
'HO2NO2    ','HONO2     ','H2O2      ','CH4       ','CO        ',  & !10
'HCHO      ','MeOOH     ','HONO      ','C2H6      ','EtOOH     ',              &
'MeCHO     ','PAN       ','C3H8      ','n-PrOOH   ','i-PrOOH   ',  & !20
'EtCHO     ','Me2CO     ','MeCOCH2OOH','PPAN      ','MeONO2    ',              &
'O3_S      ','C5H8      ','ISOOH     ','ISON      ','MACR      ',  & !30
'MACROOH   ','MPAN      ','HACET     ','MGLY      ','NALD      ',              &
'HCOOH     ','MeCO3H    ','MeCO2H    ','H2O       ','ISO2      ',  & !40
'Cl        ','ClO       ','Cl2O2     ','OClO      ','Br        ',              &
'BrO       ','BrCl      ','BrONO2    ','N2O       ','HCl       ',  & !50
'HOCl      ','HBr       ','HOBr      ','ClONO2    ','CFCl3     ',              &
'CF2Cl2    ','MeBr      ','N         ','O(3P)     ','MACRO2    ',  & !60
'C2H4      ','C4H10     ','oXYLENE   ','TBUT2ENE  ','APINENE   ',              &
'BPINENE   ','C2H2      ','BENZENE   ','CH2Br2    ','H2        ',  & !70
'DMS       ','SO2       ','H2SO4     ','MSA       ','DMSO      ',              &
'NH3       ','CS2       ','cos       ','H2S       ','H         ',  & !80
'OH        ','HO2       ','MeOO      ','EtOO      ','MeCO3     ',              &
'n-PrOO    ','i-PrOO    ','EtCO3     ','MeCOCH2OO ','MeOH      ',  & !90
'Monoterp  ','Sec_Org   ','C3H6      ','SO3       ','C4H9OOH   ',              &
'MEK       ','TOLUENE   ','NO2       ','BrO       ','HCl       ',  & !100
'Nuc_SOL_N ','Nuc_SOL_SU','Ait_SOL_N ','Ait_SOL_SU','Ait_SOL_BC',              &
'Ait_SOL_OM','Acc_SOL_N ','Acc_SOL_SU','Acc_SOL_BC','Acc_SOL_OM',  & !110
'Acc_SOL_SS','Acc_SOL_DU','Cor_SOL_N ','Cor_SOL_SU','Cor_SOL_BC',              &
'Cor_SOL_OM','Cor_SOL_SS','Cor_SOL_DU','Ait_INS_N ','Ait_INS_BC',  & !120
'Ait_INS_OM','Acc_INS_N ','Acc_INS_DU','Cor_INS_N ','Cor_INS_DU',              &
'Nuc_SOL_OM','Ait_SOL_SS','Nuc_SOL_SO','Ait_SOL_SO','Acc_SOL_SO',  & !130
'Cor_SOL_SO','Nuc_SOL_NH','Ait_SOL_NH','Acc_SOL_NH','Cor_SOL_NH',              &
'Nuc_SOL_NT','Ait_SOL_NT','Acc_SOL_NT','Cor_SOL_NT','EtOH      ',  & !140
'i-PrOH    ','n-PrOH    ','HOCH2CHO  ','HOC2H4OOH ','EtCO3H    ',  & !145
'HOCH2CO3H ','NOA       ','EtONO2    ','PASSIVE O3','AGE OF AIR',  & !150
'i-PrONO2  ','MeO2NO2   ','HOC2H4NO3 ','PHAN      ','MeSCH2OO  ',  & !155
'MeS       ','MeSO      ','MeSO2     ','MeSO3     ','MSIA      ',  & !160
'CARB14    ','CARB17    ','CARB11A   ','CARB7     ','CARB10    ',  & !165
'CARB13    ','CARB16    ','UCARB10   ','CARB3     ','CARB6     ',  & !170
'CARB9     ','CARB12    ','CARB15    ','UCARB12   ','NUCARB12  ',  & !175
'UDCARB8   ','UDCARB11  ','UDCARB14  ','TNCARB26  ','TNCARB10  ',  & !180
'RN10NO3   ','RN13NO3   ','RN16NO3   ','RN19NO3   ','RA13NO3   ',  & !185
'RA16NO3   ','RA19NO3   ','RTX24NO3  ','RN10OOH   ','RN13OOH   ',  & !190
'RN16OOH   ','RN19OOH   ','RN8OOH    ','RN11OOH   ','RN14OOH   ',  & !195
'RN17OOH   ','RU14OOH   ','RU12OOH   ','RU10OOH   ','NRU14OOH  ',  & !200
'NRU12OOH  ','RN9OOH    ','RN12OOH   ','RN15OOH   ','RN18OOH   ',  & !205
'NRN6OOH   ','NRN9OOH   ','NRN12OOH  ','RA13OOH   ','RA16OOH   ',  & !210
'RA19OOH   ','RTN28OOH  ','NRTN28OOH ','RTN26OOH  ','RTN25OOH  ',  & !215
'RTN24OOH  ','RTN23OOH  ','RTN14OOH  ','RTN10OOH  ','RTX28OOH  ',  & !220
'RTX24OOH  ','RTX22OOH  ','NRTX28OOH ','RAROH14   ','RAROH17   ',  & !225
'RU12PAN   ','RTN26PAN  ','TNCARB12  ','TNCARB11  ','RTN23NO3  ',  & !230
'CCARB12   ','TNCARB15  ','RCOOH25   ','TXCARB24  ','TXCARB22  ',  & !235
'RN9NO3    ','RN12NO3   ','RN15NO3   ','RN18NO3   ','RU14NO3   ',  & !240
'RTN28NO3  ','RTN25NO3  ','RTX28NO3  ','RTX22NO3  ','AROH14    ',  & !245
'ARNOH14   ','AROH17    ','ARNOH17   ','ANHY      ','Acc_SOL_NN',  & !250
'Cor_SOL_NN','DHPCARB9  ','HPUCARB12 ','HUCARB9   ','IEPOX     ',  & !255
'HMML      ','DHPR12OOH ','DHCARB9   ','RU12NO3   ','RU10NO3   ',  & !260
'SEC_ORG_I ','Sup_INS_N ','Sup_INS_DU','Ait_SOL_MP','Acc_SOL_MP',  & !265
'Cor_SOL_MP','Ait_INS_MP','Acc_INS_MP','Cor_INS_MP','Sup_INS_MP',  & !270
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !275
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !280
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !285
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !290
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !295
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !300
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !305
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !310
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !315
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !320
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !325
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !330
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !335
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !340
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !345
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !350
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !355
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !360
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !365
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !370
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !375
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !380
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !385
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !390
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !395
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !400
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !405
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !410
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !415
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !420
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !425
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !430
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !435
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !440
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !445
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !450
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !455
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !460
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !465
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !470
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !475
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !480
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !485
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !490
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !495
'XXX       ','XXX       ','XXX       ','XXX       '                & !499
 ]

if (l_ukca_raq .or. l_ukca_raqaero) then
  ! Overwrite some tracers for RAQ and RAQ-AERO chemistry.
  ! If MODE aerosols are used with it but their positions
  ! change in the array then the list needs to be updated.

  ! Version above has underscore O3_S
  nm_spec(26) ='O3S       '

  nm_spec(39) ='MVK       '
  nm_spec(40) ='MVKOOH    '
  nm_spec(60) ='ORGNIT    '
  nm_spec(69) ='CH3OH     '

  nm_spec(90) ='RNC2H4    '
  nm_spec(93) ='C3H6      '
  nm_spec(94) ='C4H10     '
  nm_spec(95) ='s-BuOOH   '
  nm_spec(96) ='MEK       '
  nm_spec(97) ='TOLUENE   '
  nm_spec(98) ='MEMALD    '
  nm_spec(99) ='GLY       '
  nm_spec(100)='oXYLENE   '

  ! Two tracers are inconsistent between the two RAQ schemes because they
  ! appear in a different order in ukca_chem_raq and ukca_chem_raqaero
  if (l_ukca_raq) then
    nm_spec(91) ='RNC3H6    '
    nm_spec(92) ='C2H4      '
  else ! using l_ukca_raqaero
    nm_spec(81) ='RNC3H6    '
    nm_spec(82) ='C2H4      '
  end if
end if ! l_ukca_raq or l_ukca_raqaero

! Mode components: SU: sulphate, BC: black carbon, OM: organic matter
!                  SS: sea-salt, DU: dust,         SO: organic matter 2
!                  NH: ammonium, NT: nitrate,      N:  number mixing ratio

! Make an array for the "real" nm_spec accounting for which
! tracers are on. This allows us to map from the tracer array to
! name of species and back
allocate (nm_spec_active(count(tr_ukca_a)))

! loop over all entries in nm_spec, testing whether they are on or not
! If they are on, add their names into nm_spec_active
j = 0

! If high value of print status, output on PE0 only
if (PrintStatus >= PrStatus_Oper .and. mype == 0) then
  write(umMessage,'(A)') 'NM spec active: '
  call umPrint(umMessage,src=RoutineName)
end if

do i = a_ukca_first, a_ukca_last
  if (tr_ukca_a(i)) then
    j = j + 1
    nm_spec_active(j) = nm_spec(i)

    ! If high value of print status, output on PE0 only
    if (PrintStatus >= PrStatus_Oper .and. mype == 0) then
      write(umMessage,'(I4, 1X, A)') j, nm_spec_active(j)
      call umPrint(umMessage,src=RoutineName)
    end if
  end if
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine ukca_set_nmspec


function ukca_name2index(tracer_name)

! This function is used in iau.F90 to update
! the Air Quality (AQ) Analysis increments.
! Find the index of a tracer based on the name given.
! if not found, stop with an ereport.
use ereport_mod,  only: ereport
use errormessagelength_mod, only: errormessagelength

use yomhook,              only: lhook, dr_hook
use parkind1,             only: jprb, jpim

implicit none

integer :: ukca_name2index

character(len=*), intent(in) :: tracer_name

! Local variables
character(len=*), parameter :: RoutineName='UKCA_NAME2INDEX'

integer :: i
! ErrorStatus
integer                    :: errcode
character(len=errormessagelength)   :: cmessage      ! Error return message

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check if nm_spec_active has been allocated,
! if not, call ereport
if (.not. allocated (nm_spec_active) ) then
  errcode = 1
  cmessage = 'Tried to look up contents of nm_spec_active but not allocated'
  call ereport(ModuleName//':'//RoutineName, errcode, cmessage)
end if

! Loop over all of nm_spec_active
do i = 1, size(nm_spec_active)
  ! If names match
  if (trim(adjustl(tracer_name)) == trim(adjustl(nm_spec_active(i)))) then
    ukca_name2index = i
    if (lhook) call dr_hook(ModuleName//':'//                                  &
                            RoutineName,zhook_out,zhook_handle)
    return ! Exit the function
  end if
end do

! Call ereport as we haven't found the name if we get here
errcode = 2
cmessage = 'Failed to find ' // tracer_name // ' in nm_spec_active'
call ereport(ModuleName//':'//RoutineName, errcode, cmessage)

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end function ukca_name2index


end module ukca_nmspec_mod
