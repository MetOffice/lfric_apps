! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Module to
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_UM
!
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------

module ukca_tracer_stash

implicit none

integer, parameter :: a_ukca_first   = 1    ! First UKCA tracer (STASH No)
integer, parameter :: a_ukca_last    = 499  ! Last UKCA tracer  (STASH No)

! Maximum number of UKCA tracers
integer, parameter :: a_max_ukcavars = a_ukca_last - a_ukca_first + 1

! Index to relative position.
! A_TR_INDEX(N) gives position in JTRACER for tracer number N.
! Set in SET_ATM_POINTERS.
! A_TR_INDEX(N) is the position, in the list of tracers
! actually present in D1, that tracer number N (in the list
! of all tracers selectable from the user interface) occupies,
! if it is present.
! If tracer number N is absent then A_TR_INDEX(N) is -1.
! Similarly for A_UKCA_INDEX.
! UKCA_tr_StashItem is set up in SET_ATM_POINTERS

integer :: a_ukca_index      (a_max_ukcavars)
integer :: ukca_tr_stashitem (a_max_ukcavars)

! ukca_tr_lbc_stashitem is set up in subroutine inbounda and is only
! referenced if LBC code is active.
integer :: ukca_tr_lbc_stashitem    (a_max_ukcavars)
integer :: ukca_tr_active_lbc_index (a_max_ukcavars)

end module ukca_tracer_stash

