! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module phase_change_coefs_mod

implicit none

! Module stores super-array sizes and addresses used for storing
! and accessing the coefficients in the implicit solution for
! temperature and water-vapour mixing-ratio


! Number of coefficients for process rate for a single species
integer, parameter :: n_coefs = 3

! Addresses in super-array
integer, parameter :: i_q = 1
integer, parameter :: i_t = 2
integer, parameter :: i_0 = 3


end module phase_change_coefs_mod
