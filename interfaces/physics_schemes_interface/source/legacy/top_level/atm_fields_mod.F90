! DECLARE_ATM_FIELDS_MOD
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: top_level

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module atm_fields_mod

implicit none

! 1.7a: GLOMAP_CLIM aerosol section 54
real, pointer :: gc_nd_nuc_sol(:,:,:)  ! GLOMAP_CLIM NUC (Sol) number
real, pointer :: gc_nuc_sol_su(:,:,:)  ! GLOMAP_CLIM NUC (Sol) SO4
real, pointer :: gc_nuc_sol_oc(:,:,:)  ! GLOMAP_CLIM NUC (Sol) OC
real, pointer :: gc_nd_ait_sol(:,:,:)  ! GLOMAP_CLIM AIT (Sol) number
real, pointer :: gc_ait_sol_su(:,:,:)  ! GLOMAP_CLIM AIT (Sol) SO4
real, pointer :: gc_ait_sol_bc(:,:,:)  ! GLOMAP_CLIM AIT (Sol) BC
real, pointer :: gc_ait_sol_oc(:,:,:)  ! GLOMAP_CLIM AIT (Sol) OC
real, pointer :: gc_nd_acc_sol(:,:,:)  ! GLOMAP_CLIM ACC (Sol) number
real, pointer :: gc_acc_sol_su(:,:,:)  ! GLOMAP_CLIM ACC (Sol) SO4
real, pointer :: gc_acc_sol_bc(:,:,:)  ! GLOMAP_CLIM ACC (Sol) BC
real, pointer :: gc_acc_sol_oc(:,:,:)  ! GLOMAP_CLIM ACC (Sol) OC
real, pointer :: gc_acc_sol_ss(:,:,:)  ! GLOMAP_CLIM ACC (Sol) SS
real, pointer :: gc_nd_cor_sol(:,:,:)  ! GLOMAP_CLIM COR (Sol) number
real, pointer :: gc_cor_sol_su(:,:,:)  ! GLOMAP_CLIM COR (Sol) SO4
real, pointer :: gc_cor_sol_bc(:,:,:)  ! GLOMAP_CLIM COR (Sol) BC
real, pointer :: gc_cor_sol_oc(:,:,:)  ! GLOMAP_CLIM COR (Sol) OC
real, pointer :: gc_cor_sol_ss(:,:,:)  ! GLOMAP_CLIM COR (Sol) SS
real, pointer :: gc_nd_ait_ins(:,:,:)  ! GLOMAP_CLIM AIT (Ins) number
real, pointer :: gc_ait_ins_bc(:,:,:)  ! GLOMAP_CLIM AIT (Ins) BC
real, pointer :: gc_ait_ins_oc(:,:,:)  ! GLOMAP_CLIM AIT (Ins) OC


end module atm_fields_mod
