! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module um2lfric_masked_field_adjustments_mod

use lfricinp_masked_field_adjust_type_mod, &
                                        only: lfricinp_masked_field_adjust_type

implicit none

private

type(lfricinp_masked_field_adjust_type), public  :: land_field_adjustments,    &
                                                    maritime_field_adjustments

end module um2lfric_masked_field_adjustments_mod
