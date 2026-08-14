.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

==========
RADAER API
==========

:Author: Alan J Hewitt

Description of RADAER API
=========================

Radaer will be available to the parent application via a single API module as
a minimal set of top_level subroutines. These will include an initialisation
function to correctly set up radaer for the user inputs and a runtime module.

Names of subroutines presented via the API will start ```ukca_radaer_```

All run time communication between the parent model and UKCA will be via
argument lists.

All RADAER state variables will be available to the parent model between
time steps as native FORTRAN arrays, for inspection and possible modification.
The fields in these arrays will be in a specific order and size, which will
differ depending on user configuration. The lists of fields will be determined
by field names retrieved by the parent at run time.

The RADAER interface was designed to be very lean. No significant further
developments are intended for radaer. There is scope to revisit the lean
interface in future, if significant developments are required. However,
replacing RADAER with a different package is more likely to happen in future.
