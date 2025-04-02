.. -----------------------------------------------------------------------------
    (c) Crown copyright 2025 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------
.. _psyclone_scripts:

Adding a PSyclone transformation script for a specific module
=============================================================

.. warning::
    This functionality is not yet available on LFRic Apps trunk - please see
    [LFRic Apps Ticket #624](https://code.metoffice.gov.uk/trac/lfric_apps/ticket/624)
    for updates.

Adding a transformation script to target a module involves two steps:

1. Add the script to the correct location
2. Add the name of the script/module (minus extension) to the correct variable
   in ``psyclone_transmute_file_list.mk``.


Adding a PSyclone transformation script to the correct location
---------------------------------------------------------------

Each transformation script must reside in a matching location to the target
source file as found in the **built** application.

Each transformation script must have a **matching filename minus extension**.

For example, to add a PSyclone transformation script for the ``ls_ppn.F90``
module:

* The source file is found here **before** building - **This is NOT the path
  that should be used**:

.. code-block::
    <LFRIC APPS>
    └── science/
        └── physics_schemes/
            └── source/
                └── large_scale_precipitation/
                    └── ls_ppn.F90


* In the **built** ``lfric_atm`` application, ``ls_ppn.F90`` is found here -
  **use this path**:

.. code-block::
    <lfric_atm working directory>/
    └── large_scale_precipitation/
        └──ls_ppn.F90

* This module is not written in the PSyKAl format, so requires the ``transmute``
  method of operation.

Therefore, the transformation script for this module needs to be placed here
(note the matching filename):

.. code-block::
    optimisation/
    └── <platform>/
        └── transmute/
            └── large_scale_precipitation/
                └── ls_ppn.py


Adding PSyclone transformation scripts to ``psyclone_transmute_file_list.mk``
-----------------------------------------------------------------------------

.. warning::
    This functionality is not yet available on LFRic Apps trunk - please see
    [LFRic Apps Ticket #624](https://code.metoffice.gov.uk/trac/lfric_apps/ticket/624)
    for updates.

Instead of checking every module in the built application for a matching
PSyclone transformation script, each app maintains a list of modules on which to
apply module-specific PSyclone transformations. This can be found at:

.. code-block::
    <application>
    └── psyclone_transmute_file_list.mk

Within this makefile there are environment variables that hold the names of
Fortran modules/transformation scripts that the app-specific PSyclone makefile
should target. These are:

.. code-block::
    PSYCLONE_PHSYICS_FILES_IMPORT
    PSYCLONE_PHSYICS_FILES_FCM

* If the module for which you are adding a transformation script lives in the
  LFRic Apps directory, please use ``PSYCLONE_PHSYICS_FILES_IMPORT``.

* If the module is being imported from UKCA, JULES, CASIM, UM, SOCRATES, or any
  other code base and does not have a copy in the LFRic Apps repository, please
  use the ``PSYCLONE_PHSYICS_FILES_FCM`` variable.

Module/transformation script names should be added in the following style (in
accordance with GNUMake):

.. code-block:: Make
    PSYCLONE_PHSYICS_FILES_IMPORT = \
    ls_ppn \
    lsp_taper_ndrop \
    mphys_air_density_mod
