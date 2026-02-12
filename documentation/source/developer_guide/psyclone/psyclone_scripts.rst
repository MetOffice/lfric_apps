.. -----------------------------------------------------------------------------
    (c) Crown copyright 2025 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------
.. _psyclone_scripts:

================================================================
Adding a PSyclone transformation script for a specific module
================================================================

Adding a Python transformation script to target a module involves two steps:

1. Add the script to the correct location,
2. Add the name of the script/module (minus extension) to the correct variable
   in ``psyclone_transmute_file_list.mk``.

----------------------------------------------------------------
Adding a PSyclone transformation script to the correct location
----------------------------------------------------------------

Each Python transformation script must reside in a matching location to the
target source file as found in the **built** application.

Each Fortran source file to be transformed must use either
* A global transformation scipt, such as ``global.py`` , or
* A local transformation script, such as ``local.py`` (present in the
correspondingly named folder to the target file), or
* A **matching filename minus extension**, such as ``ls_ppn.py``.

See below for how to add a PSyclone transformation script for the ``ls_ppn.F90``
module.


Storing a global transformation script, ``global.py``
+++++++++++++++++++++++++++++++++++++++++++++++++++++
All files to be optimised by transmute (see ``psyclone_transmute_file_list.mk``)
will use this script unless overwritten with a ``local.py`` or if using a matching
filename. ::

    optimisation/
    └── <platform>/
        └── transmute/
            └── global.py


Storing a local transformation script, ``local.py``
+++++++++++++++++++++++++++++++++++++++++++++++++++
All files to be optimised by transmute (see ``psyclone_transmute_file_list.mk``)
in the ``large_scale_precipitation`` directory will use this script unless overwritten
with a matching source filename.::

    optimisation/
    └── <platform>/
        └── transmute/
            └── large_scale_precipitation/
                └── local.py


For a script with matching source filename
+++++++++++++++++++++++++++++++++++++++++++

* The source file is found here **before** building - **This is NOT the path
  that should be used**::

    <LFRIC APPS>
    └── science/
        └── physics_schemes/
            └── source/
                └── large_scale_precipitation/
                    └── ls_ppn.F90


* To understand the required path for storing a matching source filename, the build directory of a built application should be consulted. In the **built** ``lfric_atm`` application, ``ls_ppn.F90`` is found here -
  **use this path**::

    <lfric_atm working directory>/
    └── large_scale_precipitation/
        └──ls_ppn.F90

* This module is not written in the PSyKAl format, so requires the ``transmute``
  method of operation.

Therefore, the transformation script for this module needs to be placed here
(note the matching filename)::

    optimisation/
    └── <platform>/
        └── transmute/
            └── large_scale_precipitation/
                └── ls_ppn.py

----------------------------------------------------------------------------------
Adding PSyclone transformation scripts to ``psyclone_transmute_file_list.mk``
----------------------------------------------------------------------------------

Instead of checking every module in the built application for a matching
PSyclone transformation script, each app maintains a list of modules on which to
apply module-specific PSyclone transformations. This can be found at::

    <application>
    └── build/
      └── psyclone_transmute_file_list.mk

Within this makefile there are several environment variables which control the
how files are passed to PSyclone Transmute. 
These are:

* ``PSYCLONE_PHYSICS_FILES``
* ``PSYCLONE_PASS_NO_SCRIPT``
* ``TRANSMUTE_INCLUDE_METHOD``
* ``PSYCLONE_DIRECTORIES``
* ``PSYCLONE_PHYSICS_EXCEPTION``


``PSYCLONE_PHYSICS_FILES`` is used with ``TRANSMUTE_INCLUDE_METHOD`` ``specify_include``.
File names (without the type), which are to be affected by PSyclone, are added to this list.
A PSyclone transformation script is applied to each respective file.
The global.py script with be used as default, overwritten with alternatives as advised above.
This is also currently the default method as we slowly transformation additional science files
across to using PSyclone.

``PSYCLONE_PASS_NO_SCRIPT`` has been implemented primarily for the ``TRANSMUTE_INCLUDE_METHOD``
``specify_include``, however can be expanded to either.
It controls whether files are to be passed to PSyclone, but not to have a script applied to them.
This can be useful to remove existing clauses from a file, but not add anything back.
files to be affected this way should only be added to this list, as they will be filtered out of
``PSYCLONE_PHYSICS_FILES`` anyway.

``TRANSMUTE_INCLUDE_METHOD`` can be set to either ``specify_include`` or ``specify_exclude``.
With ``specify_include``, ``PSYCLONE_PHYSICS_FILES`` will be used to define the list of files
to be passed to the PSyclone Transmute command, and with ``specify_exclude``, 
``PSYCLONE_DIRECTORIES`` and ``PSYCLONE_PHYSICS_EXCEPTION`` will be used.

``PSYCLONE_DIRECTORIES`` is used with ``TRANSMUTE_INCLUDE_METHOD`` ``specify_exclude``.
It provides a list of directories instead, where the build system populates the list of files to
be passed to PSyclone with all of the files in the directory.
``PSYCLONE_PHYSICS_EXCEPTION`` is used to filter the list of files found in each respective
``PSYCLONE_DIRECTORIES`` and remove them. This can often be used with .h files currently.


Module/transformation script names should be added in the following style (in
accordance with GNUMake)::

    PSYCLONE_PHYSICS_FILES = \
    ls_ppn \
    lsp_taper_ndrop \
    mphys_air_density_mod

-----------------------------------------------------------------------------
Adding PSyclone Transmute transformation functions
-----------------------------------------------------------------------------
Transformation optimisation scripts for Transmute, should be kept clean and
simple much like the existing PSyKAl scripts,
Ref: ``applications/lfric_atm/optimisation/meto-ex1a/psykal/global.py``

Functions should exist in their own Python script bucket(s) for Transmute,
like in PSyKAl in LFRic Core. Where this bucket is located is a WIP, but for
now will be here
``interfaces/physics_schemes_interface/build/psyclone_transmute``.
In the longer term, most may be held in the PSyTran repository, with
the intention to reduce code duplication 
here: `PSyTran <https://github.com/MetOffice/PSyTran>`_.

It is recommended for developers at this time to use the following process.

* Generate an LFRic Apps ticket to drive the work.
* Put script functions in a PSyclone folder in the Apps ticket here
  ``interfaces/physics_schemes_interface/build/psyclone_transmute``.
* The Python path in the ``psyclone_transmute.mk`` file may need to updated
  with this path.
* Where existing script files and functions exist, please utilise them.
* Top-level calls to these functions should occur in the transmute
  optimisation folder in the application, ``lfric_atm`` and ``ngarch``.
* Start a PSyTran issue to try and get duplicate code removed and new
  code reviewed. 
* Copy down any PSyTran code used, where the duplicated code has been
  removed in the Apps ticket so we have a local copy of the code from PSyTran,
  with a comment in the Apps ticket and PSyTran issue. 
* Currently PSyTran is not integrated into the build system.
* Link this ticket
  `Apps906 <https://code.metoffice.gov.uk/trac/lfric_apps/ticket/906#ticket>`_,
  and update it with the functions duplicated to capture technical debt.
