.. -----------------------------------------------------------------------------
    (c) Crown copyright 2025 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------
.. _psyclone_makefiles:

Science interface PSyclone makefiles
====================================

Each interface to a non-LFRic science applicaton contains a Makefile that
extends the functionality of the ``psyclone.mk`` file in LFRic Core. More
specifically, these makefiles are intended to target Fortran source files that
require a different PSyclone mode of operation than ``psykal``, which is the
default method applied when building LFRic applications.

This makefile for non-psykal modules can be found at::

    <LFRic Apps trunk>
    └── interfaces/
      └── <interface_name>/
        └── build/
          └── psyclone_transmute.mk

.. note::
    It is likely that some interfaces will not have this makefile yet. These
    will be added as they are required.

The non-psykal PSyclone makefile uses the ``OPTIMISATION_PATH`` and ``DSL``
variables to narrow the search field, as well as a list of target modules to
reduce the number of searches it must perform.

The next sections will detail each of these points and how to add a PSyclone
transformation script for a single module.


Targeting platform with OPTIMISATION_PATH
-----------------------------------------

.. note::
    As a developer, you should not have to modify this variable. It is included
    here for reference.

The ``OPTIMISATION_PATH`` variable is consistent across all applications, and is
set as::

    optimisation/<target_platform>

Where ``<target_platform>`` refers to the hardware that the application will run
on (EX machines, Archer2, etc.). This information is picked up from the
task name within rose-stem. For example, the task::

    "lfric_atm_nwp_gal9-C224_MG_ex1a_cce_production-64bit"

has specified a platform of ``ex1a``, so the OPTIMISATION_PATH variable will be
updated by rose-stem to::

    optimisation/ex1a/

.. note::
    Most platforms will be able to share optimisations and can therefore
    contain symbolic links to other platform directories containing the
    required transformation scripts.


PSyclone method with DSL
------------------------

The ``DSL`` variable directs the build to search in specified directories
within the ``OPTIMISATION_PATH``.

The purpose of ``DSL`` is to separate transformations that would otherwise
conflict. For example, in the ``lfric_atm`` application, transformation scripts
are grouped by the PSyclone method of operation they require (psykal,
transmute).

The ``DSL`` variable defaults to the LFRic-specific method of operation (psykal).
For science interface PSyclone makefiles, ``DSL`` should be set to transmute
since these makefiles are designed to target non-LFRic source (and hence need a
different method of operation).

.. note::
    As a developer, you should not have to modify this variable. It is included
    here for reference.


List of modules to target
-------------------------

Each application contains a makefile for storing the names of each module that
has a transformation script. This list is used by the science interface PSyclone
makefiles to reduce the number of searches that it must perform by only
searching for modules included in the list. These makefiles can be found in each
application at::

    <application>
    └── build/
        └── psyclone_transmute_file_list.mk

.. note::
    This is only required for modules with transformation scripts held in the
    ``transmute/`` directory (non-LFRic source or modules that are not
    targeted by the LFRic Core ``psyclone.mk`` file). Some applications will
    not have this .mk file.


Modifying the PSyclone command(s)
---------------------------------

To assist with more flexible use of the PSyclone command(s) for the two DSL methods,
PSyKAl and Transmute, a flag has been added to allow users to adjust the flags used
with the command.
The main options are still required by default, but options that could be adjustable
per site, or potentially per file are being allowed, like compiler options.
An example is the use of ``-l all``, which has caused issues with some compilers,
such as Intel or Nvidia.

* ``PSYCLONE_TRANSMUTE_EXTRAS`` - used by Transmute
* ``PSYCLONE_PSYKAL_EXTRAS`` - used by PSyKAl

This has not been expanded further, however a user should be able to override them at
a site level, and a file level.
For just a site level, this could be set in the build option in the site cylc file.
Per file, which may be needed to preserve compiler clauses, and add optimisations in, a
makefile will need to export privatised versions of these per file and imported into
the application.