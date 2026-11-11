.. -----------------------------------------------------------------------------
    (c) Crown copyright 2025 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

PSyclone Transformation Scripts
===============================

The lfric_apps repository maintains the ability to provide PSyclone with
module-specific transformation scripts. This page gives an overview of the
structure of the directories that hold the transformation scripts.

.. toctree::
   :maxdepth: 0
   :hidden:

   psyclone_scripts
   psyclone_makefiles
   psyclone_functions


Optimisation directory structure
--------------------------------

Within each application there exists an ``optimisation/`` directory that holds
all PSyclone transformation scripts. These scripts are designed to target both
LFRic and non-LFRic source code on multiple architectures, and the directory is
structured to reflect this::

      optimisation/
      └── architecture/ (minimum, cpu, gpu, etc.)
          ├── psykal
          │   ├── global.py
          │   └── sub_directory/
          └── transmute/
              ├── global.py
              └── sub_directory/

The site folder has been updated to instead represent the architecture used
by a given site instead. Given there have not been too many distinct
requirements yet, initially there is:

* minimum
* cpu
* gpu

With scope to expand this like so:
(Please use underscore rather than a dash in the names, this
causes issues with FAB):

* minimum
* cpu
* cpu_ex
* gpu
* gpu_nvidia

Unless a module-specific transformation script exists, source files are
pre-processed with the default transformation script, ``global.py``.

