.. -----------------------------------------------------------------------------
    (c) Crown copyright 2025 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------
.. _rose-stem:

PSyclone Transformation Scripts
===============================

The lfric_apps repository maintains the ability to provide PSyclone with module-
specific transformation scripts. This page gives an overview of the
structure of the directories that hold the trasnformation scripts.

.. toctree::
   :maxdepth: 0
   :hidden:

   psyclone_scripts
   psyclone_makefiles


Optimisation directory structure
--------------------------------

Within each application there exists an ``optimisation/`` directory that holds
all PSyclone transformation scripts. These scripts are designed to target both
LFRic and non-LFRic source code on multiple platforms, and the directory is
structured to reflect this::

      optimisation/
      └── platform/ (ex1a, minimum, archer2, etc.)
          ├── psykal
          │   ├── global.py
          │   └── sub_directory/
          └── transmute/
              ├── global.py
              └── sub_directory/


Unless a module-specific transformation script exists, source files are
pre-processed with the default transformation script, ``global.py``.

