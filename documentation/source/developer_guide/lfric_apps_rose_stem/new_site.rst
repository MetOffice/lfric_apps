.. -----------------------------------------------------------------------------
    (c) Crown copyright 2025 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------
.. _new_site:

Adding a New Rose Stem Site
===========================

Setting up a new site involves adding various files to your sites ``site/<SITE>/`` directory. The name of the site must match the name of the site as defined in your rose setup. This page will describe the files required, what is expected to be in them and how to override the expectations where possible. For examples of a working site see ``meto`` - this has extra files as it has multiple platforms defined.


groups.cylc
-----------

This file contains entries to the ``site_groups`` dictionary, defining groups of tasks. As shown by ``meto`` it can also contain calls to other files in order to more neatly organise the groups. For information on how to add a task to the groups file see :ref:`adding_new_test`.

kgos/
-----

This directory contains sub-directories for each application with stored checksum kgos. Each application directory then contains a directory for each platform for this site, inside which are the text file kgos.

macros
------

This directory contains ``macros_<PLATFORM>.cylc`` files for all platforms at this site. Currently these files should define 3 macros which set various directives potentially required by tasks, using the syntax specific to the platforms job runner:

* ``set_wallclock(num_minutes)`` - this sets the wallclock for the task
* ``set_task_cores(mpi_ranks)`` - this sets the number of nodes to be used. Optional arguments are ``threads, xios_nodes, ocean_nodes, river_nodes``.
* ``set_memory(mem)`` - this sets the memory for the task. The argument is a list with integer amount of memory and string for memory units.

Not all macros will be valid for a particular job, for instance set_memory may not be valid if running on a normal queue.

suite_config.cylc
-----------------

This file defines various families which are required to run on particular platform, such as environment variables, pre-scripts and default task directives. The ``meto`` version includes 2 sub-files, one for each of the platforms at the site. The exact structure of the families in these files can be decided by each site. The ``meto`` site uses a hierarchical structure with families inheriting families higher in the tree. By default the suite expects that within these files the following families are defined for each PLATFORM and COMPILER being used:

* ``<PLATFORM>_RUN_<COMPILER>``
* ``<PLATFORM>_BUILD_<COMPILER>``
* ``<PLATFORM>_MESH_<COMPILER>``
* ``<PLATFORM>_TECH_TESTS_<COMPILER>``
* ``<PLATFORM>_PLOT``
* ``<PLATFORM>_CHECK``
* ``<PLATFORM>_TECH``

However it is possible to override these by setting ``<type>_families`` variables in the task definition files. The possible types can be found in the "Default Inheritance Families" section of ``templates/default_task_definitions.cylc``. For examples of how to do this see ``site/common/lfric_coupled/tasks_lfric_coupled.cylc``.

Also in the suite_config files is where the pre-script is defined to setup any environments, etc. required at your site. To help prevent duplication the ``generate_script()`` macro is available - see the ``meto`` files for examples on how to do this.

Due to differences in their setup, lfric_coupled and lfricinputs also have their own specific suite_config files.

Task Files
----------

By default the suite will look for task configuration files in ``site/<SITE>/<APPLICATION>/tasks_<APPLICATION>.cylc``. It will first search in ``site/common/`` with the site specific files then able to add or overwrite settings on a site/platform specific basis. The task can have a definition in one or both of these places, but it must exist somewhere. The structure of the task definitions files is discussed in :ref:`adding_new_test`.

variables.cylc
--------------

Here some site specific jinja2 variables are defined and stored in the ``site_vars`` dictionary which has already been created in the ``flow.cylc``. The suite expects the following variables to be added as keys to this dictionary:
* ``known_applications`` - This is a list of all the applications defined at this site.
* ``site_platforms`` - A list of all the platforms at this site.
* ``scripts_platform`` - Optional, the platform to run scripts tasks on. Defaults to the first item in the ``site_platforms`` list.
* ``mesh_build`` - Meshes are all built with the same build settings on each platform regardless of the job settings. This is a dictionary with key as the platform and value the mesh build used, eg. ``"spice": "intel_fast-debug-64bit"``. All platforms in ``site_platforms`` need an entry here.

If your site is using lfricinputs then some additional variables are required for defining the locations of the lfricinputs kgo. See ``meto`` for details.

