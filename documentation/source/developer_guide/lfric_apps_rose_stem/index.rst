.. -----------------------------------------------------------------------------
    (c) Crown copyright 2025 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------
.. _rose-stem:

LFRic Rose Stem Test Suite
==========================

.. tip::

   Working Practices for testing changes using rose-stem can be found at https://metoffice.github.io/simulation-systems/WorkingPractices/TestSuites/lfric_apps.html

The lfric_apps repository contains a redesigned rose-stem suite written for Cylc 8. There is now a single rose-stem directory that will run on all applications within lfric_apps and the suite can be launched from anywhere within the working copy. The suite can run entire groups of tests, such as a developer group from before, or it can run any single task or subset of tasks. This page will describe how to launch the new rose-stem suite and give an overview of the structure.

.. toctree::
   :maxdepth: 0
   :hidden:

   new_site
   rose-stem_templating
   adding_new_test
   meto_groups
   perftools_rose-stem

Running the rose-stem suite
---------------------------

.. important::

   The lfric_apps rose-stem is only compatible with Cylc 8, therefore if Cylc 7 is the default on your platform you will need to run export CYLC_VERSION=8 before running the suite. Detailed Cylc 8 instructions can be found on the `Cylc User Guide <https://cylc.github.io/cylc-doc/stable/html/7-to-8/cheat-sheet.html>`_

The test suite is now launched with the ``rose stem`` command followed by ``cylc play``. For users familiar with other SSD maintained codes, there are very similar command line options available here as there. At the most basic, the commands looks like:

.. code-block:: bash

   rose stem --group=<GROUP>
   cylc play <WORKFLOW_NAME>
   cylc gui

* ``GROUP`` is the name of the groups or tasks to run. These can have multiple, comma-separated items with each item required to have an entry in your sites groups.cylc file. The groups available for the `meto` site are described at :ref:`meto_groups`.
* ``WORKFLOW_NAME`` is the name of the installed suite - with the above commands this would be the name of the working copy you are currently in. See the `​Cylc Documentation <https://cylc.github.io/cylc-doc/stable/html/index.html>`_ for more details on the naming of workflows in Cylc 8.

For example, if we were somewhere in an lfric_apps working copy with the name ``lfric_apps_wc``, then we could launch the developer group of tests as:

.. code-block:: bash

   rose stem --group=developer
   cylc play lfric_apps_wc

Other command line options are also available for the ``rose stem`` command:

* ``--source=`` This specifies the path to an lfric_core source. This can be either an fcm url or a path to a working copy. The source variable currently can '''only''' be used for the lfric_core source - others, eg. UM, will need to be set in the ``dependencies.sh`` file. To use, it, use the syntax as ``--source=. --source=fcm:lfric.xm_tr`` for the lfric_core trunk. The inclusion of the first source is required.
* ``-S <VARIABLE>`` This can set a value for any jinja2 variable in the ``[template variables]`` section of the ``rose-suite.conf``. Currently available options are below.

   * ``-S HOUSEKEEPING=<true/false>``, with the default value as ``true``. ``-S HOUSEKEEPING=false`` will turn off the Housekeeping tasks.
   * ``-S METO_USE_XCS=<false/true>``, with the default value as ``false``. ``-S METO_USE_XCS=true`` will run the test suite on Cray XCS instead of XCE/F (please note that the task names will still contain ``xc40``).
   * ``-S OVERRIDE_LOG_LEVEL=\"info\"``, which sets the logging level to "info" for all jobs. The default depends on the build of the job (e.g. fast-debug jobs are usually at "info", and full-debug are at "debug".


Rose Stem Directory Structure
-----------------------------

An example, parsed directory tree with a gungho_model app for the meto spice platform. The templating scheme is described at :ref:`rose-stem_templating`.

.. code-block::

   flow.cylc
   rose-suite.conf
   apps/
      compile/
      check_kgo/
      gungho_model/
         rose-app.conf
         opt/
            rose-app-*.conf
      plot/
   bin/
      suite_report.py
      stylist_launcher.py
   Jinja2Filters/
   site/
      meto/
         groups.cylc
         suite_config.cylc
         variables.cylc
         groups/
            groups_gungho_model.cylc
         gungho_model/
            tasks_gungho_model.cylc
         kgos/
            gungho_model/
               spice/
                  checksum_gungho_model_*.txt
         macros/
            macros_spice.cylc
   templates/
      graph/
      runtime/

