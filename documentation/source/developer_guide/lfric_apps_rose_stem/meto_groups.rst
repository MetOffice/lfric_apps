.. -----------------------------------------------------------------------------
    (c) Crown copyright 2025 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------
.. _meto_groups:

Met Office Rose Stem Groups
===========================

.. note::

   These notes only apply to the meto site - other sites may have set up their groups differently

The lfric_apps rose-stem uses groups it's tasks. Any combination of these groups or individual tasks can be selected by ``--group=group1,group2,task1,task2`` when running ``rose stem`` from the command line.

For running the test-suite for submitting a ticket, there are two main standard groups, one of which will likely need to be run:

.. list-table::
   :widths: 15 45
   :header-rows: 1

   * - Group
     - Description
   * - ``developer``
     - Suitable for most tickets, this runs a selection of groups providing a broad coverage.
   * - ``all``
     - This runs all tasks with stored kgo (including those in developer). This should be run when you have kgo changes.

In addition to these, various standard groups have also been set up to aid with more targeted testing when developing a change. These standard groups follow the naming format ``APPLICATION_PLATFORM_JOB-TYPE``. Any number of these sections can be used, but they will always be in that order. For example:

.. list-table::
   :widths: 17 45
   :header-rows: 1

   * - ``gungho_model``
     - All gungho_model tasks
   * - ``gungho_model_spice_developer``
     - The gungho_model developer tasks on spice
   * - ``xc40_unit_tests``
     - The unit tests on the xc40
   * - ``unit_tests``
     - All unit tests

The following job-types are available (they are not all necessarily available for all applications):

.. list-table::
   :widths: 20
   :header-rows: 1

   * - Job-Type
   * - ``developer``
   * - ``all``
   * - ``unit_tests``
   * - ``integration_tests``
   * - ``build``

Finally, it is anticipated that various science sections within an application will get their own groups - to view what other groups are available see the relevant applications groups file in ``rose-stem/site/meto/groups/``.
