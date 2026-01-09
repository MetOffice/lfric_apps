.. -----------------------------------------------------------------------------
    (c) Crown copyright 2025 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------
.. _rose-stem_templating:

Description of Rose Stem Templating
===================================

The LFRic Apps rose-stem suite has been set up such that a new job can be added by providing configuration settings which allow graph and runtime sections to be populated automatically by the suite. This page will briefly describe how the templating produces the cylc workflow - the main steps are all driven from the ``flow.cylc`` with various templating files handling other parts.

Initial Setup
-------------

A few initial setup tasks are completed by the suite. These include defining a ``site_vars`` dictionary and populating this with some variables specific to the site from ``site/<SITE>/variables.cylc``, defining the ``site_groups`` dictionary and populating that from ``site/<SITE>/groups.cylc`` and importing some commonly used macros from ``templates/common_macros.cylc``. It also defines the sources of a few dependencies (including LFRic Core) by searching the ``dependencies.sh`` file using a custom Jinja2 Filter.

It is assumed that rose-stem has been launched with the ``--group=`` command line option specifying some groups or individual tasks to run. The workflow now populates a list of tasks to run by expanding any groups of tasks defined in ``site_groups`` - this is done by the ``expand_task_list`` macro. It is done recursively so any groups within a group are also expanded.

Once the list of requested jobs has been acquired it is used to populate a dictionary of ``requested_tasks`` with keys as the job names and values the task configuration defined by the task files - see :ref:`adding_new_test` for instructions on how to populate these files. This step is performed by the file ``templates/populate_task_definitions.cylc`` and the build information for the job is also parsed at this stage.

Scheduling
----------

With the configuration of the tasks now available the graph section of the cylc workflow can now be defined - this is done in ``templates/graph/generate_graph.cylc``. This file loops over all tasks run and populates the graph based on the name of this task and certain configuration settings. During this section a new dictionary ``tasks_to_run`` is populated, with keys and values similar to ``requested_tasks``. ``tasks_to_run`` will provide the definitive list of tasks that need to be defined by the suite as during the graph section certain tasks are added which may not be included by the ``--groups=`` argument. For example, a run task will require a build task, a mesh task and a mesh build task for it to be able to run.

Most jobs will have their graph produced in ``templates/graph/populate_graph_sections.cylc`` but certain graphs require slightly different treatment  and therefore have their own template files. These include lfricinputs, lfric_coupled, unit/integration tests (which are covered in ``populate_graph_technical_tests.cylc``) and scripts. If a task is configured to complete continuation-runs then these graphs are set in ``populate_crun_graph.cylc``.

Once the graph has been populated a number of tasks will have been added to ``tasks_to_run`` without any configuration. The suite now once again includes ``templates/populate_task_definitions.cylc`` to provide configurations for these tasks.

Runtime
-------

All tasks which will be run are now defined in the ``tasks_to_run`` dictionary. The runtime section initially defines families with site/platform specific settings which will be inherited by the jobs by including ``site/<SITE>/suite_config.cylc`` files. It then includes ``templates/runtime/set_task_families.cylc`` which defines some families which are used to control the layout of jobs in the gui.

All jobs in ``tasks_to_run`` are then looped over and the runtime section populated by including ``templates/runtime/generate_runtime.cylc``. Depending on the type of the current task, this file will call a particular ``generate_runtime_*.cylc`` file. The options for this are:
* ``application`` - the application 'run' task for most tasks.
* ``build`` - the build task for most tasks.
* ``control`` - currently defines the export source tasks and housekeeping tasks.
* ``lfric_coupled`` - fcm_make tasks specific to the lfric_coupled application.
* ``lfricinputs`` - All lfricinputs tasks, including build, run and rose-ana tasks.
* ``mesh`` - the run task to generate the mesh.
* ``other`` - plot, check and phystest tasks.
* ``scripts`` - The runtime section for scripts. These are defined manually rather than being templated.
* ``technical_tests`` - For unit and integration tests.
