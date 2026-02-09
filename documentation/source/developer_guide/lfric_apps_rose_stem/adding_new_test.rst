.. -----------------------------------------------------------------------------
    (c) Crown copyright 2025 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------
.. _adding_new_test:

Adding a New Rose Stem Test
===========================

To add a new test requires modifying a few files:
* The ``groups.cylc`` file with the task name in the required group. The naming of tasks is described below. For the meto site, these groups are spread around multiple files in the ``groups`` directory with one per application, so add your test to the most logical place. You might want to add it to multiple groups if appropriate.
* The ``tasks_<APPLICATION>.cylc`` file with the job configuration. The configuration options are detailed below. There are 2 tasks files which are relevant for each site, one in ``site/common/<APPLICATION>/`` and one in ``site/<SITE>/<APPLICATION>/``. A configuration for a job can exist in one or both of these places, but it must be defined once. The files in ``common`` are included first with the site specific files able to then add or overwrite settings on a site specific basis. For sites with multiple platforms, there may be a tasks file for each platform, allowing customisation of the job based per platform.
* A config file located within ``rose-stem/apps/``. This is as in the old test suite, where a new app with it's unique ``rose-app.conf`` file is added, or a new ``rose-app-opt.conf`` file is added to an existing app. The default app will share the name of the application being run - this can be overwritten by setting the ``app_name`` task setting.


Task Name
---------

Tasks (with a few exceptions) follow the following naming convention. The task name you get from this is what is entered into the groups files and how the task will appear when running in Cylc.

.. code-block::

   <APPLICATION>_<CONFIG-NAME>_<PLATFORM>_<COMPILER>_<EXTRAS>

* **APPLICATION** - This is the lfric application being run by this test, eg. gungho_model, lfric_atm, etc. Make sure the application is also added to ``known_applications``, defined in ``site/<SITE>/variables.cylc``.
* | **CONFIG-NAME** - This should encapsulate the configuration of this task. Early tasks have added this as ``<opt_conf>-<resolution>``, however it can take any string (so long as it is a unique name for this particular application). As the suite grows, other naming schemes might become more suitable - this is left flexible for future developers and reviewers.
  | If the application is mesh then this **must** be set to <resolution> which must have a corresponding file at ``rose-stem/app/mesh/opt/rose-app-<resolution>.conf``.
* **PLATFORM** - Platform being run on, eg. ``azspice``, ``ex1a`` at meto. This needs to match an entry in the ``site_platforms`` list in ``site/SITE/variables.cylc``.
* **COMPILER** - The compiler being for this task. **Note**: The only naming restriction for tasks exists here, the compiler name can't contain any underscores. **Note** If it is a Perftools test, modify this section, see :ref:`perftools_rose-stem`.
* | **EXTRAS** - Extra details required by the job - primarily the build config which needs the optimisation and precision. Other information may be automatically populated here by the suite, eg. a continuation run restart number.
  | The optimisation should be one of ``fast-debug``, ``full-debug`` or ``production``.
  | The suite allows the precision to be set independently for each of ``rbl``, ``rdef``, ``rphysics``, ``rsolver`` and ``rtran`` (although not all options may be valid/runnable). The precision part of the EXTRAS string should be of the format ``<default>bit-<option>NN...``. The default value will apply to any option that isn't otherwise defined. For instance, to run just ``rbl`` at 32bit, the string would look like ``64bit-rbl32``. In order to ensure identical builds are not run with different names, the suite expects a certain build config order of ``Optimisation-MainPrecision-OtherPrecisions``. The main precision should be the most common and the Other precisions should be in alphabetical order, as above. If the build config isn't correct the suite won't launch and will write the name as it expects it to be.


Task File
---------

When adding a new test, an entry might be required in the sites ``tasks_<APPLICATION>.cylc`` file. If the config already existed and this is just adding another platform/compiler/build variant then the task entry will already exist (although it may still need modifying). If it's a completely new test then the entry will require adding. Jobs for generating the mesh don't require an entry as they can use default values  (although an entry can overwrite this).

A task definition involves defining a dictionary entry. The keys available to the dictionary are described further down. This next snippet shows an example tasks file for gungho_model on spice with a single configuration added:

.. code-block:: Jinja

    {% if task_ns.conf_name == "baroclinic-C24_MG" %}

        {% do task_dict.update({
            "opt_confs": ["baroclinic"],
            "resolution": "C24_MG",
            "DT": 3600,
            "threads": 4,
            "mpi_parts": 6,
            "tsteps": 240,
            "kgo_checks": ["checksum"],
            "plot_str": "baroclinic.py $NODAL_DATA_DIR/lfric_diag.nc $PLOT_DIR dcmip dry surface_pressure_temperature",
        }) %}

        {% if task_ns.compiler == "intel" %}
            {% do task_dict.update({"wallclock": 11}) %}
        {% endif %}

    {% elif task.startswith("build") %}

        {% do task_dict.update({
            "placeholder": true,
        }) %}

    {% endif %}

    # Set default values for Gungho Model tasks
    {% if task_dict %}
        {% do task_dict.update({"application_dir": "apps/applications/gungho_model"}) %}
        {% if "app_name" not in task_dict %}
            {% do task_dict.update({"app_name": "gungho_model"}) %}
        {% endif %}
    {% endif %}

Tasks are identified by their CONFIG-NAME as described above; for the task entered here this would be ``baroclininc-C24_MG``. The dictionary entry is then defined by the ``{% do task_dict.update({}) %}`` command. This entry also contains an example of how to change the task based on build/compiler settings - in this case increasing the wallclock to 11 minutes if using the intel compiler.

Below this the build task is defined - unless you're adding a brand new application this entry will already exist. The placeholder value is required if there are no other settings so that the dictionary has an entry in it - this is ignored later.

The final section defines some defaults for all tasks for this application. The ``application_dir`` shouldn't be set elsewhere and sets the path the the application source code in the cylc-run directory. The ``app_name`` defines the rose app that will be run - this can be overwritten if a different app is required.


Task Options
------------

The table below shows a list of possible entries for ``task_dict`` in the task definition files. Most options have the same name as the old LFRic test suite, but some have been renamed. A definitive list of default values can be found in ``rose-stem/templates/default_task_definitions.cylc``. Further options for Coupled jobs are also available and these can be found in ``rose-stem/templates/default_task_definitions_lfric_coupled.cylc``. The Required values are only relevant for Application Run tasks.

.. list-table::
   :widths: 8 8 8 25
   :header-rows: 1

   * - Key
     - Type
     - Default
     - Description
   * - **Required**
     -
     -
     -
   * - opt_confs
     - List of Strings
     - NA
     - Each string is the name of an optional config for this applications app (or for the app defined by app_name below).
   * - resolution
     - Str
     - NA
     - The mesh resolution being used. It should be an optional config of the mesh app. Not required for a canned test.
   * - **Optional**
     -
     -
     -
   * - custom_task
     - Boolean
     - False
     - If true the graph and runtime sections are populated from the custom_tasks files. The task can be defined either in the common custom_tasks files or in those for a given site. Not recommended to set to true.
   * - DT
     - Integer
     - 1
     - Time Step Size.
   * - app_name
     - Str
     - The "main" app associated with the application, usually the name of the application.
     - The name of the rose app for this task. At the moment there is only 1 app per application but this may change so this setting provides for that. For canned tests set to "canned_test" - in this case it is recommended that the conf_name contains "canned" for clarity.
   * - plot_str
     - Str
     - Empty Str
     - Acts as a boolean for whether to run a plot task. If empty then no task will be run. The string should be the call to a python plot script in app/plot/bin.
   * - phystest_str
     - Str
     - Empty Str
     - Acts as a boolean for whether to run a phystest task. If empty no task will be run. The string should be the name of a plot script in app/phystest/bin.
   * - kgo_checks
     - List
     - Empty List
     - Only valid entry is "checksum" for checksum file kgo comparison. A "fields" option for full file comparisons may be added if required.
   * - wallclock
     - Int
     - NA
     - The allowed wallclock time for this task in minutes.
   * - threads
     - Int
     - 1
     - The number of threads for this task.
   * - mpi_parts
     - Int
     - 1
     - The number of mpi parts for this task.
   * - tsteps
     - Int
     - 1
     - The number of timesteps. This is always the total number of timesteps - if crun is set the continuation runs will add up to this.
   * - crun
     - Int
     - 1
     - Sets the number of continuation runs - a value of <=1 indicates no restarts. If >1, will complete a normal run and a run with crun-1 restarts and compare their outputs. The restart run lengths will be calculated from tsteps//crun. The final run will calculate the current number of tsteps completed and subtract from the total requested, therefore not all runs are necessarily exactly the same length.
   * - crun_compare
     - Boolean
     - True if ``crun > 1`` |br| False otherwise
     - Perform a comparison between an nrun and a restarted crun. With the default settings the nrun task will compare to stored output and the final crun task will compare with the nrun. |br| If ``crun > 1`` and this is false, only the crun tasks will be performed - the nrun task will not. In this case, the output of the crun will be compared against stored kgo. |br| In all cases comparison to stored kgo is controlled by the kgo_checks setting.
   * - log_level
     - Str
     - Default based on the optimisation level: |br| fast-debug -> "info" |br| full-debug -> "trace" |br| production -> "warn"
     - Sets the log level to be used by the run task. Valid options are any valid logging levels understood by LFRic. |br| Can also run rose stem with the command line option ``-S OVERRIDE_LOG_LEVEL=\"LEVEL\"`` which would set all tasks to run with ``LEVEL`` log level. Note, the escaped quotes are important.
   * - panel_decomp
     - Str
     - "auto"
     - Panel Decomposition Settings.
   * - xproc
     - Int
     - ``mpi_parts``
     - x-axis of decomposition.
   * - yproc
     - Int
     - 1
     - y-axis of decomposition.
   * - memory
     - List of an Int and a Str, eg. ``[6, "GB"]``
     - Any default will be set in the suite_config file for the relevant platform.
     - The memory for this task. If running on a normal queue, this is not a valid option.
   * - pre_process_macros
     - Str
     - Empty Str
     -
   * - use_xios
     - Boolean
     - True
     - Controls whether to use xios. If true, sets ``use_xios_io=true`` and ``nodal_output_on_w3=false``
   * - xios_nodes
     - Int
     - 0
     - See next item.
   * - mpi_parts_xios
     - Int
     - 0
     - If ``use_xios`` is true and both this and ``xios_nodes`` are greater than 0, ``xios_server_mode`` will be true and this will set ``xios_server_ranks``. |br| Otherwise ``xios_server_mode=false`` and ``xios_server_ranks=0``.
   * - example_dir
     - Str
     - "example"
     - Directory where canned test files are stored. Path assuming starting in the application source directory.
   * - <type>_families
     - List of Strings
     - Defaults depend on type of family being used and are defined in the templates directory.
     - Each string is a family which will be inherited by the task. These families are set up by each site and if using the defaults, will be assumed to exist at this site. The valid options for <type> are run, build, mesh, plot, checksum, phystest, techtests.
   * - psyclone_transformation
     - Str
     - <SITE>-<PLATFORM>, eg. ``meto-spice``
     - This is only valid for build tasks and sets the PSYCLONE_TRANSFORMATION variable in the build tasks. Currently the only option used is "none".
   * - task_ranks_per_node
     - Int
     - NA
     - For manually setting an override for the number of MPI ranks/tasks to run per node. This is set by default on EX1A for multi-node jobs. This value further is adjusted by threads or depth and needs to be scaled appropriately. Example, under-committing XC40, we might set 18 per node to allow more memory for each rank given hardware limitations. For a threaded job of 3 threads, we would set this instead to 6 ranks per node, to maintain 18 CPUs used. Can be used in conjunction with ``task_ranks_depth_pad``.
   * - task_ranks_depth_pad
     - Int
     - NA
     - For adding a manual depth 'pad' to a task, generally used on EXs or XC40s where depth option is available in mpi launch. This is to be used in combination with ``task_ranks_per_node``, but does work on it's own. If only ``task_ranks_per_node`` is set, the ranks may bunch up together on the node, and not utilise all of the extra free memory properly (see `EX Placement <EXlink_>`_.). The 'pad' will push ranks deeper into the node, where they will utilise better cache access and this may improve runtime performance when the two options are used in tandem. If threads are set, the larger of the two values will be passed to depth, and the thread export preserved for threading purposes.

   * - **Auto-Populated**
     -
     -
     -
   * - optlevel
     - Str
     - NA
     - The optlevel being used. Populated from info in the EXTRAS part of the task name. A value must be present in the task name.
   * - precision
     - Dict
     - NA
     - Keys are the precision settings available, values the precision. Determined by the EXTRAS part of the task name. All will default to 64bit if not defined.
   * - build_conf
     - Str
     - NA
     - Set as ``<compiler>_<optlevel>_<precision>``
   * - restart
     - Logical String
     - NA
     - Controls whether the task is a restart.
   * - restart_no
     - Int
     - NA
     - The restart point for this task. Starts at 0.
   * - restart_read/write
     - Logical String
     - NA
     - Controls whether the task will read/write a restart dump.
   * - init_tstep
     - Int
     - NA
     - The initial time step for this task.
   * - last_step
     - Int
     - NA
     - The final time step for this task.
   * - prev_task
     - Str
     - NA
     - The name of the task for the previous crun point. Empty string if at 1st point.


Exceptions
----------

* Build Tasks - These follow the same format as above except for the CONFIG-NAME which is omitted as builds don't have a particular configuration. To add a build only task to a groups file, prepend with "build" such that the format is ``build_<APPLICATION>_<PLATFORM>_<COMPILER>_<EXTRAS>``.

* Script Tasks - Script tasks can be given any name, but for the suite to understand them they must be inherited by the ``scripts`` group. Examples of scripts would be ``check_style`` or ``check_config_dump``. Scripts must be defined manually as they don't follow a set template. The graph can be defined in ``rose-stem/templates/graph/populate_graph_scripts.cylc`` and the runtime section can then be defined in ``rose-stem/templates/runtime/generate_runtime_scripts.cylc``. Scripts are all run on the same platform set by the ``site_vars["script_platform"]`` variable defined in ``rose-stem/site/<SITE>/variables.cylc``. For meto, this is SPICE.

* ``lfricinputs`` tasks have different configuration options in their tasks files. See existing lfricinputs jobs for examples. They also contain their own templating system as the graph and runtime sections are both different.

* ``lfric_coupled`` tasks require additional definitions and some extra templating. These options can be seen in ``rose-stem/templates/default_task_definitions_lfric_coupled.cylc``.

* Unit/Integration Tests - These can be added by setting the CONFIG-NAME section in the task name to ``unit_tests`` or ``integration_tests``. These require a single entry to the task dictionary to work - see ``rose-stem/site/meto/gungho/tasks_gungho_xc40.cylc`` for examples. These don't require a resolution or optimisation setting.

* Canned Tests - These require ``app_name: "canned_test"`` to be set in the task dictionary. ``example_dir`` can also be set if not using the default location for the canned test files. The CONFIG-NAME can be set to anything, but ideally this should include "canned" for clarity (and only "canned" if there is just 1 example).

* Perftools can be run over certain applications in lfric_apps, see the compiler section above, and :ref:`perftools_rose-stem`

.. _EXlink: https://www-hpc/ex-user-guide/advice/placement/?h=l1+cache#nodes-diagram