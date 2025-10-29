.. -----------------------------------------------------------------------------
     (c) Crown copyright Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _lfric_atm_checkpoint:

LFRic atmosphere checkpoint/restart system
==========================================

The LFRic atmosphere ``lfric_atm`` application can be configured to generate a
checkpoint dump at the end of each model run. The checkpoint dump can be read in
by a new integration of the model allowing further timesteps to be run. The dump
is written using XIOS.

Requesting checkpoint restart
-----------------------------

Set ``checkpoint_write=.true.`` in the ``io`` namelist of the model
configuration to generate a checkpoint dump at the end of a model run. The
checkpoint dump will be named after the ``checkpoint_stem_name`` string in the
``files`` namelist appended with the number of the last timestep of the run.

Set ``checkpoint_read=.true.`` in the ``io`` namelist to restart a run from an
existing checkpoint dump. The expected start timestep will be defined by
``timestep_start`` in the ``time`` namelist. There must exist a checkpoint file
named according to the ``checkpoint_stem_name``, with the correct timestep
(which would be the ``timestep_start`` setting minus one) otherwise the run will
fail.

Some key ``lfric_atm`` configurations are regularly tested to check that a long
single run produces exactly the same results as a short initial run plus a
restarted run of the same total length. Maintaining such equivalence is
considered to be important by scientist users.

As some configurations of ``lfric_atm`` do not run all the same science every
timestep, runs that use checkpoint and restart may need to restart at a
particular frequency. In particular, aligning restart times with radiation
timesteps is important.


The checkpoint dump name and format
-----------------------------------

The name of the checkpoint dump is computed while the XIOS context is being
defined. Its name is computed from a configurable input string (that can define
a full path name plus file name prefix, or just a file name prefix in the
current working directory), and a suffix comprising the timestep number with
leading zeros allowing timesteps up to 9,999,999,999. For example, a checkpoint
written after completion of 4 timesteps may be named ``checkpoint_0000000004``.

Field data can be written to checkpoint files in one of two formats:

* "Legacy" or "stream" format. The field data in this form are written as a 1D
  stream of data, representing the whole 3D field, just as the data is held in
  memory in an LFRic field object.
* "Non-legacy" or "layer" format. The data in this form are written out much
  like it is in diagnostic output: as a series of 2D horizontal layers that are
  stacked to form the full 3D field.

In a single checkpoint file, some fields will be written in stream format, and
some will be written in layers. There are some restrictions to how fields can
be checkpointed:

1. Fields that have data points on both full- and half-levels (so anything on a
   W1 or W2 function space) can only be written directly in stream form. When
   the data is written like this, it cannot be successfully read if the parallel
   decomposition of the reading code is different between the writing and
   reading process. To ensure correct behaviour, such fields could be split into
   the full-level and half-level components which are written and read
   separately.
2. XIOS can only write discontinuous fields in layer form when the convention is
   set to "UGRID". Other field output can be written in the "CF" convention.
3. XIOS does not currently support the layer output of data on the edges of
   cells in discontinuous fields on a bi-periodic mesh. Apparently, it can't
   cope with the wrap-around cells (where the data points on the very right of
   the domain and those on the left form a single imaginary wrap-around cell).

Currently, ``lfric_atm`` outputs at least one W2 field which means
checkpoint-restart suffers the restriction described in the first point
above. Work is underway to split the field into W2H and W2V components to get
around the problem.

Eventually, all checkpointed fields should be output in layer form using the
UGRID convention - splitting W1 and W2 fields as needed. Work would still be
required to support configurations that use biperiodic meshes and need to
checkpoint such fields, though the requirement for checkpointing such
configurations is low priority.

Configuring the XIOS context
----------------------------

This section does not consider LBC prognostics or ancillary fields.

The term "prognostic fields" broadly refers to the set of model fields that must
remain in scope throughout a model timestep. The LFRic atmosphere model stores
the set of prognostic fields in a prognostic field collection held in the main
``modeldb`` data structure.

Not all prognostic fields need to be written to the checkpoint dump.
* Prognostic fields that may be written to the dump must be defined in the
  ``lfric_dictionary.xml`` file which is included in the ``iodef.xml`` file that
  configures XIOS.
* Prognostic fields that are not written to the dump are not defined in the
  ``lfric_dictionary.xml`` file

A module called ``create_physics_prognostics`` holds code that uses model
configuration to determine several things about the physics prognostics:

* Identifying which fields need to be written to the checkpoint dump at the end
  of the run.
* Identifying which of many possible prognostic fields are required for a given
  model configuration.
* Initialising those fields that are required.
* Adding initialised fields to the correct science-specific field collection.

Currently, the underpinning infrastructure cannot support doing all of these
items at the same time. A solution has been created that involves calling
``create_physics_prognostics`` twice. The procedure takes a procedure pointer as
an argument, enabling the routine to do different things on each call.

.. attention::

  There is also a ``create_gungho_prognostics`` file that manages the subset of
  prognostics used by Gungho model configurations. In part, the gungho
  prognostics are handled differently because the ``gungho_model`` application
  supports both lowest and higher order finite element order, and checkpointing
  of the latter requires a different file format. When written with XIOS, all
  gungho prognostics are output using the legacy format even though only one of
  them needs to be on a 3D mesh because it is a W2 field which cannot be split
  into levels.

Within the procedure, there is code like the following for each possible
prognostic field:

.. code-block:: fortran

  call processor%apply(make_spec('lw_down_surf_rtsi', main%radiation,    &
                       ckp=checkpoint_flag))
  call processor%apply(make_spec('qcl_at_inv_top', main%turbulence, W3,  &
                       twod=.true.))

.. attention::
   Note that in the calls above, the field that may be checkpointed does not
   provide the function space type or other information that describes the shape
   of the field, whereas the field that is not checkpointed provides the
   function space type ``W3`` and ``twod=.true`` to indicate that the field is
   on a single layer.

   Any field that `may` be checkpointed `must` have an XIOS definition (as
   pointed out above, these definitions are stored in the
   ``lfric_dictionary.xml`` file. It is possible to use this definition to infer
   the shape of the field (its function space and so forth). Therefore, there is
   no requirement to provide the argument that define these features. Indeed,
   doing so results in the definition appearing in two places and raises the
   possibility that the two definitions may diverge.

Broadly speaking, there will be logic defining whether the code is called for a
given field depending on whether the field is required; logic will set
``checkpoint_flag`` depending on whether the field needs to be checkpointed at
the end of the run (which depends on the scientific configuration, and also on
whether checkpointing has been enabled). The call also lists the field
collection to which the field will be added. In ``lfric_atm``, each field
collection is held in ``model_db`` so is accessible throughout the algorithm
layer of the model.

Key parts of the above call are:
  * ``make_spec`` returns a ``field spec`` data structure that summarises the
    requirements for the field. Optionally, this can include stating which field
    collection should store it, and what the shape of the field is in terms of
    its function space or whether it is a 2D field. Note that where the shape is
    not specified, the information will be derived from the field definition in
    the XIOS ``iodef.xml`` file.
  * ``apply`` is a procedure in the ``processor`` object that does work on the
    ``field spec``.
  * ``processor`` can be one of two different objects depending on whether this
    is the first or second call to the ``create_physics_prognostics`` code.


The first call
^^^^^^^^^^^^^^

The first call to the code in ``create_physics_prognostics`` is done while the
XIOS context is being defined. It must be done during this period because the
aim of the call is to provide a definition of the checkpoint fields to XIOS.

During this call, the ``processor`` object is the ``persistor_type`` defined in
``gungho_model_mod``. Therefore procedure ``gungho_model_mod:persistor_apply``
is applied to the ``field spec`` for each field.

If the field is to be written to or read from the checkpoint dump the ``apply``
method calls ``lfric_xios_metafile_mod:add_field``. The ``add_field`` procedure
registers an ``xios_field``, for writing to or reading from the checkpoint dump.
The ID is the name of the field prefixed with the hard-coded string literal,
``checkpoint_``. The procedure defines the checkpoint field for XIOS using XIOS
attributes copied from the definition of the original field ID (without the
``checkpoint_`` prefix) in ``lfric_dictionary.xml``.

If the field is not to be checkpointed, then nothing is done during the first
call to ``processor%apply``. In fact, the whole first call to the
``create_physics_prognostics`` code has no purpose and does nothing if
checkpoint reading or writing is not requested.

The second call
^^^^^^^^^^^^^^^
In the second call, the ``processor`` passed to ``create_prognostics_field`` is
the ``field_maker_type`` in ``field_maker_mod``.

Briefly, the ``apply`` method in ``field_maker_mod`` initialises the model field
based on the ``field_spec`` definition returned by the ``make_spec``
call. Information in the ``field_spec`` is drawn from the
``lfric_dictionary.xml`` file (if the field has a record in this file) and
information from the arguments to the ``make_spec`` call. For multidata fields
(often referred to as "tiled fields") some hardwired or configured settings are
obtained by the ``apply`` method by calling
``multidata_field_dimensions_mod:get_multidata_field_dimension``.

Problematically, the ``checkpoint_`` prefix is also used here (hard coded as a
string literal) to check if the field is valid before creating anything.

For fields that may be checkpointed, the information from
``lfric_dictionary.xml`` is extracted using the XIOS API. This can only be done
after the context definition has been closed. That is why this call is done
separately from the first call.


Simplified call tree for setting up I/O in LFRic_atm
----------------------------------------------------


.. code-block:: rst

  lfric_atm							(lfric_atm/lfric_atm.f90)
    │
    └──initialise						(gungho/driver/gungho_driver_mod.F90)
        │
        ├──initialise_infrastructure				(gungho/driver/gungho_model_mod.F90)
        │   │
        │   └──init_io						(components/driver/driver_io_mod.F90)
        │       │
        │       └──init_xios_io_context				(components/driver/driver_io_mod.F90)
        │           │
        │           ├──populate_filelist (=> init_gungho_files)	(gungho/driver/gungho_setup_io_mod.F90)
        │           │
        │           └──io_context%initialise_xios_context	(components/lfric_xios/lfric_xios_context_mod.f90)
        │               │
        │               ├──xios_context_initialise		xios
        │               │
        │               ├──xios_get_handle			xios
        │               │
        │               ├──xios_set_current_context		xios
        │               │
        │               ├──init_xios_calendar			(components/lfric_xios/lfric_xios_setup_mod.x90)
        │               │
        │               ├──init_xios_dimensions			(components/lfric_xios/lfric_xios_setup_mod.x90)
        │               │
        │               ├──setup_xios_files			(components/lfric_xios/lfric_xios_setup_mod.x90)
        │               │
        │               ├──before_close (=> before_context_close) 	(gungho/driver/gungho_model_mod.F90)
        │               │   │
        │               │   ├──persistor%init           		(gungho/driver/gungho_model_mod.F90)
        │               │   │
        │               │   ├──process_gungho_prognostics(persistor)
        │               │   │   │					(gungho/driver/create_gungho_prognostics_mod.F90)
        │               │   │   └──persistor%apply(makespec())		(gungho/driver/gungho_model_mod.F90)
        │               │   │       │
        │               │   │       └──add_field			(components/lfric-xios/lfric_xios_metafile_mod.F90)
        │               │   │           │
        │               │   │           ├──(various xios calls…) 	xios
        │               │   │           │
	│               │   │           └──handle_legacy_field 		(components/lfric-xios/lfric_xios_metafile_mod.F90)
        │               │   │
        │               │   └──process_physics_prognostics(persistor)
        │               │       │					(gungho/driver/create_physics_prognostics_mod.F90)
        │               │       └──(…)
        │               │
        │               └──xios_close_context_definitions 	xios
        │
        │
        └──create_model_data					(gungho/driver/gungho_init_fields_mod.X90)
            │
            ├──field_mapper%init				(gungho/driver/field_mapper_mod.F90)          │
            └──create_gungho_prognostics			(gungho/driver/create_gungho_prognostics_mod.F90)
                │
                ├──creator%init					(gungho/driver/field_maker_mod.F90)
                │
                └──process_gungho_prognostics(creator)		(gungho/driver/create_gungho_prognostics_mod.F90)
                    │
                    └──creator%apply(makespec())		(gungho/driver/field_maker_mod.F90)
                        │
                        └──add_real_field			(gungho/driver/field_maker_mod.F90)
                            │
                            ├──field%initialise			(infrastructure/field/field_mod.t90)
                            │
                            ├──field%set_read/write_behaviour	(infrastructure/field/field_mod.t90)
                            │
                            ├──depository%add_field		(infrastructure/field/field_collection_mod.F90)
                            │
                            └──prognostic_fields%add_field	(infrastructure/field/field_collection_mod.F90)

(For brevity, the paths are shortened. Some prefix directories and the “source”
directory have been omitted.)
