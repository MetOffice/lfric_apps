.. -----------------------------------------------------------------------------
    (c) Crown copyright 2025 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

.. _repository_contents:

Contents of the Repository
--------------------------

The repository contains the following directories:

- The ``science`` directory contains several libraries of science code and
  interface code which allows calls to external models each of which can be used
  in the development of applications. For example, the **gungho**
  library provides the dynamical core used in the |Momentum| Atmosphere
  application.
- The ``applications`` directory contains the different applications that have
  been developed using the LFRic Core and the science libraries. Applications
  include the main |Momentum| Atmosphere application as well as smaller
  applications that test subsections of different models (e.g. transport,
  solver, etc.).
- The ``rose-stem`` directory contains the system test suites for development
  testing.
- The ``build`` directory contains additional build scripts required by the
  modelling system.

Many of the directories contain directories of unit and integration
tests, directories of Rose metadata, and Makefiles for building
applications.
