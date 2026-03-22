.. -----------------------------------------------------------------------------
    (c) Crown copyright 2025 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

:html_theme.sidebar_secondary.remove: true

==========
LFRic Apps
==========

The `LFRic Apps <lfric_apps_github_>`__ repository is home to the documentation
for `LFRic Core <lfric_core_github_>`__ based science applications such as
lfric_atm and science code such as the GungHo dynamical core.

The `LFRic Core <lfric_core_github_>`__ is developed in a separate repository
and is designed to be a general modelling infrastructure for earth system
modelling. This code base contains most of the core infrastructure used by 
the applications found in the LFRic Apps repository and is documented 
`here <lfric_core_docs_>`__.

.. grid:: 2 2 4 4
    :gutter: 2

    .. grid-item-card::
        :text-align: center

        Information on getting going, from software stacks to testing

        +++
        .. button-ref:: getting_started_index
            :ref-type: ref
            :color: primary
            :outline:
            :expand:

                Getting Started

    .. grid-item-card::
        :text-align: center

        Guide on the code for users and developers

        +++
        .. button-ref:: user_guide_index
            :ref-type: ref
            :color: primary
            :outline:
            :expand:

                User Guide

    .. grid-item-card::
        :text-align: center

        Guide to the science included in the models

        +++
        .. button-ref:: science_guide_index
            :ref-type: ref
            :color: primary
            :outline:
            :expand:

                Science Guide

    .. grid-item-card::
        :text-align: center

        Guide on details of the development process

        +++
        .. button-ref:: developer_guide_index
            :ref-type: ref
            :color: primary
            :outline:
            :expand:

                Developer Guide

Development of the new LFRic atmosphere model is being done within the
|Momentum| Partnership.
The LFRic atmosphere application is accessible to Met Office partners.
Key initial aims for the |Momentum| atmosphere model are as follows:

- The model will be scientifically as good as the UM atmosphere.
- The model will scale better on future exascale platforms.
- The infrastructure will be flexible enough to support future
  evolutions of the science.

.. grid:: 2 2 4 4
    :gutter: 3 3 4 4

    .. grid-item-card::
        :text-align: center

        :material-round:`menu_book;2em`

        +++
        .. button-ref:: glossary_of_terms
            :ref-type: ref
            :color: primary
            :outline:
            :expand:

               Glossary

    .. grid-item-card::
        :text-align: center

        :material-round:`help_center;2em`

        +++
        .. button-ref:: faqs
            :ref-type: ref
            :color: primary
            :outline:
            :expand:

                FAQs

    .. grid-item-card::
        :text-align: center

        :fab:`github;fa-xl`

        +++
        .. button-link:: https://github.com/MetOffice/lfric_apps
            :ref-type: ref
            :color: primary
            :outline:
            :expand:

                GitHub

    .. grid-item-card::
        :text-align: center

        :far:`comments;fa-xl`

        +++
        .. button-link:: https://github.com/MetOffice/simulation-systems/discussions
            :ref-type: ref
            :color: primary
            :outline:
            :expand:

                Discussions


.. toctree::
    :maxdepth: 1
    :hidden:

    getting_started/index
    user_guide/index
    science_guide/index
    developer_guide/index
