# Updating the Science Guide
This README file describes how to update the Science Guide documentation for 
LFRic Apps.

## Adding new sections
To add a new section to the Science Guide, follow these steps:
1. Create a new directory under `documentation/source/science_guide/` with a 
   descriptive name for the section.
2. Inside this new directory, create an `index.rst` file. This file will 
   automatically be added to the table of contents.
3. Include a reference at the top of the `index.rst` file below the copyright 
   notice of the form ```.. _<section>_index:``` to allow cross-referencing.
4. If subsections are needed, create additional `.rst` files in the same 
   directory. Including the following table of contents tree in the sections 
   `index.rst` will pick up additional files and include them in the contents.
``` restructuredtext
.. toctree::
    :maxdepth: 1
    :glob:

    *
```
5. Follow the documentation style guide described in 
   [Sphinx Documentation Formatting](https://metoffice.github.io/lfric_core/how_to_contribute/style_guides/documentation_style_guide.html#text-formatting)

