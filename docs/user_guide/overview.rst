.. _overview:

================
Package overview
================

Here the general layout of the package, its key components and general usage
guidelines are presented.

The package consists of the following components:

    * The source code of the library, which implements the available
      functionality.
      If you are interested in contributing to the library, take a look at the
      :ref:`developer_guide`.
      The :ref:`code_reference` can be used as help when using the already
      existing modules and functions in your own code.

    * Executable python scripts that are used to drive the inflow field
      generation.
      Each generation method provided by the library
      (see :ref:`generation_methods`) has a script associated with it.
      The inflow generation methods commonly depend on a vast amount of
      parameters.
      These parameters are communicated via a configuration file, which is
      passed as a command-line argument. ::

        nameOfTheScript.py --config=configurationFileName

      The configuration file is a simple text file, with the following
      layout ::

        # This is a comment, it can describe the parameter below
        parameterOne    valueOfParameter1

        # Another comment
        parameterTwo    valueOfParameter2

      Which parameters should be present depends on the method that is being
      used.
      The parameters required by a certain method are described in the
      section of the user guide devoted to this method.

    * Utilities, which are also executable python scripts, that provide extra
      functionality, like computing statistical data or conversion between
      different formats.
      See :ref:`utilities` for more information.

    * Documentation, which includes this :ref:`user_guide`, the
      :ref:`developer_guide`, and the :ref:`code_reference`.

    * :ref:`tutorials` which provide examples on how the package can be used.

