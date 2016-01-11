================
Eddylicious
================

Eddylicious is a python package, designed for implementing methods for generating inflow boundary fields for Large Eddy Simulation (LES).

The package is divided into three main parts: readers, writers and generators.
Generators are functions that actually implement the generation of the boundary fields.
Writers are functions that implement output of the generated fields into a specific format.
Readers are functions that implement the input of either geometric data (such as point or face-center locations of the inflow plane) and/or flow field data, which can be needed for methods using a precursor-database.

Additionally, conifgurable scripts that perform the generation and output are provided.

---------------
 Implemented Generation Methods
---------------

Currently, only one generation method is implemented in the package: the rescaling procedure presented in [lund]_.

---------------
 Implemented Output Formats
---------------
Two output formats are currently supported.

* HDF5. This is a format specifically designed to be a container for large datasets.
More information can be found at the `offical homepage <https://www.hdfgroup.org/HDF5/>`_

* TimeVaryingMappedFixedValue.

[lund] T. S. Lund, X. Wu, and K. D. Squires. On the Generation of Turbulent Infow Conditions for Boundary Layer Simulations. Journal of Computational Physics, 140:233-258, 1998.