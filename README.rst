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

* HDF5. This is a file-format specifically designed to be a container for large datasets.
  More information can be found at the `offical homepage <https://www.hdfgroup.org/HDF5/>`_.
HDF5 allows to store large arrays of data in the form of datasets.

Three datasets are created by the writer.

:points: 
    A 2d array of size N-by-3, where N is the number of points. 
    The three colums are the x, y and z coordinates of the points, respectively.
    The points represent the location of the datapoints for the generated velocity fields in space.
    In the finite volume setting, it is best that this corresponds to the face centers of the faces defining the inflow boundary.
    
:time:
    A 1d array, the size corresponds to the number of time-values that the inflow boundary fields are generated for.
    Each element contains a time value.

:velocity:
   A 3d array. The first dimension corresponds to time and has the same size as the time array. 
   The second dimension is of size N and the third is of size 3, therefore for each timestep there is a N-by-3 array that contains the values of the velocity field.
   The order corresponds to the points array.

APIs for loading HDF5 datasets exist in most languages actively used in scientific computing.
For using this output format with OpenFOAM see `this repository <https://bitbucket.org/lesituu/timevaryingmappedhdf5fixedvalue>`_.

* TimeVaryingMappedFixedValue.

[lund] T. S. Lund, X. Wu, and K. D. Squires. On the Generation of Turbulent Infow Conditions for Boundary Layer Simulations. Journal of Computational Physics, 140:233-258, 1998.