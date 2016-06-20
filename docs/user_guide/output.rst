.. _output_formats:

Supported output formats
========================

Eddylicious supports several file formats for outputting the generated
inflow fields.
The choice of the file format is usually dictated by the CFD solver.

.. _hdf5_file_format:

The HDF5 file format
--------------------

HDF5 is a file format specifically developed for storing large scientific
datasets.
More details regarding HDF5 can be found `here <https://www.hdfgroup.org/>`_.

In HDF5 data can be sorted into groups.
The data itself is stored in the form of datasets.
One can think of datasets as of multidimensional arrays.

Let :math:`N_p` be the number of points at which the inflow fields are
generated and :math:`N_t` the amount of time-values for which the inflow fields
are generated.
Eddylicious creates the following datasets inside the HDF5 file.

   * ``points``, :math:`N_p \times 3` --- the points associated with the
     values of the inflow fields.
     The three columns represent the :math:`x`, :math:`y`, and :math:`z`
     coordinates respectively.

   * ``times``, :math:`N_t \times 1` --- the time values associated with the
     inflow fields.

   * ``velocity``, :math:`N_t \times N_p \times 3` --- the values of the
     velocity field.
     The first index is associated with time, the second with the available
     points and the third one with the components of the velocity field.
     Same order as in ``points`` applies.

The following parameters need to be provided in the configuration file in
order to output the velocity fields in the HDF5 file format. ::

   writer          hdf5
   writePath       path to the where the database will be stored
   hdf5FileName    name of the hdf5 file

.. _of_native_format:

OpenFOAM native format
----------------------

This is natively supported by OpenFOAM.
In order to read in boundary data from the hard-drive OpenFOAM has a special
boundary conidtion called ``timeVaryingMappedFixedValue``.
This boundary condition expects a folder called ``boundaryData/\<patchname\>``
to be located in the ``constant`` directory of the case.
Inside the folder a file named ``points`` should reside.
This file provides a list of the points where the boundary data is available.
The boundary data itself resides in folders named as the time-value associated
with the data.
The data for each available field is stored in its own file named identically
to the internal name of the field in OpenFOAM (for instance ``U`` for the
velocity field).
The format of each such file is quite similar to the :ref:`foamfile_format`
but has some additional headers.

The following parameters need to be provided in the configuration file in
order to output the velocity fields in the OpenFOAM native format ::

   writer          ofnative
   writePath       /path/to/OpenFOAM/case
   inletPatchName  name of the inlet patch
