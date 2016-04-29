Supported output formats
========================

Eddylicious supports several file formats that for outputting the generated
inflow fields.
This allows to use eddylicious with several CFD solvers.

HDF5 file format
----------------

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
      The first index is assiciated with time, the second with the available
      points and the third one with the components of the velocity field.
      Same order as in ``points`` applies.

