.. _utilities:

Utilities
=========

.. _inflowstats:

inflowStats
-----------

This utility allows to compute mean velocity and the diagonal components of the
Reynolds stress tensor of a database of inflow fields.
The database should be stored as an HDF5 file, see :ref:`hdf5_file_format` in
:ref:`output_formats`.

The utility accepts two command line arguments.

    * ``--database`` --- the HDF5 file containing the inflow fields.

    * ``--writePath`` --- the location where to write the files containing the
      computed results.

It is possible to run it in parallel using MPI.

The utility will create the following files in the location specified by
``writePath``

    * ``uMeanX``, ``uMeanY``, ``uMeanZ`` --- contain the corresponding
      component of the mean velocity field.

    * ``uPrime2MeanXX``, ``uPrime2MeanYY``, ``uPrime2MeanZZ`` --- contain the
      corresponding diagonal component of the Reynolds stress tensor.

    * ``y`` --- the locations of the datapoints.

.. _precursorstats:

precursorStats
--------------

This utility allows to compute mean velocity and the diagonal components of the
Reynolds stress tensor of a precursor database.
The database should be stored as an HDF5 file, see
:ref:`input_hdf5_file_format` in :ref:`input_formats`.

The utility accepts two command line arguments.

    * ``--database`` --- the HDF5 file containing the inflow fields.

    * ``--writePath`` --- the location where to write the files containing the
      computed results.

It is possible to run it in parallel using MPI.

The utility will create the following files in the location specified by
``writePath``

    * ``uMeanX``, ``uMeanY``, ``uMeanZ`` --- contain the corresponding
      component of the mean velocity field.

    * ``uPrime2MeanXX``, ``uPrime2MeanYY``, ``uPrime2MeanZZ`` --- contain the
      corresponding diagonal component of the Reynolds stress tensor.

    * ``y`` --- the locations of the datapoints.




.. _convertFoamFileToHDF5:

convertFoamFileToHDF5
---------------------

This utility converts a precursor database stored in the foamFile format
(see :ref:`foamfile_format`) to a database stored as a single HDF5 file,
see :ref:`hdf5_file_format`.

The utility accepts the following command line arguments.

    * ``--precursor`` --- The location of the precusor case.
      This path will be used in order to locate the folder with the sampled velocity values.

    * ``--surfaceName`` --- The name of the surface that was used for sampling
      the values.
      The name is chosen in the ``cotrolDict`` of the case.

    * ``--fileName`` --- the name of the HDF5 file that will be produced.

    * ``--uMeanFile`` --- The file containing the mean velocity profile.
      The file is assumed to have two columns, one with wall-normal coordinate,
      and the other one with the values of mean streamwise velocity.

It is possible to run the utility in parallel using MPI.