.. _installing:

======================
Installing the package
======================


Dependencies
------------
Relatively new versions of both Python 3 should work.
Python 2 is no longer supported.

Some MPI compatable with ``mpi4py`` should be installed on your system.

The package depends on the following python packages:

   * ``numpy``

   * ``scipy``

   * ``matplotlib``

   * ``mock``

   * ``Sphinx`` and ``sphinxcontrib-bibtex``

   * ``h5py`` with support for MPI I/O.

   * ``mpi4py``

Make sure these are available before you install eddylicious.


Installing from git
-------------------

Clone the following repository to a location of your choice using ``git``
https://github.com/timofeymukha/eddylicious

The catalog ``eddylicious`` is created.
Go inside the catalog.
Run ``python setup.py install`` to install the package.
As usual, run with ``sudo`` if that is necessary.
The ``--prefix`` flag can be used to install to a custom directory.

If you intend to try latest updates immediately run
``python setup.py develop`` instead.
This way, if you update the files by running ``git pull`` you don't have to
reinstall the package.

