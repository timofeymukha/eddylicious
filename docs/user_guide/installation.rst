.. _installing:

======================
Installing the package
======================


Dependencies
------------
Both python 2 and 3 should work, but the package has only been extensively
tested for python 2.

Some MPI compatable with ``mpi4py`` should be installed on your system.

The package depends on the following python packages:

   * ``numpy``

   * ``scipy``

   * ``matplotlib``

   * ``mock``

   * ``sphinxcontrib-bibtex``

   * ``h5py`` with support for MPI I/O.

   * ``mpi4py``

Using pip
---------

Run ``pip install eddylicious``.
This should do it!

Using git
---------

Clone the following repository to a location of your choice using ``git``
https://github.com/timofeymukha/eddylicious

The catalog ``eddylicious`` is created.
Go inside the catalog.
Run ``python setup.py install`` to install the package.
As usual, run with ``sudo`` if that is necessary.
The ``--prefix`` flag can be used to install to a non-default directory.

If you wish to get the latest updates immediately run
``python setup.py develop`` instead.
This way, if you update the files by running ``git pull`` you don't have to
reinstall the package.

