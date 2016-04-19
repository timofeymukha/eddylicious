from setuptools import setup

setup(name='eddylicious',
      version='0.0.1',
      description='A package for generating inflow fields for Large Eddy '
                  'Simulation',
      url='http://bitbucket.org/lesituu/eddylicious',
      download_url='https://bitbucket.org/lesituu/eddylicious/get/0.0.1.tar.gz',
      author='Timofey Mukha',
      author_email='timofey.mukha@it.uu.se',
      packages=['eddylicious'],
      scripts=['bin/convertFoamFileToHDF5.py',
               'bin/inflowStats.py',
               'bin/precursorStats.py',
               'bin/runLundRescaling.py',
               ],
      install_requires=[
                    'numpy',
                    'scipy',
                    'matplotlib',
                    'mock',
                    'sphinxcontrib-bibtex',
                       ],
      license="GNU GPL 3",
      classifiers=[
          "Development Status :: 2 - Pre-Alpha",
          "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"
      ],
      zip_safe=False)

