from setuptools import setup

setup(name='eddylicious',
      version='0.0.2',
      description='A package for generating inflow fields for LES and DNS',
      url='https://github.com/timofeymukha/eddylicious',
      download_url='https://github.com/timofeymukha/eddylicious/archive/0.0.2.tar.gz',
      author='Timofey Mukha',
      author_email='timofey.mukha@it.uu.se',
      packages=['eddylicious'],
      entry_points = {
          'console_scripts':['inflowStats=bin.inflowStats:main',
                             'precursorStats=bin.precursorStats:main',
                             'runLundRescaling=bin.runLundRescaling:main',
                             'convertFoamFileToHDF5=bin.convertFoamFileToHDF5:main']
      },
      install_requires=[
                    'numpy',
                    'scipy',
                    'matplotlib',
                    'mock',
                    'sphinxcontrib-bibtex',
                       ],
      license="GNU GPL 3",
      classifiers=[
          "Development Status :: 4 - Beta",
          "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"
      ],
      zip_safe=False)

