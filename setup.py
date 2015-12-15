from setuptools import setup

setup(name='eddylicious',
      version='0.1',
      description='A package for generating inflow velocity fields for Large Eddy Simulation',
      url='http://://bitbucket.org/lesituu/eddylicious',
      author='Timofey Mukha',
      author_email='TimofeyMukha@it.uu.se',
      packages=['eddylicious'],
      scripts=['bin/convertFoamFileToNpy.py',
               'bin/convertNpyToHDF5.py',
               'bin/runLundRescaling.py'],
      install_requires=[
                    'numpy',
                       ],
      zip_safe=False)

