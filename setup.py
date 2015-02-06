from setuptools import setup
import os

setup(name='kepler_utils',
      version='0.1',
      description='A set of utilities for analyzing and running the stellar evolution code KEPLER',
      url='https://github.com/brownjustinmichael/KEPLER-Utilities',
      author='Justin Brown',
      author_email='jumbrown@ucsc.edu',
      license='MIT',
      packages=['kepler_utils'],
      install_requires=["matplotlib","sqlalchemy","numpy","astropy","periodictable","fortranfile","pandas","celery"],
      scripts=[os.path.join (root, file) for root, subdirs, files in os.walk ("scripts") for file in files],
      zip_safe=False)