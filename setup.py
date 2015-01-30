from setuptools import setup

setup(name='kepler_utils',
      version='0.1',
      description='A set of utilities for analyzing and running the stellar evolution code KEPLER',
      url='https://github.com/brownjustinmichael/KEPLER-Utilities',
      author='Justin Brown',
      author_email='jumbrown@ucsc.edu',
      license='MIT',
      packages=['kepler_utils'],
      install_requires=["matplotlib","sqlalchemy","numpy","astropy","periodictable","fortranfile","pandas","celery"],
      zip_safe=False)