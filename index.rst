.. KEPLER Utilities documentation master file, created by
   sphinx-quickstart on Thu Jun  4 16:14:35 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

KEPLER Utilities
================
This suite of analysis routines is designed to provide a useful python interface to the 1D stellar evolution code KEPLER. Each of the major modules in this suite are designed with a particular emphasis in mind but are for the most part independent.

The :mod:`jobs <kepler_utils.jobs>` module provides an interface for generating and running KEPLER jobs either locally or remotely using celery.

The :mod:`records <kepler_utils.records>` module reads the dump and cnv output files from KEPLER. It uses the :mod:`astropy.units` module to maintain consistent units.

The :mod:`database <kepler_utils.database>` module provides a tool for creating and maintaining a database of the output files that KEPLER generates with sqlalchemy. It takes advantage of the file readers in the :mod:`records` module.

The :mod:`yields <kepler_utils.yields>` module is used to interface with KEPLER yield files. It contains a series of tools to read abundances and integrate them over a stellar population.

The :mod:`plots <kepler_utils.plots>` module provides a series of visualization tools for the various modules in the suite.

For some (hopefully) working examples, try out some of the python scripts in the "scripts" directory. Most of these should be somewhat self explanatory and provide examples on how to use the tools of this suite.

Contents:

.. toctree::
   :maxdepth: 2

   Records <kepler_utils/records/doc>
   Yields <kepler_utils/yields/doc>

..

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

