:mod:`yields` Module
********************

.. module:: kepler_utils.yields

This module provides an interface to the final nucleosynthesis produced by KEPLER. In general, the use of this module is to generate a :class:`.yields.YieldReader` object and an :class:`.integrator.Integrator` object. These two objects can then be used to generate averages of whatever yields (or other quantities) as desired.

An example of the typical use of this module:

Given a directory "yields/" containing several yield files and a single master file "yields/masses", which should look like

.. code-block:: python

    mass1 yield_file1
    mass2 yield_file2
    mass3 yield_file3
    ...

The relevant :class:`.yields.YieldReader` object and :class:`.integrator.Integrator` object can be constructed as follows, using the :meth:`.yields.YieldReader.from_directory` and :meth:`.integrator.IMFIntegrator.from_yieldreader` methods.

.. code-block:: python

	yr = YieldReader.from_directory ("yields/")
	imf = IMFIntegrator.from_yieldreader (yr)

Any desired integration can then be run with ease with :meth:`.integrator.Integrator.__call__` or :meth:`.integrator.Integrator.get_abundances`:

.. code-block:: python

	# Return the integrated yield of Fe56
	imf (yr.get_yield ("fe56"))

	# Return the average mass of the sample
	imf (yr.masses)

	# Return the total integrated abundance of the population
	imf.get_abundances (yr)

:mod:`yields.abundances` Module
-------------------------------
.. automodule:: kepler_utils.yields.abundances
	:members:
	:special-members:
	:show-inheritance:

:mod:`yields.integrator` Module
-------------------------------
.. automodule:: kepler_utils.yields.integrator
	:members:
	:special-members:
	:show-inheritance:

:mod:`yields.janka` Module
--------------------------
.. automodule:: kepler_utils.yields.janka
	:members:
	:special-members:
	:show-inheritance:

:mod:`yields.yields` Module
---------------------------
.. automodule:: kepler_utils.yields.yields
	:members:
	:special-members:
	:show-inheritance:
