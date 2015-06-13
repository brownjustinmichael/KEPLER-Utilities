:mod:`records` Module
*********************

.. module:: kepler_utils.records

This module contains the means to read output data from a KEPLER run. The basic use of this module is in the indexing methods of the :class:`.dump.DataDump` and the :class:`.cnv.CNVFile` classes. We outline their construction and use here.

To access the attributes of a KEPLER dump file, construct a :class:`.dump.DataDump` object.::

	# Load the object from the file
	dump = DataDump ("location/of/file#state")

	# Index the object for the mass coordinate
	dump ["mass coordinate"]

	# Retrieve the ncyc parameter
	dump.parameters ["ncyc"]

	# Get the hydrogen, helium, and oxygen abundances
	dump [["h1", "he4", "o16"]]

To access the attributes of a KEPLER CNV file, construct a :class:`.cnv.CNVFile` object.::

	# Load the object from the file
	cnv = CNVFile ("location/of/file.cnv")

	# Index the object for the central temperature evolution
	cnv ["tc"]

	# Retrieve the radius at the million year mark
	cnv [cnv.model_near (u.Quantity (1.0e6, u.yr))] ["rncoord"] [-1]

These can be combined with the :class:`kepler_utils.plots.abudnaces.AbundancePlot` and :class:`kepler_utils.plots.kipp.KippenhahnPlot`, respectively, for visualization.

.. toctree::
   :maxdepth: 2

   parameters
   qparameters
   cnvparameters

:mod:`records.dump` Module
--------------------------
.. automodule:: kepler_utils.records.dump
	:members:
	:special-members:
	:show-inheritance:

:mod:`records.cnv` Module
-------------------------
.. automodule:: kepler_utils.records.cnv
	:members:
	:special-members:
	:show-inheritance:
