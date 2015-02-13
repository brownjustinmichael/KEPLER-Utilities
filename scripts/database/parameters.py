#!/usr/bin/env python

from kepler_utils.database.database import Session, DumpFileEntry, SimulationEntry, Tag
import matplotlib.pyplot as plt
import matplotlib
import sqlalchemy
import numpy

session = Session ()

query = session.query (DumpFileEntry).filter (DumpFileEntry.binm10 > 21.1).join (SimulationEntry)
query = query.filter (DumpFileEntry.brumoson > 0.0).filter (DumpFileEntry.woodscon > 0.0)
query = query.filter (DumpFileEntry.scpower > 0.0)
query = query.filter (DumpFileEntry.state == 'presn')

osfactors = numpy.array ([entry.osfactor for entry in query.all ()])
scpowers = numpy.array ([entry.scpower for entry in query.all ()])

fig, axes = plt.subplots (1, 1, figsize = (18, 10))

axes.scatter (scpowers, osfactors, alpha = 0.75)

axes.set_xlabel ("SC Powers")
axes.set_ylabel ("OS Factors")

axes.set_xscale ('log')
axes.set_yscale ('log')
axes.set_ylim ((min (osfactors) / 1.2, max (osfactors) * 1.2))

plt.show ()