#!/usr/bin/env python

import database.database
import matplotlib.pyplot as plt
import matplotlib
import sqlalchemy
import numpy

session = database.database.Session ()

query = session.query (database.database.DumpFileEntry).filter (database.database.DumpFileEntry.binm10 > 14.9).filter (database.database.DumpFileEntry.binm10 < 15.1)
query = session.query (database.database.DumpFileEntry).filter (database.database.DumpFileEntry.brumoson > 0.0).filter (database.database.DumpFileEntry.woodscon > 0.0)
query = session.query (database.database.DumpFileEntry).filter (database.database.DumpFileEntry.scpower > 0.0)
query = query.filter (database.database.DumpFileEntry.state == 'presn')

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