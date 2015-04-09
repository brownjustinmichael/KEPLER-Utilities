#!/usr/bin/env python

from kepler_utils.database.database import Session, DumpFileEntry, SimulationEntry, Tag
import matplotlib.pyplot as plt
import matplotlib
import sqlalchemy
import numpy

session = Session ()

query = session.query (DumpFileEntry).filter (DumpFileEntry.binm10 < 19.0).join (SimulationEntry).filter (SimulationEntry.tags.contains (Tag.get (session, "OS/SC Grid")))
query = query.filter (SimulationEntry.tags.contains (Tag.get (session, "Low dtcp")))
query = query.filter (DumpFileEntry.state == 'presn')

entries = query.all ()

osfactors = numpy.array ([entry.osfactor for entry in entries])
scpowers = numpy.array ([entry.scpower for entry in entries])


tags = []

for entry in entries:
    found = False
    for tag in entry.simulation.tags:
        if len (tag.tag) == 1:
            tags.append (tag.tag)
            found = True
        if tag.tag == "Blue Loop":
            tags.append ("L")
            found = True
    if not found:
        tags.append (" ")

fig, axes = plt.subplots (1, 1, figsize = (18, 10))

axes.scatter (scpowers, osfactors, alpha = 0.75)

for os, sc, tag in zip (osfactors, scpowers, tags):
    axes.text (sc, os, tag)

axes.set_xlabel ("SC Powers")
axes.set_ylabel ("OS Factors")

axes.set_xscale ('log')
axes.set_yscale ('log')
axes.set_ylim ((min (osfactors) / 1.2, max (osfactors) * 1.2))

plt.show ()