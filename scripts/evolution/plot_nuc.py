#!/usr/bin/env python

import sys
import functools

import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np

import database.database as db
import database.cache as cch

if (len (sys.argv) < 2):
    raise RuntimeError ("Usage: python plot_nuc.py 'iso'")

isos = [i for i in sys.argv [1:]]

session = db.Session ()

query = db.basicQuery (session).filter (db.SimulationEntry.tags.contains (db.Tag.get (session, "OS/SC Grid"))).filter (db.SimulationEntry.tags.contains (db.Tag.get (session, "Stabilized")))
query = query.filter (db.DumpFileEntry.brumoson > 0.).filter (db.DumpFileEntry.woodscon > 0.).filter (db.DumpFileEntry.binm10 < 16.0)
query = query.filter (db.DumpFileEntry.state == 'presn').filter (db.SimulationEntry.cnvfiles.any ())

entries = [entry for sim, entry in query.all ()]

hecores = u.Quantity ([entry.cache (session, 'he_core', cch.calculate_he_core) for entry in entries])

osfactors = np.array ([entry.osfactor if (entry.brumoson > 0.0) else 1 for entry in entries])
scpowers = np.log (np.array ([entry.scpower for entry in entries]))

fig, axes = plt.subplots (len (isos), 1, sharex = True)

for ax, iso in zip (axes, isos):
    ax.set_yscale ("log")
    smass = u.Quantity ([entry.cache (session, iso + '_mass', functools.partial (cch.sum_iso, iso), loadBurn = True) for entry in entries])
    omass = u.Quantity ([entry.cache (session, 'ox16' + '_mass', functools.partial (cch.sum_iso, 'o16'), loadBurn = True) for entry in entries])
    ax.scatter (hecores, smass.to (u.solMass) / omass.to (u.solMass), scpowers * 60, c = osfactors)
    ax.set_ylabel (iso)

# ax.set_ylim ()

plt.show ()