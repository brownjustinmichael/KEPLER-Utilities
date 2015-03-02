#!/usr/bin/env python

import math

import numpy as np
import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
import matplotlib.cm as cm
import sqlalchemy

from kepler_utils.records.dump import DataDump
import kepler_utils.database.database as db
import kepler_utils.database.cache as cache
import kepler_utils.plots.kipp
from kepler_utils.plots.abundances import jTDPlot

def numToSize (num):
    return 70 * (-np.log2 (num) + 3)

session = db.Session ()

query = db.basicQuery (session).filter (db.SimulationEntry.tags.contains (db.Tag.get (session, "OS/SC Grid")))
query = query.filter (db.SimulationEntry.tags.contains (db.Tag.get (session, "Low dtcp")))
# query = query.filter (~db.SimulationEntry.tags.contains (db.Tag.get (session, "Blue Loop")))
query = query.filter (db.DumpFileEntry.brumoson > 0.).filter (db.DumpFileEntry.woodscon > 0.).filter (db.DumpFileEntry.binm10 < 19.0)#.filter (db.DumpFileEntry.binm10 > 19.0)
query = query.filter (db.DumpFileEntry.state == 'presn').filter (db.SimulationEntry.cnvfiles.any ())

entries = [entry for sim, entry in query.all ()]
sims = [sim for sim, entry in query.all ()]
cnvs = [sim.cnvfiles [0] for sim in sims]

if len (entries) == 0:
    raise ValueError ("No entries in query")

hecores = u.Quantity ([entry.cache (session, 'he_core', cache.calculate_he_core) for entry in entries])
prehecores = u.Quantity ([sim.getStateDump ("hdep").cache (session, 'pre_he_core', cache.calculate_he_core) for sim in sims])
earlyhecores = u.Quantity ([sim.getStateDump ("heign").cache (session, 'early_he_core', cache.calculate_he_core) for sim in sims])
prehcontain = u.Quantity ([sim.getStateDump ("hdep").cache (session, 'early_h_contained', cache.calculate_h_contained) for sim in sims])
earlyhcontain = u.Quantity ([sim.getStateDump ("heign").cache (session, 'early_h_contained', cache.calculate_h_contained) for sim in sims])
cocores = u.Quantity ([entry.cache (session, 'co_core', cache.calculate_co_core) for entry in entries])
coratio = u.Quantity ([entry.cache (session, 'co_ratio', cache.calculate_co_ratio_in_core) for entry in entries])
compactness = u.Quantity ([entry.cache (session, 'compactness', cache.calculate_compactness_2_5) for entry in entries])
fecores = u.Quantity ([entry.cache (session, 'fe_core', cache.calculate_fe_core) for entry in entries])

lines = [0 if db.Tag.get (session, "Low dtcp") in sim.tags else 3 for sim in sims]

osfactors = np.array ([entry.osfactor if (entry.brumoson > 0.0) else 1 for entry in entries])
scpowers = np.array ([entry.scpower for entry in entries])

fig, axes = plt.subplots (2, 2, sharex = False, figsize = (18, 10))

scs = []

x1 = hecores.to (u.solMass)
x2 = hecores.to (u.solMass)
colors = osfactors
sizes = numToSize (scpowers)

# Plot 1
ax = axes [0] [0]
ax.set_ylabel ("C/O")
scs.append (ax.scatter (x1, cocores.to (u.solMass), c = colors, s = sizes, linewidths = lines, picker = True, alpha = 0.75))

# Plot 2
ax = axes [0] [1]
ax.set_ylabel ("C/O in C/O Core")
scs.append (ax.scatter (x2, coratio, c = colors, s = sizes, linewidths = lines, picker = True, alpha = 0.75))

# Plot 3
ax = axes [1] [0]
ax.set_ylabel ("Compactness")

scs.append (ax.scatter (x1, compactness, c = colors, s = sizes, linewidths = lines, picker = True, alpha = 0.75))

ax.set_xlabel ("Core Mass")

# Plot 4
ax = axes [1] [1]
ax.set_ylabel ("Fe Core")
ax.plot (x2, x2)
scs.append (ax.scatter (x2, fecores.to (u.solMass), c = colors, s = sizes, linewidths = lines, picker = True, alpha = 0.75))
# ax.set_yscale ("log")

ax.set_xlabel ("Core Mass")

figs = {}
cnv_lines = {}

def onabunclose (event):
    line = cnv_lines.pop (event.canvas)
    l = line.pop (0)
    fig = l.get_figure ()
    l.remove ()
    fig.canvas.draw ()

def onclose (event):
    plt.figure (fig.number)
    i = figs.pop (event.canvas)
    lines [i] = 0
    for sc in scs:
        sc.set_lw (lines)
    plt.draw ()

def onclick (event):
    time = event.xdata
    i = figs [event.canvas]

    if time is None:
        return

    for entry in session.query (db.DumpFileEntry).filter (db.DumpFileEntry.simulation == sims [i]).order_by (db.DumpFileEntry.toffset, db.DumpFileEntry.time).all ():
        if u.Quantity (entry.time + entry.toffset, 's') > u.Quantity (time, 'year'):
            newfig, ax = jTDPlot (entry.file)
            record = DataDump (entry.file, False)
            ax.plot (record ["mass coordinate"].to (u.solMass) [1:], (np.log (record ["tn"] [1:].value) - np.log (record ["tn"] [:-1].value)) / (np.log (record ["pn"] [1:].value) - np.log (record ["pn"] [:-1].value)), label = "Grad")
            ax.plot (record ["mass coordinate"].to (u.solMass), (record ["xln"] / 1.e6 / u.solLum).to (1), label = "L/Lsun")
            ax.plot (record ["mass coordinate"].to (u.solMass), (record ["dn"] / 1.e3), label = "D")
            ax.plot (record ["mass coordinate"].to (u.solMass), (record ["pn"] / 1.e19), label = "P")
            ax.plot (record ["mass coordinate"].to (u.solMass), (record ["tn"] / 1.e9), label = "T")
            ax.plot (record ["mass coordinate"].to (u.solMass), (record ["rn"] / 1.e9), label = "r")
            ax.plot (record ["mass coordinate"].to (u.solMass) [:-1], -np.diff (record ["tn"]) / np.diff (record ["rn"]), label = "dT")
            ax.plot (record ["mass coordinate"].to (u.solMass) [:-1], -np.diff (record ["pn"]) / np.diff (record ["rn"]) / 1.e9, label = "dp")
            # ax.plot (record ["mass coordinate"].to (u.solMass), (record ["xkn"]), label = "O")
            ax.legend (bbox_to_anchor=(1.05, 1.2), loc=2, borderaxespad=1)
            newfig.canvas.mpl_connect('close_event', onabunclose)
            cnv_lines [newfig.canvas] = event.inaxes.plot (u.Quantity ([entry.time + entry.toffset, entry.time + entry.toffset], 's').to (u.year), event.inaxes.get_ylim (), color = 'black', lw = 5, alpha = 0.5)
            event.canvas.draw ()
            break

    plt.show ()

def onpick (event):
    ind = event.ind

    for i in ind:
        # plt.figure (fig.number)
        lines [i] = 3
        event.artist.set_lw (lines)
        # plt.draw ()
        event.canvas.draw ()

        # extent = range (cnvs [i].simulation.getStateDump ("hdep").ncyc, cnvs [i].simulation.getStateDump ("heign").ncyc)

        newfig, ax = kepler_utils.plots.kipp.jTDPlot (cnvs [i].get_data (), logspace = True)#, extent = extent)

        ax.set_title (cnvs [i].name)

        figs [newfig.canvas] = i
        newfig.canvas.mpl_connect('close_event', onclose)
        newfig.canvas.mpl_connect('button_press_event', onclick)

    plt.show ()

fig.subplots_adjust (right = 0.8, hspace = 0)
cbar_ax = fig.add_axes ([0.85, 0.15, 0.05, 0.7])
cb = fig.colorbar (scs [0], cax = cbar_ax)
cb.set_label ("Overshoot Factor")

# ls = []
# powers = []
# labels = []
# for i in scpowers:
#     if i not in powers:
#         powers.append (i)
# powers.sort ()
# for i in powers:
#     ls.append (plt.scatter ([], [], s = numToSize (i), linewidth = 0))
#     labels.append (str (i))
# axes [1].legend (ls, labels, scatterpoints = 1, title = "SC Power", loc = 'lower right')
# ax.set_yscale ('log')
fig.canvas.mpl_connect('pick_event', onpick)

print ("Showing...")

plt.show ()