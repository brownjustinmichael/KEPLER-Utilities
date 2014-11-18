import math

import numpy as np
import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
import sqlalchemy

import records.cnv
import database.database as db
import database.cache
import plots.kipp
import plots.abundances

def numToSize (num):
    return 100 * (-np.log2 (num) + 3.5)

session = db.Session ()

query = db.basicQuery (session).filter (db.SimulationEntry.tags.contains (db.Tag.get (session, "OS/SC Grid"))).filter (db.SimulationEntry.binm10 > 16.0)
query = query.filter (db.DumpFileEntry.brumoson > 0.).filter (db.DumpFileEntry.woodscon > 0.)
query = query.filter (db.DumpFileEntry.scpower >= 0.9).filter (db.DumpFileEntry.scpower <= 1.1)
query = query.filter (db.DumpFileEntry.state == 'presn').filter (db.SimulationEntry.cnvfiles.any ())

allquery = db.basicQuery (session).filter (db.SimulationEntry.tags.contains (db.Tag.get (session, "OS/SC Grid"))).filter (db.SimulationEntry.binm10 > 16.0)
allquery = allquery.filter (db.DumpFileEntry.brumoson > 0.).filter (db.DumpFileEntry.woodscon > 0.)
allquery = allquery.filter (db.DumpFileEntry.osfactor > 0.0).filter (db.DumpFileEntry.scpower > 0.0)
allquery = allquery.filter (db.DumpFileEntry.state == 'presn').filter (db.SimulationEntry.cnvfiles.any ())

entries = [entry for sim, entry in query.all ()]
sims = [sim for sim, entry in query.all ()]
cnvs = [sim.cnvfiles [0] for sim in sims]
allentries = [entry for sim, entry in allquery.all ()]
allsims = [sim for sim, entry in allquery.all ()]

if len (entries) == 0:
    raise ValueError ("No entries in query")
    
hecores = u.Quantity ([entry.cache (session, 'he_core', database.cache.calculate_he_core) for entry in entries])
earlyhecores = u.Quantity ([sim.get_state_dump ("hdep").cache (session, 'he_core', database.cache.calculate_he_core) for sim in sims])
cocores = u.Quantity ([entry.cache (session, 'co_core', database.cache.calculate_co_core) for entry in entries])

allhecores = u.Quantity ([entry.cache (session, 'he_core', database.cache.calculate_he_core) for entry in allentries])
allearlyhecores = u.Quantity ([sim.get_state_dump ("hdep").cache (session, 'he_core', database.cache.calculate_he_core) for sim in allsims])
allcocores = u.Quantity ([entry.cache (session, 'co_core', database.cache.calculate_co_core) for entry in allentries])

tcores = np.genfromtxt ("tcores.dat", names = True)
lines = np.zeros (len (entries))

osfactors = np.array ([entry.osfactor if (entry.brumoson > 0.0) else 1 for entry in entries])
scpowers = np.array ([entry.scpower for entry in entries])

fig, axes = plt.subplots (1, 1, sharex = True)

scs = []

# Plot 1
ax = axes
ax.set_ylabel ("C/O Core Mass (solar masses)")
scs.append (ax.scatter (allhecores.to (u.solMass), allcocores.to (u.solMass), 15.0, c = 'black', alpha = 0.5))
scs.append (ax.scatter (hecores.to (u.solMass), cocores.to (u.solMass), 70.0, c = osfactors, picker = True, label = "Brown & Woosley (in prep)"))
# scs.append (ax.scatter (tcores ["HeCore"], tcores ["COCore"], 70.0, marker = "s", c = tcores ["qr"], cmap = plt.get_cmap ("gray"), norm = matplotlib.colors.LogNorm (vmin = 0.001, vmax = 0.1), label = "Sukhbold & Woosley (2013)"))
ylim = ax.get_ylim ()
xlim = ax.get_xlim ()
x = np.arange (5.6, 8.5, 0.1)
y = x - 1.7
ax.plot (x, y, "--", c = 'black')
ax.set_ylim (ylim)
ax.set_xlim (xlim)

# # Plot 2
# ax = axes [1]
# ax.set_ylabel ("Early He")
# scs.append (ax.scatter (hecores.to (u.solMass), earlyhecores.to (u.solMass), c = scpowers, linewidths = lines, picker = True, alpha = 0.75))
# ylim = ax.get_ylim ()
# xlim = ax.get_xlim ()
# ax.plot ((0,10), (0,10), color = 'black')
# ax.set_ylim (ylim)
# ax.set_xlim (xlim)

axes.set_xlabel ("He Core Mass (solar masses)")

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
            newfig = plots.abundances.jTDPlot (entry.file)
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
        
        newfig, ax = plots.kipp.jTDPlot (cnvs [i].get_data ())
        
        ax.set_title (cnvs [i].name)
        
        figs [newfig.canvas] = i
        newfig.canvas.mpl_connect('close_event', onclose)
        newfig.canvas.mpl_connect('button_press_event', onclick)

    plt.show ()

fig.subplots_adjust (right = 0.80, hspace = 0)
cbar_ax = fig.add_axes ([0.85, 0.15, 0.05, 0.7])
cb = fig.colorbar (scs [1], cax = cbar_ax)
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
axes.legend (loc = 'upper left')
fig.canvas.mpl_connect('pick_event', onpick)

plt.show ()