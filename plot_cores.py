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

query = db.basicQuery (session).filter (db.SimulationEntry.tags.contains (db.Tag.get (session, "OS/SC Grid"))).filter (db.SimulationEntry.tags.contains (db.Tag.get (session, "Stabilized")))
query = query.filter (db.DumpFileEntry.brumoson > 0.).filter (db.DumpFileEntry.woodscon > 0.).filter (db.DumpFileEntry.binm10 < 16.0)
query = query.filter (db.DumpFileEntry.state == 'presn').filter (db.SimulationEntry.cnvfiles.any ())

entries = [entry for sim, entry in query.all ()]
sims = [sim for sim, entry in query.all ()]
cnvs = [sim.cnvfiles [0] for sim in sims]

if len (entries) == 0:
    raise ValueError ("No entries in query")
    
hecores = u.Quantity ([entry.cache (session, 'he_core', database.cache.calculate_he_core) for entry in entries])
earlyhecores = u.Quantity ([sim.getStateDump ("hdep").cache (session, 'he_core', database.cache.calculate_he_core) for sim in sims])
cocores = u.Quantity ([entry.cache (session, 'co_core', database.cache.calculate_co_core) for entry in entries])
sicores = u.Quantity ([entry.cache (session, 'si_core', database.cache.calculate_si_core) for entry in entries])
fecores = u.Quantity ([entry.cache (session, 'fe_core', database.cache.calculate_fe_core) for entry in entries])
tasbsg = u.Quantity ([cnv.cache (session, 'tasbsg', database.cache.calculate_tasbsg) for cnv in cnvs])

lines = np.zeros (len (entries))

for i in range (len (sims)):
    if db.Tag.get (session, "Stabilized") in sims [i].tags:
        lines [i] = 2.0

osfactors = np.array ([entry.osfactor if (entry.brumoson > 0.0) else 1 for entry in entries])
scpowers = np.array ([entry.scpower for entry in entries])

fig, axes = plt.subplots (2, 2, sharex = True, figsize = (18, 10))

scs = []

# Plot 1
ax = axes [0] [0]
ax.set_ylabel ("C/O")
scs.append (ax.scatter (hecores.to (u.solMass), cocores.to (u.solMass), numToSize (scpowers), c = osfactors, linewidths = lines, picker = True, alpha = 0.75))

# Plot 2
ax = axes [0] [1]
ax.set_ylabel ("Early He")
scs.append (ax.scatter (hecores.to (u.solMass), earlyhecores.to (u.solMass), numToSize (scpowers), c = osfactors, linewidths = lines, picker = True, alpha = 0.75))
ylim = ax.get_ylim ()
xlim = ax.get_xlim ()
ax.plot ((0,10), (0,10), color = 'black')
ax.set_ylim (ylim)
ax.set_xlim (xlim)

# Plot 3
ax = axes [1] [0]
ax.set_ylabel ("Si")
scs.append (ax.scatter (hecores.to (u.solMass), sicores.to (u.solMass), numToSize (scpowers), c = osfactors, linewidths = lines, picker = True, alpha = 0.75))

# Plot 4
ax = axes [1] [1]
ax.set_ylabel ("Time as BSG")
scs.append (ax.scatter (hecores.to (u.solMass), tasbsg.to (u.solMass), numToSize (scpowers), c = osfactors, linewidths = lines, picker = True, alpha = 0.75))

# axes.flat [1].set_xlabel ("He Core Mass")

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
        
        # extent = range (cnvs [i].simulation.getStateDump ("hdep").ncyc, cnvs [i].simulation.getStateDump ("heign").ncyc)
        
        newfig, ax = plots.kipp.jTDPlot (cnvs [i].get_data (), logspace = True)#, extent = extent)
        
        ax.set_title (cnvs [i].name)
        
        figs [newfig.canvas] = i
        newfig.canvas.mpl_connect('close_event', onclose)
        newfig.canvas.mpl_connect('button_press_event', onclick)

    plt.show ()

fig.subplots_adjust (right = 0.8, hspace = 0)
cbar_ax = fig.add_axes ([0.85, 0.15, 0.05, 0.7])
cb = fig.colorbar (scs [0], cax = cbar_ax)
cb.set_label ("Overshoot Factor")

ls = []
powers = []
labels = []
for i in scpowers: 
    if i not in powers:
        powers.append (i)
powers.sort ()
for i in powers:
    ls.append (plt.scatter ([], [], s = numToSize (i), linewidth = 0))
    labels.append (str (i))
# axes [1].legend (ls, labels, scatterpoints = 1, title = "SC Power", loc = 'lower right')
# ax.set_yscale ('log')
fig.canvas.mpl_connect('pick_event', onpick)

plt.show ()