import math

import numpy as np
import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
import sqlalchemy

import records.cnv
import database.database
import database.cache
import plots.kipp
import plots.abundances

def numToSize (num):
    return 100 * (-np.log2 (num) + 3.5)

session = database.database.Session ()

query = session.query (database.database.DumpFileEntry).filter (database.database.DumpFileEntry.binm10 > 14.9).filter (database.database.DumpFileEntry.binm10 < 15.1)
query = query.filter (database.database.DumpFileEntry.brumoson > 0.).filter (database.database.DumpFileEntry.woodscon > 0.)
query = query.filter (database.database.DumpFileEntry.osfactor >= 0.09).filter (database.database.DumpFileEntry.osfactor <= 1.0)
query = query.filter (database.database.DumpFileEntry.scpower >= 0.1).filter (database.database.DumpFileEntry.scpower <= 8.1)
query = query.filter (database.database.DumpFileEntry.state == 'presn')

entries = [entry for entry in query.all ()]
# data = [entry.get_data () for entry in entries]
sims = [entry.simulation for entry in entries]
cnvs = [sim.cnvfile for sim in sims]

hecores = u.Quantity ([entry.cache (session, 'he_core', database.cache.calculate_he_core) for entry in entries])
earlyhecores = u.Quantity ([sim.get_state_dump ("hdep").cache (session, 'he_core', database.cache.calculate_he_core) for sim in sims])
lum = u.Quantity ([sim.get_state_dump ("heign").xlum * u.erg / u.s for sim in sims])
mass = u.Quantity ([sim.get_state_dump ("heign").totm * u.g for sim in sims])
massbetween = u.Quantity ([sim.get_state_dump ("hign").totm * u.g - sim.get_state_dump ("hdep").totm * u.g for sim in sims])

cocores = u.Quantity ([entry.cache (session, 'co_core', database.cache.calculate_co_core) for entry in entries])
necores = u.Quantity ([entry.cache (session, 'ne_core', database.cache.calculate_ne_core) for entry in entries])
sicores = u.Quantity ([entry.cache (session, 'si_core', database.cache.calculate_si_core) for entry in entries])
fecores = u.Quantity ([entry.cache (session, 'fe_core', database.cache.calculate_fe_core) for entry in entries])

tasbsg = u.Quantity ([cnv.cache (session, 'tasbsg', database.cache.calculate_tasbsg) for cnv in cnvs])

lines = np.zeros (len (entries))

osfactors = np.array ([entry.osfactor if (entry.brumoson > 0.0) else 1 for entry in entries])
scpowers = np.array ([entry.scpower for entry in entries])

fig, axes = plt.subplots (3, 2, sharex = False, figsize = (18, 10))

scs = []

# Plot 1
ax = axes [0] [0]
ax.set_ylabel ("C/O")
scs.append (ax.scatter (hecores.to (u.solMass), cocores.to (u.solMass), numToSize (scpowers), c = osfactors, linewidths = lines, picker = True, norm = clrs.LogNorm (), alpha = 0.75))

# Plot 2
ax = axes [1] [0]
ax.set_ylabel ("Early He")
scs.append (ax.scatter (hecores.to (u.solMass), earlyhecores.to (u.solMass), numToSize (scpowers), c = osfactors, linewidths = lines, picker = True, norm = clrs.LogNorm (), alpha = 0.75))
ylim = ax.get_ylim ()
xlim = ax.get_xlim ()
ax.fill_between ((0,10), (0,10), 100, color = 'black', alpha = 0.25)
ax.set_ylim (ylim)
ax.set_xlim (xlim)

# Plot 3
ax = axes [2] [0]
ax.set_ylabel ("Time as BSG (yr)")
scs.append (ax.scatter (hecores.to (u.solMass), tasbsg.to (u.year), numToSize (scpowers), c = osfactors, linewidths = lines, picker = True, norm = clrs.LogNorm (), alpha = 0.75))
ax.set_yscale ('log')

ax = axes [0] [1]
ax.set_ylabel ("C/O")
scs.append (ax.scatter (massbetween.to (u.solMass), cocores.to (u.solMass), numToSize (scpowers), c = osfactors, linewidths = lines, picker = True, norm = clrs.LogNorm (), alpha = 0.75))
# ax.set_xscale ('log')

ax = axes [1] [1]
ax.set_ylabel ("Early He")
scs.append (ax.scatter (massbetween.to (u.solMass), earlyhecores.to (u.solMass), numToSize (scpowers), c = osfactors, linewidths = lines, picker = True, norm = clrs.LogNorm (), alpha = 0.75))
# ax.set_xscale ('log')

# Plot 6
ax = axes [2] [1]
ax.set_ylabel ("Time as BSG (yr)")
scs.append (ax.scatter (massbetween.to (u.solMass), tasbsg.to (u.year), numToSize (scpowers), c = osfactors, linewidths = lines, picker = True, norm = clrs.LogNorm (), alpha = 0.75))
ax.set_yscale ('log')
# ax.set_xscale ('log')


ax.set_yscale ('log')

axes.flat [5].set_xlabel ("Mass Lost on Main Sequence")
axes.flat [4].set_xlabel ("He Core Mass")

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
    
    for entry in session.query (database.database.DumpFileEntry).filter (database.database.DumpFileEntry.simulation == sims [i]).order_by (database.database.DumpFileEntry.toffset, database.database.DumpFileEntry.time).all ():
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
        
        newfig = plots.kipp.jTDPlot (cnvs [i].get_data ())
        plots.abundances.jTDPlot (sims [i].get_state_dump ("heign").get_data ())
        plots.abundances.jTDPlot (sims [i].get_state_dump ("hdep").get_data ())

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
axes [0] [1].legend (ls, labels, scatterpoints = 1, title = "SC Power")
# ax.set_yscale ('log')
fig.canvas.mpl_connect('pick_event', onpick)

plt.show ()