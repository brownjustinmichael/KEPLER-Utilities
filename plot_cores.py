import cache
import database
import abundances
import numpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
import matplotlib
import astropy.units as u
import kipp
import sqlalchemy
import cnv

def numToSize (num):
    return 100 * (-num + 3.5)

session = database.Session ()

query = session.query (database.DumpFileEntry).filter (database.DumpFileEntry.binm10 > 14.9).filter (database.DumpFileEntry.binm10 < 15.1)
query = query.filter (database.DumpFileEntry.brumoson > 0.).filter (database.DumpFileEntry.woodscon > 0.)
query = query.filter (database.DumpFileEntry.osfactor >= 0.1).filter (database.DumpFileEntry.osfactor <= 1.0)
query = query.filter (database.DumpFileEntry.scpower >= 0.9).filter (database.DumpFileEntry.scpower <= 2.1)
query = query.filter (database.DumpFileEntry.state == 'presn')

entries = [entry for entry in query.all ()]
# data = [entry.get_data () for entry in entries]
sims = [entry.simulation for entry in entries]

hecores = u.Quantity ([entry.cache (session, 'he_core', cache.calculate_he_core) for entry in entries])
cocores = u.Quantity ([entry.cache (session, 'co_core', cache.calculate_co_core) for entry in entries])
necores = u.Quantity ([entry.cache (session, 'ne_core', cache.calculate_ne_core) for entry in entries])
sicores = u.Quantity ([entry.cache (session, 'si_core', cache.calculate_si_core) for entry in entries])
fecores = u.Quantity ([entry.cache (session, 'fe_core', cache.calculate_fe_core) for entry in entries])

tasbsg = u.Quantity ([sim.cnvfile.cache (session, 'tasbsg', cache.calculate_tasbsg) for sim in sims])

lines = numpy.zeros (len (entries))

osfactors = numpy.array ([entry.osfactor if (entry.brumoson > 0.0) else 1 for entry in entries])
scpowers = numpy.array ([entry.scpower for entry in entries])

fig, axes = plt.subplots (2, 2, sharex = True, figsize = (18, 10))

scs = []

# Plot 1
ax = axes [0] [0]
ax.set_ylabel ("C/O")
scs.append (ax.scatter (hecores.to (u.solMass), cocores.to (u.solMass), numToSize (scpowers), c = osfactors, linewidths = lines, picker = True, norm = clrs.LogNorm (), alpha = 0.75))

# Plot 2
ax = axes [0] [1]
ax.set_ylabel ("Ne")
scs.append (ax.scatter (hecores.to (u.solMass), necores.to (u.solMass), numToSize (scpowers), c = osfactors, linewidths = lines, picker = True, norm = clrs.LogNorm (), alpha = 0.75))

# Plot 3
ax = axes [1] [0]
ax.set_ylabel ("Si")
scs.append (ax.scatter (hecores.to (u.solMass), sicores.to (u.solMass), numToSize (scpowers), c = osfactors, linewidths = lines, picker = True, norm = clrs.LogNorm (), alpha = 0.75))

# Plot 4
ax = axes [1] [1]
ax.set_ylabel ("Time as BSG (yr)")
scs.append (ax.scatter (hecores.to (u.solMass), tasbsg.to (u.year), numToSize (scpowers), c = osfactors, linewidths = lines, picker = True, norm = clrs.LogNorm (), alpha = 0.75))

ax.set_yscale ('log')

# for corearray in cores [1:]:
#     ax = next (axiter)
#     ax.set_ylabel (next (nameiter) + " Core Mass")
#     print (corearray)
#     scs.append (ax.scatter (hecores, corearray.to (u.solMass), numToSize (scpowers), c = osfactors, linewidths = lines, picker = True, norm = clrs.LogNorm (), alpha = 0.75))

axes.flat [2].set_xlabel ("He Core Mass")
axes.flat [3].set_xlabel ("He Core Mass")

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
    
    for entry in session.query (database.DumpFileEntry).filter (database.DumpFileEntry.simulation == sims [i]).order_by (database.DumpFileEntry.toffset, database.DumpFileEntry.time).all ():
        if u.Quantity (entry.time + entry.toffset, 's') > u.Quantity (time, 'year'):
            newfig = abundances.jTDPlot (entry.file)
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
        
        newfig = kipp.jTDPlot (sims [i].cnvfile.file)
        # newfig = abundances.jTDPlot (entries [i].get_data ())

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