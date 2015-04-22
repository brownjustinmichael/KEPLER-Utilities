#!/opt/local/bin/python

import numpy as np
import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
import database.cache

import database.database
import plots.shiftlog
import plots.abundances

session = database.database.Session ()

query = session.query (database.database.SimulationEntry).filter (database.database.SimulationEntry.binm10 > 14.9).filter (database.database.SimulationEntry.binm10 < 15.1)
query = query.filter (database.database.SimulationEntry.osfactor >= 0.1)
query = query.filter (database.database.SimulationEntry.osfactor <= 1.0)
query = query.filter (database.database.SimulationEntry.scpower == 1.0)
# query = query.filter (database.DumpFileEntry.state == 'sidep')

sims = [entry for entry in query.all ()]

ordered = [session.query (database.database.DumpFileEntry).filter_by (simulation = sim).filter (database.database.DumpFileEntry.ncyc < 15000).order_by (database.database.DumpFileEntry.ncyc).all () for sim in sims]

fig, ax = plt.subplots (1, 1, sharex = True)

datas = []
namess = []
timess = []
radiis = []
hecores = []
liness = []
colorss = []
for sim in ordered:
    # data = [dumpentry.get_data () for dumpentry in sim]
    if len (sim) > 0:
        datas.extend (list (sim))

        # names = [dump.namep for dump in sim]
        # namess.extend (list (names))
        times = u.Quantity ([dumpentry.time + dumpentry.toffset for dumpentry in sim], u.s)
        timess.extend (list (times))
        
        radii = u.Quantity ([dumpentry.radius for dumpentry in sim], u.cm)
        radiis.extend (list (radii))
        
        hecore = u.Quantity ([dumpentry.cache (session, "hecore", database.cache.calculate_he_core) for dumpentry in sim])
        hecores.extend (list (hecore))
        lines = np.zeros (len (hecores))
        liness.extend (list (lines))

        colors = np.array ([dump.osfactor for dump in sim])
        colorss.extend (list (colors))
    
        ax.plot (hecore.to (u.solMass), radii.to (u.solRad), c = 'black')

figs = {}

sc = ax.scatter (u.Quantity (hecores).to (u.solMass), u.Quantity (radiis).to (u.solRad), linewidth = liness, s = 100,  picker = True, c = colorss, norm = clrs.LogNorm (), alpha = 0.75)

ax.set_yscale ('log')

def onclose (event):
    plt.figure (fig.number)
    i = figs.pop (event.canvas)
    liness [i] = 0
    sc.set_lw (liness)
    plt.draw ()

def onpick (event):
    ind = event.ind
    
    for i in ind:
        plt.figure (fig.number)
        liness [i] = 3
        event.artist.set_lw (liness)
        plt.draw ()
        
        record = datas [i].get_data ()

        newfig = plt.figure ()
        ax = newfig.add_axes ([0.1, 0.1, 0.6, 0.75])
        newfig.suptitle ("Overshoot Factor = %.3f, %s" % (colorss [i], "ON" if record.parameters ['brumoson'] > 0 else "OFF"))

        abun = plots.abundances.AbundancePlot (ax, record)
        newplots = abun.plotAll (ymin = 10.**-4)

        ax.legend (bbox_to_anchor=(1.05, 1.2), loc=2, borderaxespad=1)

        figs [newfig.canvas] = i
        newfig.canvas.mpl_connect('close_event', onclose)

    plt.show ()

fig.subplots_adjust (right = 0.8, hspace = 0)
cbar_ax = fig.add_axes ([0.85, 0.15, 0.05, 0.7])
fig.colorbar (sc, cax = cbar_ax)
ax.set_yscale ('log')
fig.canvas.mpl_connect('pick_event', onpick)

plt.show ()