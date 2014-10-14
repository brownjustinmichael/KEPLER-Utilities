import numpy as np
import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as clrs

import database.database
import plots.shiftlog
import plots.abundances

session = database.database.Session ()

query = session.query (database.database.SimulationEntry).filter (database.database.SimulationEntry.binm10 > 79.9).filter (database.database.SimulationEntry.binm10 < 80.1)
query = query.filter (database.database.SimulationEntry.osfactor >= 0.1)
query = query.filter (database.database.SimulationEntry.osfactor <= 1.0)
# query = query.filter (database.DumpFileEntry.state == 'sidep')

sims = [entry for entry in query.all ()]

ordered = [session.query (database.database.DumpFileEntry).filter_by (simulation = sim).filter (database.database.DumpFileEntry.ncyc % 50000 == 0).order_by (database.database.DumpFileEntry.ncyc).all () for sim in sims]

fig, ax = plt.subplots (1, 1, sharex = True)

datas = []
namess = []
timess = []
nabunds = []
liness = []
colorss = []
for sim in ordered:
    data = [dumpentry.get_data () for dumpentry in sim]
    if len (data) > 0:
        datas.extend (list (data))

        names = [dump.namep for dump in data]
        namess.extend (list (names))
        times = u.Quantity ([dump.parameters ['time'] + dump.parameters ['toffset'] for dump in data])
        timess.extend (list (times))

        nabund = u.Quantity ([np.sum (dump ['n14'] * dump ['xm']) for dump in data])
        nabunds.extend (list (nabund))
        lines = np.zeros (len (nabund))
        liness.extend (list (lines))

        colors = np.array ([dump.parameters ['osfactor'] for dump in data])
        colorss.extend (list (colors))
    
        ax.plot (times.to (u.year), nabund.to (u.solMass), c = 'black')

figs = {}

sc = ax.scatter (u.Quantity (timess).to (u.year), u.Quantity (nabunds).to (u.solMass), linewidth = liness, s = 100,  picker = True, c = colorss, norm = clrs.LogNorm (), alpha = 0.75)

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
        
        record = datas [i]

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