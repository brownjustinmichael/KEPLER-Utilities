from itertools import cycle

import matplotlib.pyplot as plt
from sqlalchemy.sql import func, or_

from kepler_utils.database.database import Session, DumpFileEntry, SimulationEntry, cache
import kepler_utils.database.cache
from kepler_utils.plots.cores import CorePlot
from kepler_utils.plots.parser import PlotArgumentParser

parser = PlotArgumentParser (inputFile = False)
parser.add_argument ("--tags", nargs = "*", default = [])
parser.add_argument ("--states", nargs = "*", default = [])
parser.add_argument ("--xkeys", nargs = "*", default = ["he_core"])
parser.add_argument ("--ykeys", nargs = "*", default = ["co_core"])
parser.add_argument ("--skey", default = "scpower")
parser.add_argument ("--ckey", default = "osfactor")
parser.add_argument ("--checkcache", dest = "checkcache", action = 'store_true')
parser.add_argument ("--exclude", nargs = "*", default = [])

namespace = parser.parse_args ()

session = Session ()

q = session.query ().select_from (DumpFileEntry).join (SimulationEntry)
q = q.filter (DumpFileEntry.binm10 < 16.)
q = q.filter (SimulationEntry.contains (session, [tag for tag in namespace.tags]))
if (len (namespace.exclude) > 0):
    q = q.filter (~SimulationEntry.contains (session, [tag for tag in namespace.exclude]))

fig, axes = plt.subplots (len (namespace.ykeys), len (namespace.xkeys))
try:
    axes [0]
except TypeError:
    axes = [axes]
    
try:
    axes [0] [0]
except TypeError:
    axes = [axes]

if namespace.checkcache:
    caches = []
    for key in namespace.xkeys + namespace.ykeys + [namespace.ckey, namespace.skey]:
        try:
            getattr (DumpFileEntry, key)
        except AttributeError:
            caches.append (key)

    cache (session, q.add_entity (SimulationEntry).filter (DumpFileEntry.complete).all (), {cache: getattr (kepler_utils.database.cache, cache) for cache in caches}, states = [state for state in namespace.states])

for i, ykey, axisList in zip (range (len (namespace.ykeys)), namespace.ykeys, axes):
    axisList [0].set_ylabel (ykey)
    for j, xkey, axis in zip (range (len (namespace.xkeys)), namespace.xkeys, axisList):
        colors = cycle (["blue", "green", "red"])
        
        axis.set_sharex = axes [0] [j]
        axis.set_sharey = axisList [0]

        for state in namespace.states:
            qs = q.filter (DumpFileEntry.state == state)
            cign = CorePlot (qs, axis)
            
            s = cign.plotScatter (xkey, ykey, skey = namespace.skey, ckey = namespace.ckey, label = state, resize = True)

for xkey, axis in zip (namespace.xkeys, axes [-1]):
    axis.set_xlabel (xkey)

axis.legend (loc = "upper left")

parser.exit (fig)