from itertools import cycle

import matplotlib.pyplot as plt
from sqlalchemy.sql import func

from kepler_utils.database.database import Session, DumpFileEntry, SimulationEntry, cache
import kepler_utils.database.cache
from kepler_utils.plots.cores import CorePlot
from kepler_utils.plots.parser import PlotArgumentParser

parser = PlotArgumentParser (inputFile = False)
parser.add_argument ("--caches", nargs = "*", default = ["he_core"])
parser.add_argument ("--tags", nargs = "*", default = [])
parser.add_argument ("--states", nargs = "*", default = [])
parser.add_argument ("--binkey", default = "scpower")
parser.add_argument ("--binlog", dest = "binlog", action = 'store_true')
parser.add_argument ("--checkcache", dest = "checkcache", action = 'store_true')

parser.set_defaults (elements = True)
namespace = parser.parse_args ()

session = Session ()

q = session.query ().select_from (DumpFileEntry).join (SimulationEntry)
q = q.filter (DumpFileEntry.binm10 < 16.)
q = q.filter (SimulationEntry.contains (session, [tag for tag in namespace.tags]))

fig, axes = plt.subplots (len (namespace.caches), 1, sharex = True)
try:
    axes [0]
except TypeError:
    axes = [axes]

if namespace.checkcache:
    cache (session, q.add_entity (SimulationEntry).filter (DumpFileEntry.state == "presn").all (), {cache: getattr (kepler_utils.database.cache, cache) for cache in namespace.caches}, states = [state for state in namespace.states])

for cache, axis in zip (namespace.caches, axes):
    colors = cycle (["blue", "green", "red"])

    for state in namespace.states:
        qs = q.filter (DumpFileEntry.state == state)
        cign = CorePlot (axis, qs)
    
        binkey = getattr (DumpFileEntry, namespace.binkey)
        if namespace.binlog:
            binkey = func.log (binkey)
        
        s, e = cign.plotScatter (getattr (DumpFileEntry, namespace.binkey), cache, binkey = binkey, errorbars = True, color = next (colors), label = state)

    axis.set_ylabel (cache)

if namespace.binlog:
    axis.set_xscale ("log")
axis.set_xlabel (namespace.binkey)
axis.legend (loc = "upper left")