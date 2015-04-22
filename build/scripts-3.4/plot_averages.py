#!/opt/local/bin/python

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

def filter_ossc (query, os = None, sc = None):
    if sc is not None:
        query = query.filter (db.DumpFileEntry.scpower >= sc / 1.01, db.DumpFileEntry.scpower <= sc * 1.01)
    if os is not None:
        query = query.filter (db.DumpFileEntry.osfactor >= os - 0.05, db.DumpFileEntry.osfactor <= os + 0.05)
    return query

def numToSize (num):
    return 70 * (-np.log2 (num) + 6)

session = db.Session ()

osfactors = np.arange (0.1, 1.1, 0.1)
scpowers = 2.0 ** np.arange (-3, 6, 1)

query = db.basicQuery (session).filter (db.SimulationEntry.tags.contains (db.Tag.get (session, "OS/SC Grid")))
query = query.filter (~db.SimulationEntry.tags.contains (db.Tag.get (session, "Blue Loop")))
query = query.filter (db.SimulationEntry.tags.contains (db.Tag.get (session, "Low dtcp")))
query = query.filter (db.DumpFileEntry.binm10 < 19.0)#.filter (db.DumpFileEntry.binm10 > 19.0)
query = query.filter (db.DumpFileEntry.state == 'presn').filter (db.SimulationEntry.cnvfiles.any ())
results = [filter_ossc (query, osfactor).order_by (db.DumpFileEntry.scpower).all () for osfactor in osfactors]

result = []
for i in range (len (osfactors)):
    result.append ([])
    kmin = 0
    for j in range (len (scpowers)):
        for sim, dump in results [i] [kmin:]:
            if sim.scpower < scpowers [j] / 1.01:
                kmin += 1
                continue
            if sim.scpower > scpowers [j] * 1.01:
                result [-1].append ((None, None))
                break
            result [-1].append ((sim, dump))
            break
    while len (result [-1]) < len (scpowers):
        result [-1].append ((None, None))

entries = [[entry for sim, entry in array] for array in result]
sims = [[sim for sim, entry in array] for array in result]

hecores = np.array ([[entry.cache (session, 'he_core', cache.calculate_he_core).to (u.solMass).value if entry is not None else np.nan for entry in array] for array in entries])
prehecores = np.array ([[sim.getStateDump ("hdep").cache (session, 'pre_he_core', cache.calculate_he_core).to (u.solMass).value if sim is not None else np.nan for sim in array] for array in sims])
# earlyhecores = u.Quantity ([sim.getStateDump ("heign").cache (session, 'early_he_core', cache.calculate_he_core) for sim in sims])
prehcontain = np.array ([[sim.getStateDump ("hdep").cache (session, 'pre_h_contained', cache.calculate_h_contained).to (u.solMass).value if sim is not None else np.nan for sim in array] for array in sims])
earlyhcontain = np.array ([[sim.getStateDump ("heign").cache (session, 'early_h_contained', cache.calculate_h_contained).to (u.solMass).value if sim is not None else np.nan for sim in array] for array in sims])
cocores = np.array ([[entry.cache (session, 'co_core', cache.calculate_co_core).to (u.solMass).value if entry is not None else np.nan for entry in array] for array in entries])
coratio = np.array ([[entry.cache (session, 'co_ratio', cache.calculate_co_ratio_in_core) if entry is not None else np.nan for entry in array] for array in entries])
compactness = np.array ([[entry.cache (session, 'compactness', cache.calculate_compactness_2_5) if entry is not None else np.nan for entry in array] for array in entries])
fecores = np.array ([[entry.cache (session, 'fe_core', cache.calculate_fe_core).to (u.solMass).value if entry is not None else np.nan for entry in array] for array in entries])

axis_0 = 0
axis_1 = 1
data = {}

for d, label in [(hecores, "he"), (prehecores, "phe"), (earlyhcontain, "eh"), (prehcontain, "ph"), (cocores, "co"), (coratio, "c/o"), (compactness, "com"), (fecores, "fe")]:
    temp = np.atleast_2d (np.nanmean (d, axis_0))
    if axis_0 == 1:
        temp = temp.T
    data [label] = np.array ((np.nanmean (d - temp, axis_1), np.nanstd (d - temp, axis_1)))

results = [query.all () for scpower in scpowers]

osfactors = np.array ([[osfactor for scpower in scpowers] for osfactor in osfactors])

scpowers = np.array ([[scpower for scpower in scpowers] for osfactor in osfactors])

fig, axes = plt.subplots (2, 2, sharex = False, figsize = (18, 10))

scs = []

if axis_1 == 0:
    x = np.mean (scpowers, axis_1)
else:
    x = np.mean (osfactors, axis_1)
for ax, d, label in zip (axes.flat, [data ["co"], data ["c/o"], data ["com"], data ["fe"]], ["C/O Core", "C/O", "Compactness", "Fe Core"]):
    ax.set_ylabel (label)
    scs.append (ax.scatter (x, d [0,:], s = numToSize (np.mean (scpowers, axis_1)), c = np.mean (osfactors, axis_1), picker = True, alpha = 0.75))
    ax.errorbar (x, d [0,:], yerr = d [1,:])
    if axis_1 == 0:
        ax.set_xscale ("log")

for ax in axes [-1]:
    if axis_1 == 0:
        ax.set_xlabel ("SC Power")
    else:
        ax.set_xlabel ("OS Factor")

# fig.subplots_adjust (right = 0.8, hspace = 0)
# cbar_ax = fig.add_axes ([0.85, 0.15, 0.05, 0.7])
# cb = fig.colorbar (scs [0], cax = cbar_ax)
# cb.set_label ("Overshoot Factor")

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

print ("Showing...")

plt.show ()