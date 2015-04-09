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
import kepler_utils.plots.kipp
from kepler_utils.plots.abundances import jTDPlot

from scipy.optimize import curve_fit

def numToSize (num):
    return 70 * (-np.log2 (num) + 6)

def calculate_T_hshell (datadump):
    return datadump ["tn"] [np.argmax (datadump ["h1"] > 0.01)]

def calculate_d_hshell (datadump):
    return datadump ["dn"] [np.argmax (datadump ["h1"] > 0.01)]

def calculate_L_hshell (datadump):
    return datadump ["xln"] [np.argmax (datadump ["h1"] > 0.01)]

def calculate_T_center (datadump):
    return datadump ["tn"] [0]

def calculate_d_center (datadump):
    return datadump ["dn"] [0]

def calculate_he_core (datadump):
    return datadump ['mass coordinate'] [np.argmax (datadump ['h1'] > 0.1)]
    
def calculate_he_abun (datadump):
    return datadump ['he4'] [0]
    
def fitting (xdata, a, b, c, d, f, g):
    return (d * np.arctan (a * xdata [:,0] + b * xdata [:,1] + c * xdata [:,2] + g) + f)
    
session = db.Session ()

query = db.basicQuery (session).filter (db.SimulationEntry.tags.contains (db.Tag.get (session, "OS/SC Grid")))
query = query.filter (db.SimulationEntry.tags.contains (db.Tag.get (session, "Low dtcp")))
# query = query.filter (db.SimulationEntry.tags.contains (db.Tag.get (session, "A")))
query = query.filter (~db.SimulationEntry.tags.contains (db.Tag.get (session, "C")))
query = query.filter (~db.SimulationEntry.tags.contains (db.Tag.get (session, "Blue Loop")))
query = query.filter (db.DumpFileEntry.brumoson > 0.).filter (db.DumpFileEntry.woodscon > 0.).filter (db.DumpFileEntry.binm10 < 19.0)#.filter (db.DumpFileEntry.binm10 > 19.0)
query = query.filter (db.DumpFileEntry.state == 'presn').filter (db.SimulationEntry.cnvfiles.any ())

entries = [entry for sim, entry in query.all ()]
sims = [sim for sim, entry in query.all ()]
cnvs = [sim.cnvfiles [0] for sim in sims]

if len (entries) == 0:
    raise ValueError ("No entries in query")

states = ["%d" % i for i in range (14000, 32000, 2000)]
# states = [["heign", 'T_hshell_ign'], ["heburn", 'T_hshell_burn'], ["hedep", 'T_hshell_dep']]

tcenter = np.array ([(sim.getStateDump (state).cache (session, "T_center_" + state, calculate_T_center)).value for state in states for sim in sims])
tcenter /= max (tcenter)
dcenter = np.array ([(sim.getStateDump (state).cache (session, "d_center_" + state, calculate_d_center)).value for state in states for sim in sims])
dcenter /= max (dcenter)
lum = np.array ([((sim.getStateDump (state).xlum) * u.erg / u.s).to (u.solLum).value for state in states for sim in sims])
abun = np.array ([sim.getStateDump (state).cache (session, "he_center_" + state, calculate_he_abun).value for state in states for sim in sims])
cores = np.array ([sim.getStateDump (state).cache (session, "hecore_" + state, calculate_he_core).value for state in states for sim in sims])
cores /= max (cores)
radius = np.array ([sim.getStateDump (state).radius for state in states for sim in sims])
radius /= max (radius)

popt, pcov = curve_fit (fitting, np.array ([tcenter, dcenter, abun]).T, radius)

print (popt)

osfactors = np.array ([[entry.osfactor] * len (states) if (entry.brumoson > 0.0) else 1 for entry in entries])
scpowers = np.array ([[entry.scpower] * len (states) for entry in entries])

fig, axes = plt.subplots (1, 1, sharex = False, figsize = (18, 10))

ax = axes
ax.set_xlabel ("Temperature (at H shell)")
ax.set_ylabel ("Radius (solar radii)")
# ax.set_yscale ("log")
# ax.set_xscale ("log")
array = tcenter
yarray = radius
# plot = ax.plot (array.T, yarray.T, c = "black")
sc = ax.scatter (array, yarray, c = osfactors, s = numToSize (scpowers), alpha = 0.75, vmin = 0.1, vmax = 1.0)
ax.plot (tcenter, fitting (np.array ([tcenter, dcenter, abun]).T, *popt))

fig.subplots_adjust (right = 0.8, hspace = 0)
cbar_ax = fig.add_axes ([0.85, 0.15, 0.05, 0.7])
cb = fig.colorbar (sc, cax = cbar_ax)
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

print ("Showing...")

plt.show ()