#!/usr/bin/env python

import math

import numpy as np
import astropy.units as u
import astropy.constants as const
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
import matplotlib.cm as cm
import sqlalchemy

from kepler_utils.records.dump import DataDump
import kepler_utils.database.database as db
import kepler_utils.plots.kipp
from kepler_utils.plots.abundances import jTDPlot

from pandas import DataFrame
from pandas.tools.plotting import scatter_matrix

def numToSize (num):
    return 70 * (-np.log2 (num) + 6)

def calculate_T_hshell (datadump):
    return datadump ["tn"] [np.argmax (datadump ["h1"] > 0.001)]

def calculate_r_hshell (datadump):
    return datadump ["rn"] [np.argmax (datadump ["h1"] > 0.001)]

def calculate_d_hshell (datadump):
    return datadump ["dn"] [np.argmax (datadump ["h1"] > 0.001)]

def calculate_L_hshell (datadump):
    return datadump ["xln"] [np.argmax (datadump ["h1"] > 0.001)]

def calculate_L_h (datadump):
    return datadump ["xln"] [-np.argmax ((datadump ["snn"] > 0.01 * np.max (datadump ["snn"])) [::-1])]

def calculate_T_center (datadump):
    return datadump ["tn"] [0]

def calculate_P_center (datadump):
    return datadump ["pn"] [0]

def calculate_efold_T_center (datadump):
    return datadump ["mass coordinate"] [np.argmax (datadump ["tn"] < datadump ["tn"] [0] / math.e)]

def calculate_efold_T_hshell (datadump):
    return datadump ["mass coordinate"] [np.argmax (datadump ["tn"] < calculate_T_hshell (datadump) / math.e)] - datadump ["mass coordinate"] [np.argmax (datadump ["h1"] > 0.001)]

def calculate_d_center (datadump):
    return datadump ["dn"] [0]

def calculate_p_center (datadump):
    return datadump ["pn"] [0]

def calculate_he_core (datadump):
    return datadump ['mass coordinate'] [np.argmax (datadump ['h1'] > 0.001)]
    
def calculate_he_abun (datadump):
    return datadump ['he4'] [0]
    
def calculate_U_core (datadump):
    argmax = np.argmax (datadump ["h1"] > 0.001)
    return np.sum (const.G * datadump ["xm"] [:argmax] * datadump ["mass coordinate"] [:argmax] / (datadump ["rn"] [:argmax]))
    
def calculate_U (datadump):
    return np.sum (const.G * datadump ["xm"] * datadump ["mass coordinate"] / (datadump ["rn"]))
    
def calculate_h_10fold (datadump):
    return datadump ["mass coordinate"] [np.argmax (datadump ["h1"] > 0.2)] - datadump ["mass coordinate"] [np.argmax (datadump ["h1"] > 0.001)]
    
def calculate_U_shell (datadump):
    argmax = np.argmax (datadump ["h1"] > 0.001)
    return const.G * datadump ["mass coordinate"] [argmax] / datadump ["rn"] [argmax]
    
def calculate_L_edd (datadump):
    ledd = (4 * math.pi * (const.G * datadump ["mass coordinate"] [1:] + const.k_B / 0.62 / const.m_p * datadump ["rn"] [1:] ** 2 * (np.diff (datadump ["tn"]) + datadump ["tn"] [1:] * np.diff (datadump ["dn"]) / datadump ["dn"] [1:]) / np.diff (datadump ["rn"])) * const.c / datadump ["xkn"] [1:]).to (u.solLum)

    return (np.mean (ledd [np.argmax (ledd) - 10:np.argmax (ledd) + 10]))

def calculate_hshell_abun (datadump):
    return (1.0 * u.solMass).to (u.g) / (datadump ["mass coordinate"] [np.argmax (np.cumsum (datadump ["xm"] * datadump ["h1"]) > 1.0 * u.solMass)] - datadump ["mass coordinate"] [np.argmax (datadump ["h1"] > 0.001)])
    
def conv_extent (datadump):
    return datadump ["mass coordinate"] [-10 - np.argmin (datadump ["icon"] [:-10:-1] == "conv")]
    
def calculate_K_core (datadump):
    argmax = np.argmax (datadump ["h1"] > 0.001)
    return np.sum (const.k_B * datadump ["xm"] [:argmax] * datadump ["tn"] [:argmax])
    
def calculate_K (datadump):
    return np.sum (const.k_B * datadump ["xm"] * datadump ["tn"])

def calculate_kappa (datadump):
    return datadump ["xkn"] [-1]
    
def calculate_p_support (datadump):
    return 4 * math.pi * const.k_B / 0.62 / const.m_p * datadump ["rn"] [-1] ** 2 * (datadump ["tn"] [-1] - datadump ["tn"] [-2] + datadump ["tn"] [-1] * (datadump ["dn"] [-1] - datadump ["dn"] [-2]) / datadump ["dn"] [-1]) / (datadump ["rn"] [-1] - datadump ["rn"] [-2]) * const.c / datadump ["xkn"] [-1]
    
def calculate_L_ratio (datadump):
    return np.max (datadump ["xln"]) / datadump ["xln"] [-1]

session = db.Session ()

# for tag in ["D", "A", "B"]:
query = db.basicQuery (session).filter (db.SimulationEntry.tags.contains (db.Tag.get (session, "OS/SC Grid")))
query = query.filter (db.SimulationEntry.tags.contains (db.Tag.get (session, "Low dtcp")))
# query = query.filter (db.SimulationEntry.name == "s15t.44")
# query = query.filter (db.SimulationEntry.tags.contains (db.Tag.get (session, tag)))
# query = query.filter (~db.SimulationEntry.tags.contains (db.Tag.get (session, "A")))
# query = query.filter (~db.SimulationEntry.tags.contains (db.Tag.get (session, "B")))
# query = query.filter (~db.SimulationEntry.tags.contains (db.Tag.get (session, "D")))
# query = query.filter (~db.SimulationEntry.tags.contains (db.Tag.get (session, "Blue Loop")))
# query = query.filter (~db.SimulationEntry.tags.contains (db.Tag.get (session, "C")))
query = query.filter (db.DumpFileEntry.brumoson > 0.).filter (db.DumpFileEntry.woodscon > 0.).filter (db.DumpFileEntry.binm10 < 19.0)#.filter (db.DumpFileEntry.binm10 > 19.0)
query = query.filter (db.DumpFileEntry.state == 'presn').filter (db.SimulationEntry.cnvfiles.any ())

sims = [sim for sim, entry in query.all ()]

def getTag (tags):
    for tag in ["A", "B", "C", "D", "Blue Loop"]:
        if db.Tag.get (session, tag) in tags:
            return tag
    return None

colors = {"A": "red", "B": "green", "D": "red", "Blue Loop": "blue", "C": "Green", None: "white"}

tags = [getTag (sim.tags) for sim in sims]

states = ["heign"] + ["%d" % i for i in range (18000, 32000, 2000)] + ["hedep"]

results = db.cache (session, sims, {"L_he": calculate_L_hshell, "M_core": calculate_he_core, "U_core": calculate_U_core, "U_total": calculate_U, "h_10fold": calculate_h_10fold, "L_edd": calculate_L_edd, "L_h": calculate_L_h, "r_hshell": calculate_r_hshell, "habun": calculate_hshell_abun, "heabun": calculate_he_abun, "T_center": calculate_T_center, "P_center": calculate_P_center}, states)

results ["L_edd"] = results ["L_edd"].to (u.erg / u.s)
results ["L_h"] -= results ["L_he"]

results.pop ("U_core")
results.pop ("U_total")
# results ["U_env"] = results ["U_total"] - results ["U_core"]

# results.pop ("h_10fold")
# results.pop ("L_h")
# results.pop ("L_edd")

results ["radius"] = np.array ([[sim.getStateDump (state).radius for state in states] for sim in sims]) * u.cm
results ["L"] = np.array ([[sim.getStateDump (state).xlum for state in states] for sim in sims]) * u.erg / u.s
results ["mass"] = np.array ([[sim.getStateDump (state).totm for state in states] for sim in sims]) * u.g
results ["mass loss"] = np.array ([[sim.getStateDump (state).xmlossr for state in states] for sim in sims]) * u.g / u.s

results ["osfactors"] = u.Quantity ([u.Quantity ([sim.osfactor] * len (states)) for sim in sims])
results ["scpowers"] = u.Quantity ([u.Quantity ([sim.scpower] * len (states)) for sim in sims])

# dR = u.Quantity (np.concatenate ([np.diff (results ["radius"]), np.array ([[0]] * len (results ["radius"]))], 1).value, results ["radius"].unit)
# dL_h = u.Quantity (np.concatenate ([np.diff (results ["L_h"]), np.array ([[0]] * len (results ["L_h"]))], 1).value, results ["L_h"].unit)
# dL_he = u.Quantity (np.concatenate ([np.diff (results ["L_he"]), np.array ([[0]] * len (results ["L_he"]))], 1).value, results ["L_he"].unit)

# results ["dL/dR"] = (dL_h + dL_he) / dR
# results ["dL_h/dR"] = (dL_h) / dR
# results ["dL_he/dR"] = (dL_he) / dR

dfresults = {}
for name in results:
    dfresults [name] = np.log10 (results [name].flatten ().value)

# dfresults ["dL_h/dR"] = (results ["dL_h"] / results ["dR"]).flatten ()
# dfresults ["dL_he/dR"] = (results ["dL_he"] / results ["dR"]).flatten ()

# dfresults ["dL/dR"] = results ["dL/dR"].flatten ()
# dfresults ["dL_h/dR"] = results ["dL_h/dR"].flatten ()
# dfresults ["dL_he/dR"] = results ["dL_he/dR"].flatten ()

# dfresults ["??"] = results ["xlum"].flatten () - math.pi * 4 * const.G * results ["mass"].flatten () * const.c / u.Quantity (2.2, u.cm ** 2 / u.g)

dfresults ["told"] = np.array ([[(sim.getStateDump (state).told - sim.getStateDump ("heign").told) / (sim.getStateDump ("hedep").told - sim.getStateDump ("heign").told) for state in states] for sim in sims]).flatten ()

mask = dfresults ["told"] > 1.0
# mask = np.logical_or (mask, np.logical_or (np.abs (((results ["L_h"] + results ["L_he"]) / results ["L"]).flatten () - 1) > 0.2), np.logical_or (dfresults ["L_h"] > 39, dfresults ["L_edd"] > 39)))
# mask = np.logical_or (mask, np.abs (dfresults ["dL/dR"]) > 2e25 * u.erg / u.s / u.cm)

dfresults.pop ("told")
dfresults.pop ("scpowers")
dfresults.pop ("osfactors")

for name in dfresults:
    dfresults [name] [mask] = np.nan

df = DataFrame.from_dict (dfresults)

scatter_matrix (df)

# fig, axes = plt.subplots (1, 1, sharex = False, figsize = (18, 10))
#
# ax = axes
# ax.set_xlabel ("Temperature (at H shell)")
# ax.set_ylabel ("Radius (solar radii)")
# # ax.set_yscale ("log")
# # ax.set_xscale ("log")
# array = results ["mass loss"]
# yarray = results ["radius"]
# plot = ax.plot (array.T, yarray.T, c = "black")
# sc = ax.scatter (array, yarray, c = results ["osfactors"], s = numToSize (results ["scpowers"]), alpha = 0.75, vmin = 0.1, vmax = 1.0)
#
#
# fig.subplots_adjust (right = 0.8, hspace = 0)
# cbar_ax = fig.add_axes ([0.85, 0.15, 0.05, 0.7])
# cb = fig.colorbar (sc, cax = cbar_ax)
# cb.set_label ("Overshoot Factor")
#
# # ls = []
# # powers = []
# # labels = []
# # for i in scpowers:
# #     if i not in powers:
# #         powers.append (i)
# # powers.sort ()
# # for i in powers:
# #     ls.append (plt.scatter ([], [], s = numToSize (i), linewidth = 0))
# #     labels.append (str (i))
# # axes [1].legend (ls, labels, scatterpoints = 1, title = "SC Power", loc = 'lower right')
# # ax.set_yscale ('log')

print ("Showing...")

plt.show ()