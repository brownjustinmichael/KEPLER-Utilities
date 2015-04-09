#!/usr/bin/env python

import math

import numpy as np
import astropy.units as u
import astropy.constants as const
# import matplotlib
# import matplotlib.pyplot as plt
# import matplotlib.colors as clrs
# import matplotlib.cm as cm
import sqlalchemy

from kepler_utils.records.dump import DataDump
import kepler_utils.database.database as db
import kepler_utils.plots.kipp
from kepler_utils.plots.abundances import jTDPlot

# from pandas import DataFrame
# from pandas.tools.plotting import scatter_matrix

def numToSize (num):
    return 70 * (-np.log2 (num) + 6)

def calculate_T_hshell (datadump):
    return datadump ["tn"] [np.argmax (datadump ["h1"] > 0.01)]

def calculate_r_hshell (datadump):
    return datadump ["rn"] [np.argmax (datadump ["h1"] > 0.01)]

def calculate_d_hshell (datadump):
    return datadump ["dn"] [np.argmax (datadump ["h1"] > 0.01)]

def calculate_L_hshell (datadump):
    return datadump ["xln"] [np.argmax (datadump ["h1"] > 0.01)]

def calculate_T_center (datadump):
    return datadump ["tn"] [0]

def calculate_efold_T_center (datadump):
    return datadump ["mass coordinate"] [np.argmax (datadump ["tn"] < datadump ["tn"] [0] / math.e)]

def calculate_efold_T_hshell (datadump):
    return datadump ["mass coordinate"] [np.argmax (datadump ["tn"] < calculate_T_hshell (datadump) / math.e)] - datadump ["mass coordinate"] [np.argmax (datadump ["h1"] > 0.01)]

def calculate_d_center (datadump):
    return datadump ["dn"] [0]

def calculate_p_center (datadump):
    return datadump ["pn"] [0]

def calculate_he_core (datadump):
    return datadump ['mass coordinate'] [np.argmax (datadump ['h1'] > 0.01)]
    
def calculate_he_abun (datadump):
    return datadump ['he4'] [0]
    
def calculate_lum (datadump):
    return datadump ['xln'] [-100]
    
def calculate_U_core (datadump):
    argmax = np.argmax (datadump ["h1"] > 0.01)
    print (argmax)
    return np.sum (const.G * datadump ["xm"] [:argmax] * datadump ["mass coordinate"] [:argmax] / (datadump ["rn"] [:argmax]) ** 2)
    
def calculate_U (datadump):
    return np.sum (const.G * datadump ["xm"] * datadump ["mass coordinate"] / (datadump ["rn"]) ** 2)
    
def calculate_h_10fold (datadump):
    return datadump ["mass coordinate"] [np.argmax (datadump ["h1"] > 0.2)] - datadump ["mass coordinate"] [np.argmax (datadump ["h1"] > 0.01)]
    
session = db.Session ()

# for tag in ["D", "A", "B"]:
query = db.basicQuery (session).filter (db.SimulationEntry.tags.contains (db.Tag.get (session, "OS/SC Grid")))
query = query.filter (db.SimulationEntry.tags.contains (db.Tag.get (session, "Low dtcp")))
# query = query.filter (db.SimulationEntry.tags.contains (db.Tag.get (session, tag)))
# query = query.filter (~db.SimulationEntry.tags.contains (db.Tag.get (session, "A")))
# query = query.filter (~db.SimulationEntry.tags.contains (db.Tag.get (session, "B")))
# query = query.filter (~db.SimulationEntry.tags.contains (db.Tag.get (session, "D")))
# query = query.filter (~db.SimulationEntry.tags.contains (db.Tag.get (session, "Blue Loop")))
# query = query.filter (~db.SimulationEntry.tags.contains (db.Tag.get (session, "C")))
query = query.filter (db.DumpFileEntry.brumoson > 0.).filter (db.DumpFileEntry.woodscon > 0.).filter (db.DumpFileEntry.binm10 < 19.0)#.filter (db.DumpFileEntry.binm10 > 19.0)
query = query.filter (db.DumpFileEntry.state == 'presn').filter (db.SimulationEntry.cnvfiles.any ())

sims = [sim for sim, entry in query.all ()]

states = ["heign"] + ["%d" % i for i in range (18000, 32000, 2000)] + ["hedep"]

results = db.cache (session, sims, {"T_center": calculate_T_center, "L_hshell": calculate_L_hshell, "r_hshell": calculate_r_hshell, "he_core": calculate_he_core, "xlum": calculate_lum, "he_abun": calculate_he_abun, "T_efold_center": calculate_efold_T_center, "T_efold_hshell": calculate_efold_T_hshell, "U_core": calculate_U_core, "U_total": calculate_U, "h_10fold": calculate_h_10fold}, states)

results ["U_env"] = results ["U_total"] - results ["U_core"]
results ["L-L_he"] = results ["xlum"] - results ["L_hshell"]
results.pop ("xlum")
results.pop ("U_total")

results.pop ("r_hshell")
results.pop ("U_core")

#Try radius of 2.5 solar masses, grad_ad



# as h shell burns outward, t_hshell goes down (why?)
# Star needs so much radiation (why?) to prevent collapse, so as Hshell luminosity falls, He core luminosity grows
# collapsing core leads to expanding cool envelope which then convects
# Halted by degenerate core?

results ["radius"] = np.array ([[sim.getStateDump (state).radius for state in states] for sim in sims]) * u.cm

# results ["osfactors"] = u.Quantity ([u.Quantity ([sim.osfactor] * len (states)) for sim in sims])
# results ["scpowers"] = u.Quantity ([u.Quantity ([sim.scpower] * len (states)) for sim in sims])

# dfresults = {}
# for name in results:
#     dfresults [name] = np.log10 (results [name].flatten ().value)
#
# dfresults ["told"] = np.array ([[(sim.getStateDump (state).told - sim.getStateDump ("heign").told) / (sim.getStateDump ("hedep").told - sim.getStateDump ("heign").told) for state in states] for sim in sims]).flatten ()
#
# mask = dfresults ["told"] > 1.0
# for name in dfresults:
#     dfresults [name] [mask] = np.nan
#
# df = DataFrame.from_dict (dfresults)
#
# scatter_matrix (df)
#
# # fig, axes = plt.subplots (1, 1, sharex = False, figsize = (18, 10))
# #
# # ax = axes
# # ax.set_xlabel ("Temperature (at H shell)")
# # ax.set_ylabel ("Radius (solar radii)")
# # # ax.set_yscale ("log")
# # # ax.set_xscale ("log")
# # array = results ["T_hshell"]
# # yarray = results ["L_hshell"]
# # plot = ax.plot (array.T, yarray.T, c = "black")
# # sc = ax.scatter (array, yarray, c = results ["osfactors"], s = numToSize (results ["scpowers"]), alpha = 0.75, vmin = 0.1, vmax = 1.0)
# #
# #
# # fig.subplots_adjust (right = 0.8, hspace = 0)
# # cbar_ax = fig.add_axes ([0.85, 0.15, 0.05, 0.7])
# # cb = fig.colorbar (sc, cax = cbar_ax)
# # cb.set_label ("Overshoot Factor")
# #
# # # ls = []
# # # powers = []
# # # labels = []
# # # for i in scpowers:
# # #     if i not in powers:
# # #         powers.append (i)
# # # powers.sort ()
# # # for i in powers:
# # #     ls.append (plt.scatter ([], [], s = numToSize (i), linewidth = 0))
# # #     labels.append (str (i))
# # # axes [1].legend (ls, labels, scatterpoints = 1, title = "SC Power", loc = 'lower right')
# # # ax.set_yscale ('log')
# #
# # print ("Showing...")
#
# plt.show ()