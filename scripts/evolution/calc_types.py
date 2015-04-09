#!/usr/bin/env python

import math

import numpy as np
import astropy.units as u
import sqlalchemy

from kepler_utils.records.dump import DataDump
from kepler_utils.records.cnv import CNVFile
import kepler_utils.database.database as db

def numToSize (num):
    return 70 * (-np.log2 (num) + 6)

def calculate_t_rsg (datacnv):
    results = []
    for model in datacnv.models:
        previous = None
        for i in range (len (model ["yzip"]) - 1, -1, -1):
            if model ["yzip"] [i] == "C":
                if i == 0:
                    results.append (model ["xmcoord"] [0])
                    break
                results.append (model ["xmcoord"] [model ["iconv"] [i - 1]])
                break
    smoothed = moving_average (np.diff (results), 50) [14000:32000]
    argmin = np.argmin (smoothed)
    return datacnv.models [14000 + argmin] ["timesec"] * u.s

def calculate_conv_min (datacnv):
    results = []
    for model in datacnv.models:
        previous = None
        for i in range (len (model ["yzip"]) - 1, -1, -1):
            if model ["yzip"] [i] == "C":
                if i == 0:
                    results.append (model ["xmcoord"] [0])
                    break
                results.append (model ["xmcoord"] [model ["iconv"] [i - 1]])
                break
    smoothed = moving_average (np.diff (results), 50) [14000:32000]
    argmin = np.argmin (smoothed)
    return (datacnv.models [argmin] ["xmcoord"] [-1] * u.g - results [argmin] * u.g)

def calculate_t_bsg (datacnv):
    cmax = calculate_conv_max (datacnv)
    if cmax != cmax or cmax > (1.0 * u.solMass):
        return np.nan * u.s
    results = []
    for model in datacnv.models:
        previous = None
        for i in range (len (model ["yzip"]) - 1, -1, -1):
            if model ["yzip"] [i] == "C":
                if i == 0:
                    results.append (model ["xmcoord"] [0])
                    break
                results.append (model ["xmcoord"] [model ["iconv"] [i - 1]])
                break
    smoothed = moving_average (np.diff (results), 50) [14000:32000]
    argmax = np.argmax (smoothed)
    return (datacnv.models [14000 + argmax] ["timesec"] * u.s)

def calculate_t_conv_max (datacnv):
    results = []
    for model in datacnv.models:
        previous = None
        for i in range (len (model ["yzip"]) - 1, -1, -1):
            if model ["yzip"] [i] == "C":
                if i == 0:
                    results.append (model ["xmcoord"] [-1])
                    break
                results.append (model ["xmcoord"] [model ["iconv"] [i - 1]])
                break
    argmin = 32000
    end = 32000
    while argmin > end - 100 and argmin > 20000:
        end -= 1000
        smoothed = moving_average (results, 50) [14000:end]
        argmin = 14000 + np.argmin (smoothed)
    argmax = argmin + np.argmax (results [argmin + np.argmax (smoothed [argmin - 14000:] > 1.0e25):])
    if argmax < argmin + 100:
        return np.nan * u.s
    else:
        return (datacnv.models [argmax] ["timesec"] * u.s)

def calculate_conv_max (datacnv):
    results = []
    for model in datacnv.models:
        previous = None
        for i in range (len (model ["yzip"]) - 1, -1, -1):
            if model ["yzip"] [i] == "C":
                if i == 0:
                    results.append (model ["xmcoord"] [-1])
                    break
                results.append (model ["xmcoord"] [model ["iconv"] [i - 1]])
                break
    argmin = 32000
    end = 32000
    while argmin > end - 100 and argmin > 20000:
        end -= 1000
        smoothed = moving_average (results, 50) [14000:end]
        argmin = 14000 + np.argmin (smoothed)
    argmax = argmin + np.argmax (results [argmin + np.argmax (smoothed [argmin - 14000:] > 1.0e25):])
    if argmax < argmin + 100:
        return np.nan * u.g
    else:
        return (datacnv.models [argmax] ["xmcoord"] [-1] * u.g - results [argmax] * u.g)
    
def calculate_t_shell (datacnv):
    results = []
    heabun = []
    for model in datacnv.models:
        first = False
        found = False
        for i in range (len (model ["yzip"])):
            if model ["yzip"] [i] == "C":
                if first:
                    results.append (model ["xmcoord"] [model ["iconv"] [i]])
                    found = True
                    break
                else:
                    first = True
        if not found:
            results.append (model ["xmcoord"] [0])
        heabun.append (model ["abun_cnv"] [4])
    argmax = np.argmax (heabun)
    argmax = argmax + np.argmax (np.array (heabun [argmax:]) < 0.999 * heabun [argmax])
    return (datacnv.models [argmax + np.argmax (results [argmax:32000])] ["timesec"] * u.s)
    
def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n
    
session = db.Session ()

query = db.basicQuery (session).filter (db.SimulationEntry.tags.contains (db.Tag.get (session, "OS/SC Grid")))
query = query.filter (db.SimulationEntry.tags.contains (db.Tag.get (session, "Low dtcp")))
# query = query.filter (db.SimulationEntry.tags.contains (db.Tag.get (session, "B")))
# query = query.filter (~db.SimulationEntry.tags.contains (db.Tag.get (session, "C")))
# query = query.filter (~db.SimulationEntry.tags.contains (db.Tag.get (session, "Blue Loop")))
query = query.filter (db.DumpFileEntry.binm10 < 19.0)
query = query.filter (db.DumpFileEntry.state == 'presn').filter (db.SimulationEntry.cnvfiles.any ())

sims = [sim for sim, dump in query.all ()]

results = db.cnv_cache (session, sims, {"t_rsg": calculate_t_rsg, "t_bsg": calculate_t_bsg, "conv_max": calculate_conv_max, "t_conv_max": calculate_t_conv_max, "t_shell": calculate_t_shell, "conv_min": calculate_conv_min})

# Add depth of first envelope drop

results ["heign"] = [sim.getStateDump ("heign").told for sim in sims] * u.s
results ["hedep"] = [sim.getStateDump ("hedep").told for sim in sims] * u.s

osfactors = np.array ([entry.osfactor if (entry.brumoson > 0.0) else 1 for dump, entry in query.all ()])
scpowers = np.array ([entry.scpower for dump, entry in query.all ()])
