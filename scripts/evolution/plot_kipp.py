#!/usr/bin/env python

import sys
from kepler_utils.plots.parser import PlotArgumentParser

import numpy as np
import astropy.units as u
import matplotlib
import matplotlib.colors
import matplotlib.patches
import matplotlib.pyplot as plt

import kepler_utils.records.cnv
import kepler_utils.plots.kipp
import kepler_utils.plots.shiftlog
import kepler_utils.database.database as db

# matplotlib.rc ('text', usetex = True)

parser = PlotArgumentParser ()
parser.add_argument ('--useModels', dest = 'models', action = 'store_true')
parser.set_defaults (models = False)
parser.add_argument ('--logSpace', default = None)
parser.add_argument ('--points', default = 400)
parser.add_argument ('--min', default = None)
parser.add_argument ('--max', default = None)

namespace = parser.parse_args ()

if namespace.logSpace is None:
    namespace.logSpace = not namespace.models

try:
    cnv_record = kepler_utils.records.cnv.CNVFile (namespace.input_file)
except FileNotFoundError:
    session = db.Session ()
    cnv_record = session.query (db.CNVFileEntry).filter (db.CNVFileEntry.file.contains (namespace.input_file)).first ().get_data ()

if namespace.min == None and namespace.max == None:
    extent = None
else:
    if namespace.min is not None:
        mini = cnv_record.modelNear (cnv_record [-1] ["timesec"] * u.s - u.yr * 10.0 ** float (namespace.min))
    else:
        mini = cnv_record.modelNear (0.0 * u.s)

    if namespace.min is not None:
        maxi = cnv_record.modelNear (cnv_record [-1] ["timesec"] * u.s - u.yr * 10.0 ** float (namespace.max))
    else:
        maxi = cnv_record.modelNear ((cnv_record [-1] ["timesec"] - 1.0e-5) * u.s)
    print (mini, maxi)
    extent = range (mini, maxi)

# Generate a matplotlib figure with one subplot
fig = plt.figure (figsize = (18,10))
ax = plt.subplot (111)

# Initialize the KippenhahnPlot object with the given axis and record file
kippplot = kepler_utils.plots.kipp.KippenhahnPlot (ax, cnv_record, useModels = namespace.models)

energy, = kippplot.plotEnergy (logspace = namespace.logSpace, points = namespace.points, extent = extent)
cb = kippplot.addEnergyColorBar (energy)

kippplot.plotConvection (logspace = namespace.logSpace, points = namespace.points, extent = extent)

# Generate the outer edge of the star
mass, = kippplot.plotMax ('xmcoord', 1.0 / kepler_utils.plots.kipp.msun, color = 'black', label = "Total Mass")

# Set the axis to be logarithmically scaled around the end of the star's life, with exponents increasing backward
if namespace.logSpace:
    ax.set_xscale ('shiftlog', base = 10, zero = kippplot.tend, sign = -1.0)

# Set the plot limits and show the grid and legend
ax.set_xlim (kippplot.xmin if namespace.max == None else kippplot.xmax - 10.0 ** float (namespace.min), kippplot.xmax - 10.**-5 if namespace.max == None else (kippplot.xmax - 10. ** float (namespace.max)))
ax.set_ylim (kippplot.mmin, kippplot.mmax)
plt.grid ()
plt.legend ()

# This can equivalently be run with
# plots.kipp.jTDPlot (sys.argv [1], logspace = False)
