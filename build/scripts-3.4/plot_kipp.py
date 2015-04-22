#!/opt/local/bin/python

import sys
import argparse

import numpy as np
import astropy.units as u
import matplotlib
import matplotlib.colors
import matplotlib.patches
import matplotlib.pyplot as plt

import kepler_utils.records.cnv
import kepler_utils.plots.kipp
import kepler_utils.plots.shiftlog

# matplotlib.rc ('text', usetex = True)

parser = argparse.ArgumentParser ()
parser.add_argument ('input_file', default = None)
parser.add_argument ('--output', default = None)
parser.add_argument ('--useModels', default = False)
parser.add_argument ('--logSpace', default = None)
parser.add_argument ('--points', default = 400)

namespace = parser.parse_args ()

if namespace.logSpace is None:
    namespace.logSpace = not namespace.useModels

cnv_record = kepler_utils.records.cnv.CNVFile (namespace.input_file)

# Generate a matplotlib figure with one subplot
fig = plt.figure (figsize = (18,10))
ax = plt.subplot (111)

# Initialize the KippenhahnPlot object with the given axis and record file
kippplot = kepler_utils.plots.kipp.KippenhahnPlot (ax, cnv_record, useModels = namespace.useModels)

energy, = kippplot.plotEnergy (logspace = namespace.logSpace, points = namespace.points)
cb = kippplot.addEnergyColorBar (energy)

kippplot.plotConvection (logspace = namespace.logSpace, points = namespace.points)

# Generate the outer edge of the star
mass, = kippplot.plotMax ('xmcoord', 1.0 / kepler_utils.plots.kipp.msun, color = 'black', label = "Total Mass")

# Set the axis to be logarithmically scaled around the end of the star's life, with exponents increasing backward
if namespace.logSpace:
    ax.set_xscale ('shiftlog', base = 10, zero = kippplot.tend, sign = -1.0)

# Set the plot limits and show the grid and legend
ax.set_xlim (kippplot.xmin, kippplot.xmax - 10.**-5)
ax.set_ylim (kippplot.mmin, kippplot.mmax)
plt.grid ()
plt.legend ()

# This can equivalently be run with
# plots.kipp.jTDPlot (sys.argv [1], logspace = False)

if (namespace.output is not  None):
    plt.savefig (namespace.output)
else:
    # Plot the result
    plt.show ()
