import cnv
import kipp
import matplotlib
import matplotlib.colors
import matplotlib.pyplot as plt
import shiftlog
import matplotlib.patches
import astropy.units as u
import numpy as np

import sys

# sys.argv.append ("/Users/justinbrown/Codes/kepler/run/s15/s15o0s3.cnv")

# matplotlib.rc ('text', usetex = True)

# Read in the KEPLER cnv output file into cnv_record
if (len (sys.argv) < 2):
    print ("Usage: plot_kipp cnv_file (output_file)")
    exit ()

cnv_record = cnv.CNVFile (sys.argv [1])

# Generate a matplotlib figure with one subplot
fig = plt.figure (figsize = (18,10))
ax = plt.subplot (111)

# Initialize the KippenhahnPlot object with the given axis and record file
kippplot = kipp.KippenhahnPlot (ax, cnv_record)

energy, = kippplot.plotEnergy ()
cb = kippplot.addEnergyColorBar (energy)

kippplot.plotConvection ()

# Generate the outer edge of the star
mass, = kippplot.plotMax ('xmcoord', 1.0 / kipp.msun, color = 'black', label = "Total Mass")

# Set the axis to be logarithmically scaled around the end of the star's life, with exponents increasing backward
ax.set_xscale ('shiftlog', base = 10, zero = kippplot.tend, sign = -1.0)

# Set the plot limits and show the grid and legend
ax.set_xlim (kippplot.tmin, kippplot.tmax - 10.**-5)
ax.set_ylim (kippplot.mmin, kippplot.mmax)
plt.grid ()
plt.legend ()

# This can equivalently be run with
# kipp.jTDPlot (sys.argv [1])

if (len (sys.argv) > 2):
    plt.savefig (sys.argv [2])
else:
    # Plot the result
    plt.show ()
