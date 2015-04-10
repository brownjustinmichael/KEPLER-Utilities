#!/usr/bin/env python

import matplotlib.pyplot as plt
import kepler_utils.plots.abundances
import kepler_utils.records.dump

import sys
import math
import numpy as np
import astropy.units as u
import astropy.constants as const

if (len (sys.argv) < 2):
    print ("Usage: plot_abundances dump_file (output_file)")
    exit ()

record = kepler_utils.records.dump.DataDump (sys.argv [1], False)

fig = plt.figure ()
ax = fig.add_axes ([0.1, 0.1, 0.6, 0.75])

abun = kepler_utils.plots.abundances.AbundancePlot (ax, record)
plots = abun.plotAll (10.**-4)

ax.legend (bbox_to_anchor=(1.05, 1.2), loc=2, borderaxespad=1)

# Alternately, this could be done equivalently with 
# abundances.jTDPlot (sys.argv [1], 10.**-4)
    
if (len (sys.argv) > 2):
    plt.savefig (sys.argv [2])
else:
    # Plot the result
    plt.show ()