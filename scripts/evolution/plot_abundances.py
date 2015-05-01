#!/usr/bin/env python

from kepler_utils.plots.parser import PlotArgumentParser

import matplotlib.pyplot as plt
import kepler_utils.plots.abundances
import kepler_utils.records.dump

import sys
import math
import numpy as np
import astropy.units as u
import astropy.constants as const

parser = PlotArgumentParser ()
namespace = parser.parse_args ()

record = kepler_utils.records.dump.DataDump (namespace.input_file, False)

fig = plt.figure ()
ax = fig.add_axes ([0.2, 0.2, 0.6, 0.75])

abun = kepler_utils.plots.abundances.AbundancePlot (ax, record)
plots = abun.plotAll (10.**-4)

ax.legend (bbox_to_anchor=(1.0, 1.0), loc=2, borderaxespad=1)

# Alternately, this could be done equivalently with 
# abundances.jTDPlot (sys.argv [1], 10.**-4)
