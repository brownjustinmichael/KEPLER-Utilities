#!/opt/local/bin/python

import matplotlib.pyplot as plt
import kepler_utils.plots.abundances
import kepler_utils.records.dump

import sys
import math
import argparse
import numpy as np
import astropy.units as u
import astropy.constants as const

parser = argparse.ArgumentParser ()
parser.add_argument ('input_file', default = None)
parser.add_argument ('--output', default = None)
parser.add_argument ('--useCells', default = False)

namespace = parser.parse_args ()

record = kepler_utils.records.dump.DataDump (namespace.input_file, False)

fig, axes = plt.subplots (4, 1, sharex = True, figsize = (6, 10))
if namespace.useCells:
    x = np.arange (len (record ["mass coordinate"].to (u.solMass)))
else:
    x = record ["mass coordinate"].to (u.solMass)

for ax, name in zip (axes, ["rn", "dn", "tn", "xln"]):
    ax.plot (x, (record [name]), label = name)
    ax.set_ylabel (name)
    if np.any (record [name] > 0):
        ax.set_yscale ("log")

if namespace.useCells:
    ax.set_xlabel ("cells")
else:
    ax.set_xlabel ("m")

# ax.legend (bbox_to_anchor=(1.05, 1.2), loc=2, borderaxespad=1)

plt.tight_layout ()

if (namespace.output is not None):
    plt.savefig (namespace.output)
else:
    # Plot the result
    plt.show ()