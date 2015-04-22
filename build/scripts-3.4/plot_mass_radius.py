#!/opt/local/bin/python

import sys
from math import pi

import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.constants import sigma_sb

import kepler_utils.records.cnv as cnv

# Read in the KEPLER cnv output file into cnv_record
if (len (sys.argv) < 2):
    print ("Usage: plot_hr cnv_file (output_file)")
    exit ()

records = [cnv.CNVFile (arg) for arg in sys.argv [1:]]

fig, axes = plt.subplots (2, 1, sharex = True)

for cnv_record in records:
    radii = u.Quantity ([model ['rncoord'] [-1] for model in cnv_record], u.cm)
    mass = u.Quantity ([model ['xmcoord'] [-1] for model in cnv_record], u.g)
    hecoremass = u.Quantity ([model ["xmcoord"] [int (model ["inuc"] [np.argmax (model ["nuc"] [3:]) + 3])] if len (model ["nuc"]) > 3 else 0 for model in cnv_record], u.g)
    points = []
    times = cnv_record ['timesec']
    # timezero = u.Quantity (0.0, u.yr)
    # timediff = u.Quantity (1e4, u.yr)
    # for i in range (len (times)):
    #     if times [i] > timezero + timediff:
    #         timezero = times [i]
    #         points.append (i)
    axes [0].plot (times.to (u.yr) [times > 1.1e7 * u.yr], (mass).to (u.solMass) [times > 1.1e7 * u.yr])
    axes [0].plot (times.to (u.yr) [times > 1.1e7 * u.yr], (hecoremass).to (u.solMass) [times > 1.1e7 * u.yr])
    axes [0].plot (times.to (u.yr) [times > 1.1e7 * u.yr], (mass - hecoremass).to (u.solMass) [times > 1.1e7 * u.yr])
    axes [1].plot (times.to (u.yr) [times > 1.1e7 * u.yr], radii.to (u.solRad) [times > 1.1e7 * u.yr])
# ax.set_xscale ('log')
# ax.set_yscale ('log')

axes [1].set_xlabel ("Time (Years)")
axes [0].set_ylabel ("Mass (Solar Masses)")
axes [1].set_ylabel ("Radius (Solar Radii)")

# if (len (sys.argv) > 2):
    # plt.savefig (sys.argv [2])
# else:
    # Plot the result
plt.show ()
