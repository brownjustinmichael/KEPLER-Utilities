import sys
from math import pi

import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.constants import sigma_sb

import records.cnv as cnv

# sys.argv.append ("/Users/justinbrown/Codes/kepler/run/87a/87a.cnv")
sys.argv.append ("/Users/justinbrown/Codes/kepler/run/87a/87h1..cnv")
sys.argv.append ("/Users/justinbrown/Codes/kepler/run/87a/87h.9.cnv")
sys.argv.append ("/Users/justinbrown/Codes/kepler/run/87a/87h.8.cnv")
sys.argv.append ("/Users/justinbrown/Codes/kepler/run/87a/87h.7.cnv")
sys.argv.append ("/Users/justinbrown/Codes/kepler/run/87a/87h.6.cnv")
sys.argv.append ("/Users/justinbrown/Codes/kepler/run/87a/87h.5.cnv")
sys.argv.append ("/Users/justinbrown/Codes/kepler/run/87a/87h.4.cnv")
sys.argv.append ("/Users/justinbrown/Codes/kepler/run/87a/87h.3.cnv")
sys.argv.append ("/Users/justinbrown/Codes/kepler/run/87a/87h.2.cnv")
sys.argv.append ("/Users/justinbrown/Codes/kepler/run/87a/87h.1.cnv")
sys.argv.append ("none")

# Read in the KEPLER cnv output file into cnv_record
if (len (sys.argv) < 2):
    print ("Usage: plot_hr cnv_file (output_file)")
    exit ()

records = [cnv.CNVFile (arg) for arg in sys.argv [1:-1]]

fig = plt.figure ()
ax = plt.subplot (111)

for cnv_record in records:
    radii = u.Quantity ([model [-1] for model in cnv_record ['rncoord']])
    lum = u.Quantity (cnv_record ['xlum_cnv'])
    temp = (lum / 4 / pi / radii ** 2 / sigma_sb) ** (0.25)
    mass = u.Quantity ([model [-1] for model in cnv_record ['xmcoord']])
    points = []
    times = cnv_record ['timesec']
    timezero = u.Quantity (0.0, u.yr)
    timediff = u.Quantity (1e4, u.yr)
    for i in range (len (times)):
        if times [i] > timezero + timediff:
            timezero = times [i]
            points.append (i)
    ax.plot (temp [points [0]:].to (u.K), lum.to (u.solLum) [points [0]:])
    colors = (-times.to (u.year) [points] + times.to (u.year) [-1]) / u.year
    colors [colors < 0.1] = 0.1
    ax.scatter (temp [points].to (u.K), lum.to (u.solLum) [points], c = np.log (colors), s = 40, lw = 0.0)
plt.gca().invert_xaxis()
ax.set_xlim (right = 0.)
ax.set_ylim ((np.min (lum.to (u.solLum) [points]).value * 0.5, np.max (lum.to (u.solLum) [points]).value / 0.5))
ax.set_yscale ('log')

ax.set_xlabel ("$T_\mathrm{eff}$ (K)")
ax.set_ylabel ("$L$ (solar luminosities)")

# if (len (sys.argv) > 2):
    # plt.savefig (sys.argv [2])
# else:
    # Plot the result
plt.show ()
