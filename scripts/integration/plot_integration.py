#!/usr/bin/env python

from kepler_utils.yields.yields import YieldReader
from kepler_utils.yields.integrator import IMFIntegrator, Integrator
import astropy.units as u
import numpy as np

from kepler_utils.plots.yields import YieldPlot

import matplotlib.pyplot as plt

wyr = YieldReader (explosions = False)
wimf = IMFIntegrator (wyr.get_masses ())

yr = YieldReader (directory = "yields/yields_w18/", masses = np.concatenate ([np.arange (9.5, 12.0, 0.5), np.arange (12.0, 30.0, 0.1), np.array ([30., 31., 32., 33., 35., 40., 45., 50., 55., 60., 70., 80., 100., 120.])]) * u.solMass, winds = False)
imf = IMFIntegrator (yr.get_masses ())

typeiayr = YieldReader (directory = "yields/1a/")
typeiaimf = Integrator ([1.0])

fig, ax2 = plt.subplots (1, 1, figsize = (7, 6))

isos = yr.isotopes

fyr = yr + wyr
yields = (fyr).get_yield (isos [0])
for iso in isos [1:]:
    yields += (fyr).get_yield (iso)

print ((imf + wimf) (yields))
print ((0.5 * (imf + wimf)) ((fyr).get_yield ("o16")))

first = True
for b, ls, marker, color in [(1.0, "-", "o", "g"), (0.5, "--", "s", "r"), (0.0, "-.", "p", "b")]:
    # if b != 1.0:
    #     continue
    
    final_imf = imf + wimf
    final_yr = yr + b * wyr
    
    a, final_imf, final_yr = final_imf.makeSolar ("o16", "fe56", other = typeiaimf, yr = final_yr, other_yr = typeiayr)
    
    print ("Need a fraction of " + str (a) + " to come from massive stars.")

    yp = YieldPlot (final_yr, imfIntegrator = final_imf)
    # lines1 = yp.plot (ax1, removeIsotopes = ["al26", "k40", "fe60"], record = "ia+no_winds.dat", ls = ls, names = first, marker = marker, ms = 10, alpha = 0.75, relativeIso = "o16")
    #
    # ax1.set_ylim ((0.1,10))
    # ax1.set_xlim ((10, 51))
    #
    # ax1.set_xlabel ("")

    lines2 = yp.plot (ax2, removeIsotopes = ["al26", "k40", "fe60"], record = "ia+no_winds.dat", ls = ls, names = first, marker = marker, ms = 10, alpha = 0.75, relativeIso = "o16", label = "Exp+Ia+%s*Winds" % b)

    ax2.set_ylim ((0.1,10))
    ax2.set_xlim ((37, 90))
    
    ax2.legend (loc = "upper right")
    
    first = False

# plt.savefig ("full_ia.pdf")

plt.show ()
