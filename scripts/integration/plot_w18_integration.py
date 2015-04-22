#!/usr/bin/env python

from kepler_utils.yields.yields import YieldReader
from kepler_utils.yields.integrator import IMFIntegrator, Integrator
import astropy.units as u
import numpy as np

from kepler_utils.plots.yields import YieldPlot

import matplotlib.pyplot as plt

import argparse

parser = argparse.ArgumentParser ()
parser.add_argument ('--isotopes', dest = 'elements', action = 'store_false')
parser.set_defaults (elements = True)
parser.add_argument ('--ia', default = False, type = bool)
parser.add_argument ('--upper', default = None)
parser.add_argument ('--lower', default = None)
parser.add_argument ('--output', default = None)
parser.add_argument ('yields', default = "yields/yields_w18/")
parser.add_argument ('plots', nargs = "*")

namespace = parser.parse_args ()

# Grab the wind yields and IMF from the WH07 models; do not include the explosions from this set
wyr = YieldReader (explosions = False)
wimf = IMFIntegrator (wyr.get_masses ())

# Grab the yields and IMF from T's W18 calibrated runs, for which the masses must be stated explicitly because we're missing models; do not include the winds from these sets
# yr = YieldReader (directory = namespace.yields, masses = np.concatenate ([np.arange (9.5, 12.0, 0.5), np.arange (12.0, 30.0, 0.1), np.array ([30., 31., 32., 33., 35., 40., 45., 50., 55., 60., 70., 80., 100., 120.])]) * u.solMass, winds = False)
yr = YieldReader (directory = "yields/y_data_W18_special_13x6_14x9/", masses = np.arange (13.6, 14.9, 0.1) * u.solMass, winds = False)
yr += YieldReader (directory = "yields/y_data_W18_special_15x0_30x0/", masses = np.arange (15.0, 30.0, 0.1) * u.solMass, winds = False)
yr += YieldReader (directory = "yields/y_data_W18_special_31_120/", masses = np.array ([31, 32, 33, 35, 40, 45, 50, 55, 60, 70, 80, 100, 120]) * u.solMass, winds = False)
# yr = YieldReader (directory = "yields/y_data_N20_special_13x6_120/", masses = np.array (list (np.arange (13.6, 30.0, 0.1)) + [31, 32, 33, 35, 40, 45, 50, 55, 60, 70, 80, 100, 120]) * u.solMass, winds = False)
# yr = YieldReader (winds = False)
imf = IMFIntegrator (yr.get_masses ())

# Grab the yields from our Ia model of choice; if using the w7 model, we'll need to be more explicit since it isn't stored in a kepler yield format
typeiayr = YieldReader (directory = "yields/1a/")
# typeiayr = YieldReader (directory = "yields/w7/", keplerYield = False, totalYieldName = "W7")
typeiaimf = Integrator ([1.0])

# For more information, we can print the average production of o16 with full and half wind contributions
print ((imf + wimf) ((yr + 0.0 * wyr).get_yield ("o16")))
print (yr.get_yield ("o16"))
print ((imf + wimf) ((yr + 0.5 * wyr).get_yield ("o16")))

# Specify the upper and lower limits of integration; beware that non-IMF objects (such as Type Ias) do not have masses in the system and thus will be included regardless of these limits
imfUpperLimit = None if namespace.upper is None else (float (namespace.upper) * u.solMass)
imfLowerLimit = None if namespace.lower is None else (float (namespace.lower) * u.solMass)

plots = namespace.plots if len (namespace.plots) != 0 else (0, 1, 2)
elements = namespace.elements
# Initialize the matplotlib figure
fig, axes = plt.subplots (len (plots), 1, figsize = (7, 6 * len (plots)), sharey = True)
if len (plots) == 1:
    axes = (axes,)

if elements:
    xlim = ((0, 22), (20, 42), (40, 62))
else:
    xlim = ((10, 42), (40, 92), (90, 172))


# Iterate over the wind contributions
first = True
for b, ls, marker, color in [(1.0, "-", "o", "g"), (0.5, "--", "s", "r"), (0.0, "-.", "^", "b")]:
    if not elements and b != 1.0:
        continue
    
    final_imf = imf + wimf
    final_yr = yr + b * wyr
    
    if not namespace.ia:
        a = 0.0
    else:
        a, final_imf, final_yr = final_imf.makeSolar ("o16", "fe56", other = typeiaimf, yr = final_yr, other_yr = typeiayr)

    print ("Need a fraction of " + str (a) + " to come from Type Ias.")
    yp = YieldPlot (final_yr, imfIntegrator = final_imf)

    for ax, plot in zip (axes, plots):
        lines1 = yp.plot (ax, removeIsotopes = ["al26", "k40", "fe60"], ls = ls, names = first, marker = marker, ms = 10, alpha = 0.75, relativeIso = "o16", imfUpperLimit = imfUpperLimit, imfLowerLimit = imfLowerLimit, elements = elements, label = "Exp+" + ("Ia+" if a != 0.0 else "") + "%s*Winds" % b)

        ax.set_ylim ((0.03,10))
        ax.set_xlim (xlim [plot])

        xlabel = ax.get_xlabel ()
        ax.set_xlabel ("")

    ax.set_xlabel (xlabel)
    if elements:
        ax.legend (loc = "upper right")

    first = False

if namespace.output:
    plt.savefig (namespace.output)
else:
    plt.show ()
