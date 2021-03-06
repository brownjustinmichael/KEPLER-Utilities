#!/usr/bin/env python

from kepler_utils.yields.yields import YieldReader
from kepler_utils.yields.integrator import IMFIntegrator, Integrator
import astropy.units as u
import numpy as np

from kepler_utils.plots.yields import YieldPlot

import matplotlib.pyplot as plt

from kepler_utils.plots.parser import PlotArgumentParser

parser = PlotArgumentParser (False)
parser.add_argument ('--isotopes', dest = 'elements', action = 'store_false')
parser.set_defaults (elements = True)
parser.add_argument ('--ia', default = None)
parser.add_argument ('--upper', default = None)
parser.add_argument ('--lower', default = None)
parser.add_argument ('yields', default = "yields/yields_w18/")
parser.add_argument ('plots', nargs = "*", type = int)

namespace = parser.parse_args ()

# Grab the wind yields and IMF from the WH07 models; do not include the explosions from this set
wyr = YieldReader.from_directory (explosions = False)
wimf = IMFIntegrator.from_yieldreader (wyr)

# Grab the yields and IMF from T's W18 calibrated runs, for which the masses must be stated explicitly because we're missing models; do not include the winds from these sets
yr = YieldReader.from_directory (directory = namespace.yields, winds = False)
imf = IMFIntegrator.from_yieldreader (yr)

# Grab the yields from our Ia model of choice; if using the w7 model, we'll need to be more explicit since it isn't stored in a kepler yield format
typeiayr = None
if namespace.ia == "ia":
    typeiayr = YieldReader.from_file ("yields/1a/1pi1a.yield", 0.0 * u.solMass)
if namespace.ia == "w7":
    typeiayr = YieldReader.from_file ("yields/w7/iwamoto.dat", keplerYield = False, totalYieldName = "W7")
    
print (typeiayr)

if typeiayr is not None:
    typeiaimf = Integrator (u.Quantity ([0.0 * u.solMass]), [1.0])

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
    
    if typeiayr is None:
        a = 1.0
    else:
        a, final_imf, final_yr = final_imf.makeSolar ("o16", "fe56", other = typeiaimf, yr = final_yr, other_yr = typeiayr)

    print ("Need a fraction of " + str (1 - a) + " to come from Type Ias.")
    yp = YieldPlot (final_yr, imfIntegrator = final_imf)

    for ax, plot in zip (axes, plots):
        lines1 = yp.plot (ax, removeIsotopes = ["al26", "k40", "fe60"], ls = ls, names = first, marker = marker, ms = 10, alpha = 0.75, relativeIso = "o16", imfUpperLimit = imfUpperLimit, imfLowerLimit = imfLowerLimit, elements = elements, label = "Exp+" + ("Ia+" if a != 0.0 else "") + "%s*Winds" % b, moveLabels = {"Na": (0.5, 0.8), "Ca": (1, 1.0), "B": (0.5, 0.0), "Ar": (1, 0.0), "Sc": (-1.0, -10.0), "Mn": (0.0, -10.0), "Fe": (-1.5, 0.0), "Co": (-0.0, -10.0), "Ni": (-0.5, 0.0), "Ga": (-1.0, 0.0), "Se": (1.0, 0.0), "Rb": (0.5, 0.0), "As": (-0.5, 0.0)})

        ax.set_ylim ((0.03,10))
        ax.set_xlim (xlim [plot])

        xlabel = ax.get_xlabel ()
        ax.set_xlabel ("")

    ax.set_xlabel (xlabel)
    if elements:
        ax.legend (loc = "upper right")

    first = False

parser.exit (fig)
