from kepler_utils.yields.yields import YieldReader
from kepler_utils.yields.integrator import IMFIntegrator, Integrator
import astropy.units as u
import numpy as np

from kepler_utils.plots.yields import YieldPlot

import matplotlib.pyplot as plt

wyr = YieldReader (explosions = False)
wimf = IMFIntegrator (wyr.get_masses ())

yr = YieldReader (directory = "yields/y_data/", masses = np.arange (15.0, 30.0, 0.1) * u.solMass, winds = False)
imf = IMFIntegrator (yr.get_masses ())

typeiayr = YieldReader (directory = "yields/1a/")
typeiaimf = Integrator ([1.0])

print (typeiaimf (typeiayr.get_yield ("fe56")))

fig = plt.figure (figsize = (15,5))
ax = plt.subplot (111)

print (imf (yr.get_yield ("o16")) + wimf (wyr.get_yield ("o16")))

first = True
for b, ls, marker in [(0.0, "-", "o"), (1.0, "--", "v")]:
    star_imf = imf + wimf
    star_yr = yr + b * wyr
    
    a, final_imf, final_yr = star_imf.makeSolar ("o16", "fe56", other = typeiaimf, yr = star_yr, other_yr = typeiayr)
    # final_imf = star_imf
    # final_yr = star_yr
    
    print ("Need a fraction of " + str (a) + " to come from massive stars.")

    yp = YieldPlot (final_yr, imfIntegrator = final_imf)
    lines = yp.plot (ax, removeIsotopes = ["al26", "k40", "fe60"], record = "ia+no_winds.dat", ls = ls, names = first, marker = marker)#, imfUpperLimit = 120 * u.solMass, imfLowerLimit = 15 * u.solMass)
    first = False
    ax.set_color_cycle (None)

ax.set_ylim ((0.1,100))
ax.set_xlim ((0.0, 100))

plt.show ()
