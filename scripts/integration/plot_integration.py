import yields.yields as yld
from yields.integrator import IMFIntegrator
import astropy.units as u
import numpy as np

from kepler_utils.plots.yields import YieldPlot
import kepler_utils.yields.abundances as ab

import matplotlib.pyplot as plt

wyr = yld.YieldReader (windType = True)
wimf = IMFIntegrator (wyr.get_masses ())

masses = np.arange (12.0, 120.0, 0.1) * u.solMass
yr = yld.YieldReader (directory = "yields/y_data/", masses = masses)
imf = IMFIntegrator (masses)

typeiayr = yld.YieldReader (directory = "yields/1a/", masses = [0.0 * u.solMass])
typeiaimf = IMFIntegrator ([0.0 * u.solMass], None, [1.0])

fig = plt.figure (figsize = (15,5))
ax = plt.subplot (111)

b = 0.0

isotopes = yr.isotopes
yields = yr.get_yield (isotopes [0])
wyields = wyr.get_wind (isotopes [0])

for isotope in isotopes [1:]:
    yields += yr.get_yield (isotope)
    wyields += wyr.get_wind (isotope)

print (wyields)
print (imf (yields) + 0.5 * wimf (wyields))

print (imf (yr.get_yield ("o16")) + 0.5 * wimf (wyr.get_wind ("o16")))

first = True
for b, ls, marker in [(0.5, "-", "o"), (1.0, "--", "v")]:
    
    ab1 = (imf + b * wimf).getAbundances (yr + wyr)
    # ab1 = (imf).getAbundances (yr)
    ab2 = typeiaimf.getAbundances (typeiayr)

    x = ab1 ['o16'] / ab.solar ['o16'] * ab1.total
    z = ab1 ['fe56'] / ab.solar ['fe56'] * ab1.total

    y = ab2 ['o16'] / ab.solar ['o16'] * ab2.total
    w = ab2 ['fe56'] / ab.solar ['fe56'] * ab2.total

    a = float ((w - y) / (x + w - y - z))

    print (a)

    a = 1.0

    yp = YieldPlot (yr + wyr + typeiayr, imfIntegrator = a * (imf + b * wimf) + (1. - a) * typeiaimf)
    lines = yp.plot (ax, removeIsotopes = ["al26", "k40", "fe60"], record = "ia+no_winds.dat", ls = ls, names = first, marker = marker, imfUpperLimit = 120 * u.solMass, imfLowerLimit = 30 * u.solMass)
    first = False
    ax.set_color_cycle (None)

ax.set_ylim ((0.1,100))
ax.set_xlim ((0.0, 100))

plt.show ()
