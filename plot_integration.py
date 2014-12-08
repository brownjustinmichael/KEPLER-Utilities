import yields.yields as yld
import astropy.units as u

from plots.yields import YieldPlot

import matplotlib.pyplot as plt

yr = yld.YieldReader ()

# masses = np.arange (15.0, 30.0, 0.1) * u.solMass
# yr = yld.YieldReader (directory = "yields/y_data/", masses = masses)

fig = plt.figure (figsize = (15,5))
ax = plt.subplot (111)

yp = YieldPlot (yr)
yp.plot (ax, windMultiplier = 1.0, removeIsotopes = ["al26", "k40", "fe60"], imfUpperLimit = 30 * u.solMass, imfLowerLimit = 15 * u.solMass)

plt.show ()