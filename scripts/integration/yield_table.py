from sys import argv

from kepler_utils.yields.yields import YieldReader, Isotope
from kepler_utils.yields.integrator import IMFIntegrator, Integrator
import astropy.units as u
import numpy as np

from kepler_utils.plots.yields import YieldPlot

if len (argv) > 1:
    masses = argv [1:]
else:
    masses = [25]

w18yr = YieldReader (directory = "yields/y_data_W18_special_15x0_30x0/", masses = np.arange (15.0, 30.0, 0.1) * u.solMass, winds = False)
whyr = YieldReader (explosions = False)

table = None
indices = []
for mass in masses:
    index = np.argmax (w18yr.masses >= float (mass) * u.solMass)
    ejecta = w18yr.yields.iloc [index : index + 1].reset_index () [[str (iso) for iso in w18yr.isotopes]]
    index = np.argmax (whyr.masses >= float (mass) * u.solMass)
    winds = whyr.yields.iloc [index : index + 1].reset_index () [[str (iso) for iso in whyr.isotopes]]
    
    for df in [ejecta, winds]:
        if table is None:
            table = df
        else:
            table = table.append (df, ignore_index = True)
    indices += [mass + " \Msun Ejecta", mass + " \Msun Winds"]

isos = [Isotope (col) for col in table.columns]
isos.sort ()
table = table [[str (iso) for iso in isos]]
table.columns = [iso.getLabel () for iso in isos]
table = table.fillna (0.0)
table.index = indices

print (table.transpose ().to_latex (escape = False))